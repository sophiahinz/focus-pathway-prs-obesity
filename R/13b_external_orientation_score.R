# ==============================================================================
# 13b_external_orientation_score.R
# Externally oriented scores, step 2: pathway SNP-sum-scores with SNP orientation
# taken from the external BMI GWAS lookup of 13a, for all 54 pathways
# (Results).
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================

#=================================================================================
# Variants computed:
#   - 01_Main_LDpruned_All, score_mode = "oriented"
#   - 05_Sensitivity_BetaWeighted_LDpruned, score_mode = "beta_weighted"
#     (external counterpart of the internal beta-weighted score)
# Not computed:
#   - "raw" and ridge/elastic-net-weighted score modes (no external analogue)
#   - 07_Sensitivity_KEGGonly_LDpruned (database axis, orientation-independent)
#   - SKAT (orientation-invariant by construction)
#   - 06_Sensitivity_LoF_VEP_LDpruned: the LoF-burden orientation is derived
#     from the VEP-annotated LoF allele only (create_lof_orientation_table())
#     and does not depend on the univariate regression, so it is identical
#     under external orientation.
#
# Architecture: no PLINK/QC/HWE rerun. The scoring/model functions and the
# orchestration loops are taken verbatim from 07_pathway_pipeline.R and read
# the cached, QC'd, LD-pruned per-pathway datasets. The only substituted
# input is the per-pathway "univariate regression results" table
# (SNP / Category / Estimate_Add / OR_Add / Pvalue_Add) read by
# run_snp_sumscore(): instead of the internal univariate logistic
# regression, it is derived from the external lookup (13a) via
# create_external_orientation_table(), gated by the pre-specified coverage
# threshold (>=50% AND >=10 SNPs).
#
# The model is fitted on the full sample: external orientation carries no
# FoCus-internal information, so there is no leakage risk (see section 1b in
# the main loop).
#
# Outputs are written to <FOCUS_RESULTS_DIR>/OutExt, separate from the
# internal-orientation outputs. Cached datasets are read from
# <FOCUS_RESULTS_DIR>/OutPW; only the small Global_TrainTest_Split.rds file
# is copied into OutExt, because the reused functions read it via
# output_file().
#
# Optional: per-pathway results and session info are written to Appendix B of
# the methods/results log, if that file exists. Pathways already present in
# the output accumulator CSV with the current SCRIPT_VERSION are skipped on a
# re-run.
#=================================================================================

SCRIPT_VERSION <- "Phase6_v1.3_Stage1"   # version label stored in the result files (read by 13c)

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------

sensitivity_base <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Sensitivity_External_SNP_Orientation")
main_project_base <- Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus")

# Output location (writes) - the reused functions below call output_file()
# for every new result file, so this is the only output_dir in scope.
# Short folder name because of the Windows MAX_PATH limit (260 characters).
results_dir <- Sys.getenv("FOCUS_RESULTS_DIR", unset = "path/to/results")
output_dir <- file.path(results_dir, "OutExt")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_file <- function(...) {
  file.path(output_dir, paste0(...))
}

# Internal-orientation pipeline output location (reads only) - cached
# QC'd/LD-pruned datasets and the per-pathway Global_TrainTest_Split.rds.
pipeline_output_dir <- file.path(results_dir, "OutPW")

pipeline_output_file <- function(...) {
  file.path(pipeline_output_dir, paste0(...))
}

# Genome-wide population-stratification PCs (06_global_PCA.R).
genomewide_pcs_file <- file.path(main_project_base, "Global_PCA", "Global_Genomewide_PCs.rds")

# ---------------------------------------------------------------------------
# LIBRARIES
# ---------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(pROC)
library(caret)
library(ggplot2)
library(broom)
library(car)
library(vcd)
library(margins)
library(pscl)
library(ggeffects)
library(grid)
library(stringr)
# (glmnet is not loaded: the ridge/elastic-net branch of run_snp_sumscore()
# is never executed here, because analysis_plan does not use
# score_mode %in% c("ridge_weighted", "elasticnet_weighted").)

# ---------------------------------------------------------------------------
# EXTERNAL ORIENTATION LOOKUP (output of 13a) + coverage threshold +
# create_external_orientation_table() - copied verbatim from
# 13a_external_orientation_prepare.R so that this script is self-contained.
# ---------------------------------------------------------------------------

external_lookup_rds <- file.path(
  sensitivity_base, "External_SNP_Data",
  "External_Orientation_GIANT_Yengo2018_Locke2015.rds"
)

if (!file.exists(external_lookup_rds)) {
  stop(
    "External orientation lookup not found: ", external_lookup_rds,
    ". Run 13a_external_orientation_prepare.R first."
  )
}

external_lookup <- readRDS(external_lookup_rds)

# Dual coverage threshold (percentage AND absolute-count floor)
MIN_COVERAGE_PCT <- 0.50
MIN_SNP_COUNT    <- 10

pathway_passes_coverage_threshold <- function(n_matched, n_total_main_model) {
  (n_matched / n_total_main_model) >= MIN_COVERAGE_PCT & n_matched >= MIN_SNP_COUNT
}

# Genuine-rsID filter (verbatim from 13a)
is_genuine_rsid <- function(x) grepl("^rs[0-9]+$", x)

# NO-OP STUB: run_multivariate_sumscore_analysis() (reused verbatim below)
# unconditionally calls plot_analysis_flow_summary() to draw/save a flowchart
# plot. Plotting functions are not included in this script; the stub skips
# the plot. Its return value (analysis_flow_plot) is not used downstream.
plot_analysis_flow_summary <- function(
    summary_df, pathway_name, outfile_prefix, score_mode, analysis_config, show_plot = TRUE) {
  cat("  (plot skipped: plot_analysis_flow_summary is a stub)\n")
  invisible(NULL)
}

# Verbatim from 13a - maps external Effect_Allele/Beta onto the
# pathway's own matrix-coded allele (Matrix_Allele1/Matrix_Allele2, from that
# pathway's .bim file), flags palindromic SNPs, and derives whether the
# external Beta needs sign-inversion to be expressed "per copy of the
# matrix-counted allele" - i.e. on the same footing as the internal
# univariate regression's Estimate_Add.
create_external_orientation_table <- function(external_lookup, bim_table, matrix_snps) {

  bim_table_rsid_only <- bim_table %>%
    mutate(Is_Genuine_RsID = is_genuine_rsid(SNP))

  n_total       <- nrow(bim_table_rsid_only)
  n_non_rsid    <- sum(!bim_table_rsid_only$Is_Genuine_RsID)
  cat(sprintf(
    "  SNP-ID QC: %d/%d (%.1f%%) are non-rsID pseudo-IDs (e.g. GSA-rs...) and are excluded from external matching by construction.\n",
    n_non_rsid, n_total, 100 * n_non_rsid / n_total
  ))

  bim_table_rsid_only %>%
    filter(Is_Genuine_RsID) %>%
    left_join(external_lookup, by = "SNP") %>%
    filter(SNP %in% matrix_snps) %>%
    mutate(
      Matched = !is.na(Effect_Allele),

      Is_Palindromic = Matched & (
        (Effect_Allele == "A" & Other_Allele == "T") |
        (Effect_Allele == "T" & Other_Allele == "A") |
        (Effect_Allele == "C" & Other_Allele == "G") |
        (Effect_Allele == "G" & Other_Allele == "C")
      ),

      Allele_Match = case_when(
        !Matched ~ "No_external_match",
        Is_Palindromic ~ "Palindromic_excluded",
        Effect_Allele == Matrix_Allele2 ~ "Effect_is_matrix_counted_allele",
        Effect_Allele == Matrix_Allele1 ~ "Effect_is_matrix_other_allele",
        TRUE ~ "No_allele_match"
      ),

      Invert_For_External_Score = case_when(
        Allele_Match == "Effect_is_matrix_counted_allele" ~ FALSE,
        Allele_Match == "Effect_is_matrix_other_allele" ~ TRUE,
        TRUE ~ NA
      )
    )
}

# ---------------------------------------------------------------------------
# Build the external-orientation equivalent of the internal univariate
# regression results table. run_snp_sumscore() (below, reused verbatim) does
# not know or care HOW this table was produced - it just needs the same
# columns (SNP, Category, Estimate_Add, OR_Add, Pvalue_Add) at the same file
# path it always reads (via make_output_name(..., analysis_configs$main_ld)).
#
# NOTE ON SCALE (flagged for the Methods write-up, not a bug): Estimate_Add
# here is the external GIANT/Yengo beta (continuous-BMI-SD scale), whereas
# the internal Estimate_Add is a log-OR from obesity ~ genotype (binary
# outcome). For "oriented" scoring only the SIGN matters, so this is
# unaffected. For "beta_weighted" scoring the two are not on the same
# numeric scale - this is expected/desired (the external-beta-weighted score
# is deliberately weighted by published, non-FoCus-derived effect sizes) but
# means its magnitude/interpretation is not directly comparable to the
# internal beta-weighted score and should be described as such in Results.
# ---------------------------------------------------------------------------

build_external_regression_table <- function(pathway_name, main_ld_bim_path, matrix_snps) {

  if (!file.exists(main_ld_bim_path)) {
    stop("Main-model .bim file not found: ", main_ld_bim_path)
  }

  bim_table <- read.table(main_ld_bim_path, header = FALSE, stringsAsFactors = FALSE)
  colnames(bim_table) <- c("CHR", "SNP", "CM", "BP", "Allele1", "Allele2")

  bim_table <- bim_table %>%
    transmute(
      SNP = trimws(as.character(SNP)),
      Matrix_Allele1 = toupper(trimws(as.character(Allele1))),
      Matrix_Allele2 = toupper(trimws(as.character(Allele2)))
    ) %>%
    distinct(SNP, .keep_all = TRUE)

  ext_table <- create_external_orientation_table(external_lookup, bim_table, matrix_snps)

  usable <- ext_table %>%
    filter(Allele_Match %in% c("Effect_is_matrix_counted_allele", "Effect_is_matrix_other_allele")) %>%
    mutate(
      Estimate_Add = ifelse(Invert_For_External_Score, -Beta, Beta),
      OR_Add       = exp(Estimate_Add),
      Pvalue_Add   = P,
      Category     = ifelse(Estimate_Add > 0, "Risk", "Protective")
    )

  list(
    results = usable %>% dplyr::select(SNP, Category, Estimate_Add, OR_Add, Pvalue_Add),
    n_total_matrix_snps = length(matrix_snps),
    n_matched = nrow(usable),
    n_no_external_match = sum(ext_table$Allele_Match == "No_external_match"),
    n_palindromic_excluded = sum(ext_table$Allele_Match == "Palindromic_excluded"),
    n_no_allele_match = sum(ext_table$Allele_Match == "No_allele_match")
  )
}

# ---------------------------------------------------------------------------
# AUTO-LOGGING into Methods_Results_Log_External_SNP_Orientation.md,
# Appendix B - idempotent upsert-by-pathway (same convention as 11 / 13a).
# ---------------------------------------------------------------------------

methods_log_path <- file.path(results_dir, "Methods_Results_Log.md")

upsert_marker_table <- function(md_path, begin_marker, end_marker, key_col, new_rows_df) {

  if (!file.exists(md_path)) {
    warning("Methods log not found, skipping auto-log: ", md_path)
    return(invisible(NULL))
  }

  lines <- readLines(md_path, warn = FALSE, encoding = "UTF-8")
  begin_idx <- grep(begin_marker, lines, fixed = TRUE)
  end_idx   <- grep(end_marker, lines, fixed = TRUE)

  if (length(begin_idx) != 1 || length(end_idx) != 1 || end_idx <= begin_idx) {
    warning("Could not find marker pair ", begin_marker, " / ", end_marker, " - skipping auto-log.")
    return(invisible(NULL))
  }

  block <- lines[(begin_idx + 1):(end_idx - 1)]
  block <- block[nchar(trimws(block)) > 0]

  header_line <- NULL
  sep_line <- NULL
  existing_rows <- character(0)

  for (ln in block) {
    if (is.null(header_line) && startsWith(trimws(ln), "|") && grepl(key_col, ln, fixed = TRUE) == FALSE && grepl("---", ln, fixed = TRUE) == FALSE && is.null(header_line)) {
      # first non-separator pipe-line encountered before we've set header is ambiguous;
      # instead identify header/sep by position (first two pipe-lines).
    }
  }

  pipe_lines <- block[startsWith(trimws(block), "|")]
  if (length(pipe_lines) >= 2) {
    header_line <- pipe_lines[1]
    sep_line <- pipe_lines[2]
    existing_rows <- if (length(pipe_lines) > 2) pipe_lines[3:length(pipe_lines)] else character(0)
  } else {
    stop("Could not parse existing table header/separator inside marker block.")
  }

  col_names <- trimws(strsplit(gsub("^\\||\\|$", "", header_line), "\\|")[[1]])
  key_pos <- which(col_names == key_col)
  if (length(key_pos) != 1) stop("Key column '", key_col, "' not found in table header.")

  existing_keys <- sapply(existing_rows, function(r) {
    cells <- trimws(strsplit(gsub("^\\||\\|$", "", r), "\\|")[[1]])
    if (length(cells) >= key_pos) cells[key_pos] else NA_character_
  })

  new_row_strings <- apply(new_rows_df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  new_keys <- new_rows_df[[key_col]]

  keep_existing <- !(existing_keys %in% new_keys)
  merged_rows <- c(existing_rows[keep_existing], new_row_strings)

  new_block <- c(header_line, sep_line, merged_rows)

  new_lines <- c(
    lines[1:begin_idx],
    new_block,
    lines[end_idx:length(lines)]
  )

  writeLines(new_lines, md_path, useBytes = TRUE)
  invisible(NULL)
}

overwrite_marker_block <- function(md_path, begin_marker, end_marker, content_lines) {

  if (!file.exists(md_path)) {
    warning("Methods log not found, skipping auto-log: ", md_path)
    return(invisible(NULL))
  }

  lines <- readLines(md_path, warn = FALSE, encoding = "UTF-8")
  begin_idx <- grep(begin_marker, lines, fixed = TRUE)
  end_idx   <- grep(end_marker, lines, fixed = TRUE)

  if (length(begin_idx) != 1 || length(end_idx) != 1 || end_idx <= begin_idx) {
    warning("Could not find marker pair ", begin_marker, " / ", end_marker, " - skipping auto-log.")
    return(invisible(NULL))
  }

  new_lines <- c(
    lines[1:begin_idx],
    content_lines,
    lines[end_idx:length(lines)]
  )

  writeLines(new_lines, md_path, useBytes = TRUE)
  invisible(NULL)
}

log_phase6_sessioninfo <- function() {
  session_lines <- c(
    "",
    paste0("- **Last run:** ", format(Sys.time(), "%Y-%m-%d %H:%M")),
    paste0("- **R version:** ", R.version.string),
    paste0("- **Platform:** ", R.version$platform),
    paste0("- **Package versions:** data.table ", as.character(packageVersion("data.table")),
           ", dplyr ", as.character(packageVersion("dplyr")),
           ", pROC ", as.character(packageVersion("pROC")),
           ", caret ", as.character(packageVersion("caret"))),
    paste0("- **Script version:** ", SCRIPT_VERSION),
    ""
  )
  overwrite_marker_block(
    methods_log_path,
    "<!-- AUTOLOG-EXTORIENT-PHASE6-SESSIONINFO:BEGIN -->",
    "<!-- AUTOLOG-EXTORIENT-PHASE6-SESSIONINFO:END -->",
    session_lines
  )
}

# ---------------------------------------------------------------------------
# MANUAL COVARIATE DECISION (copied verbatim from 07_pathway_pipeline.R).
# Fixed global object consumed by run_multivariate_sumscore_analysis() for
# every pathway/configuration; defined after review of the covariate-precheck
# decision tables of the internal-orientation run.
# ---------------------------------------------------------------------------

manual_covariate_decisions <- list(
  default = list(
    keep = c(
      "gender_cat",
      "smoking_cat",
      "school_cert_cat",
      "graduation_cat",
      "sample_age",
      "sports_activity_min_week",
      "everyday_activity_min_week",
      "energy_kcal",
      "gPC1",
      "gPC2",
      "gPC3"
    ),
    remove = c(
      "Batch"
    ),
    comment = paste0(
      "Batch only as sensitivity analysis; ",
      "school_cert_cat and graduation_cat retained on content grounds ",
      "despite Cramer's V REVIEW."
    )
  )
)

# ---------------------------------------------------------------------------
# PATHWAY LIST - all 54 screened pathways (not restricted to the robust
# candidates). Read from the tracking table, "Pathway" column; the
# Robust_Candidate flag is ignored here.
# ---------------------------------------------------------------------------

pathway_tracking_file <- pipeline_output_file("AllPathways_Results_Tracking_Table.csv")
pathway_tracking <- fread(pathway_tracking_file, sep = ";")
pathway_list <- unique(trimws(as.character(pathway_tracking$Pathway)))
pathway_list <- pathway_list[nchar(pathway_list) > 0]

cat(sprintf("\nPathway list loaded: %d pathways.\n", length(pathway_list)))
print(pathway_list)


# ---------------------------------------------------------------------------
# REUSED (verbatim, unmodified) FUNCTIONS from
# 07_pathway_pipeline.R
# Only the input they are fed (analysis_configs$main_ld's regression_file,
# via build_external_regression_table() above) differs from the internal run.
# ---------------------------------------------------------------------------

# ---- clean_pathway_display_name / sanitize_pathway_filename ----
clean_pathway_display_name <- function(name) {
  name <- gsub(" - Homo sapiens \\(human\\)", "", name)
  name <- gsub(" \\(KEGG erweitert\\)", "", name)
  name
}

# ---- ----
sanitize_pathway_filename <- function(name) {
  name <- gsub("/", "-", name, fixed = TRUE)
  name <- gsub('[\\\\:*?"<>|]', "-", name)
  name <- gsub("\\s+", "_", trimws(name))
  name <- gsub("_+", "_", name)
  name
}

# ---- make_analysis_suffix / make_output_name ----
make_analysis_suffix <- function(config, score_mode = NULL) {
  
  suffix <- config$label
  
  if (!is.null(score_mode)) {
    suffix <- paste(
      suffix,
      tools::toTitleCase(score_mode),
      sep = "_"
    )
  }
  
  suffix
}

make_output_name <- function(
    pathway_name,
    analysis_type,
    analysis_name,
    result_name,
    config,
    score_mode = NULL,
    ext = "csv") {
  
  output_file(
    "Pathway_",
    pathway_name,
    "_",
    analysis_type,
    "_",
    analysis_name,
    "_",
    result_name,
    "_",
    make_analysis_suffix(config, score_mode),
    ".",
    ext
  )
}


# ---------------------------------------------------------------------------
# PATH-SAFE SHORT NAMES (Windows MAX_PATH = 260 characters) - ONLY for output
# file names. pathway_name itself remains unchanged everywhere else (data access under
# Outputs_Pathway/, Global_TrainTest_Split.rds, content of the result CSVs,
# console output) - the abbreviation is used exclusively when
# building the output file path via make_output_name(), for the
# one pathway whose full name in combination with the longest
# result/configuration suffix (Demographics_Covariates_CategoricalSummary,
# sensitivity_betaweighted / beta_weighted) exceeds the 260-character path limit under
# OutExt/ by 1 character (261 characters).
# ---------------------------------------------------------------------------

path_safe_pathway_names <- list(
  "Mitochondrial_Iron_Sulfur_Cluster_Biogenesis" = "Mito_FeS_Cluster_Biogenesis"
)

path_safe_name <- function(pathway_name) {
  if (pathway_name %in% names(path_safe_pathway_names)) {
    return(path_safe_pathway_names[[pathway_name]])
  }
  pathway_name
}

# Wrap make_output_name() (above, template-verbatim): identical logic,
# but with path-safe short name instead of the full pathway_name when building
# the file path.
make_output_name_original <- make_output_name

make_output_name <- function(
    pathway_name,
    analysis_type,
    analysis_name,
    result_name,
    config,
    score_mode = NULL,
    ext = "csv") {

  make_output_name_original(
    pathway_name = path_safe_name(pathway_name),
    analysis_type = analysis_type,
    analysis_name = analysis_name,
    result_name = result_name,
    config = config,
    score_mode = score_mode,
    ext = ext
  )
}

# ---- score-cutoff/metric helper functions ----
library(dplyr)
library(caret)
library(pROC)

#===================================================================================
# HELPER FUNCTIONS
#===================================================================================

calc_metrics <- function(pred, obs, cutoff) {
  
  pred <- as.numeric(pred)
  obs  <- as.numeric(obs)
  
  valid <- complete.cases(pred, obs)
  
  pred <- pred[valid]
  obs  <- obs[valid]
  
  if (
    length(pred) == 0 ||
    !is.finite(cutoff)
  ) {
    
    return(
      data.frame(
        TN = NA_integer_,
        FP = NA_integer_,
        FN = NA_integer_,
        TP = NA_integer_,
        Sensitivity = NA_real_,
        Specificity = NA_real_,
        PPV = NA_real_,
        NPV = NA_real_,
        Accuracy = NA_real_
      )
    )
  }
  
  pred_class <- ifelse(
    pred >= cutoff,
    1,
    0
  )
  
  cm <- table(
    Predicted = factor(
      pred_class,
      levels = c(0, 1)
    ),
    Observed = factor(
      obs,
      levels = c(0, 1)
    )
  )
  
  TN <- unname(cm[1, 1])
  FP <- unname(cm[2, 1])
  FN <- unname(cm[1, 2])
  TP <- unname(cm[2, 2])
  
  data.frame(
    TN = TN,
    FP = FP,
    FN = FN,
    TP = TP,
    
    Sensitivity = ifelse(
      (TP + FN) > 0,
      TP / (TP + FN),
      NA_real_
    ),
    
    Specificity = ifelse(
      (TN + FP) > 0,
      TN / (TN + FP),
      NA_real_
    ),
    
    PPV = ifelse(
      (TP + FP) > 0,
      TP / (TP + FP),
      NA_real_
    ),
    
    NPV = ifelse(
      (TN + FN) > 0,
      TN / (TN + FN),
      NA_real_
    ),
    
    Accuracy = ifelse(
      sum(cm) > 0,
      (TP + TN) / sum(cm),
      NA_real_
    )
  )
}


calc_score_cutoff <- function(
    prob_cutoff,
    beta0,
    beta1) {
  
  if (
    !is.finite(prob_cutoff) ||
    prob_cutoff <= 0 ||
    prob_cutoff >= 1 ||
    !is.finite(beta0) ||
    !is.finite(beta1) ||
    beta1 == 0
  ) {
    return(NA_real_)
  }
  
  (
    qlogis(prob_cutoff) - beta0
  ) / beta1
}


get_score_cutoff_direction <- function(beta1) {
  
  if (!is.finite(beta1) || beta1 == 0) {
    return(NA_character_)
  }
  
  if (beta1 > 0) {
    return("Obesity predicted when score >= cutoff")
  }
  
  "Obesity predicted when score <= cutoff"
}


safe_roc <- function(obs, pred) {
  
  valid <- complete.cases(obs, pred)
  
  obs  <- obs[valid]
  pred <- pred[valid]
  
  if (
    length(obs) == 0 ||
    length(unique(obs)) < 2 ||
    length(unique(pred)) < 2
  ) {
    return(NULL)
  }
  
  tryCatch(
    pROC::roc(
      response = obs,
      predictor = pred,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    ),
    error = function(e) NULL
  )
}


get_n_score_snps <- function(df) {
  
  n_score_snps <- attr(
    df,
    "n_score_snps"
  )
  
  if (
    is.null(n_score_snps) ||
    length(n_score_snps) == 0
  ) {
    return(NA_integer_)
  }
  
  as.integer(n_score_snps[1])
}

# ---- run_snp_sumscore() ----
run_snp_sumscore <- function(config_name, config, score_mode) {
  
  cat("\n====================================================\n")
  cat("Analysis:", config_name, "| Score mode:", score_mode, "\n")
  cat("====================================================\n")
  
  if (!score_mode %in% c(
    "oriented",
    "raw",
    "beta_weighted",
    "ridge_weighted",
    "elasticnet_weighted",
    "lof_burden"
  )) {
    stop(
      "Unknown score mode: ",
      score_mode
    )
  }
  
  if (
    score_mode %in% c("ridge_weighted", "elasticnet_weighted") &&
    config$type != "beta_weighted"
  ) {
    stop(
      "score_mode = 'ridge_weighted'/'elasticnet_weighted' may only be used with ",
      "config$type = 'beta_weighted' (same dataset ",
      "as the beta-weighted sensitivity analysis)."
    )
  }
  
  if (
    score_mode == "lof_burden" &&
    config$type != "lof_vep"
  ) {
    stop(
      "score_mode = 'lof_burden' may only be used with ",
      "config$type = 'lof_vep'."
    )
  }
  
  #------------------------------------------------------------
  # Output files
  #------------------------------------------------------------
  
  individuals_file <- make_output_name(
    pathway_name,
    "Univariate",
    "SNP_SumScore",
    "Individuals",
    config,
    score_mode
  )
  
  percentiles_file <- make_output_name(
    pathway_name,
    "Univariate",
    "SNP_SumScore",
    "Percentiles",
    config,
    score_mode
  )
  
  weights_file <- make_output_name(
    pathway_name,
    "Univariate",
    "SNP_SumScore",
    "Weights",
    config,
    score_mode
  )
  
  #------------------------------------------------------------
  # Load data
  #------------------------------------------------------------
  pathway_data <- readRDS(
    config$dataset
  )
  
  geno_filtered <- pathway_data$genotypes
  
  pheno_filtered <- as.data.frame(
    pathway_data$pheno,
    stringsAsFactors = FALSE
  )
  
  if (is.null(rownames(geno_filtered))) {
    stop(
      "The genotype matrix has no row names/IDs."
    )
  }
  
  if (!"id" %in% colnames(pheno_filtered)) {
    stop(
      "The phenotype data contain no column 'id'."
    )
  }
  
  rownames(geno_filtered) <- as.character(
    rownames(geno_filtered)
  )
  
  pheno_filtered$id <- as.character(
    pheno_filtered$id
  )
  
  #------------------------------------------------------------
  # Effect sizes from the univariate LD-pruned main analysis
  #------------------------------------------------------------
  
  regression_file <- make_output_name(
    pathway_name = pathway_name,
    analysis_type = "Univariate",
    analysis_name = "SNP_Logistic_Regression",
    result_name = "Results",
    config = analysis_configs$main_ld
  )
  
  if (!file.exists(regression_file)) {
    stop(
      "Univariate regression results not found:\n",
      regression_file
    )
  }
  
  results <- read.csv2(
    regression_file,
    stringsAsFactors = FALSE
  )
  
  required_result_cols <- c(
    "SNP",
    "Category",
    "Estimate_Add",
    "OR_Add",
    "Pvalue_Add"
  )
  
  missing_result_cols <- setdiff(
    required_result_cols,
    colnames(results)
  )
  
  if (length(missing_result_cols) > 0) {
    stop(
      "The following columns are missing in the regression results: ",
      paste(missing_result_cols, collapse = ", ")
    )
  }
  
  results <- results %>%
    mutate(
      SNP = trimws(as.character(SNP)),
      Estimate_Add = as.numeric(Estimate_Add),
      OR_Add = as.numeric(OR_Add),
      Pvalue_Add = as.numeric(Pvalue_Add)
    )
  
  # Table for equally weighted scores
  res_cat <- results %>%
    filter(
      !is.na(SNP),
      SNP != "",
      !is.na(Category),
      Category != ""
    ) %>%
    dplyr::select(
      SNP,
      Category
    ) %>%
    distinct(
      SNP,
      .keep_all = TRUE
    )
  
  # Table for β-weighted score
  res_beta <- results %>%
    filter(
      !is.na(SNP),
      SNP != "",
      is.finite(Estimate_Add)
    ) %>%
    dplyr::select(
      SNP,
      Beta = Estimate_Add,
      OR = OR_Add,
      Pvalue = Pvalue_Add,
      Category
    ) %>%
    distinct(
      SNP,
      .keep_all = TRUE
    )
  
  #------------------------------------------------------------
  # Match IDs robustly and sort identically
  #------------------------------------------------------------
  
  geno_ids <- as.character(
    rownames(geno_filtered)
  )
  
  pheno_ids <- as.character(
    pheno_filtered$id
  )
  
  if (anyDuplicated(geno_ids) > 0) {
    stop(
      "Duplicate IDs in the genotype matrix for ",
      config_name,
      "."
    )
  }
  
  if (anyDuplicated(pheno_ids) > 0) {
    stop(
      "Duplicate IDs in the phenotype data for ",
      config_name,
      "."
    )
  }
  
  # Keep the order of the genotype matrix
  common_ids <- geno_ids[
    geno_ids %in% pheno_ids
  ]
  
  if (length(common_ids) == 0) {
    
    warning(
      "No shared genotype/phenotype IDs for ",
      config_name
    )
    
    return(NULL)
  }
  
  geno_df <- geno_filtered[
    common_ids,
    ,
    drop = FALSE
  ] |>
    as.data.frame(
      stringsAsFactors = FALSE
    )
  
  pheno_match <- match(
    common_ids,
    pheno_filtered$id
  )
  
  if (anyNA(pheno_match)) {
    stop(
      "At least one genotype ID could not be assigned to the ",
      "phenotype data."
    )
  }
  
  pheno_df <- pheno_filtered[
    pheno_match,
    ,
    drop = FALSE
  ]
  
  pheno_df$id <- as.character(
    pheno_df$id
  )
  
  rownames(geno_df) <- common_ids
  rownames(pheno_df) <- common_ids
  
  if (!identical(
    as.character(rownames(geno_df)),
    as.character(pheno_df$id)
  )) {
    
    stop(
      "ID order after genotype/phenotype matching ",
      "does not match."
    )
  }
  
  #------------------------------------------------------------
  # Base set of SNPs depending on score mode
  #------------------------------------------------------------
  
  if (score_mode %in% c("beta_weighted", "ridge_weighted", "elasticnet_weighted")) {
    
    # Only SNPs with a valid, training-based estimated beta
    # (ridge_weighted/elasticnet_weighted use the same SNP base as
    # beta_weighted, so that all three weighting variants are directly comparable)
    valid_snps <- intersect(
      res_beta$SNP,
      colnames(geno_df)
    )
    
    #----------------------------------------------------------
    # Ridge / elastic net: shrunken weights instead of raw log(OR)
    #----------------------------------------------------------
    # Replaces res_beta$Beta with jointly (not SNP-by-SNP univariately)
    # estimated, regularized coefficients - exclusively on the
    # training data, analogous to the Estimate_Add estimation of the beta_weighted mode.
    if (score_mode %in% c("ridge_weighted", "elasticnet_weighted")) {
      
      library(glmnet)
      
      global_split_glmnet <- readRDS(
        output_file("Pathway_", pathway_name, "_Global_TrainTest_Split.rds")
      )
      
      train_ids_glmnet <- intersect(common_ids, global_split_glmnet$train_ids)
      
      if (length(train_ids_glmnet) < 10) {
        stop(
          "Too few training individuals for ridge/elastic-net weighting: ",
          config_name, " / ", score_mode
        )
      }
      
      x_train_raw <- as.matrix(
        geno_df[train_ids_glmnet, valid_snps, drop = FALSE]
      )
      storage.mode(x_train_raw) <- "numeric"
      
      # Mean imputation of missing genotypes per SNP (training data only),
      # consistent with the imputation in the univariate regression.
      x_train <- apply(x_train_raw, 2, function(col) {
        col[is.na(col)] <- mean(col, na.rm = TRUE)
        col
      })
      
      y_train <- as.numeric(
        as.character(pheno_df[train_ids_glmnet, "obesity"])
      )
      
      glmnet_alpha <- if (score_mode == "ridge_weighted") 0 else 0.5
      
      set.seed(123)
      cv_fit <- cv.glmnet(
        x = x_train,
        y = y_train,
        family = "binomial",
        alpha = glmnet_alpha
      )
      
      glmnet_coefs <- as.matrix(coef(cv_fit, s = "lambda.1se"))
      glmnet_beta  <- glmnet_coefs[valid_snps, 1]
      names(glmnet_beta) <- valid_snps
      
      # Replace the univariate beta column (log(OR)) with shrunken coefficients;
      # category/OR/p-value from the univariate analysis are kept
      # in the weights file for classification/documentation.
      res_beta <- res_beta %>%
        filter(SNP %in% valid_snps) %>%
        mutate(Beta = glmnet_beta[SNP])
      
      cat(
        "Ridge/elastic-net weights computed (", score_mode,
        ", alpha = ", glmnet_alpha,
        ", lambda.1se = ", signif(cv_fit$lambda.1se, 4),
        ", n_train = ", length(train_ids_glmnet), ")\n",
        sep = ""
      )
    }
    
  } else if (score_mode == "lof_burden") {
    
    # LoF SNPs are defined exclusively biologically via VEP.
    # No dependence on obesity-based risk/protective categories.
    valid_snps <- colnames(geno_df)
    
  } else if (score_mode %in% c("oriented", "raw")) {
    
    # Use the same SNP base for oriented and raw,
    # so that both analyses differ only by orientation.
    valid_snps <- intersect(
      res_cat$SNP,
      colnames(geno_df)
    )
    
  } else {
    
    stop(
      "Unsupported score mode in SNP selection: ",
      score_mode
    )
  }
  
  #------------------------------------------------------------
  # SNP selection depending on analysis configuration
  #------------------------------------------------------------
  lof_orientation_table <- NULL
  
  if (config$type == "functional_vep") {
    
    functional_snps <- load_functional_vep_snps(
      config$vep_file
    )
    
    valid_snps <- intersect(
      valid_snps,
      functional_snps
    )
  }
  
  if (config$type == "lof_vep") {
    
    if (score_mode != "lof_burden") {
      stop(
        "The configuration 'lof_vep' may only be run with ",
        "score_mode = 'lof_burden'."
      )
    }
    
    lof_orientation_table <- create_lof_orientation_table(
      vep_file = config$vep_file,
      bim_table = lof_bim,
      matrix_snps = colnames(geno_df)
    )
    
    if (nrow(lof_orientation_table) == 0) {
      warning(
        "No LoF SNPs could be matched between VEP, BIM and ",
        "genotype matrix."
      )
      return(NULL)
    }
    
    no_match_snps <- lof_orientation_table %>%
      filter(
        Allele_Match == "No_allele_match"
      ) %>%
      pull(SNP) %>%
      unique()
    
    if (length(no_match_snps) > 0) {
      warning(
        "The following VEP LoF alleles match no matrix allele ",
        "and are excluded: ",
        paste(no_match_snps, collapse = ", ")
      )
    }
    
    lof_orientation_table <- lof_orientation_table %>%
      filter(
        Allele_Match != "No_allele_match",
        !is.na(Invert_For_LoF_Score)
      ) %>%
      distinct(
        SNP,
        .keep_all = TRUE
      )
    
    valid_snps <- intersect(
      valid_snps,
      lof_orientation_table$SNP
    )
  }
  
  if (config$type == "database_kegg_only") {
    
    snp_gene_map <- read.csv2(
      config$snp_gene_map,
      stringsAsFactors = FALSE
    )
    
    kegg_only_snps <- snp_gene_map %>%
      filter(
        GENE %in% config$kegg_genes
      ) %>%
      pull(SNP) %>%
      unique()
    
    valid_snps <- intersect(
      valid_snps,
      kegg_only_snps
    )
  }
  
  # Align order with the genotype matrix and remove duplicates
  valid_snps <- colnames(geno_df)[
    colnames(geno_df) %in% unique(valid_snps)
  ]
  
  if (length(valid_snps) == 0) {
    
    warning(
      "No valid SNPs for ",
      config_name,
      " / ",
      score_mode,
      " found."
    )
    
    return(NULL)
  }
  
 
  
  #------------------------------------------------------------
  # Prepare genotype matrix
  #------------------------------------------------------------
  
  geno_score <- geno_df[
    ,
    valid_snps,
    drop = FALSE
  ]
  
  geno_score[] <- lapply(
    geno_score,
    function(x) {
      
      if (is.list(x)) {
        x <- unlist(x)
      }
      
      as.numeric(x)
    }
  )
  
  geno_matrix_original <- as.matrix(geno_score)
  storage.mode(geno_matrix_original) <- "numeric"
  
  qc_range <- range(
    geno_matrix_original,
    na.rm = TRUE
  )
  
  if (
    !all(is.finite(qc_range)) ||
    qc_range[1] < 0 ||
    qc_range[2] > 2
  ) {
    stop(
      "Genotype values outside the expected range 0–2: ",
      paste(qc_range, collapse = " to ")
    )
  }
  
  #------------------------------------------------------------
  # Exclude individuals with more than 10 % missing SNPs
  #------------------------------------------------------------
  
  indiv_missing <- rowMeans(
    is.na(geno_matrix_original)
  )
  
  keep <- indiv_missing <= 0.10
  
  ids_before_missingness <- as.character(
    rownames(geno_score)
  )
  
  ids_after_missingness <- ids_before_missingness[
    keep
  ]
  
  geno_score <- geno_score[
    ids_after_missingness,
    ,
    drop = FALSE
  ]
  
  pheno_match_missingness <- match(
    ids_after_missingness,
    pheno_df$id
  )
  
  if (anyNA(pheno_match_missingness)) {
    stop(
      "At least one ID could not be assigned to the phenotype data after the missingness filter ",
      ""
    )
  }
  
  pheno_df <- pheno_df[
    pheno_match_missingness,
    ,
    drop = FALSE
  ]
  
  pheno_df$id <- as.character(
    pheno_df$id
  )
  
  rownames(pheno_df) <- pheno_df$id
  
  if (!identical(
    as.character(rownames(geno_score)),
    as.character(pheno_df$id)
  )) {
    
    stop(
      "ID order after the individual missingness filter ",
      "does not match."
    )
  }
  
  #------------------------------------------------------------
  # Training-based genotype imputation for ALL score modes
  #------------------------------------------------------------
  
  global_split <- readRDS(
    output_file(
      "Pathway_",
      pathway_name,
      "_Global_TrainTest_Split.rds"
    )
  )
  
  train_rows <- rownames(geno_score) %in%
    global_split$train_ids
  
  if (!any(train_rows)) {
    stop(
      "No training individuals available for genotype imputation."
    )
  }
  
  # Numeric original matrix after exclusion of individuals
  # with > 10 % missing genotypes
  geno_matrix <- as.matrix(geno_score)
  storage.mode(geno_matrix) <- "numeric"
  
  # Determine SNP means exclusively in the training dataset
  train_snp_means <- colMeans(
    geno_matrix[
      train_rows,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )
  
  # SNPs without a single observed genotype in the training set
  invalid_imputation_snps <- names(
    train_snp_means
  )[
    !is.finite(train_snp_means)
  ]
  
  if (length(invalid_imputation_snps) > 0) {
    
    warning(
      "The following SNPs have no observed genotype values in the training set ",
      "and are excluded: ",
      paste(
        invalid_imputation_snps,
        collapse = ", "
      )
    )
    
    valid_imputation_snps <- setdiff(
      colnames(geno_matrix),
      invalid_imputation_snps
    )
    
    geno_matrix <- geno_matrix[
      ,
      valid_imputation_snps,
      drop = FALSE
    ]
    
    train_snp_means <- train_snp_means[
      valid_imputation_snps
    ]
    
    valid_snps <- intersect(
      valid_snps,
      valid_imputation_snps
    )
    
    if (length(valid_snps) == 0) {
      warning(
        "After excluding non-imputable SNPs, no SNPs remain."
      )
      return(NULL)
    }
  }
  
  # Replace missing values in train AND test exclusively with the
  # SNP means computed in the training set
  n_imputed_total <- 0L
  
  n_imputed_by_snp <- setNames(
    integer(ncol(geno_matrix)),
    colnames(geno_matrix)
  )
  
  for (j in seq_len(ncol(geno_matrix))) {
    
    missing_j <- is.na(
      geno_matrix[, j]
    )
    
    n_imputed_by_snp[j] <- sum(
      missing_j
    )
    
    if (any(missing_j)) {
      
      geno_matrix[
        missing_j,
        j
      ] <- train_snp_means[j]
      
      n_imputed_total <- n_imputed_total +
        sum(missing_j)
    }
  }
  
  if (anyNA(geno_matrix)) {
    stop(
      "After the training-based imputation, missing ",
      "genotype values remain."
    )
  }
  
  # Documentation of the imputation values used
  imputation_table <- data.frame(
    SNP = colnames(geno_matrix),
    Training_Mean_Genotype = as.numeric(
      train_snp_means[
        colnames(geno_matrix)
      ]
    ),
    N_Imputed_All = as.integer(
      n_imputed_by_snp[
        colnames(geno_matrix)
      ]
    ),
    Imputation_Source = "Training_Set",
    stringsAsFactors = FALSE
  )
  
  imputation_file <- make_output_name(
    pathway_name,
    "Univariate",
    "SNP_SumScore",
    "Imputation_Values",
    config,
    score_mode
  )
  
  write.csv2(
    imputation_table,
    imputation_file,
    row.names = FALSE
  )
  
  #------------------------------------------------------------
  # Compute score
  #------------------------------------------------------------
  
  protective <- character(0)
  beta_table_used <- NULL
  prs_raw <- NULL
  snp_sum <- NULL
  
  if (score_mode == "oriented") {
    
    protective <- res_cat %>%
      filter(
        tolower(Category) == "protective"
      ) %>%
      pull(SNP) %>%
      intersect(
        colnames(geno_matrix)
      )
    
    # Re-orient protective SNPs only after imputation
    if (length(protective) > 0) {
      
      protective_idx <- match(
        protective,
        colnames(geno_matrix)
      )
      
      geno_matrix[
        ,
        protective_idx
      ] <- 2 - geno_matrix[
        ,
        protective_idx,
        drop = FALSE
      ]
    }
    
    snp_sum <- rowSums(
      geno_matrix
    )
  }
  
  if (score_mode == "raw") {
    
    snp_sum <- rowSums(
      geno_matrix
    )
  }
  
  if (score_mode == "lof_burden") {
    
    if (is.null(lof_orientation_table)) {
      stop(
        "For the LoF burden no ",
        "allele orientation table was created."
      )
    }
    
    # Create copy:
    # The original geno_matrix remains unchanged.
    geno_matrix_lof <- geno_matrix
    
    # Align the orientation table exactly with the SNP order
    # of the matrix actually used
    lof_orientation_used <- lof_orientation_table[
      match(
        colnames(geno_matrix_lof),
        lof_orientation_table$SNP
      ),
      ,
      drop = FALSE
    ]
    
    if (
      any(is.na(lof_orientation_used$SNP)) ||
      !all(
        lof_orientation_used$SNP ==
        colnames(geno_matrix_lof)
      )
    ) {
      stop(
        "An unambiguous allele orientation could not be ",
        "assigned to all LoF SNPs used."
      )
    }
    
    # Identify SNPs for which the original matrix
    # counts the counter-allele of the VEP LoF allele
    invert_snps <- lof_orientation_used$SNP[
      lof_orientation_used$Invert_For_LoF_Score %in% TRUE
    ]
    
    # Re-orient only the temporary LoF matrix
    if (length(invert_snps) > 0) {
      
      invert_idx <- match(
        invert_snps,
        colnames(geno_matrix_lof)
      )
      
      geno_matrix_lof[
        ,
        invert_idx
      ] <- 2 - geno_matrix_lof[
        ,
        invert_idx,
        drop = FALSE
      ]
    }
    
    #----------------------------------------------------------
    # Document LoF allele orientation used
    #----------------------------------------------------------
    
    lof_orientation_file <- make_output_name(
      pathway_name,
      "Univariate",
      "SNP_SumScore",
      "LoF_Allele_Orientation",
      config,
      score_mode
    )
    
    lof_orientation_used$Applied_Transformation <- ifelse(
      lof_orientation_used$Invert_For_LoF_Score %in% TRUE,
      "2_minus_original_dosage",
      "original_dosage"
    )
    
    write.csv2(
      lof_orientation_used,
      lof_orientation_file,
      row.names = FALSE
    )
    
    #----------------------------------------------------------
    # Compute LoF burden
    #----------------------------------------------------------
    
    snp_sum <- rowSums(
      geno_matrix_lof
    )
  }
  
  if (score_mode %in% c("beta_weighted", "ridge_weighted", "elasticnet_weighted")) {
    
    #----------------------------------------------------------
    # Align β weights with the order of the imputed matrix
    #----------------------------------------------------------
    
    beta_table_used <- res_beta[
      match(
        colnames(geno_matrix),
        res_beta$SNP
      ),
      ,
      drop = FALSE
    ]
    
    if (
      any(is.na(beta_table_used$SNP)) ||
      any(!is.finite(beta_table_used$Beta))
    ) {
      stop(
        "A valid beta value could not be assigned to all ",
        "selected SNPs."
      )
    }
    
    stopifnot(
      all(
        beta_table_used$SNP ==
          colnames(geno_matrix)
      )
    )
    
    beta_vector <- beta_table_used$Beta
    names(beta_vector) <- beta_table_used$SNP
    
    # Additionally document unweighted dosage
    prs_raw <- rowSums(
      geno_matrix
    )
    
    # β-weighted score:
    # Sum of genotype dosage × internally estimated log(OR)
    snp_sum <- as.numeric(
      geno_matrix %*% beta_vector
    )
    
    beta_table_used$Training_Mean_Genotype <-
      train_snp_means[
        beta_table_used$SNP
      ]
    
    beta_table_used$N_Imputed_All <-
      n_imputed_by_snp[
        beta_table_used$SNP
      ]
    
    beta_table_used$Weight_Source <- switch(
      score_mode,
      "ridge_weighted" = "Training_Set_Ridge_CV_Lambda1se",
      "elasticnet_weighted" = "Training_Set_ElasticNet_Alpha05_CV_Lambda1se",
      "Training_Set_Univariate_Logistic_Regression"
    )
    
    beta_table_used$Imputation_Source <-
      "Training_Set"
    
    write.csv2(
      beta_table_used,
      weights_file,
      row.names = FALSE
    )
  }
  
  if (is.null(snp_sum)) {
    stop(
      "For score mode ",
      score_mode,
      " no score was computed."
    )
  }
  
  #------------------------------------------------------------
  # LoF burden information
  #------------------------------------------------------------
  if (score_mode == "lof_burden") {
    
    cat(
      "LoF burden score: sum of biologically annotated ",
      "LoF allele dosages.\n",
      sep = ""
    )
    
    cat(
      "No orientation based on the obesity association.\n"
    )
    
    cat(
      "SNPs inverted for LoF burden: ",
      length(invert_snps),
      "\n"
    )
    
    cat(
      "LoF orientation file: ",
      lof_orientation_file,
      "\n"
    )
  }
  
  if (
    length(snp_sum) != nrow(pheno_df)
  ) {
    stop(
      "The length of the computed score does not match the ",
      "number of individuals."
    )
  }
  
  if (any(!is.finite(snp_sum))) {
    stop(
      "The computed SNP sum score contains invalid values."
    )
  }
  
  #------------------------------------------------------------
  # Result dataset
  #------------------------------------------------------------
  
  df <- pheno_df %>%
    mutate(
      SNP_SumScore = as.numeric(snp_sum),
      Score_Mode = score_mode,
      Analysis_Config = config_name,
      N_Score_SNPs = ncol(geno_matrix),
      N_Imputed_Genotypes = n_imputed_total,
      Imputation_Method = "Training_Mean_Genotype"
    )
  
  if (score_mode %in% c("beta_weighted", "ridge_weighted", "elasticnet_weighted")) {
    
    # Unweighted dosage sum as additional information
    df$PRS_Raw <- as.numeric(prs_raw)
    
  } else {
    
    df$PRS_Raw <- NA_real_
  }
  
  #------------------------------------------------------------
  # Export individuals
  #------------------------------------------------------------
  
  score_table <- data.frame(
    id = as.character(df$id),
    SNP_SumScore = df$SNP_SumScore,
    PRS_Raw = df$PRS_Raw,
    Score_Mode = score_mode,
    Analysis_Config = config_name,
    N_Score_SNPs = ncol(geno_matrix),
    N_Imputed_Genotypes = n_imputed_total,
    Imputation_Method = "Training_Mean_Genotype",
    stringsAsFactors = FALSE
  )
  
  write.csv2(
    score_table,
    individuals_file,
    row.names = FALSE
  )
  
  #------------------------------------------------------------
  # Percentiles based exclusively on the training set
  #------------------------------------------------------------
  
  global_split <- readRDS(
    output_file(
      "Pathway_",
      pathway_name,
      "_Global_TrainTest_Split.rds"
    )
  )
  
  df_train <- df %>%
    filter(id %in% global_split$train_ids)
  
  if (nrow(df_train) == 0) {
    stop("No training individuals available for percentile calculation.")
  }
  
  percentile_table <- data.frame(
    Percentile = seq(0, 100, by = 1),
    
    SNP_SumScore = as.numeric(
      quantile(
        df_train$SNP_SumScore,
        probs = seq(0, 1, by = 0.01),
        na.rm = TRUE
      )
    ),
    
    Score_Mode = score_mode,
    Analysis_Config = config_name,
    Based_On = "Training_Set",
    stringsAsFactors = FALSE
  )
  
  write.csv2(
    percentile_table,
    percentiles_file,
    row.names = FALSE
  )
  
  #------------------------------------------------------------
  # Console output and attributes
  #------------------------------------------------------------
  
  cat(
    "SNPs used:",
    ncol(geno_matrix),
    "\n"
  )
  
  cat(
    "Genotype imputation for all score modes: ",
    "SNP means exclusively from the training set.\n",
    sep = ""
  )
  
  cat(
    "Imputed genotype values (total):",
    n_imputed_total,
    "\n"
  )
  
  cat(
    "Imputation file:",
    imputation_file,
    "\n"
  )
  
  if (score_mode == "oriented") {
    
    cat(
      "Inverted protective SNPs:",
      length(protective),
      "\n"
    )
    
    cat(
      "No β weighting.\n"
    )
  }
  
  if (score_mode == "raw") {
    
    cat(
      "No SNP inversion and no weighting.\n"
    )
  }
  
  if (score_mode == "beta_weighted") {
    
    cat(
      "Weighting: β = Estimate_Add = log(OR)\n"
    )
    
    cat(
      "No SNP inversion.\n"
    )
    
    cat(
      "Weights table:",
      weights_file,
      "\n"
    )
  }
  
  if (score_mode %in% c("ridge_weighted", "elasticnet_weighted")) {
    
    cat(
      "Weighting: β = ",
      ifelse(score_mode == "ridge_weighted", "Ridge", "Elastic Net (alpha=0.5)"),
      " coefficient (cv.glmnet, lambda.1se), training data\n",
      sep = ""
    )
    
    cat(
      "No SNP inversion.\n"
    )
    
    cat(
      "Weights table:",
      weights_file,
      "\n"
    )
  }
  
  if (score_mode == "lof_burden") {
    
    cat(
      "LoF burden score: sum of biologically annotated ",
      "LoF allele dosages.\n",
      sep = ""
    )
    
    cat(
      "No orientation based on the obesity association.\n"
    )
  }
  
  cat(
    "Score file:",
    individuals_file,
    "\n"
  )
  
  #------------------------------------------------------------
  # Attributes for subsequent model analysis
  #------------------------------------------------------------
  
  attr(
    df,
    "n_score_snps"
  ) <- ncol(geno_matrix)
  
  attr(
    df,
    "score_mode"
  ) <- score_mode
  
  attr(
    df,
    "analysis_config"
  ) <- config_name
  
  attr(
    df,
    "imputation_method"
  ) <- "Training_Mean_Genotype"
  
  attr(
    df,
    "n_imputed_genotypes"
  ) <- n_imputed_total
  
  attr(
    df,
    "imputation_values"
  ) <- imputation_table
  
  if (score_mode %in% c("beta_weighted", "ridge_weighted", "elasticnet_weighted")) {
    
    attr(
      df,
      "beta_weights"
    ) <- beta_table_used
  }
  
  return(df)
}

# ---- run_snp_sumscore_model() ----
run_snp_sumscore_model <- function(
    df,
    result_name) {
  
  cat("\n====================================================\n")
  cat("SNP sum score model:", result_name, "\n")
  cat("====================================================\n")
  
  if (is.null(df)) {
    warning("No dataset for ", result_name)
    return(NULL)
  }
  
  required_columns <- c(
    "id",
    "obesity",
    "SNP_SumScore",
    "Score_Mode",
    "Analysis_Config"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(df)
  )
  
  if (length(missing_columns) > 0) {
    
    warning(
      "The following columns are missing for ",
      result_name,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
    
    return(NULL)
  }
  
  analysis_configs_present <- unique(
    as.character(df$Analysis_Config)
  )
  
  score_modes_present <- unique(
    as.character(df$Score_Mode)
  )
  
  analysis_configs_present <- analysis_configs_present[
    !is.na(analysis_configs_present)
  ]
  
  score_modes_present <- score_modes_present[
    !is.na(score_modes_present)
  ]
  
  if (length(analysis_configs_present) != 1) {
    
    warning(
      "Not exactly one Analysis_Config in ",
      result_name
    )
    
    return(NULL)
  }
  
  if (length(score_modes_present) != 1) {
    
    warning(
      "Not exactly one Score_Mode in ",
      result_name
    )
    
    return(NULL)
  }
  
  analysis_config <- analysis_configs_present[1]
  score_mode      <- score_modes_present[1]
  
  if (!analysis_config %in% names(analysis_configs)) {
    
    warning(
      "Unknown analysis_config: ",
      analysis_config
    )
    
    return(NULL)
  }
  
  config <- analysis_configs[[analysis_config]]
  
  n_score_snps <- get_n_score_snps(df)
  
  #-------------------------------------------------------------------------
  # Consistency check for the LoF analysis
  #-------------------------------------------------------------------------
  
  if (
    score_mode == "lof_burden" &&
    config$type != "lof_vep"
  ) {
    warning(
      "score_mode = 'lof_burden' is only permitted for ",
      "config$type = 'lof_vep': ",
      result_name
    )
    return(NULL)
  }
  
  if (
    config$type == "lof_vep" &&
    score_mode != "lof_burden"
  ) {
    warning(
      "A LoF VEP configuration may only be used with ",
      "score_mode = 'lof_burden': ",
      result_name
    )
    return(NULL)
  }
  
  #-------------------------------------------------------------------------
  # Complete cases
  #-------------------------------------------------------------------------
  
  df_clean <- df %>%
    mutate(
      id = as.character(id),
      obesity = as.numeric(as.character(obesity)),
      SNP_SumScore = as.numeric(SNP_SumScore)
    ) %>%
    filter(
      !is.na(id),
      id != "",
      !is.na(obesity),
      is.finite(SNP_SumScore)
    ) %>%
    distinct(
      id,
      .keep_all = TRUE
    )
  
  if (nrow(df_clean) < 20) {
    warning("Too few complete cases for ", result_name)
    return(NULL)
  }
  
  df_clean$obesity_bin <- ifelse(
    df_clean$obesity == 1,
    1,
    0
  )
  
  if (length(unique(df_clean$obesity_bin)) < 2) {
    warning("Only one outcome class available for ", result_name)
    return(NULL)
  }
  
  #-------------------------------------------------------------------------
  # Apply global train/test split
  #-------------------------------------------------------------------------
  
  global_split <- readRDS(
    output_file(
      "Pathway_",
      pathway_name,
      "_Global_TrainTest_Split.rds"
    )
  )
  
  train_set <- df_clean %>%
    filter(id %in% global_split$train_ids)
  
  test_set <- df_clean %>%
    filter(id %in% global_split$test_ids)
  
  if (
    nrow(train_set) < 10 ||
    nrow(test_set) < 5
  ) {
    
    warning(
      "Train or test dataset is too small for ",
      result_name
    )
    
    return(NULL)
  }
  
  if (
    length(unique(train_set$obesity_bin)) < 2 ||
    length(unique(test_set$obesity_bin)) < 2
  ) {
    
    warning(
      "Train or test contains only one outcome class for ",
      result_name
    )
    
    return(NULL)
  }
  
  if (
    length(unique(train_set$SNP_SumScore)) < 2
  ) {
    
    warning(
      "The SNP sum score has no variation in the training dataset for ",
      result_name
    )
    
    return(NULL)
  }
  
  #-------------------------------------------------------------------------
  # Score distribution in the training dataset
  #-------------------------------------------------------------------------
  
  score_mean_train <- mean(
    train_set$SNP_SumScore,
    na.rm = TRUE
  )
  
  score_sd_train <- sd(
    train_set$SNP_SumScore,
    na.rm = TRUE
  )
  
  if (
    !is.finite(score_sd_train) ||
    score_sd_train <= 0
  ) {
    
    warning(
      "The standard deviation of the score is invalid for ",
      result_name
    )
    
    return(NULL)
  }
  
  # Standardized version:
  # exclusively with mean and SD of the training dataset
  train_set$SNP_SumScore_SD <- (
    train_set$SNP_SumScore - score_mean_train
  ) / score_sd_train
  
  test_set$SNP_SumScore_SD <- (
    test_set$SNP_SumScore - score_mean_train
  ) / score_sd_train
  
  #-------------------------------------------------------------------------
  # Model with score in original units
  #-------------------------------------------------------------------------
  
  model_train <- tryCatch(
    glm(
      obesity_bin ~ SNP_SumScore,
      data = train_set,
      family = binomial
    ),
    error = function(e) NULL
  )
  
  if (is.null(model_train)) {
    warning("Model could not be estimated for ", result_name)
    return(NULL)
  }
  
  coef_table <- summary(model_train)$coefficients
  
  if (!"SNP_SumScore" %in% rownames(coef_table)) {
    
    warning(
      "No SNP sum score coefficient available for ",
      result_name
    )
    
    return(NULL)
  }
  
  beta0 <- unname(
    coef(model_train)["(Intercept)"]
  )
  
  beta1 <- unname(
    coef_table[
      "SNP_SumScore",
      "Estimate"
    ]
  )
  
  se1 <- unname(
    coef_table[
      "SNP_SumScore",
      "Std. Error"
    ]
  )
  
  wald_pvalue <- unname(
    coef_table[
      "SNP_SumScore",
      "Pr(>|z|)"
    ]
  )
  
  OR <- exp(beta1)
  
  ci <- tryCatch(
    exp(
      confint.default(model_train)[
        "SNP_SumScore",
      ]
    ),
    error = function(e) c(NA_real_, NA_real_)
  )
  
  #-------------------------------------------------------------------------
  # Effect per one standard deviation
  #-------------------------------------------------------------------------
  
  model_train_sd <- tryCatch(
    glm(
      obesity_bin ~ SNP_SumScore_SD,
      data = train_set,
      family = binomial
    ),
    error = function(e) NULL
  )
  
  beta_per_sd <- NA_real_
  se_per_sd   <- NA_real_
  OR_per_sd   <- NA_real_
  ci_per_sd   <- c(NA_real_, NA_real_)
  p_per_sd    <- NA_real_
  
  if (!is.null(model_train_sd)) {
    
    coef_sd <- summary(model_train_sd)$coefficients
    
    if ("SNP_SumScore_SD" %in% rownames(coef_sd)) {
      
      beta_per_sd <- coef_sd[
        "SNP_SumScore_SD",
        "Estimate"
      ]
      
      se_per_sd <- coef_sd[
        "SNP_SumScore_SD",
        "Std. Error"
      ]
      
      p_per_sd <- coef_sd[
        "SNP_SumScore_SD",
        "Pr(>|z|)"
      ]
      
      OR_per_sd <- exp(beta_per_sd)
      
      ci_per_sd <- tryCatch(
        exp(
          confint.default(model_train_sd)[
            "SNP_SumScore_SD",
          ]
        ),
        error = function(e) c(
          NA_real_,
          NA_real_
        )
      )
    }
  }
  
  #-------------------------------------------------------------------------
  # Null model and likelihood-ratio test
  #-------------------------------------------------------------------------
  
  model_null <- glm(
    obesity_bin ~ 1,
    data = train_set,
    family = binomial
  )
  
  lrt_pvalue <- tryCatch(
    anova(
      model_null,
      model_train,
      test = "LRT"
    )$`Pr(>Chi)`[2],
    error = function(e) NA_real_
  )
  
  #-------------------------------------------------------------------------
  # Predictions
  #-------------------------------------------------------------------------
  
  train_set$pred_prob <- predict(
    model_train,
    newdata = train_set,
    type = "response"
  )
  
  test_set$pred_prob <- predict(
    model_train,
    newdata = test_set,
    type = "response"
  )
  
  #-------------------------------------------------------------------------
  # ROC
  #-------------------------------------------------------------------------
  
  roc_train <- safe_roc(
    obs = train_set$obesity_bin,
    pred = train_set$pred_prob
  )
  
  roc_obj <- safe_roc(
    obs = test_set$obesity_bin,
    pred = test_set$pred_prob
  )
  
  if (
    is.null(roc_train) ||
    is.null(roc_obj)
  ) {
    
    warning(
      "ROC could not be computed for ",
      result_name
    )
    
    return(NULL)
  }
  
  auc_val <- as.numeric(
    auc(roc_obj)
  )
  
  auc_ci <- tryCatch(
    as.numeric(
      ci.auc(roc_obj)
    ),
    error = function(e) c(
      NA_real_,
      NA_real_,
      NA_real_
    )
  )
  
  #-------------------------------------------------------------------------
  # Cutoffs
  #-------------------------------------------------------------------------
  
  cutoff_0.5 <- 0.5
  
  cutoff_youden <- tryCatch(
    {
      
      cutoff_result <- coords(
        roc_train,
        x = "best",
        ret = "threshold",
        best.method = "youden",
        transpose = FALSE
      )
      
      as.numeric(
        cutoff_result[1]
      )
    },
    error = function(e) NA_real_
  )
  
  score_cutoff_0.5 <- calc_score_cutoff(
    prob_cutoff = cutoff_0.5,
    beta0 = beta0,
    beta1 = beta1
  )
  
  score_cutoff_youden <- calc_score_cutoff(
    prob_cutoff = cutoff_youden,
    beta0 = beta0,
    beta1 = beta1
  )
  
  cutoff_direction <- get_score_cutoff_direction(
    beta1
  )
  
  metrics_0.5 <- calc_metrics(
    pred = test_set$pred_prob,
    obs = test_set$obesity_bin,
    cutoff = cutoff_0.5
  )
  
  metrics_youden <- calc_metrics(
    pred = test_set$pred_prob,
    obs = test_set$obesity_bin,
    cutoff = cutoff_youden
  )
  
  #-------------------------------------------------------------------------
  # Batch sensitivity analysis
  #
  # Main analysis only, not sensitivity_betaweighted.
  #-------------------------------------------------------------------------
  
  batch_sensitivity_results <- NULL
  
  if (
    analysis_config == "main_ld" &&
    "Batch" %in% colnames(train_set)
  ) {
    
    train_set$Batch <- factor(
      train_set$Batch
    )
    
    test_set$Batch <- factor(
      test_set$Batch,
      levels = levels(train_set$Batch)
    )
    
    if (
      nlevels(train_set$Batch) >= 2
    ) {
      
      model_batch <- tryCatch(
        glm(
          obesity_bin ~ SNP_SumScore + Batch,
          data = train_set,
          family = binomial
        ),
        error = function(e) NULL
      )
      
      if (!is.null(model_batch)) {
        
        batch_coef <- summary(
          model_batch
        )$coefficients
        
        if (
          "SNP_SumScore" %in%
          rownames(batch_coef)
        ) {
          
          batch_ci <- tryCatch(
            exp(
              confint.default(model_batch)[
                "SNP_SumScore",
              ]
            ),
            error = function(e) c(
              NA_real_,
              NA_real_
            )
          )
          
          batch_beta <- batch_coef[
            "SNP_SumScore",
            "Estimate"
          ]
          
          batch_sensitivity_results <- data.frame(
            Pathway = pathway_name,
            Analysis = result_name,
            Score_Mode = score_mode,
            Analysis_Config = analysis_config,
            Sensitivity_Type = "Batch_Adjusted",
            N_Train = nobs(model_batch),
            
            Estimate = batch_beta,
            StdError = batch_coef[
              "SNP_SumScore",
              "Std. Error"
            ],
            
            OR = exp(batch_beta),
            CI_lower = batch_ci[1],
            CI_upper = batch_ci[2],
            
            Pvalue = batch_coef[
              "SNP_SumScore",
              "Pr(>|z|)"
            ],
            
            OR_per_SD = exp(
              batch_beta * score_sd_train
            ),
            
            stringsAsFactors = FALSE
          )
          
          batch_file <- make_output_name(
            pathway_name,
            "Univariate",
            "SNP_SumScore",
            "Batch_Sensitivity_Results",
            config,
            score_mode
          )
          
          write.csv2(
            batch_sensitivity_results,
            batch_file,
            row.names = FALSE
          )
          
          cat(
            "Batch sensitivity analysis saved.\n"
          )
        }
      }
    }
  }
  
  #-------------------------------------------------------------------------
  # Main results
  #-------------------------------------------------------------------------
  
  results_out <- data.frame(
    
    Pathway = rep(
      pathway_name,
      2
    ),
    
    Analysis = rep(
      result_name,
      2
    ),
    
    Score_Mode = rep(
      score_mode,
      2
    ),
    
    Analysis_Config = rep(
      analysis_config,
      2
    ),
    
    Score_Type = rep(
      dplyr::case_when(
        score_mode == "beta_weighted" ~
          "Beta-weighted PRS",
        
        score_mode == "ridge_weighted" ~
          "Ridge-weighted PRS",
        
        score_mode == "elasticnet_weighted" ~
          "Elastic-Net-weighted PRS",
        
        score_mode == "lof_burden" ~
          "LoF allele burden",
        
        score_mode == "oriented" ~
          "Obesity-oriented SNP SumScore",
        
        score_mode == "raw" ~
          "Unweighted raw SNP SumScore",
        
        TRUE ~
          "Other SNP Score"
      ),
      2
    ),
    
    Score_Unit = rep(
      dplyr::case_when(
        score_mode == "lof_burden" ~
          "One additional LoF allele copy",
        
        score_mode == "beta_weighted" ~
          "One unit of beta-weighted score",
        
        score_mode == "ridge_weighted" ~
          "One unit of ridge-weighted score",
        
        score_mode == "elasticnet_weighted" ~
          "One unit of elastic-net-weighted score",
        
        score_mode == "oriented" ~
          "One additional obesity-oriented allele",
        
        score_mode == "raw" ~
          "One additional counted allele",
        
        TRUE ~
          "One score unit"
      ),
      2
    ),
    
    N_Score_SNPs = rep(
      n_score_snps,
      2
    ),
    
    N_Total = rep(
      nrow(df_clean),
      2
    ),
    
    N_Train = rep(
      nrow(train_set),
      2
    ),
    
    N_Test = rep(
      nrow(test_set),
      2
    ),
    
    N_Train_NonObese = rep(
      sum(train_set$obesity_bin == 0),
      2
    ),
    
    N_Train_Obese = rep(
      sum(train_set$obesity_bin == 1),
      2
    ),
    
    N_Test_NonObese = rep(
      sum(test_set$obesity_bin == 0),
      2
    ),
    
    N_Test_Obese = rep(
      sum(test_set$obesity_bin == 1),
      2
    ),
    
    Score_Mean_Train = rep(
      score_mean_train,
      2
    ),
    
    Score_SD_Train = rep(
      score_sd_train,
      2
    ),
    
    Estimate_per_Unit = rep(
      beta1,
      2
    ),
    
    StdError_per_Unit = rep(
      se1,
      2
    ),
    
    OR_per_Unit = rep(
      OR,
      2
    ),
    
    CI_lower_per_Unit = rep(
      ci[1],
      2
    ),
    
    CI_upper_per_Unit = rep(
      ci[2],
      2
    ),
    
    Wald_pvalue = rep(
      wald_pvalue,
      2
    ),
    
    Estimate_per_SD = rep(
      beta_per_sd,
      2
    ),
    
    StdError_per_SD = rep(
      se_per_sd,
      2
    ),
    
    OR_per_SD = rep(
      OR_per_sd,
      2
    ),
    
    CI_lower_per_SD = rep(
      ci_per_sd[1],
      2
    ),
    
    CI_upper_per_SD = rep(
      ci_per_sd[2],
      2
    ),
    
    Pvalue_per_SD = rep(
      p_per_sd,
      2
    ),
    
    LRT_pvalue = rep(
      lrt_pvalue,
      2
    ),
    
    AIC = rep(
      AIC(model_train),
      2
    ),
    
    AUC = rep(
      auc_val,
      2
    ),
    
    AUC_CI_lower = rep(
      auc_ci[1],
      2
    ),
    
    AUC_CI_median = rep(
      auc_ci[2],
      2
    ),
    
    AUC_CI_upper = rep(
      auc_ci[3],
      2
    ),
    
    Cutoff_Type = c(
      "Probability_0.5",
      "Youden_Training"
    ),
    
    Cutoff_Prob = c(
      cutoff_0.5,
      cutoff_youden
    ),
    
    SNP_SumScore_Cutoff = c(
      score_cutoff_0.5,
      score_cutoff_youden
    ),
    
    Score_Cutoff_Direction = rep(
      cutoff_direction,
      2
    ),
    
    TN = c(
      metrics_0.5$TN,
      metrics_youden$TN
    ),
    
    FP = c(
      metrics_0.5$FP,
      metrics_youden$FP
    ),
    
    FN = c(
      metrics_0.5$FN,
      metrics_youden$FN
    ),
    
    TP = c(
      metrics_0.5$TP,
      metrics_youden$TP
    ),
    
    Sensitivity = c(
      metrics_0.5$Sensitivity,
      metrics_youden$Sensitivity
    ),
    
    Specificity = c(
      metrics_0.5$Specificity,
      metrics_youden$Specificity
    ),
    
    PPV = c(
      metrics_0.5$PPV,
      metrics_youden$PPV
    ),
    
    NPV = c(
      metrics_0.5$NPV,
      metrics_youden$NPV
    ),
    
    Accuracy = c(
      metrics_0.5$Accuracy,
      metrics_youden$Accuracy
    ),
    
    stringsAsFactors = FALSE
  )
  
  #-------------------------------------------------------------------------
  # Save regression results
  #-------------------------------------------------------------------------
  
  regression_file <- make_output_name(
    pathway_name,
    "Univariate",
    "SNP_SumScore",
    "Regression_Results",
    config,
    score_mode
  )
  
  write.csv2(
    results_out,
    regression_file,
    row.names = FALSE
  )
  
  #-------------------------------------------------------------------------
  # Save test predictions
  #-------------------------------------------------------------------------
  
  prediction_columns <- c(
    "id",
    "obesity",
    "obesity_bin",
    "SNP_SumScore",
    "SNP_SumScore_SD",
    "pred_prob",
    "Score_Mode",
    "Analysis_Config"
  )
  
  # Additionally keep PRS_Raw for the β-weighted score
  if ("PRS_Raw" %in% colnames(test_set)) {
    prediction_columns <- append(
      prediction_columns,
      "PRS_Raw",
      after = 4
    )
  }
  
  prediction_table <- test_set %>%
    dplyr::select(
      all_of(prediction_columns)
    )
  
  prediction_table$Predicted_Class_0.5 <- ifelse(
    prediction_table$pred_prob >= cutoff_0.5,
    1,
    0
  )
  
  prediction_table$Predicted_Class_Youden <- ifelse(
    prediction_table$pred_prob >= cutoff_youden,
    1,
    0
  )
  
  predictions_file <- make_output_name(
    pathway_name,
    "Univariate",
    "SNP_SumScore",
    "Predictions_TestSet",
    config,
    score_mode
  )
  
  write.csv2(
    prediction_table,
    predictions_file,
    row.names = FALSE
  )
  
  #-------------------------------------------------------------------------
  # Save ROC data
  #-------------------------------------------------------------------------
  
  roc_df <- data.frame(
    Threshold = roc_obj$thresholds,
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    Specificity = roc_obj$specificities,
    Sensitivity = roc_obj$sensitivities,
    Analysis = result_name,
    Score_Mode = score_mode,
    Analysis_Config = analysis_config,
    stringsAsFactors = FALSE
  )
  
  roc_file <- make_output_name(
    pathway_name,
    "Univariate",
    "SNP_SumScore",
    "ROC_Data",
    config,
    score_mode
  )
  
  write.csv2(
    roc_df,
    roc_file,
    row.names = FALSE
  )
  
  #-------------------------------------------------------------------------
  # Console output
  #-------------------------------------------------------------------------
  if (score_mode == "lof_burden") {
    
    cat(
      "LoF model: biologically oriented allele burden.\n"
    )
    
    cat(
      "Interpretation OR per unit: OR per additional ",
      "LoF-Allelkopie.\n",
      sep = ""
    )
  }
  cat("Score mode:", score_mode, "\n")
  cat("Analysis configuration:", analysis_config, "\n")
  cat("Score-SNPs:", n_score_snps, "\n")
  cat("AUC Test:", round(auc_val, 3), "\n")
  cat("OR per score unit:", round(OR, 3), "\n")
  cat("OR per SD:", round(OR_per_sd, 3), "\n")
  cat("Saved:", regression_file, "\n")
  
  return(
    list(
      results = results_out,
      model = model_train,
      model_standardized = model_train_sd,
      roc_train = roc_train,
      roc = roc_obj,
      predictions = prediction_table,
      roc_data = roc_df,
      batch_sensitivity = batch_sensitivity_results,
      score_mean_train = score_mean_train,
      score_sd_train = score_sd_train,
      n_score_snps = n_score_snps
    )
  )
}

# ---- run_covariate_precheck() ----
library(margins)
library(pscl)
library(ggeffects)
library(ggplot2)

# COVARIATE PRE-CHECK
run_covariate_precheck <- function(df, cov_categorical, cov_numeric, config, score_mode) {
  
  cat("Starting covariate pre-analyses...\n")
  
  #============================================================================
  # Consistency check of the analysis configuration
  #============================================================================
  
  if (
    score_mode == "lof_burden" &&
    config$type != "lof_vep"
  ) {
    stop(
      "score_mode = 'lof_burden' may only be used in the covariate pre-check ",
      "with config$type = 'lof_vep'."
    )
  }
  
  if (
    config$type == "lof_vep" &&
    score_mode != "lof_burden"
  ) {
    stop(
      "A LoF VEP configuration may only be used in the covariate pre-check ",
      "with score_mode = 'lof_burden'."
    )
  }
  
  #============================================================================
  # 1) Remove outliers of numeric covariates directly
  #============================================================================
  
  outlier_log <- list()
  outlier_rows_all <- c()
  
  n_before_outlier_removal <- nrow(df)
  
  for (v in cov_numeric) {
    
    x <- df[[v]]
    
    Q1 <- quantile(x, 0.25, na.rm = TRUE)
    Q3 <- quantile(x, 0.75, na.rm = TRUE)
    IQRv <- IQR(x, na.rm = TRUE)
    
    lower <- Q1 - 1.5 * IQRv
    upper <- Q3 + 1.5 * IQRv
    
    outlier_idx <- which(x < lower | x > upper)
    
    outlier_log[[v]] <- data.frame(
      Variable = v,
      n_outliers = length(outlier_idx),
      Lower_Cutoff = lower,
      Upper_Cutoff = upper
    )
    
    outlier_rows_all <- unique(c(outlier_rows_all, outlier_idx))
  }
  
  if (length(outlier_rows_all) > 0) {
    df <- df[-outlier_rows_all, , drop = FALSE]
  }
  
  n_after_outlier_removal <- nrow(df)
  
  outlier_df <- bind_rows(outlier_log) %>%
    mutate(
      n_before_outlier = n_before_outlier_removal,
      n_after_outlier = n_after_outlier_removal,
      n_removed = n_before_outlier_removal - n_after_outlier_removal
    )
  
  write.csv2(
    outlier_df,
    make_output_name(
      pathway_name,
      "Multivariate",
      "Covariate_Precheck",
      "Outliers",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  cat("\n========================================\n")
  cat("Outlier removal numeric covariates\n")
  cat("========================================\n")
  if (nrow(outlier_df) > 0) {
    
    for (i in seq_len(nrow(outlier_df))) {
      
      cat(
        outlier_df$Variable[i],
        ": ",
        outlier_df$n_outliers[i],
        " Outlier\n",
        sep = ""
      )
    }
  }
  cat("Before outlier removal: ", n_before_outlier_removal, "\n")
  cat("After outlier removal:", n_after_outlier_removal, "\n")
  cat("Removed:               ", n_before_outlier_removal - n_after_outlier_removal, "\n")
  
  #============================================================================
  # 2) Dataset for pre-analyses
  #============================================================================
  
  df_clean <- df %>%
    select(
      id,
      obesity,
      SNP_SumScore,
      all_of(cov_categorical),
      all_of(cov_numeric)
    )
  
  if ("Batch" %in% colnames(df)) {
    
    df_clean$Batch <- factor(
      df$Batch
    )
    
    cov_categorical_all <- unique(
      c(
        cov_categorical,
        "Batch"
      )
    )
    
  } else {
    
    cov_categorical_all <- cov_categorical
  }
  
  #============================================================================
  # 3) Numeric correlation
  #============================================================================
  
  cov_numeric_present <- cov_numeric[
    cov_numeric %in% colnames(df_clean)
  ]
  
  if (length(cov_numeric_present) < 2) {
    
    warning(
      "Fewer than two numeric covariates available. ",
      "Correlation analysis is skipped."
    )
    
    cor_mat <- matrix(
      numeric(0),
      nrow = 0,
      ncol = 0
    )
    
    cor_summary <- data.frame(
      Variable = character(),
      MaxCorrelation = numeric(),
      MaxCor_with = character(),
      HighCorrelation = character(),
      stringsAsFactors = FALSE
    )
    
  } else {
    
    cor_mat <- cor(
      df_clean[
        ,
        cov_numeric_present,
        drop = FALSE
      ],
      use = "pairwise.complete.obs",
      method = "pearson"
    )
    
    write.csv2(
      cor_mat,
      make_output_name(
        pathway_name,
        "Multivariate",
        "Covariate_Precheck",
        "Correlation_Numeric",
        config,
        score_mode
      )
    )
    
    cor_long <- as.data.frame(as.table(cor_mat)) %>%
      filter(Var1 != Var2) %>%
      mutate(abs_cor = abs(Freq))
    
    cor_summary <- cor_long %>%
      group_by(Var1) %>%
      slice_max(abs_cor, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(
        Variable = as.character(Var1),
        MaxCorrelation = Freq,
        MaxCor_with = as.character(Var2),
        HighCorrelation = case_when(
          abs(MaxCorrelation) > 0.90 ~ "CRITICAL",
          abs(MaxCorrelation) > 0.80 ~ "PROBLEMATIC",
          abs(MaxCorrelation) >= 0.50 ~ "OBSERVE",
          TRUE ~ "OK"
        )
      )
  }
  
  #============================================================================
  # 4) Cramer's V
  #============================================================================
  
  cramer_results <- data.frame(
    Var1 = character(),
    Var2 = character(),
    CramersV = numeric(),
    stringsAsFactors = FALSE
  )
  
  if (length(cov_categorical_all) >= 2) {
    for (i in 1:(length(cov_categorical_all) - 1)) {
      for (j in (i + 1):length(cov_categorical_all)) {
        
        tbl <- table(
          df_clean[[cov_categorical_all[i]]],
          df_clean[[cov_categorical_all[j]]]
        )
        
        cv <- tryCatch(
          assocstats(tbl)$cramer,
          error = function(e) NA
        )
        
        cramer_results <- rbind(
          cramer_results,
          data.frame(
            Var1 = cov_categorical_all[i],
            Var2 = cov_categorical_all[j],
            CramersV = cv
          )
        )
      }
    }
  }
  
  write.csv2(
    cramer_results,
    make_output_name(
      pathway_name,
      "Multivariate",
      "Covariate_Precheck",
      "CramersV",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  if (nrow(cramer_results) > 0) {
    
    cramer_summary <- bind_rows(
      
      cramer_results %>%
        transmute(
          Variable = Var1,
          MaxCramersV_with = Var2,
          CramersV = CramersV
        ),
      
      cramer_results %>%
        transmute(
          Variable = Var2,
          MaxCramersV_with = Var1,
          CramersV = CramersV
        )
      
    ) %>%
      filter(
        is.finite(CramersV)
      ) %>%
      group_by(
        Variable
      ) %>%
      slice_max(
        CramersV,
        n = 1,
        with_ties = FALSE
      ) %>%
      ungroup() %>%
      transmute(
        Variable,
        MaxCramersV = CramersV,
        MaxCramersV_with,
        HighCramerV = case_when(
          MaxCramersV > 0.80 ~ "PROBLEMATIC",
          MaxCramersV > 0.50 ~ "CRITICAL_CHECK",
          MaxCramersV >= 0.30 ~ "OBSERVE",
          TRUE ~ "OK"
        )
      )
    
  } else {
    
    cramer_summary <- data.frame(
      Variable = character(),
      MaxCramersV = numeric(),
      MaxCramersV_with = character(),
      HighCramerV = character(),
      stringsAsFactors = FALSE
    )
  }
  
  #============================================================================
  # 5) VIF
  #============================================================================
  
  vif_covariates <- c(cov_categorical, cov_numeric)
  
  if ("Batch" %in% colnames(df_clean)) {
    vif_covariates <- unique(
      c(
        vif_covariates,
        "Batch"
      )
    )
  }
  
  formula_vif <- as.formula(
    paste("obesity ~", paste(vif_covariates, collapse = " + "))
  )
  
  model_vif <- glm(
    formula_vif,
    data = df_clean,
    family = binomial
  )
  
  vif_values <- car::vif(model_vif)
  
  if (is.matrix(vif_values)) {
    vif_df <- data.frame(
      Variable = rownames(vif_values),
      VIF = vif_values[, "GVIF"]^(1 / (2 * vif_values[, "Df"]))
    )
  } else {
    vif_df <- data.frame(
      Variable = names(vif_values),
      VIF = as.numeric(vif_values)
    )
  }
  
  vif_df <- vif_df %>%
    mutate(
      VIF_flag = case_when(
        VIF > 10 ~ "CRITICAL",
        VIF >= 5 ~ "PROBLEMATIC",
        VIF >= 3 ~ "MODERATE",
        TRUE ~ "OK"
      )
    )
  
  write.csv2(
    vif_df,
    make_output_name(
      pathway_name,
      "Multivariate",
      "Covariate_Precheck",
      "VIF",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  #============================================================================
  # 6) Decision Table
  #============================================================================
  
  decision_table <- data.frame(
    Variable = unique(vif_covariates)
  ) %>%
    left_join(vif_df, by = "Variable") %>%
    left_join(cor_summary, by = "Variable") %>%
    left_join(cramer_summary, by = "Variable") %>%
    mutate(
      Recommendation = case_when(
        VIF_flag %in% c("CRITICAL", "PROBLEMATIC") ~ "REVIEW",
        HighCorrelation %in% c("CRITICAL", "PROBLEMATIC") ~ "REVIEW",
        HighCramerV %in% c("CRITICAL_CHECK", "PROBLEMATIC") ~ "REVIEW",
        TRUE ~ "KEEP"
      ),
      Manual_Decision = NA_character_,
      Comment = NA_character_
    )
  
  print(decision_table)
  
  write.csv2(
    decision_table,
    make_output_name(
      pathway_name,
      "Multivariate",
      "Covariate_Precheck",
      "DecisionTable",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  cat("Covariate pre-analyses completed.\n")
  
  return(
    list(
      df = df,
      outlier_df = outlier_df,
      cor_summary = cor_summary,
      cramer_results = cramer_results,
      cramer_summary = cramer_summary,
      vif_df = vif_df,
      decision_table = decision_table,
      cov_categorical_all = cov_categorical_all
    )
  )
}

# ---- run_multivariate_sumscore_analysis() ----
run_multivariate_sumscore_analysis <- function(df_score, result_name) {
  
  cat("\n====================================================\n")
  cat("Multivariable analysis:", result_name, "\n")
  cat("====================================================\n")
  
  if (is.null(df_score)) {
    warning("No SNP sum score dataset for ", result_name)
    return(NULL)
  }
  
  score_mode <- unique(
    as.character(df_score$Score_Mode)
  )
  
  score_mode <- score_mode[
    !is.na(score_mode)
  ]
  
  analysis_config <- unique(
    as.character(df_score$Analysis_Config)
  )
  
  analysis_config <- analysis_config[
    !is.na(analysis_config)
  ]
  
  if (
    length(score_mode) != 1 ||
    length(analysis_config) != 1
  ) {
    
    warning(
      "Score_Mode or Analysis_Config is not unique for ",
      result_name,
      ""
    )
    
    return(NULL)
  }
  
  score_mode <- score_mode[1]
  analysis_config <- analysis_config[1]
  
  if (!analysis_config %in% names(analysis_configs)) {
    warning(
      "Unknown analysis configuration: ",
      analysis_config
    )
    return(NULL)
  }
  
  config <- analysis_configs[[analysis_config]]
  
  #============================================================================
  # Consistency check of analysis and score configuration
  #============================================================================
  
  if (
    score_mode == "lof_burden" &&
    config$type != "lof_vep"
  ) {
    warning(
      "Invalid combination for ",
      result_name,
      ": score_mode = 'lof_burden' is only permitted with ",
      "config$type = 'lof_vep'."
    )
    return(NULL)
  }
  
  if (
    config$type == "lof_vep" &&
    score_mode != "lof_burden"
  ) {
    warning(
      "Invalid combination for ",
      result_name,
      ": a LoF VEP configuration may only be used with ",
      "score_mode = 'lof_burden'."
    )
    return(NULL)
  }
  
  if (score_mode == "lof_burden") {
    cat(
      "LoF burden model: the burden already oriented to the VEP LoF allele ",
      "is analysed in the multivariable model.\n",
      sep = ""
    )
  }
  
  #============================================================================
  # 1) Base data
  #============================================================================
  
  df <- df_score
  
  stopifnot("id" %in% colnames(df))
  stopifnot("obesity" %in% colnames(df))
  stopifnot("SNP_SumScore" %in% colnames(df))
  
  #============================================================================
  # Check score variation
  #============================================================================
  
  score_values <- suppressWarnings(
    as.numeric(
      df$SNP_SumScore
    )
  )
  
  finite_score_values <- score_values[
    is.finite(score_values)
  ]
  
  if (
    length(finite_score_values) == 0 ||
    length(unique(finite_score_values)) < 2
  ) {
    warning(
      "The genetic score has insufficient variation for ",
      result_name,
      "."
    )
    return(NULL)
  }
  
  if (score_mode == "lof_burden") {
    
    cat(
      "LoF burden 0 for ",
      sum(finite_score_values == 0),
      " of ",
      length(finite_score_values),
      " individuals.\n",
      sep = ""
    )
  }
  
  #============================================================================
  # 2) Recoding
  #============================================================================
  
  df$gender_num <- as.numeric(as.character(df$gender))
  df$gender_num[df$gender_num %in% c(0, 7777777, 9999999)] <- NA
  
  df$gender_cat <- factor(
    df$gender_num,
    levels = c(1, 2),
    labels = c("male", "female")
  )
  
  df$smoking_habits_3m_num <- as.numeric(as.character(df$smoking_habits_3m))
  df$smoking_habits_3m_num[
    df$smoking_habits_3m_num %in% c(
      0,
      7777777,
      9999999
    )
  ] <- NA
  
  df$smoking_cat <- factor(
    dplyr::case_when(
      df$smoking_habits_3m_num %in% c(1, 2) ~ "never_smoker",
      df$smoking_habits_3m_num == 3 ~ "former_smoker",
      df$smoking_habits_3m_num == 4 ~ "current_smoker",
      TRUE ~ NA_character_
    ),
    levels = c(
      "never_smoker",
      "former_smoker",
      "current_smoker"
    )
  )
  
  df$school_certificate_A_num <- as.numeric(as.character(df$school_certificate_A))
  df$school_certificate_A_num[!df$school_certificate_A_num %in% 1:5] <- NA
  
  df$school_cert_cat <- factor(
    ifelse(df$school_certificate_A_num %in% c(1, 2, 5), "risk", "reference"),
    levels = c("reference", "risk")
  )
  
  df$graduation_num <- as.numeric(as.character(df$graduation))
  df$graduation_num[df$graduation_num %in% c(6, 7777777, 9999999)] <- NA
  
  df$graduation_cat <- factor(
    ifelse(df$graduation_num %in% c(1, 2, 5), "risk", "reference"),
    levels = c("reference", "risk")
  )
  
  df$sample_age <- as.numeric(as.character(df$sample_age))
  df$sports_activity_min_week <- as.numeric(as.character(df$sports_activity_min_week))
  df$everyday_activity_min_week <- as.numeric(as.character(df$everyday_activity_min_week))
  df$energy_kcal <- as.numeric(gsub(",", ".", as.character(df$energy_kcal)))
  
  #============================================================================
  # 3) Batch
  #============================================================================
  
  batch_col <- NULL
  
  if ("Batch" %in% colnames(df)) {
    
    batch_col <- "Batch"
    
  } else if ("batch" %in% colnames(df)) {
    
    batch_col <- "batch"
  }
  
  if (!is.null(batch_col)) {
    
    df$Batch <- factor(
      as.character(df[[batch_col]]),
      levels = c(
        "Ambu",
        "MD"
      )
    )
    
  } else {
    
    df$Batch <- NULL
  }
  
  #============================================================================
  # 4a) Load genome-wide PCs (population stratification)
  #============================================================================
  # Not computed per pathway from the pathway SNPs, but
  # generated once in 06_global_PCA.R from the complete GSA array (QC + LD pruning +
  # KING-cutoff sensitivity check + PCA) and loaded here unchanged
  # (see genomewide_pcs_file above).
  
  df$id <- as.character(df$id)
  
  global_gPCs <- readRDS(genomewide_pcs_file)
  global_gPCs$id <- as.character(global_gPCs$id)
  
  n_before_gpc_join <- nrow(df)
  
  df <- merge(
    df,
    global_gPCs,
    by = "id",
    all.x = FALSE
  )
  
  if (nrow(df) < n_before_gpc_join) {
    warning(
      n_before_gpc_join - nrow(df),
      " individuals in ", result_name,
      " have no genome-wide PCs and were removed."
    )
  }
  
  if (nrow(df) == 0) {
    stop(
      "No individuals with genome-wide PCs for ",
      result_name,
      "."
    )
  }
  
  #============================================================================
  # 4b) Take over outliers from the covariate pre-check
  #============================================================================
  
  # Number after genotype/phenotype matching and PCA
  n_after_id_matching <- nrow(df)
  
  # Retrieve matching pre-check result
  precheck <- all_covariate_prechecks[[result_name]]
  
  if (is.null(precheck)) {
    stop(
      "No covariate pre-check for ",
      result_name,
      " available."
    )
  }
  
  if (is.null(precheck$df) || !"id" %in% colnames(precheck$df)) {
    stop(
      "The covariate pre-check contains no cleaned dataset ",
      "with ID column for ",
      result_name,
      "."
    )
  }
  
  # Keep IDs that passed the outlier pre-check
  ids_after_outlier_removal <- unique(
    as.character(precheck$df$id)
  )
  
  df <- df %>%
    filter(
      as.character(id) %in% ids_after_outlier_removal
    )
  
  # Number after outlier removal
  n_after_outlier_removal <- nrow(df)
  
  cat("\n========================================\n")
  cat("Outlier removal applied\n")
  cat("========================================\n")
  cat("After ID matching/PCA:   ", n_after_id_matching, "\n")
  cat("After outlier removal:", n_after_outlier_removal, "\n")
  cat("Removed:               ",
      n_after_id_matching - n_after_outlier_removal,
      "\n")
  
  #============================================================================
  # 5) Covariates
  #============================================================================
  
  covariates_keep <- manual_covariate_decisions$default$keep
  covariates_removed <- manual_covariate_decisions$default$remove
  
  # Use only covariates actually present in df
  covariates_present <- intersect(
    covariates_keep,
    colnames(df)
  )
  
  # Document missing covariates
  missing_covariates <- setdiff(
    covariates_keep,
    colnames(df)
  )
  
  if (length(missing_covariates) > 0) {
    warning(
      "The following selected covariates are missing in df and are not used: ",
      paste(missing_covariates, collapse = ", ")
    )
  }
  
  cov_categorical <- intersect(
    covariates_present,
    c(
      "gender_cat",
      "smoking_cat",
      "school_cert_cat",
      "graduation_cat",
      "Batch"
    )
  )
  
  cov_numeric <- setdiff(
    covariates_present,
    cov_categorical
  )
  
  # data.table-conform selection of factor variables
  factor_vars <- c(
    covariates_present[
      vapply(
        df[
          ,
          covariates_present,
          drop = FALSE
        ],
        is.factor,
        logical(1)
      )
    ],
    "obesity"
  )
  
  factor_vars <- unique(
    factor_vars
  )
  
  #============================================================================
  # 6) Covariate pre-analyses
  #============================================================================
  cat("Covariate selection based on VIF, correlation, Cramer's V and manual content-based decision.\n")
  cat("\n")
  cat("Final covariate selection\n")
  
  cat("KEEP:\n")
  cat(paste(covariates_keep, collapse = ", "), "\n\n")
  
  cat("REMOVE:\n")
  cat(paste(covariates_removed, collapse = ", "), "\n\n")
  
  covariate_selection <- data.frame(
    Variable = c(covariates_keep, covariates_removed),
    Manual_Decision = c(
      rep("KEEP", length(covariates_keep)),
      rep("REMOVE_OR_SENSITIVITY_ONLY", length(covariates_removed))
    ),
    Comment = rep(
      manual_covariate_decisions$default$comment,
      length(c(covariates_keep, covariates_removed))
    )
  )
  
  write.csv2(
    covariate_selection,
    make_output_name(
      pathway_name,
      "Multivariate",
      "Covariate_Precheck",
      "Covariate_Selection",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  #============================================================================
  # 7) Factor-level check
  #============================================================================
  
  # Ensure that only existing variables are checked
  factor_vars <- unique(
    factor_vars[
      factor_vars %in% colnames(df)
    ]
  )
  
  if (length(factor_vars) == 0) {
    
    warning("No factor variables available for the factor-level check.")
    
    factor_level_table <- data.frame(
      Variable = character(),
      Level = character(),
      Count = numeric(),
      stringsAsFactors = FALSE
    )
    
    factor_check_summary <- data.frame(
      Variable = character(),
      IsFactor = logical(),
      N_Levels = integer(),
      Min_Level_Count = numeric(),
      Has_At_Least_2_Levels = logical(),
      Flag = character(),
      stringsAsFactors = FALSE
    )
    
  } else {
    
    factor_level_table <- dplyr::bind_rows(
      lapply(factor_vars, function(v) {
        
        tab <- table(
          df[[v]],
          useNA = "ifany"
        )
        
        data.frame(
          Variable = v,
          Level = names(tab),
          Count = as.numeric(tab),
          stringsAsFactors = FALSE
        )
      })
    )
    
    factor_check_summary <- dplyr::bind_rows(
      lapply(factor_vars, function(v) {
        
        tab <- table(
          df[[v]],
          useNA = "no"
        )
        
        n_levels <- length(tab)
        min_level_count <- if (n_levels > 0) {
          min(as.numeric(tab))
        } else {
          NA_real_
        }
        
        data.frame(
          Variable = v,
          IsFactor = is.factor(df[[v]]),
          N_Levels = n_levels,
          Min_Level_Count = min_level_count,
          Has_At_Least_2_Levels = n_levels >= 2,
          Flag = dplyr::case_when(
            !is.factor(df[[v]]) && v != "obesity" ~ "NOT_FACTOR",
            n_levels < 2 ~ "ONLY_ONE_LEVEL",
            !is.na(min_level_count) && min_level_count < 10 ~ "LOW_LEVEL_COUNT",
            TRUE ~ "OK"
          ),
          stringsAsFactors = FALSE
        )
      })
    )
  }
  
  write.csv2(
    factor_level_table,
    make_output_name(
      pathway_name,
      "Multivariate",
      "Covariate_Precheck",
      "FactorLevels",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  write.csv2(
    factor_check_summary,
    make_output_name(
      pathway_name,
      "Multivariate",
      "Covariate_Precheck",
      "FactorLevel_Check",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  #============================================================================
  # 8) Final analysis dataset
  #============================================================================
  
  required_vars <- c(
    "id",
    "obesity",
    "SNP_SumScore",
    covariates_present
  )
  
  n_before_complete_case <- nrow(df)
  
  df_final <- df %>%
    filter(if_all(all_of(required_vars), ~ !is.na(.)))
  
  n_after_complete_case <- nrow(df_final)
  
  if (nrow(df_final) < 30) {
    warning("Too few complete cases for ", result_name)
    return(NULL)
  }
  
  if (length(unique(df_final$obesity)) < 2) {
    warning("Only one outcome class for ", result_name)
    return(NULL)
  }
  
  write.csv2(
    df_final,
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "Individuals",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  saveRDS(
    df_final,
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "Individuals",
      config,
      score_mode,
      ext = "rds"
    )
  )
  
  #============================================================================
  # Analysis Flow Summary
  #============================================================================
  
  analysis_snps <- attr(df_score, "n_score_snps")
  
  if (is.null(analysis_snps) || is.na(analysis_snps)) {
    analysis_snps <- NA_integer_
  }
  
  analysis_flow <- data.frame(
    Pathway = pathway_name,
    Analysis = result_name,
    Score_Mode = score_mode,
    Analysis_Config = analysis_config,
    Step = c(
      "Final Analysis Dataset",
      "Genomewide PCs Joined",
      "Outlier Removed",
      "Covariate Selection",
      "Complete Cases (Covariates)",
      "Final Analysis Population",
      "Multivariate Regression",
      "SNP Sum Score"
    ),
    SNPs = rep(analysis_snps, 8),
    Individuals = c(
      nrow(df_score),
      n_after_id_matching,
      n_after_outlier_removal,
      n_after_outlier_removal,
      n_after_complete_case,
      nrow(df_final),
      nrow(df_final),
      nrow(df_final)
    ),
    stringsAsFactors = FALSE
  )
  
  print(analysis_flow)
  
  write.csv2(
    analysis_flow,
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "Analysis_Flow_Summary",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  flowchart_prefix <- sub(
    "\\.csv$",
    "",
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "Analysis_Flowchart",
      config,
      score_mode
    )
  )
  
  analysis_flow_plot <- plot_analysis_flow_summary(
    summary_df = analysis_flow,
    pathway_name = pathway_name,
    outfile_prefix = flowchart_prefix,
    score_mode = score_mode,
    analysis_config = analysis_config,
    show_plot = TRUE
  )
  
  numeric_summary <- bind_rows(
    lapply(c("SNP_SumScore", cov_numeric), function(v) {
      
      data.frame(
        Variable = v,
        N = sum(!is.na(df_final[[v]])),
        Mean = mean(df_final[[v]], na.rm = TRUE),
        SD = sd(df_final[[v]], na.rm = TRUE),
        Median = median(df_final[[v]], na.rm = TRUE),
        Q1 = quantile(df_final[[v]], 0.25, na.rm = TRUE),
        Q3 = quantile(df_final[[v]], 0.75, na.rm = TRUE),
        Min = min(df_final[[v]], na.rm = TRUE),
        Max = max(df_final[[v]], na.rm = TRUE),
        Type = "numeric"
      )
    })
  )
  
  categorical_summary <- bind_rows(
    lapply(factor_vars, function(v) {
      
      tab <- table(df_final[[v]])
      
      data.frame(
        Variable = v,
        Level = names(tab),
        N = as.numeric(tab),
        Percent = round(100 * as.numeric(tab) / sum(tab), 2),
        Type = "categorical"
      )
    })
  )
  
  write.csv2(
    numeric_summary,
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "Demographics_Covariates_NumericSummary",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  write.csv2(
    categorical_summary,
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "Demographics_Covariates_CategoricalSummary",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  #============================================================================
  # 9) Multivariable model train/test
  #============================================================================
  
  covariates_formula <- paste(
    covariates_present,
    collapse = " + "
  )
  
  model_formula <- reformulate(
    termlabels = c(
      "SNP_SumScore",
      covariates_present
    ),
    response = "obesity"
  )
  
  model_formula_standardized <- reformulate(
    termlabels = c(
      "SNP_SumScore_SD",
      covariates_present
    ),
    response = "obesity"
  )
  
  global_split <- readRDS(
    output_file(
      "Pathway_",
      pathway_name,
      "_Global_TrainTest_Split.rds"
    )
  )
  
  global_split$train_ids <- as.character(
    global_split$train_ids
  )
  
  global_split$test_ids <- as.character(
    global_split$test_ids
  )
  
  df_final$id <- as.character(
    df_final$id
  )
  
  train_set <- df_final %>%
    filter(id %in% global_split$train_ids)
  
  test_set <- df_final %>%
    filter(id %in% global_split$test_ids)
  
  if (
    nrow(train_set) < 20 ||
    nrow(test_set) < 10
  ) {
    
    warning(
      "Too few individuals in the train or test dataset for ",
      result_name,
      ". Train: ",
      nrow(train_set),
      ", Test: ",
      nrow(test_set)
    )
    
    return(NULL)
  }
  
  if (
    length(unique(train_set$obesity)) < 2 ||
    length(unique(test_set$obesity)) < 2
  ) {
    
    warning(
      "Train or test dataset contains only one outcome class for ",
      result_name
    )
    
    return(NULL)
  }
  
  #------------------------------------------------------------
  # Standardize SNP score using the training set
  #------------------------------------------------------------
  
  score_mean_train <- mean(
    train_set$SNP_SumScore,
    na.rm = TRUE
  )
  
  score_sd_train <- sd(
    train_set$SNP_SumScore,
    na.rm = TRUE
  )
  
  if (
    !is.finite(score_sd_train) ||
    score_sd_train <= 0
  ) {
    
    warning(
      "SNP sum score has insufficient variation for ",
      result_name
    )
    
    return(NULL)
  }
  
  train_set$SNP_SumScore_SD <-
    (
      train_set$SNP_SumScore -
        score_mean_train
    ) /
    score_sd_train
  
  test_set$SNP_SumScore_SD <-
    (
      test_set$SNP_SumScore -
        score_mean_train
    ) /
    score_sd_train
  
  #------------------------------------------------------------
  # Fit models
  #------------------------------------------------------------
  model <- glm(
    model_formula,
    data = train_set,
    family = binomial
  )
  
  model_standardized <- glm(
    model_formula_standardized,
    data = train_set,
    family = binomial
  )
  
  coef_sd <- summary(
    model_standardized
  )$coefficients[
    "SNP_SumScore_SD",
    ,
    drop = FALSE
  ]
  
  ci_sd <- tryCatch(
    exp(
      confint.default(
        model_standardized
      )[
        "SNP_SumScore_SD",
        ,
        drop = TRUE
      ]
    ),
    error = function(e) {
      c(
        NA_real_,
        NA_real_
      )
    }
  )
  
  OR_per_SD <- exp(
    coef_sd[
      "SNP_SumScore_SD",
      "Estimate"
    ]
  )
  
  #------------------------------------------------------------
  # Regression results of the model with original score
  #------------------------------------------------------------
  
  res <- broom::tidy(
    model,
    conf.int = TRUE,
    conf.level = 0.95
  ) %>%
    rename(
      CI_lower = conf.low,
      CI_upper = conf.high
    ) %>%
    mutate(
      OR = exp(estimate),
      OR_CI_lower = exp(CI_lower),
      OR_CI_upper = exp(CI_upper),
      Pathway = pathway_name,
      Analysis = result_name,
      Score_Mode = score_mode,
      Analysis_Config = analysis_config,
      Score_Type = dplyr::case_when(
        score_mode == "oriented" ~
          "Obesity-oriented SNP SumScore",
        
        score_mode == "raw" ~
          "Unweighted raw SNP SumScore",
        
        score_mode == "beta_weighted" ~
          "Beta-weighted genetic score",
        
        score_mode == "ridge_weighted" ~
          "Ridge-weighted genetic score",
        
        score_mode == "elasticnet_weighted" ~
          "Elastic-Net-weighted genetic score",
        
        score_mode == "lof_burden" ~
          "Loss-of-function allele burden",
        
        TRUE ~
          "Genetic score"
      ),
      
      Score_Unit = dplyr::case_when(
        score_mode == "lof_burden" ~
          "One additional LoF allele copy",
        
        score_mode == "oriented" ~
          "One additional obesity-oriented allele",
        
        score_mode == "raw" ~
          "One additional counted allele",
        
        score_mode == "beta_weighted" ~
          "One beta-weighted score unit",
        
        score_mode == "ridge_weighted" ~
          "One ridge-weighted score unit",
        
        score_mode == "elasticnet_weighted" ~
          "One elastic-net-weighted score unit",
        
        TRUE ~
          "One score unit"
      ),
      N_total = nrow(df_final),
      N_train = nrow(train_set),
      N_test = nrow(test_set),
      AIC = AIC(model),
      BIC = BIC(model),
      McFadden_R2 = as.numeric(
        pscl::pR2(model)["McFadden"]
      ),
      Score_Mean_Train = score_mean_train,
      Score_SD_Train = score_sd_train,
      OR_Score_Per_SD = OR_per_SD,
      OR_Score_Per_SD_CI_lower = ci_sd[1],
      OR_Score_Per_SD_CI_upper = ci_sd[2]
    )
  
  write.csv2(
    res,
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "Regression_Results",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  #============================================================================
  # 10) Marginal Effects
  #============================================================================
  
  marg_effects <- tryCatch({
    margins_obj <- margins::margins(model)
    as.data.frame(summary(margins_obj))
  }, error = function(e) NULL)
  
  if (!is.null(marg_effects)) {
    write.csv2(
      marg_effects,
      make_output_name(
        pathway_name,
        "Multivariate",
        "SNP_SumScore",
        "MarginalEffects",
        config,
        score_mode
      ),
      row.names = FALSE
    )
  }
  
  #============================================================================
  # 11) ROC + cutoffs test set
  #============================================================================
  
  train_set$pred_prob <- predict(
    model,
    newdata = train_set,
    type = "response"
  )
  
  test_set$pred_prob <- predict(
    model,
    newdata = test_set,
    type = "response"
  )
  
  roc_train <- roc(
    train_set$obesity,
    train_set$pred_prob,
    quiet = TRUE
  )
  
  roc_obj <- roc(
    test_set$obesity,
    test_set$pred_prob,
    quiet = TRUE
  )
  
  auc_val <- as.numeric(auc(roc_obj))
  
  auc_ci <- tryCatch(
    as.numeric(ci.auc(roc_obj)),
    error = function(e) {
      warning(
        "ci.auc() failed: ",
        conditionMessage(e),
        " -- AUC confidence interval is set to NA."
      )
      c(NA_real_, auc_val, NA_real_)
    }
  )
  
  youden_cutoff <- tryCatch(
    as.numeric(
      coords(
        roc_train,
        x = "best",
        best.method = "youden",
        ret = "threshold"
      )
    ),
    error = function(e) NA_real_
  )
  
  calc_metrics <- function(obs, pred_prob, cutoff) {
    
    pred_class <- ifelse(pred_prob >= cutoff, 1, 0)
    
    TP <- sum(obs == 1 & pred_class == 1)
    TN <- sum(obs == 0 & pred_class == 0)
    FP <- sum(obs == 0 & pred_class == 1)
    FN <- sum(obs == 1 & pred_class == 0)
    
    data.frame(
      Sensitivity = ifelse((TP + FN) > 0, TP / (TP + FN), NA),
      Specificity = ifelse((TN + FP) > 0, TN / (TN + FP), NA),
      PPV = ifelse((TP + FP) > 0, TP / (TP + FP), NA),
      NPV = ifelse((TN + FN) > 0, TN / (TN + FN), NA),
      Accuracy = (TP + TN) / length(obs)
    )
  }
  
  metrics_05 <- calc_metrics(test_set$obesity, test_set$pred_prob, 0.5)
  metrics_youden <- calc_metrics(test_set$obesity, test_set$pred_prob, youden_cutoff)
  
  beta0 <- coef(model)["(Intercept)"]
  beta_snp <- coef(model)["SNP_SumScore"]
  
  calc_score_cutoff <- function(prob_cutoff, beta0, beta_snp) {
    if (is.na(beta_snp) || beta_snp == 0) return(NA)
    (log(prob_cutoff / (1 - prob_cutoff)) - beta0) / beta_snp
  }
  
  score_cutoff_05 <- NA_real_
  score_cutoff_youden <- NA_real_
  
  formula_no_snp <- reformulate(
    termlabels = covariates_present,
    response = "obesity"
  )
  
  model_no_snp <- glm(
    formula_no_snp,
    data = train_set,
    family = binomial
  )
  
  lrt_pvalue <- anova(model_no_snp, model, test = "LRT")$`Pr(>Chi)`[2]
  
  snp_coef <- summary(model)$coefficients["SNP_SumScore", ]
  snp_ci <- exp(confint.default(model)["SNP_SumScore", ])
  
  base_row <- data.frame(
    Pathway = pathway_name,
    Analysis = result_name,
    Score_Mode = unique(df_final$Score_Mode)[1],
    Analysis_Config = unique(df_final$Analysis_Config)[1],
    Score_Type = dplyr::case_when(
      score_mode == "oriented" ~
        "Obesity-oriented SNP SumScore",
      
      score_mode == "raw" ~
        "Unweighted raw SNP SumScore",
      
      score_mode == "beta_weighted" ~
        "Beta-weighted genetic score",
      
      score_mode == "ridge_weighted" ~
        "Ridge-weighted genetic score",
      
      score_mode == "elasticnet_weighted" ~
        "Elastic-Net-weighted genetic score",
      
      score_mode == "lof_burden" ~
        "Loss-of-function allele burden",
      
      TRUE ~
        "Genetic score"
    ),
    
    Score_Unit = dplyr::case_when(
      score_mode == "lof_burden" ~
        "One additional LoF allele copy",
      
      score_mode == "oriented" ~
        "One additional obesity-oriented allele",
      
      score_mode == "raw" ~
        "One additional counted allele",
      
      score_mode == "beta_weighted" ~
        "One beta-weighted score unit",
      
      score_mode == "ridge_weighted" ~
        "One ridge-weighted score unit",
      
      score_mode == "elasticnet_weighted" ~
        "One elastic-net-weighted score unit",
      
      TRUE ~
        "One score unit"
    ),
    Estimate = snp_coef["Estimate"],
    StdError = snp_coef["Std. Error"],
    OR = exp(snp_coef["Estimate"]),
    CI_lower = snp_ci[1],
    CI_upper = snp_ci[2],
    Pvalue = snp_coef["Pr(>|z|)"],
    LRT_pvalue = lrt_pvalue,
    AUC = auc_val,
    AUC_CI_lower = auc_ci[1],
    AUC_CI_median = auc_ci[2],
    AUC_CI_upper = auc_ci[3],
    N_train = nrow(train_set),
    N_test = nrow(test_set),
    Score_Cutoff_Interpretation =
      "Not uniquely defined in multivariable model"
  )
  
  cutoff_table <- bind_rows(
    cbind(
      base_row,
      Cutoff_Type = "Fixed_0.5",
      Cutoff_Prob = 0.5,
      SNP_SumScore_Cutoff = score_cutoff_05,
      metrics_05
    ),
    cbind(
      base_row,
      Cutoff_Type = "Youden",
      Cutoff_Prob = youden_cutoff,
      SNP_SumScore_Cutoff = score_cutoff_youden,
      metrics_youden
    )
  )
  
  write.csv2(
    cutoff_table,
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "Cutoff_Results",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  #============================================================================
  # 12) Batch sensitivity optional – main analysis only
  #============================================================================
  
  batch_results <- NULL
  
  if (
    grepl("^main_ld_", result_name) &&
    "Batch" %in% colnames(train_set) &&
    length(unique(na.omit(train_set$Batch))) > 1
  ) {
    
    
    # Do not include batch twice in the formula
    covariates_without_batch <- setdiff(
      covariates_present,
      "Batch"
    )
    
    formula_no_batch <- reformulate(
      termlabels = c(
        "SNP_SumScore",
        covariates_without_batch
      ),
      response = "obesity"
    )
    
    formula_with_batch <- reformulate(
      termlabels = c(
        "SNP_SumScore",
        "Batch",
        covariates_without_batch
      ),
      response = "obesity"
    )
    
    model_no_batch <- glm(
      formula_no_batch,
      data = train_set,
      family = binomial
    )
    
    model_with_batch <- glm(
      formula_with_batch,
      data = train_set,
      family = binomial
    )
    
    coef_no <- summary(model_no_batch)$coefficients["SNP_SumScore", ]
    coef_ba <- summary(model_with_batch)$coefficients["SNP_SumScore", ]
    
    OR_no <- exp(coef_no["Estimate"])
    OR_ba <- exp(coef_ba["Estimate"])
    
    CI_no <- exp(confint.default(model_no_batch)["SNP_SumScore", ])
    CI_ba <- exp(confint.default(model_with_batch)["SNP_SumScore", ])
    
    batch_results <- data.frame(
      Pathway = pathway_name,
      Analysis = result_name,
      Sensitivity_Type = "Batch_Adjusted_Main_Model_Only",
      OR_no_batch = OR_no,
      CI_low_no_batch = CI_no[1],
      CI_high_no_batch = CI_no[2],
      P_no_batch = coef_no["Pr(>|z|)"],
      OR_with_batch = OR_ba,
      CI_low_with_batch = CI_ba[1],
      CI_high_with_batch = CI_ba[2],
      P_with_batch = coef_ba["Pr(>|z|)"],
      Delta_OR = OR_ba - OR_no,
      Delta_percent = 100 * (OR_ba - OR_no) / OR_no,
      N = nrow(train_set)
    )
    
    write.csv2(
      batch_results,
      make_output_name(
        pathway_name,
        "Multivariate",
        "SNP_SumScore",
        "Batch_Sensitivity",
        config,
        score_mode
      ),
      row.names = FALSE
    )
  }
  
  #============================================================================
  # 13) Plots
  #============================================================================
  
  roc_df <- data.frame(
    TPR = rev(roc_obj$sensitivities),
    FPR = rev(1 - roc_obj$specificities)
  )
  
  roc_plot <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_line(linewidth = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    theme_minimal(base_size = 13) +
    labs(
      title = paste0("Pathway_", pathway_name, " – Multivariate Logistic Regression"),
      subtitle = paste("AUC =", round(auc_val, 3)),
      x = "1 - Specificity",
      y = "Sensitivity"
    )
  
  ggsave(
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "ROC",
      config,
      score_mode,
      ext = "pdf"
    ),
    roc_plot,
    width = 7,
    height = 5
  )
  
  ggsave(
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "ROC",
      config,
      score_mode,
      ext = "png"
    ),
    roc_plot,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  res_standardized <- broom::tidy(
    model_standardized,
    conf.int = TRUE,
    conf.level = 0.95
  ) %>%
    rename(
      CI_lower = conf.low,
      CI_upper = conf.high
    ) %>%
    mutate(
      OR = exp(estimate),
      OR_CI_lower = exp(CI_lower),
      OR_CI_upper = exp(CI_upper)
    )
  
  res_plot <- res_standardized %>%
    filter(term != "(Intercept)") %>%
    filter(
      !is.na(OR),
      !is.na(OR_CI_lower),
      !is.na(OR_CI_upper)
    ) %>%
    mutate(
      term_friendly = case_when(
        term == "SNP_SumScore_SD" ~ case_when(
          score_mode == "oriented" ~
            "Oriented SNP SumScore (per 1 SD)",
          
          score_mode == "raw" ~
            "Raw SNP SumScore (per 1 SD)",
          
          score_mode == "beta_weighted" ~
            "Beta-weighted Genetic Score (per 1 SD)",
          
          score_mode == "ridge_weighted" ~
            "Ridge-weighted Genetic Score (per 1 SD)",
          
          score_mode == "elasticnet_weighted" ~
            "Elastic-Net-weighted Genetic Score (per 1 SD)",
          
          score_mode == "lof_burden" ~
            "Loss-of-Function Allele Burden (per 1 SD)",
          
          TRUE ~
            "Genetic Score (per 1 SD)"
        ),
        term == "gender_catfemale" ~ "Gender female",
        term == "smoking_catformer_smoker" ~ "Former smoker",
        term == "smoking_catcurrent_smoker" ~ "Current smoker",
        term == "school_cert_catrisk" ~ "School certificate risk",
        term == "graduation_catrisk" ~ "Graduation risk",
        term == "sample_age" ~ "Age",
        term == "sports_activity_min_week" ~ "Sports activity",
        term == "everyday_activity_min_week" ~ "Everyday activity",
        term == "energy_kcal" ~ "Energy intake",
        term == "gPC1" ~ "Genomic PC1",
        term == "gPC2" ~ "Genomic PC2",
        term == "gPC3" ~ "Genomic PC3",
        TRUE ~ term
      ),
      is_pc = grepl("^gPC", term),
      direction = case_when(
        OR > 1 ~ "Risk",
        OR < 1 ~ "Protective",
        TRUE ~ "Unclear"
      )
    )

  # Derive sensible axis limits from the non-PC covariates (incl. genetic
  # score), so that the partly extreme PC confidence intervals do not distort
  # the scale and the direction/magnitude of the other covariates remains
  # readable. coord_flip(ylim = ...) only zooms the view - no
  # data points are removed; error bars of PCs extending beyond the
  # edge are visibly clipped at the plot margin.
  res_plot_non_pc <- res_plot %>% filter(!is_pc)
  if (nrow(res_plot_non_pc) > 0) {
    or_lo <- min(res_plot_non_pc$OR_CI_lower, na.rm = TRUE)
    or_hi <- max(res_plot_non_pc$OR_CI_upper, na.rm = TRUE)
    or_pad <- (log(or_hi) - log(or_lo)) * 0.2
    if (!is.finite(or_pad) || or_pad <= 0) or_pad <- 0.2
    or_lo <- exp(log(or_lo) - or_pad)
    or_hi <- exp(log(or_hi) + or_pad)
  } else {
    or_lo <- min(res_plot$OR_CI_lower, na.rm = TRUE)
    or_hi <- max(res_plot$OR_CI_upper, na.rm = TRUE)
  }

  forest_plot <- ggplot(res_plot, aes(x = reorder(term_friendly, OR), y = OR, color = direction)) +
    geom_point(aes(shape = is_pc), size = 3) +
    geom_errorbar(
      aes(ymin = OR_CI_lower, ymax = OR_CI_upper),
      width = 0.2
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    scale_y_log10() +
    scale_color_manual(
      values = c("Risk" = "firebrick", "Protective" = "steelblue", "Unclear" = "gray50"),
      guide = "none"
    ) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    coord_flip(ylim = c(or_lo, or_hi)) +
    theme_minimal(base_size = 13) +
    labs(
      title = paste0(pathway_name, " – Multivariate SNP Sum Model"),
      x = "",
      y = "Odds ratio, log scale, 95% CI (axis adapted to non-PC covariates)"
    )

  ggsave(
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "ForestPlot",
      config,
      score_mode,
      ext = "pdf"
    ),
    forest_plot,
    width = 7,
    height = 5
  )
  
  ggsave(
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "ForestPlot",
      config,
      score_mode,
      ext = "png"
    ),
    forest_plot,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  score_label <- switch(
    score_mode,
    
    oriented =
      "Oriented SNP SumScore",
    
    raw =
      "Raw SNP SumScore",
    
    beta_weighted =
      "Beta-weighted Genetic Score",
    
    ridge_weighted =
      "Ridge-weighted Genetic Score",
    
    elasticnet_weighted =
      "Elastic-Net-weighted Genetic Score",
    
    lof_burden =
      "Loss-of-Function Allele Burden",
    
    "Genetic Score"
  )
  
  partial_plot <- tryCatch({
    
    ggeff <- ggpredict(model, terms = "SNP_SumScore [all]")
    
    ggplot(ggeff, aes(x = x, y = predicted)) +
      geom_line(linewidth = 1) +
      geom_ribbon(
        aes(ymin = conf.low, ymax = conf.high),
        alpha = 0.2
      ) +
      theme_minimal(base_size = 13) +
      labs(
        title = paste0(
          pathway_name,
          " – Adjusted ",
          score_label,
          " Effect"
        ),
        x = score_label,
        y = "Predicted probability of obesity"
      )
    
  }, error = function(e) NULL)
  
  if (!is.null(partial_plot)) {
    
    ggsave(
      make_output_name(
        pathway_name,
        "Multivariate",
        "SNP_SumScore",
        "PartialEffect",
        config,
        score_mode,
        ext = "pdf"
      ),
      partial_plot,
      width = 7,
      height = 5
    )
    
    ggsave(
      make_output_name(
        pathway_name,
        "Multivariate",
        "SNP_SumScore",
        "PartialEffect",
        config,
        score_mode,
        ext = "png"
      ),
      partial_plot,
      width = 7,
      height = 5,
      dpi = 300
    )
  }
  
  if (score_mode == "lof_burden") {
    
    lof_coef <- summary(model)$coefficients[
      "SNP_SumScore",
      ,
      drop = FALSE
    ]
    
    cat(
      "LoF burden OR per additional LoF allele copy:",
      round(
        exp(lof_coef[1, "Estimate"]),
        3
      ),
      "\n"
    )
    
    cat(
      "LoF burden OR per 1 SD:",
      round(
        OR_per_SD,
        3
      ),
      "\n"
    )
  }
  
  cat("Done:", result_name, "\n")
  
  return(
    list(
      regression = res,
      cutoff = cutoff_table,
      batch = batch_results,
      model = model,
      model_standardized = model_standardized,
      roc = roc_obj,
      roc_data = data.frame(
        FPR = 1 - roc_obj$specificities,
        TPR = roc_obj$sensitivities
      ),
      df_final = df_final,
      train_set = train_set,
      test_set = test_set,
      score_mean_train = score_mean_train,
      score_sd_train = score_sd_train,
      covariates_keep = covariates_present
    )
  )
}


# ---------------------------------------------------------------------------
# RESULTS ACCUMULATOR + CHECKPOINTING
# ---------------------------------------------------------------------------

results_accumulator_file <- file.path(
  output_dir,
  "AllPathways_ExternalOrientation_Score_Results_Stage1.csv"
)

if (file.exists(results_accumulator_file)) {
  all_pathway_results <- fread(results_accumulator_file, sep = ";")
  all_pathway_results <- as.data.frame(all_pathway_results, stringsAsFactors = FALSE)
} else {
  all_pathway_results <- data.frame(
    Pathway = character(0), Script_Version = character(0), Status = character(0),
    N_Main_Model_SNPs = integer(0), N_Externally_Matched = integer(0),
    Coverage_Pct = numeric(0), Coverage_Threshold_Met = logical(0),
    N_No_External_Match = integer(0), N_Palindromic_Excluded = integer(0),
    N_No_Allele_Match = integer(0),
    Oriented_OR = numeric(0), Oriented_CI_lower = numeric(0), Oriented_CI_upper = numeric(0),
    Oriented_Pvalue = numeric(0),
    BetaWeighted_OR = numeric(0), BetaWeighted_CI_lower = numeric(0),
    BetaWeighted_CI_upper = numeric(0), BetaWeighted_Pvalue = numeric(0),
    Run_Date = character(0), Error_Message = character(0),
    stringsAsFactors = FALSE
  )
}

already_done <- all_pathway_results$Pathway[
  all_pathway_results$Script_Version == SCRIPT_VERSION &
  all_pathway_results$Status == "OK"
]

cat(sprintf(
  "\n%d of %d pathways already completed under %s - these will be skipped.\n",
  sum(pathway_list %in% already_done), length(pathway_list), SCRIPT_VERSION
))

# ---------------------------------------------------------------------------
# MAIN PER-PATHWAY LOOP
# ---------------------------------------------------------------------------

for (pathway_name in pathway_list) {

  if (pathway_name %in% already_done) {
    cat("Skipping (already completed):", pathway_name, "\n")
    next
  }

  cat("\n\n#####################################################\n")
  cat("PATHWAY:", pathway_name, "\n")
  cat("#####################################################\n")

  skip_reason <- NA_character_
  coverage_stats <- NULL
  multivar_summary <- NULL

  run_outcome <- tryCatch({

    #-----------------------------------------------------------------------
    # 1) Point at the ORIGINAL cached, QC'd, LD-pruned dataset + copy the
    #    one cached artifact (Global_TrainTest_Split.rds) that the reused
    #    functions read via output_file() (i.e. from output_dir).
    #-----------------------------------------------------------------------

    main_ld_dataset  <- pipeline_output_file("Pathway_", pathway_name, "_Combined_QC_LDpruned.rds")
    main_ld_bim_path <- file.path(
      results_dir, "Working_Plink", "plink_pathway_outputs",
      paste0("Working_Plink_Pathway_", pathway_name, "_after_R_QC_LDpruned.bim")
    )

    if (!file.exists(main_ld_dataset)) stop("Cached dataset not found: ", main_ld_dataset)
    if (!file.exists(main_ld_bim_path)) stop(".bim file not found: ", main_ld_bim_path)

    split_src <- pipeline_output_file("Pathway_", pathway_name, "_Global_TrainTest_Split.rds")
    split_dst <- output_file("Pathway_", pathway_name, "_Global_TrainTest_Split.rds")
    if (!file.exists(split_dst) && file.exists(split_src)) {
      file.copy(split_src, split_dst, overwrite = FALSE)
    }
    if (!file.exists(split_src)) stop("Global_TrainTest_Split.rds not found: ", split_src)

    #-----------------------------------------------------------------------
    # 1b) Full-sample fit: external orientation carries no FoCus-internal
    #     information, so there is no leakage risk in fitting on every
    #     available individual. The reused template functions filter on
    #     global_split$train_ids / $test_ids internally; the copied split
    #     object is therefore overridden so that both "train" and "test"
    #     resolve to the full pathway-specific sample. Idempotent.
    #     Test-set-specific outputs (ROC/cutoff on "test_set") become
    #     identical to the training ones and are not reported for this score.
    #-----------------------------------------------------------------------
    split_obj <- readRDS(split_dst)
    full_ids <- union(split_obj$train_ids, split_obj$test_ids)
    split_obj$train_ids <- full_ids
    split_obj$test_ids  <- full_ids
    saveRDS(split_obj, split_dst)

    #-----------------------------------------------------------------------
    # 2) Build the external-orientation regression-equivalent table and
    #    write it to EXACTLY the path run_snp_sumscore() will read
    #    (analysis_configs$main_ld's regression_file), gated by the
    #    pre-specified coverage threshold.
    #-----------------------------------------------------------------------

    pathway_data_check <- readRDS(main_ld_dataset)
    matrix_snps <- colnames(pathway_data_check$genotypes)
    rm(pathway_data_check)

    ext_table_result <- build_external_regression_table(
      pathway_name = pathway_name,
      main_ld_bim_path = main_ld_bim_path,
      matrix_snps = matrix_snps
    )

    coverage_stats <- ext_table_result
    coverage_pct <- ext_table_result$n_matched / ext_table_result$n_total_matrix_snps

    cat(sprintf(
      "Coverage: %d / %d SNPs (%.1f%%) externally orientable.\n",
      ext_table_result$n_matched, ext_table_result$n_total_matrix_snps, 100 * coverage_pct
    ))

    if (!pathway_passes_coverage_threshold(ext_table_result$n_matched, ext_table_result$n_total_matrix_snps)) {
      skip_reason <- "coverage_below_threshold"
      stop("SKIP: coverage threshold not met (>=50% AND >=10 SNPs).")
    }

    #-----------------------------------------------------------------------
    # 3) analysis_configs / analysis_plan - this script: main_ld (oriented) +
    #    sensitivity_betaweighted (beta_weighted) only. Labels match the
    #    internal pipeline's exactly, so output filenames/suffixes line up
    #    for easy side-by-side comparison.
    #-----------------------------------------------------------------------

    analysis_configs <- list(
      main_ld = list(
        dataset = main_ld_dataset,
        label = "01_Main_LDpruned_All",
        type = "all_ld"
      ),
      sensitivity_betaweighted = list(
        dataset = main_ld_dataset,
        label = "05_Sensitivity_BetaWeighted_LDpruned",
        type = "beta_weighted"
      )
    )

    analysis_plan <- data.frame(
      config_name = c("main_ld", "sensitivity_betaweighted"),
      score_mode  = c("oriented", "beta_weighted"),
      stringsAsFactors = FALSE
    )

    results_ld_file <- make_output_name(
      pathway_name, "Univariate", "SNP_Logistic_Regression", "Results",
      analysis_configs$main_ld
    )
    write.csv2(ext_table_result$results, results_ld_file, row.names = FALSE)

    #-----------------------------------------------------------------------
    # 4) Reused (verbatim) orchestration - unchanged from the template. Each
    #    block re-initializes its own accumulator list at its start, so
    #    running this per pathway inside the loop is safe.
    #-----------------------------------------------------------------------

    # ---- run_snp_sumscore() calling loop ----
#===================================================================================
# EXECUTION
#===================================================================================

all_snp_sumscore_results <- list()

for (i in seq_len(nrow(analysis_plan))) {
  
  config_name <- analysis_plan$config_name[i]
  score_mode  <- analysis_plan$score_mode[i]
  
  cat(
    "\n----------------------------------------------------\n",
    "Starting analysis ",
    i,
    " of ",
    nrow(analysis_plan),
    ": ",
    config_name,
    " / ",
    score_mode,
    "\n",
    "----------------------------------------------------\n",
    sep = ""
  )
  
  #------------------------------------------------------------
  # Check configuration
  #------------------------------------------------------------
  
  if (!config_name %in% names(analysis_configs)) {
    
    warning(
      "Analysis configuration not found: ",
      config_name
    )
    
    next
  }
  
  config <- analysis_configs[[config_name]]
  
  #------------------------------------------------------------
  # Safeguard permitted LoF combination once more
  #------------------------------------------------------------
  
  if (
    score_mode == "lof_burden" &&
    config$type != "lof_vep"
  ) {
    
    warning(
      "Invalid combination skipped: ",
      config_name,
      " / ",
      score_mode,
      ". 'lof_burden' is only permitted for type = 'lof_vep'."
    )
    
    next
  }
  
  if (
    config$type == "lof_vep" &&
    score_mode != "lof_burden"
  ) {
    
    warning(
      "Invalid combination skipped: ",
      config_name,
      " / ",
      score_mode,
      ". A LoF configuration may only be used with ",
      "score_mode = 'lof_burden'."
    )
    
    next
  }
  
  #------------------------------------------------------------
  # Result name
  #------------------------------------------------------------
  
  result_name <- paste(
    config_name,
    score_mode,
    sep = "_"
  )
  
  #------------------------------------------------------------
  # Run analysis
  #------------------------------------------------------------
  
  result <- tryCatch(
    
    run_snp_sumscore(
      config_name = config_name,
      config = config,
      score_mode = score_mode
    ),
    
    error = function(e) {
      
      warning(
        "Error for ",
        result_name,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  #------------------------------------------------------------
  # Save result
  #------------------------------------------------------------
  
  if (is.null(result)) {
    
    warning(
      "For ",
      result_name,
      " no result was produced."
    )
    
    next
  }
  
  all_snp_sumscore_results[[result_name]] <- result
  
  cat(
    "Analysis completed: ",
    result_name,
    "\n",
    sep = ""
  )
}

    # ---- run_snp_sumscore_model() calling loop ----
#===================================================================================
# LOOP OVER ALL SNP SUM SCORE RESULTS
#===================================================================================

all_snp_sumscore_model_results <- list()

for (result_name in names(all_snp_sumscore_results)) {
  
  df_current <- all_snp_sumscore_results[[result_name]]
  
  cat(
    "\n----------------------------------------------------\n",
    "Starting SNP sum-score model: ",
    result_name,
    "\n",
    "----------------------------------------------------\n",
    sep = ""
  )
  
  if (is.null(df_current)) {
    
    warning(
      "No SNP sum score dataset for ",
      result_name,
      ". Analysis is skipped."
    )
    
    next
  }
  
  model_result <- tryCatch(
    
    run_snp_sumscore_model(
      df = df_current,
      result_name = result_name
    ),
    
    error = function(e) {
      
      warning(
        "Error in ",
        result_name,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  if (is.null(model_result)) {
    
    warning(
      "For ",
      result_name,
      " no model result was produced."
    )
    
    next
  }
  
  all_snp_sumscore_model_results[[result_name]] <-
    model_result
  
  cat(
    "Model analysis completed: ",
    result_name,
    "\n",
    sep = ""
  )
}

#===================================================================================
# EXPORT OVERALL SUMMARY
#===================================================================================

snp_sumscore_model_summary <- bind_rows(
  lapply(
    all_snp_sumscore_model_results,
    function(x) {
      
      if (
        is.null(x) ||
        is.null(x$results)
      ) {
        return(NULL)
      }
      
      x$results
    }
  )
)

if (nrow(snp_sumscore_model_summary) > 0) {
  
  summary_file <- output_file(
    "Pathway_",
    pathway_name,
    "_Univariate_SNP_SumScore_Model_Results_All_Configurations.csv"
  )
  
  write.csv2(
    snp_sumscore_model_summary,
    summary_file,
    row.names = FALSE
  )
  
  cat("\n====================================================\n")
  cat("Overall summary saved:\n")
  cat(summary_file, "\n")
  cat("Successful model analyses:\n")
  
  summary_columns <- c(
    "Analysis_Config",
    "Score_Mode",
    "Score_Type"
  )
  
  # If Score_Unit was added to results_out,
  # it is also shown in the summary.
  if (
    "Score_Unit" %in%
    colnames(snp_sumscore_model_summary)
  ) {
    summary_columns <- c(
      summary_columns,
      "Score_Unit"
    )
  }
  
  print(
    unique(
      snp_sumscore_model_summary[
        ,
        summary_columns,
        drop = FALSE
      ]
    )
  )
  
  #------------------------------------------------------------
  # Additional check of the LoF analysis
  #------------------------------------------------------------
  
  lof_summary <- snp_sumscore_model_summary %>%
    filter(
      Score_Mode == "lof_burden"
    )
  
  if (nrow(lof_summary) > 0) {
    
    cat("\nLoF burden analysis successfully included:\n")
    
    print(
      lof_summary[
        ,
        intersect(
          c(
            "Analysis",
            "Analysis_Config",
            "Score_Mode",
            "Score_Type",
            "Score_Unit",
            "N_Score_SNPs",
            "N_Train",
            "N_Test",
            "OR_per_Unit",
            "CI_lower_per_Unit",
            "CI_upper_per_Unit",
            "Wald_pvalue",
            "OR_per_SD",
            "Pvalue_per_SD",
            "AUC"
          ),
          colnames(lof_summary)
        ),
        drop = FALSE
      ]
    )
    
  } else {
    
    warning(
      "The overall summary contains no ",
      "LoF-burden analysis."
    )
  }
  
  cat(
    "\nNumber of successful analyses: ",
    length(all_snp_sumscore_model_results),
    "\n",
    sep = ""
  )
  
} else {
  
  warning(
    "No SNP sum score model results could ",
    "be merged."
  )
}

#===================================================================================
# OPTIONAL SEPARATE SUMMARY: ONE RESULT PER ANALYSIS
# Youden row as primary performance row
#===================================================================================

if (
  exists("snp_sumscore_model_summary") &&
  nrow(snp_sumscore_model_summary) > 0
) {
  
  snp_sumscore_model_primary_summary <-
    snp_sumscore_model_summary %>%
    filter(
      Cutoff_Type == "Youden_Training"
    )
  
  if (
    nrow(
      snp_sumscore_model_primary_summary
    ) > 0
  ) {
    
    write.csv2(
      snp_sumscore_model_primary_summary,
      output_file(
        "Pathway_",
        pathway_name,
        "_Univariate_SNP_SumScore_Model_Primary_Results_All_Configurations.csv"
      ),
      row.names = FALSE
    )
  }
}


    # ---- covariate-precheck calling loop ----
# COVARIATE PRE-CHECK LOOP - DECISION TABLES ONLY
all_covariate_prechecks <- list()

for (result_name in names(all_snp_sumscore_results)) {
  
  df_score <- all_snp_sumscore_results[[result_name]]
  
  if (is.null(df_score)) next
  
  cat("\n====================================================\n")
  cat("Covariate pre-check:", result_name, "\n")
  cat("====================================================\n")
  
  df <- df_score
  
  # Recoding as later in the analysis
  df$gender_num <- as.numeric(as.character(df$gender))
  df$gender_num[df$gender_num %in% c(0, 7777777, 9999999)] <- NA
  
  df$gender_cat <- factor(
    df$gender_num,
    levels = c(1, 2),
    labels = c("male", "female")
  )
  
  df$smoking_habits_3m_num <- as.numeric(as.character(df$smoking_habits_3m))
  
  df$smoking_habits_3m_num[
    df$smoking_habits_3m_num %in% c(0, 7777777, 9999999)
  ] <- NA
  
  df$smoking_cat <- factor(
    case_when(
      df$smoking_habits_3m_num %in% c(1, 2) ~ "never_smoker",
      df$smoking_habits_3m_num == 3 ~ "former_smoker",
      df$smoking_habits_3m_num == 4 ~ "current_smoker",
      TRUE ~ NA_character_
    ),
    levels = c("never_smoker", "former_smoker", "current_smoker")
  )
  
  df$school_certificate_A_num <- as.numeric(as.character(df$school_certificate_A))
  df$school_certificate_A_num[!df$school_certificate_A_num %in% 1:5] <- NA
  
  df$school_cert_cat <- factor(
    ifelse(df$school_certificate_A_num %in% c(1, 2, 5), "risk", "reference"),
    levels = c("reference", "risk")
  )
  
  df$graduation_num <- as.numeric(as.character(df$graduation))
  df$graduation_num[df$graduation_num %in% c(6, 7777777, 9999999)] <- NA
  
  df$graduation_cat <- factor(
    ifelse(df$graduation_num %in% c(1, 2, 5), "risk", "reference"),
    levels = c("reference", "risk")
  )
  
  df$sample_age <- as.numeric(as.character(df$sample_age))
  df$sports_activity_min_week <- as.numeric(as.character(df$sports_activity_min_week))
  df$everyday_activity_min_week <- as.numeric(as.character(df$everyday_activity_min_week))
  df$energy_kcal <- as.numeric(gsub(",", ".", as.character(df$energy_kcal)))
  
  #===========================================================================
  # Batch
  #===========================================================================
  
  batch_col <- NULL
  
  # Prefers the spelling "Batch" used in the main dataset.
  # The lower-case variant is recognised for compatibility only.
  if ("Batch" %in% colnames(df)) {
    
    batch_col <- "Batch"
    
  } else if ("batch" %in% colnames(df)) {
    
    batch_col <- "batch"
  }
  
  if (!is.null(batch_col)) {
    
    df$Batch <- factor(
      as.character(df[[batch_col]]),
      levels = c(
        "Ambu",
        "MD"
      )
    )
    
  } else {
    
    df$Batch <- NULL
  }
  
  #===========================================================================
  # Determine score and analysis configuration unambiguously
  #===========================================================================
  
  score_mode <- unique(
    as.character(df$Score_Mode)
  )
  
  score_mode <- score_mode[
    !is.na(score_mode)
  ]
  
  analysis_config <- unique(
    as.character(df$Analysis_Config)
  )
  
  analysis_config <- analysis_config[
    !is.na(analysis_config)
  ]
  
  if (
    length(score_mode) != 1 ||
    length(analysis_config) != 1
  ) {
    
    warning(
      "Score_Mode or Analysis_Config is not unique for ",
      result_name,
      ""
    )
    
    next
  }
  
  score_mode <- score_mode[1]
  analysis_config <- analysis_config[1]
  
  if (!analysis_config %in% names(analysis_configs)) {
    warning(
      "Analysis configuration not found: ",
      analysis_config,
      " for ",
      result_name
    )
    next
  }
  
  config <- analysis_configs[[analysis_config]]
  
  #===========================================================================
  # Consistency check for the LoF configuration
  #===========================================================================
  
  if (
    score_mode == "lof_burden" &&
    config$type != "lof_vep"
  ) {
    warning(
      "Invalid combination in ",
      result_name,
      ": score_mode = 'lof_burden' is only permitted with ",
      "config$type = 'lof_vep'."
    )
    next
  }
  
  if (
    config$type == "lof_vep" &&
    score_mode != "lof_burden"
  ) {
    warning(
      "Invalid combination in ",
      result_name,
      ": a LoF VEP configuration may only be used with ",
      "score_mode = 'lof_burden'."
    )
    next
  }
  
  if (score_mode == "lof_burden") {
    cat(
      "LoF burden analysis: the already biologically oriented ",
      "LoF allele burden is used in the covariate pre-check.\n",
      sep = ""
    )
  }
  
  #===========================================================================
  # Check score for sufficient variation
  #===========================================================================
  
  score_values <- suppressWarnings(
    as.numeric(
      df$SNP_SumScore
    )
  )
  
  score_values <- score_values[
    is.finite(score_values)
  ]
  
  if (
    length(score_values) == 0 ||
    length(unique(score_values)) < 2
  ) {
    warning(
      "The SNP sum score has insufficient variation for ",
      result_name,
      "."
    )
    next
  }
  
  if (score_mode == "lof_burden") {
    
    n_zero_burden <- sum(
      score_values == 0
    )
    
    cat(
      "LoF burden distribution: ",
      n_zero_burden,
      " of ",
      length(score_values),
      " individuals have a burden of 0.\n",
      sep = ""
    )
  }
  
  #===========================================================================
  # Load genome-wide PCs (population stratification)
  #===========================================================================
  # Not computed per pathway from the pathway SNPs, but
  # generated once in 06_global_PCA.R from the complete GSA array (QC + LD pruning +
  # KING-cutoff sensitivity check + PCA) and loaded here unchanged
  # (see genomewide_pcs_file above).
  
  df$id <- as.character(df$id)
  
  global_gPCs <- readRDS(genomewide_pcs_file)
  global_gPCs$id <- as.character(global_gPCs$id)
  
  n_before_gpc_join <- nrow(df)
  
  df <- merge(
    df,
    global_gPCs,
    by = "id",
    all.x = FALSE
  )
  
  if (nrow(df) < n_before_gpc_join) {
    warning(
      n_before_gpc_join - nrow(df),
      " individuals in ", result_name,
      " have no genome-wide PCs and were removed."
    )
  }
  
  if (nrow(df) == 0) {
    warning(
      "No individuals with genome-wide PCs for ",
      result_name
    )
    next
  }
  
  cov_categorical <- c(
    "gender_cat",
    "smoking_cat",
    "school_cert_cat",
    "graduation_cat"
  )
  
  cov_numeric <- c(
    "sample_age",
    "sports_activity_min_week",
    "everyday_activity_min_week",
    "energy_kcal",
    "gPC1",
    "gPC2",
    "gPC3"
  )
  
  precheck <- run_covariate_precheck(
    df = df,
    cov_categorical = cov_categorical,
    cov_numeric = cov_numeric,
    config = config,
    score_mode = score_mode
  )
  
  all_covariate_prechecks[[result_name]] <- precheck
  
  cat("\nAUTOMATIC COVARIATE ASSESSMENT FOR:", result_name, "\n")
  print(precheck$decision_table)
}

    # ---- run_multivariate_sumscore_analysis() calling loop ----

options(nwarnings = 10000)

all_multivar_results <- list()

for (result_name in names(all_snp_sumscore_results)) {
  
  df_score <- all_snp_sumscore_results[[result_name]]
  
  if (is.null(df_score)) {
    warning(
      "No score dataset for ",
      result_name,
      "."
    )
    next
  }
  
  analysis_result <- tryCatch(
    
    run_multivariate_sumscore_analysis(
      df_score = df_score,
      result_name = result_name
    ),
    
    error = function(e) {
      
      warning(
        "Error in the multivariable analysis ",
        result_name,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  if (!is.null(analysis_result)) {
    all_multivar_results[[result_name]] <- analysis_result
  }
}

"main_ld_oriented" %in% names(all_multivar_results)

w <- names(warnings())
w[grepl("main_ld_oriented", w)]

names(all_multivar_results)

multivar_summary <- bind_rows(
  lapply(
    all_multivar_results,
    function(x) {
      if (is.null(x)) return(NULL)
      x$cutoff
    }
  )
)

multivar_summary_file <- output_file(
  "Pathway_",
  pathway_name,
  "_Multivariate_SNP_SumScore_Cutoff_Results_All_Configurations.csv"
)

write.csv2(
  multivar_summary,
  multivar_summary_file,
  row.names = FALSE
)




    #-----------------------------------------------------------------------
    # 5) multivar_summary (built by the reused orchestration above) now
    #    holds one row per (Analysis_Config, Score_Mode, Cutoff_Type) -
    #    OR/CI/Pvalue are identical across Cutoff_Type for the same
    #    config/score_mode (same underlying model), so any one row suffices.
    #-----------------------------------------------------------------------

    "OK"

  }, error = function(e) {
    conditionMessage(e)
  })

  run_date <- format(Sys.Date(), "%Y-%m-%d")

  if (identical(run_outcome, "OK")) {

    get_stat <- function(df, config_name, score_mode_name, col) {
      if (is.null(df) || nrow(df) == 0) return(NA_real_)
      if (!all(c("Analysis_Config", "Score_Mode", col) %in% colnames(df))) return(NA_real_)
      hit <- df[df$Analysis_Config == config_name & df$Score_Mode == score_mode_name, , drop = FALSE]
      if (nrow(hit) == 0) return(NA_real_)
      as.numeric(hit[[col]][1])
    }

    new_row <- data.frame(
      Pathway = pathway_name,
      Script_Version = SCRIPT_VERSION,
      Status = "OK",
      N_Main_Model_SNPs = coverage_stats$n_total_matrix_snps,
      N_Externally_Matched = coverage_stats$n_matched,
      Coverage_Pct = round(100 * coverage_stats$n_matched / coverage_stats$n_total_matrix_snps, 1),
      Coverage_Threshold_Met = TRUE,
      N_No_External_Match = coverage_stats$n_no_external_match,
      N_Palindromic_Excluded = coverage_stats$n_palindromic_excluded,
      N_No_Allele_Match = coverage_stats$n_no_allele_match,
      Oriented_OR = get_stat(multivar_summary, "main_ld", "oriented", "OR"),
      Oriented_CI_lower = get_stat(multivar_summary, "main_ld", "oriented", "CI_lower"),
      Oriented_CI_upper = get_stat(multivar_summary, "main_ld", "oriented", "CI_upper"),
      Oriented_Pvalue = get_stat(multivar_summary, "main_ld", "oriented", "Pvalue"),
      BetaWeighted_OR = get_stat(multivar_summary, "sensitivity_betaweighted", "beta_weighted", "OR"),
      BetaWeighted_CI_lower = get_stat(multivar_summary, "sensitivity_betaweighted", "beta_weighted", "CI_lower"),
      BetaWeighted_CI_upper = get_stat(multivar_summary, "sensitivity_betaweighted", "beta_weighted", "CI_upper"),
      BetaWeighted_Pvalue = get_stat(multivar_summary, "sensitivity_betaweighted", "beta_weighted", "Pvalue"),
      Run_Date = run_date,
      Error_Message = NA_character_,
      stringsAsFactors = FALSE
    )

  } else if (!is.na(skip_reason)) {

    new_row <- data.frame(
      Pathway = pathway_name, Script_Version = SCRIPT_VERSION, Status = "SKIPPED_COVERAGE",
      N_Main_Model_SNPs = if (!is.null(coverage_stats)) coverage_stats$n_total_matrix_snps else NA_integer_,
      N_Externally_Matched = if (!is.null(coverage_stats)) coverage_stats$n_matched else NA_integer_,
      Coverage_Pct = if (!is.null(coverage_stats)) round(100 * coverage_stats$n_matched / coverage_stats$n_total_matrix_snps, 1) else NA_real_,
      Coverage_Threshold_Met = FALSE,
      N_No_External_Match = if (!is.null(coverage_stats)) coverage_stats$n_no_external_match else NA_integer_,
      N_Palindromic_Excluded = if (!is.null(coverage_stats)) coverage_stats$n_palindromic_excluded else NA_integer_,
      N_No_Allele_Match = if (!is.null(coverage_stats)) coverage_stats$n_no_allele_match else NA_integer_,
      Oriented_OR = NA_real_, Oriented_CI_lower = NA_real_, Oriented_CI_upper = NA_real_, Oriented_Pvalue = NA_real_,
      BetaWeighted_OR = NA_real_, BetaWeighted_CI_lower = NA_real_, BetaWeighted_CI_upper = NA_real_, BetaWeighted_Pvalue = NA_real_,
      Run_Date = run_date, Error_Message = "Coverage threshold not met",
      stringsAsFactors = FALSE
    )

  } else {

    warning("ERROR for ", pathway_name, ": ", run_outcome)

    new_row <- data.frame(
      Pathway = pathway_name, Script_Version = SCRIPT_VERSION, Status = "ERROR",
      N_Main_Model_SNPs = NA_integer_, N_Externally_Matched = NA_integer_,
      Coverage_Pct = NA_real_, Coverage_Threshold_Met = NA,
      N_No_External_Match = NA_integer_, N_Palindromic_Excluded = NA_integer_, N_No_Allele_Match = NA_integer_,
      Oriented_OR = NA_real_, Oriented_CI_lower = NA_real_, Oriented_CI_upper = NA_real_, Oriented_Pvalue = NA_real_,
      BetaWeighted_OR = NA_real_, BetaWeighted_CI_lower = NA_real_, BetaWeighted_CI_upper = NA_real_, BetaWeighted_Pvalue = NA_real_,
      Run_Date = run_date, Error_Message = substr(run_outcome, 1, 500),
      stringsAsFactors = FALSE
    )
  }

  all_pathway_results <- all_pathway_results[all_pathway_results$Pathway != pathway_name, , drop = FALSE]
  all_pathway_results <- rbind(all_pathway_results, new_row)

  # Checkpoint after EVERY pathway (as in 08) so an interrupted
  # run only loses the pathway not yet completed.
  write.csv2(all_pathway_results, results_accumulator_file, row.names = FALSE)
}

cat("\n\n=====================================================\n")
cat("Batch complete.\n")
print(table(all_pathway_results$Status))
cat("=====================================================\n")

# ---------------------------------------------------------------------------
# AUTO-LOG to Methods_Results_Log_External_SNP_Orientation.md, Appendix B
# ---------------------------------------------------------------------------

log_rows <- all_pathway_results[all_pathway_results$Script_Version == SCRIPT_VERSION, ]
log_rows_ok <- log_rows[log_rows$Status == "OK", ]

if (nrow(log_rows_ok) > 0) {

  fmt_or_ci <- function(or, lo, hi, p) {
    if (is.na(or)) return("n/a")
    sprintf("%.2f [%.2f-%.2f], p=%.3g", or, lo, hi, p)
  }

  md_rows <- data.frame(
    Pathway = log_rows_ok$Pathway,
    Coverage = paste0(log_rows_ok$N_Externally_Matched, "/", log_rows_ok$N_Main_Model_SNPs,
                       " (", log_rows_ok$Coverage_Pct, "%)"),
    Oriented_OR_CI_p = mapply(fmt_or_ci, log_rows_ok$Oriented_OR, log_rows_ok$Oriented_CI_lower,
                               log_rows_ok$Oriented_CI_upper, log_rows_ok$Oriented_Pvalue),
    BetaWeighted_OR_CI_p = mapply(fmt_or_ci, log_rows_ok$BetaWeighted_OR, log_rows_ok$BetaWeighted_CI_lower,
                                   log_rows_ok$BetaWeighted_CI_upper, log_rows_ok$BetaWeighted_Pvalue),
    Run_Date = log_rows_ok$Run_Date,
    stringsAsFactors = FALSE
  )

  upsert_marker_table(
    methods_log_path,
    "<!-- AUTOLOG-EXTORIENT-PHASE6-STAGE1:BEGIN -->",
    "<!-- AUTOLOG-EXTORIENT-PHASE6-STAGE1:END -->",
    key_col = "Pathway",
    new_rows_df = md_rows
  )
}

log_phase6_sessioninfo()

cat("\nAuto-logged to:", methods_log_path, "\n")

