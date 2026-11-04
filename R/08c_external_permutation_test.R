# Pathway-restricted polygenic scores for energy-substrate metabolism and obesity (FoCus cohort)
# Script: 08c_external_permutation_test.R | Source file: Protocol_SNP_Analysis_FoCus_Pathway_Phase7_External_Permutation_Test_REPAIR.R
# Local paths replaced by the environment variable FOCUS_BASE_DIR (project root); analysis code unchanged.
# Individual-level FoCus data are not included (see README).

# ==============================================================================
# REPAIR VERSION (2026-10-09): copy of Protocol_SNP_Analysis_FoCus_Pathway_Phase7_External_Permutation_Test.R
# Re-run after correction of the defective phenotype CSVs. Changes ONLY to paths
# (lines marked "# REPAIR"); analysis code unchanged. Writes only to FoCus/Repair.
# ==============================================================================

# ==============================================================================
# Phase7_External_Permutation_Test.R  (v1.0)
#
# Purpose: external counterpart to Phase4_Permutation_Test.R (Q1, competitive
# permutation test). Tests whether the internal Q1 null finding (none of the
# robust candidates shows enrichment relative to a random,
# size-matched SNP set) is independent of the reorientation
# bias explanation -- i.e. whether it reproduces under EXTERNAL (non-circular)
# orientation (manuscript section [Outlook]: "applying the
# competitive permutation test to the externally-oriented scores should, if
# the account developed in this Discussion is correct, again yield a
# predominantly null result").
#
# Target population (default, see ROBUST_ONLY below): the robust
# candidates (Robust_Candidate == TRUE in AllPathways_Results_Tracking_Table.
# csv) that also have status "OK" in the external stage-1 run (Phase6_v1.3_Stage1)
# (i.e. meet the coverage threshold, methods log section
# 4.4). ROBUST_ONLY <- FALSE instead tests all externally testable pathways
# pathways -- easy to switch if desired.
#
# Design (agreed in the project, 2026-09-05):
#   - Draw size per permutation = FULL LD-pruned pathway SNP number before
#     any external filtering (option A) -- identical to the number that also
#     is drawn internally (phase 4). The drawn SNPs are then subjected to
#     EXACTLY THE SAME external orientation/filter logic as in the
#     real pathway (Phase6/Phase5: rsID filter, GIANT/Yengo match,
#     palindrome exclusion). The usable SNP number per draw therefore varies
#     naturally (mirrors exactly what also happens for the real pathway)
#     -- NO oversampling, no rejection of draws with low
#     coverage (that would distort the null distribution by a selection
#     that does not exist for the real pathway).
#   - Model is fitted on the FULL sample (no train/test split),
#     identical to the corrected full-sample fix in Phase6 v1.3 (methods log
#     section 7): the split file is read from Outputs_Ext, where
#     train_ids/test_ids are already set to the full sample.
#   - B_FINAL = 1000, identical deterministic seed base
#     (pathway_seed_base()) as internally -- the random draws are
#     thus aligned 1:1 with the internal permutations (same column
#     indices per (pathway, b)), only the orientation differs.
#   - Genome-wide sampling pool (QC-pruned, minus whole pathway genes) and the
#     associated cache/imputation logic are taken over UNCHANGED from phase 4
#     (same .bed/.rds cache files, no recomputation).
#
# Checkpointing/resume and error handling per pathway: identical to phase 4
# (status "OK" + same B in the checkpoint CSV is skipped;
# FORCE_RERUN <- TRUE forces a re-run of all pathways; an error in one
# pathway does not abort the batch).
# ==============================================================================

library(bigsnpr)
library(bigstatsr)
library(dplyr)

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
main_project_base <- Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus")

orig_output_dir <- file.path(main_project_base, "Repair", "OutPW")   # REPAIR   # genotypes/ranges (reused unchanged)
output_dir      <- file.path(main_project_base, "Repair", "OutExt")   # REPAIR       # external split file + new results
ranges_dir      <- file.path(main_project_base, "Working_Plink")
genomewide_pcs_file <- file.path(main_project_base, "Global_PCA", "Global_Genomewide_PCs.rds")
pool_bed_file       <- file.path(main_project_base, "Global_PCA", "FO_GSA_Merged_QC_Pruned.bed")

sensitivity_base <- file.path(main_project_base, "Sensitivity_External_SNP_Orientation")
methods_log_path <- file.path(main_project_base, "Repair", "Methods_Results_Log_Repair.md")   # REPAIR
tracking_table_file <- file.path(orig_output_dir, "AllPathways_Results_Tracking_Table.csv")
external_stage1_file <- file.path(output_dir, "AllPathways_ExternalOrientation_Score_Results_Stage1.csv")

external_lookup_rds <- file.path(
  sensitivity_base, "External_SNP_Data",
  "External_Orientation_GIANT_Yengo2018_Locke2015.rds"
)

orig_output_file <- function(...) file.path(orig_output_dir, paste0(...))
output_file      <- function(...) file.path(output_dir, paste0(...))
ranges_file <- function(pathway_name) {
  file.path(main_project_base, "Repair", "Working_Plink", "plink_pathway_outputs",   # REPAIR
            paste0("Working_Plink_Pathway_", pathway_name, "_after_R_QC_LDpruned.bim")
  )
}
gene_ranges_file <- function(pathway_name) {
  file.path(ranges_dir, paste0("pathway_", pathway_name, "_gene_ranges.txt"))
}

checkpoint_file <- output_file("AllPathways_External_Permutation_Test_Results.csv")

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
B_FINAL <- 1000
FORCE_RERUN <- FALSE     # TRUE = recompute all pathways, including those already completed
ROBUST_ONLY <- TRUE      # TRUE (default) = only robust candidates with external status "OK"
                         # FALSE = all externally testable (status "OK") pathways

REQUIRED_EXTERNAL_SCRIPT_VERSION <- "Phase6_v1.3_Stage1"   # full-sample-corrected version

cov_categorical <- c("gender_cat", "smoking_cat", "school_cert_cat", "graduation_cat")
cov_numeric <- c(
  "sample_age", "sports_activity_min_week", "everyday_activity_min_week",
  "energy_kcal", "gPC1", "gPC2", "gPC3"
)

pathway_seed_base <- function(pathway_name) {
  sum(utf8ToInt(pathway_name)) * 100000L
}

# Section 4.1 (Phase5, verbatim): genuine-rsID-Filter
is_genuine_rsid <- function(x) grepl("^rs[0-9]+$", x)

# Section 4.4 (Phase5/Phase6, verbatim): dual coverage threshold
MIN_COVERAGE_PCT <- 0.50
MIN_SNP_COUNT    <- 10
pathway_passes_coverage_threshold <- function(n_matched, n_total_main_model) {
  (n_matched / n_total_main_model) >= MIN_COVERAGE_PCT & n_matched >= MIN_SNP_COUNT
}

# Section 7 (Phase5, verbatim): map external beta/effect allele to the
# matrix-coded alleles of the respective SNP set, flag palindromes,
# derive inversion need. Reused here identically to Phase6
# -- only the input table (bim_table/matrix_snps) differs:
# in Phase6 it is the respective pathway, here the per-permutation
# randomly drawn SNPs.
create_external_orientation_table <- function(external_lookup, bim_table, matrix_snps) {

  bim_table_rsid_only <- bim_table %>%
    mutate(Is_Genuine_RsID = is_genuine_rsid(SNP))

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

# ------------------------------------------------------------------
# Load robust candidates + external status dynamically (single source of
# truth, identical to the principle in phase 4)
# ------------------------------------------------------------------
if (!file.exists(tracking_table_file)) stop("Tracking table not found: ", tracking_table_file)
if (!file.exists(external_stage1_file)) stop("External stage-1 results not found: ", external_stage1_file)

tracking_table <- read.csv2(tracking_table_file, stringsAsFactors = FALSE)
if (!"Robust_Candidate" %in% colnames(tracking_table)) {
  stop("Column 'Robust_Candidate' not in ", tracking_table_file, " found.")
}
robust_pathways <- tracking_table$Pathway[
  tracking_table$Robust_Candidate %in% c("TRUE", "True", "1", "WAHR", TRUE)
]

ext_stage1 <- read.csv2(external_stage1_file, stringsAsFactors = FALSE)
ok_ext_pathways <- ext_stage1$Pathway[
  ext_stage1$Status == "OK" & ext_stage1$Script_Version == REQUIRED_EXTERNAL_SCRIPT_VERSION
]

if (length(ok_ext_pathways) == 0) {
  stop(
    "No pathways with status OK and Script_Version == '",
    REQUIRED_EXTERNAL_SCRIPT_VERSION, "' in ", external_stage1_file,
    " found -- run/check Phase6 (full-sample version) first."
  )
}

PATHWAYS_FINAL <- if (ROBUST_ONLY) {
  intersect(robust_pathways, ok_ext_pathways)
} else {
  ok_ext_pathways
}

cat(
  "Externer Permutationstest: ", length(PATHWAYS_FINAL), " Pathway(s) geladen ",
  "(ROBUST_ONLY=", ROBUST_ONLY, "; robust candidates total: ", length(robust_pathways),
  "; externally OK total: ", length(ok_ext_pathways), ").\n",
  "  ", paste(PATHWAYS_FINAL, collapse = ", "), "\n",
  sep = ""
)

#==================================================================================
# 1) Load + impute the genome-wide QC-pruned pool ONCE (with cache) --
#    taken over unchanged from phase 4 (identical cache files, no
#    recomputation).
#==================================================================================
pool_rds_file <- sub("\\.bed$", ".rds", pool_bed_file)
# REPAIR: pool caches in Global_PCA read only - never regenerate/overwrite
stopifnot(file.exists(pool_rds_file), file.exists(sub("\\.rds$", "_imputed_meta.rds", pool_rds_file)))   # REPAIR

if (!file.exists(pool_rds_file)) {
  cat("Creating backing file from", pool_bed_file, "(one-time, takes a while)...\n")
  snp_readBed(pool_bed_file)
}

pool <- snp_attach(pool_rds_file)
map_pool <- pool$map
fam_pool <- pool$fam

stopifnot(all(c("chromosome", "physical.pos", "marker.ID") %in% colnames(map_pool)))
map_pool$chromosome <- as.character(map_pool$chromosome)

# PLEASE CHECK (analogous to Phase4/Phase6: confirm-don't-assume): bigsnpr
# names the two allele columns of a loaded .bed/.bim pair
# by default "allele1"/"allele2" (columns 5/6 of the .bim, in the same
# order as in the pathway-specific .bim files that Phase5/6
# use for Matrix_Allele1/Matrix_Allele2) -- if this is named differently in this
# bigsnpr version, stop here clearly with column names instead of
# to silently compute wrongly.
allele_col_candidates <- list(c("allele1", "allele2"), c("a0", "a1"))
allele_cols <- NULL
for (cand in allele_col_candidates) {
  if (all(cand %in% colnames(map_pool))) { allele_cols <- cand; break }
}
if (is.null(allele_cols)) {
  stop(
    "Allele columns not found in the genome-wide pool (expected allele1/allele2 or a0/a1). ",
    "Vorhandene Spalten: ", paste(colnames(map_pool), collapse = ", ")
  )
}
map_pool$Matrix_Allele1 <- toupper(as.character(map_pool[[allele_cols[1]]]))
map_pool$Matrix_Allele2 <- toupper(as.character(map_pool[[allele_cols[2]]]))

id_col_fam <- if ("sample.ID" %in% colnames(fam_pool)) "sample.ID" else "IID"
pool_sample_ids <- as.character(fam_pool[[id_col_fam]])

pool_imputed_meta_file <- sub("\\.rds$", "_imputed_meta.rds", pool_rds_file)

if (file.exists(pool_imputed_meta_file)) {
  cat("Loading cached imputed genome-wide pool:", pool_imputed_meta_file, "\n")
  pool_imputed <- readRDS(pool_imputed_meta_file)
  G_pool <- pool_imputed$genotypes
} else {
  cat(
    "Impute genomweiten Pool einmalig (mean2, ", ncol(pool$genotypes),
    " SNPs) -- may take a few minutes depending on the machine...\n", sep = ""
  )
  G_pool_imputed <- bigsnpr::snp_fastImputeSimple(
    pool$genotypes, method = "mean2", ncores = 1
  )
  pool_imputed <- list(genotypes = G_pool_imputed, map = map_pool, fam = fam_pool)
  saveRDS(pool_imputed, pool_imputed_meta_file)
  G_pool <- G_pool_imputed
  cat("Imputed pool cached at:", pool_imputed_meta_file, "\n")
}

cat(
  "Genomweiter Pool bereit:", nrow(map_pool), "SNPs x", length(pool_sample_ids),
  "Personen.\n"
)

#==================================================================================
# 2) Load external orientation lookup table (Phase5 output) -- identical
#    to Phase6.
#==================================================================================
if (!file.exists(external_lookup_rds)) {
  stop(
    "External orientation lookup not found: ", external_lookup_rds,
    ". Run Phase5_External_Orientation_Prepare.R first."
  )
}
external_lookup <- readRDS(external_lookup_rds)

#==================================================================================
# 3) Helper: reconstruct training/analysis covariate frame per pathway --
#    identical to phase 4's build_train_covariate_frame(), BUT: reads the
#    split file (already corrected by the full-sample fix) from Outputs_Ext,
#    so that train_ids == test_ids == full pathway-specific sample
#    is -- identical to the fitting population of the real external score in
#    Phase6 v1.3.
#==================================================================================
build_analysis_covariate_frame <- function(pathway_name) {

  final_dataset_ld_file <- orig_output_file(
    "Pathway_", pathway_name, "_Combined_QC_LDpruned.rds"
  )
  if (!file.exists(final_dataset_ld_file)) {
    stop("File not found: ", final_dataset_ld_file)
  }
  pathway_data <- readRDS(final_dataset_ld_file)
  pheno_full <- as.data.frame(pathway_data$pheno, stringsAsFactors = FALSE)
  pheno_full$id <- as.character(pheno_full$id)

  split_file <- output_file(
    "Pathway_", pathway_name, "_Global_TrainTest_Split.rds"
  )
  if (!file.exists(split_file)) {
    stop(
      "File not found: ", split_file,
      " -- Phase6 (full-sample fix) must first have run successfully for this pathway ",
      "(copies/corrects the split file to Outputs_Ext)."
    )
  }
  global_split <- readRDS(split_file)
  analysis_ids <- as.character(union(global_split$train_ids, global_split$test_ids))

  df <- pheno_full

  df$gender_num <- as.numeric(as.character(df$gender))
  df$gender_num[df$gender_num %in% c(0, 7777777, 9999999)] <- NA
  df$gender_cat <- factor(df$gender_num, levels = c(1, 2), labels = c("male", "female"))

  df$smoking_habits_3m_num <- as.numeric(as.character(df$smoking_habits_3m))
  df$smoking_habits_3m_num[df$smoking_habits_3m_num %in% c(0, 7777777, 9999999)] <- NA
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

  global_gPCs <- readRDS(genomewide_pcs_file)
  global_gPCs$id <- as.character(global_gPCs$id)
  df <- merge(df, global_gPCs, by = "id", all.x = FALSE)

  outlier_rows_all <- c()
  for (v in cov_numeric) {
    x <- df[[v]]
    Q1 <- quantile(x, 0.25, na.rm = TRUE)
    Q3 <- quantile(x, 0.75, na.rm = TRUE)
    IQRv <- IQR(x, na.rm = TRUE)
    outlier_idx <- which(x < Q1 - 1.5 * IQRv | x > Q3 + 1.5 * IQRv)
    outlier_rows_all <- unique(c(outlier_rows_all, outlier_idx))
  }
  if (length(outlier_rows_all) > 0) df <- df[-outlier_rows_all, , drop = FALSE]

  df_analysis <- df[df$id %in% analysis_ids, c("id", "obesity", cov_categorical, cov_numeric)]
  df_analysis <- df_analysis[complete.cases(df_analysis), ]

  if (nrow(df_analysis) < 10) {
    stop("Too few complete cases for ", pathway_name)
  }

  df_analysis
}

#==================================================================================
# 4) Helper: observed external test statistic from the Phase6 v1.3 output
#    (main_ld_oriented, term SNP_SumScore) -- identical principle to phase 4's
#    get_observed_statistic(), only from Outputs_Ext instead of Outputs_Pathway.
#==================================================================================
get_observed_statistic_external <- function(pathway_name) {
  regression_results_file <- output_file(
    "Pathway_", pathway_name,
    "_Multivariate_SNP_SumScore_Regression_Results_01_Main_LDpruned_All_Oriented.csv"
  )
  if (!file.exists(regression_results_file)) {
    stop("File not found: ", regression_results_file)
  }
  reg <- read.csv2(regression_results_file, stringsAsFactors = FALSE)
  reg_ref <- reg[reg$term == "SNP_SumScore", ]
  if (nrow(reg_ref) != 1) stop("SNP_SumScore row not unique in ", regression_results_file)

  train_or_ref <- as.numeric(reg_ref$OR_Score_Per_SD[1])
  list(train_or_ref = train_or_ref, observed_stat = abs(log(train_or_ref)))
}

#==================================================================================
# 5) Helper: permissible sampling pool (genome-wide pool MINUS whole genes of the
#    pathway) -- taken over unchanged from phase 4.
#==================================================================================
get_eligible_pool_indices <- function(pathway_name) {

  rf <- gene_ranges_file(pathway_name)
  if (!file.exists(rf)) stop("Range file not found: ", rf)

  ranges <- read.table(
    rf, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
    col.names = c("CHR", "START", "END", "GENE")
  )
  ranges$CHR <- as.character(ranges$CHR)

  excluded <- rep(FALSE, nrow(map_pool))

  for (chr in unique(ranges$CHR)) {
    chr_ranges <- ranges[ranges$CHR == chr, ]
    chr_idx <- which(map_pool$chromosome == chr)
    if (length(chr_idx) == 0) next
    pos <- map_pool$physical.pos[chr_idx]

    in_any_gene <- rep(FALSE, length(pos))
    for (g in seq_len(nrow(chr_ranges))) {
      in_any_gene <- in_any_gene |
        (pos >= chr_ranges$START[g] & pos <= chr_ranges$END[g])
    }
    excluded[chr_idx[in_any_gene]] <- TRUE
  }

  which(!excluded)
}

#==================================================================================
# 6) Helper: number of SNPs in the actual (LD-pruned) pathway -- the full
#    draw size BEFORE external filtering (option A, 2026-09-05).
#    Identical to phase 4's get_n_snps_pathway().
#==================================================================================
get_n_snps_pathway <- function(pathway_name) {
  f <- orig_output_file("Pathway_", pathway_name, "_Combined_QC_LDpruned.rds")
  if (!file.exists(f)) stop("File not found: ", f)
  d <- readRDS(f)
  ncol(d$genotypes)
}

#==================================================================================
# 7) Core: ONE external permutation for a pathway.
#
#    Difference from phase 4's run_one_permutation(): the internal, inline
#    computed orientation (big_univLogReg() on the training data -> sign
#    -> "protective" flip) is replaced by a lookup against the external
#    GIANT/Yengo table (identical logic as for the real
#    pathway in Phase6: rsID filter, allele match, palindrome exclusion,
#    inversion flag). Only SNPs with a valid external allele match
#    enter the score -- the usable number varies per draw
#    (no oversampling, no coverage rejection, see script header).
#==================================================================================
run_one_permutation_external <- function(
  pathway_name, df_analysis, eligible_idx, n_snps_draw, row_idx_analysis_in_pool, seed
) {

  set.seed(seed)
  drawn_cols <- sample(eligible_idx, size = n_snps_draw, replace = FALSE)
  drawn_snp_ids <- map_pool$marker.ID[drawn_cols]

  bim_table_drawn <- data.frame(
    SNP = drawn_snp_ids,
    Matrix_Allele1 = map_pool$Matrix_Allele1[drawn_cols],
    Matrix_Allele2 = map_pool$Matrix_Allele2[drawn_cols],
    stringsAsFactors = FALSE
  )

  ext_table <- create_external_orientation_table(external_lookup, bim_table_drawn, drawn_snp_ids)

  usable <- ext_table %>%
    filter(Allele_Match %in% c("Effect_is_matrix_counted_allele", "Effect_is_matrix_other_allele"))

  n_usable <- nrow(usable)
  if (n_usable == 0) return(list(stat = NA_real_, n_usable = 0L))

  usable_cols <- drawn_cols[match(usable$SNP, drawn_snp_ids)]
  invert_idx <- which(usable$Invert_For_External_Score)

  geno <- G_pool[row_idx_analysis_in_pool, usable_cols, drop = FALSE]

  if (anyNA(geno)) {
    col_means <- colMeans(geno, na.rm = TRUE)
    for (j in seq_len(ncol(geno))) {
      na_idx <- is.na(geno[, j])
      if (any(na_idx)) geno[na_idx, j] <- col_means[j]
    }
  }

  if (length(invert_idx) > 0) {
    geno[, invert_idx] <- 2 - geno[, invert_idx, drop = FALSE]
  }

  sum_score <- rowSums(geno)
  score_sd <- sd(sum_score)
  if (!is.finite(score_sd) || score_sd == 0) return(list(stat = NA_real_, n_usable = n_usable))
  sum_score_sd <- (sum_score - mean(sum_score)) / score_sd

  df_perm <- df_analysis
  df_perm$SNP_SumScore_SD <- sum_score_sd

  model_formula <- reformulate(
    termlabels = c("SNP_SumScore_SD", cov_categorical, cov_numeric),
    response = "obesity"
  )

  fit <- tryCatch(
    glm(model_formula, data = df_perm, family = binomial),
    error = function(e) NULL
  )
  if (is.null(fit)) return(list(stat = NA_real_, n_usable = n_usable))

  coefs <- summary(fit)$coefficients
  if (!"SNP_SumScore_SD" %in% rownames(coefs)) return(list(stat = NA_real_, n_usable = n_usable))

  list(stat = abs(coefs["SNP_SumScore_SD", "Estimate"]), n_usable = n_usable)
}

#==================================================================================
# 8) Full run per pathway: B_FINAL permutations, with error handling --
#    Structure identical to phase 4's run_permtest_for_pathway(), additionally
#    mean/sd of the usable SNP number per draw as diagnostic columns.
#==================================================================================
run_permtest_for_pathway_external <- function(pathway_name, B) {

  result_row <- data.frame(
    pathway = pathway_name, Status = NA_character_,
    n_snps_pathway_drawn = NA_integer_, n_eligible_pool = NA_integer_,
    B = NA_integer_, n_valid = NA_integer_,
    mean_n_usable_per_draw = NA_real_, sd_n_usable_per_draw = NA_real_,
    observed_stat = NA_real_, external_or_ref = NA_real_,
    null_mean = NA_real_, null_sd = NA_real_, p_value = NA_real_,
    median_time_per_perm_s = NA_real_, mean_time_per_perm_s = NA_real_,
    total_time_s = NA_real_,
    stringsAsFactors = FALSE
  )

  tryCatch({

    df_analysis <- build_analysis_covariate_frame(pathway_name)
    row_idx_analysis_in_pool <- match(df_analysis$id, pool_sample_ids)
    if (anyNA(row_idx_analysis_in_pool)) {
      stop(
        "Not all analysis IDs found in the genome-wide pool (",
        sum(is.na(row_idx_analysis_in_pool)), " missing)."
      )
    }

    eligible_idx <- get_eligible_pool_indices(pathway_name)
    n_snps_draw <- get_n_snps_pathway(pathway_name)
    observed <- get_observed_statistic_external(pathway_name)

    if (n_snps_draw > length(eligible_idx)) {
      stop(
        "Pathway hat mehr SNPs (", n_snps_draw, ") than the permissible ",
        "sampling pool contains elements (", length(eligible_idx), ")."
      )
    }

    seed_base <- pathway_seed_base(pathway_name)
    null_stats <- numeric(B)
    n_usable_vec <- numeric(B)
    timings <- numeric(B)

    for (b in seq_len(B)) {
      t0 <- Sys.time()
      res_b <- run_one_permutation_external(
        pathway_name, df_analysis, eligible_idx, n_snps_draw,
        row_idx_analysis_in_pool, seed = seed_base + b
      )
      null_stats[b] <- res_b$stat
      n_usable_vec[b] <- res_b$n_usable
      timings[b] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    }

    n_valid <- sum(!is.na(null_stats))
    if (n_valid == 0) stop("All ", B, " permutations failed (NA).")

    p_value <- (1 + sum(null_stats >= observed$observed_stat, na.rm = TRUE)) /
      (1 + n_valid)

    result_row$Status <- "OK"
    result_row$n_snps_pathway_drawn <- n_snps_draw
    result_row$n_eligible_pool <- length(eligible_idx)
    result_row$B <- B
    result_row$n_valid <- n_valid
    result_row$mean_n_usable_per_draw <- mean(n_usable_vec)
    result_row$sd_n_usable_per_draw <- sd(n_usable_vec)
    result_row$observed_stat <- observed$observed_stat
    result_row$external_or_ref <- observed$train_or_ref
    result_row$null_mean <- mean(null_stats, na.rm = TRUE)
    result_row$null_sd <- sd(null_stats, na.rm = TRUE)
    result_row$p_value <- p_value
    result_row$median_time_per_perm_s <- median(timings)
    result_row$mean_time_per_perm_s <- mean(timings)
    result_row$total_time_s <- sum(timings)

  }, error = function(e) {
    result_row$Status <<- paste0("FEHLER: ", conditionMessage(e))
  })

  result_row
}

#==================================================================================
# 9) MD-Log: Appendix C (external permutation) -- Upsert + Session-Info,
#    identical to the pattern of phase 4 (Appendix E) / Phase6 (Appendix B), new
#    markers AUTOLOG-PERM-EXTERNAL(-SESSIONINFO), so that the internal Appendix E
#    (AUTOLOG-PERM-PILOT) remains untouched.
#==================================================================================
log_external_permutation_to_md <- function(md_path, results_table) {

  if (!file.exists(md_path)) {
    warning("Methods/results log not found: ", md_path, " - not written.")
    return(invisible(NULL))
  }

  log_lines <- readLines(md_path, warn = FALSE, encoding = "UTF-8")
  begin_marker <- which(trimws(log_lines) == "<!-- AUTOLOG-PERM-EXTERNAL:BEGIN -->")
  end_marker <- which(trimws(log_lines) == "<!-- AUTOLOG-PERM-EXTERNAL:END -->")

  if (length(begin_marker) != 1 || length(end_marker) != 1 || end_marker <= begin_marker) {
    warning("Marker for Appendix C (external) not found in ", md_path, " - not written.")
    return(invisible(NULL))
  }

  header_end <- begin_marker + 2
  data_start <- header_end + 1
  data_end <- end_marker - 1
  existing_data_rows <- if (data_start <= data_end) log_lines[data_start:data_end] else character(0)

  safe_cell <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x) || x == "") return("")
    gsub("\\|", "/", gsub("[\r\n]+", " ", as.character(x)))
  }
  extract_pathway_from_row <- function(row) {
    cells <- strsplit(row, "\\|")[[1]]
    if (length(cells) < 2) return(NA_character_)
    trimws(cells[2])
  }

  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M")
  new_pathways <- trimws(as.character(results_table$pathway))

  new_rows <- vapply(seq_len(nrow(results_table)), function(i) {
    r <- results_table[i, ]
    if (identical(r$Status, "OK")) {
      sprintf(
        "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
        safe_cell(r$pathway), timestamp,
        safe_cell(r$n_snps_pathway_drawn), safe_cell(r$B),
        safe_cell(round(r$mean_n_usable_per_draw, 1)),
        safe_cell(round(r$observed_stat, 4)), safe_cell(round(r$null_mean, 4)),
        safe_cell(round(r$null_sd, 4)), safe_cell(round(r$p_value, 4)),
        safe_cell(round(r$total_time_s, 1))
      )
    } else {
      sprintf(
        "| %s | %s | | | | | | | %s | |",
        safe_cell(r$pathway), timestamp, safe_cell(r$Status)
      )
    }
  }, character(1))

  if (length(existing_data_rows) > 0) {
    existing_pathways <- vapply(existing_data_rows, extract_pathway_from_row, character(1))
    kept_existing_rows <- existing_data_rows[
      is.na(existing_pathways) | !(existing_pathways %in% new_pathways)
    ]
  } else {
    kept_existing_rows <- character(0)
  }

  final_data_rows <- c(kept_existing_rows, new_rows)
  updated_lines <- c(log_lines[1:header_end], final_data_rows, log_lines[end_marker:length(log_lines)])
  writeLines(updated_lines, md_path, useBytes = TRUE)

  invisible(NULL)
}

log_external_permutation_sessioninfo_to_md <- function(md_path) {

  if (!file.exists(md_path)) return(invisible(NULL))

  log_lines <- readLines(md_path, warn = FALSE, encoding = "UTF-8")
  begin_marker <- which(trimws(log_lines) == "<!-- AUTOLOG-PERM-EXTERNAL-SESSIONINFO:BEGIN -->")
  end_marker <- which(trimws(log_lines) == "<!-- AUTOLOG-PERM-EXTERNAL-SESSIONINFO:END -->")

  if (length(begin_marker) != 1 || length(end_marker) != 1 || end_marker <= begin_marker) {
    return(invisible(NULL))
  }

  pkg_version <- function(pkg) {
    tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "nicht installiert")
  }

  info_lines <- c(
    "",
    paste0("- **Letzter Vollauf:** ", format(Sys.time(), "%Y-%m-%d %H:%M")),
    paste0("- **R-Version:** ", R.version.string),
    paste0("- **Plattform:** ", R.version$platform, " (", Sys.info()[["sysname"]], ")"),
    paste0(
      "- **Paketversionen:** bigsnpr ", pkg_version("bigsnpr"),
      ", bigstatsr ", pkg_version("bigstatsr"),
      ", dplyr ", pkg_version("dplyr")
    ),
    paste0("- **Skriptversion:** Phase7_External_Permutation_Test.R v1.0 (B=", B_FINAL,
           ", ROBUST_ONLY=", ROBUST_ONLY, ")"),
    ""
  )

  updated_lines <- c(log_lines[1:begin_marker], info_lines, log_lines[end_marker:length(log_lines)])
  writeLines(updated_lines, md_path, useBytes = TRUE)
  invisible(NULL)
}

#==================================================================================
# 10) Batch run with checkpointing/resume -- identical to the pattern of phase 4.
#==================================================================================
existing_results <- if (file.exists(checkpoint_file) && !FORCE_RERUN) {
  read.csv2(checkpoint_file, stringsAsFactors = FALSE)
} else {
  NULL
}

already_done <- character(0)
if (!is.null(existing_results) && nrow(existing_results) > 0) {
  already_done <- existing_results$pathway[
    existing_results$Status == "OK" & existing_results$B == B_FINAL
  ]
}

if (length(already_done) > 0) {
  cat(
    length(already_done), "of", length(PATHWAYS_FINAL),
    "pathways already with B=", B_FINAL,
    "successfully completed (checkpoint) -- will be skipped:\n  ",
    paste(already_done, collapse = ", "), "\n", sep = ""
  )
}

pathways_to_run <- setdiff(PATHWAYS_FINAL, already_done)

all_results <- if (!is.null(existing_results)) {
  existing_results[!(existing_results$pathway %in% pathways_to_run), ]
} else {
  data.frame()
}

n_total <- length(PATHWAYS_FINAL)
run_start_time <- Sys.time()

for (i in seq_along(pathways_to_run)) {

  pathway_name <- pathways_to_run[i]

  cat("\n==============================\n")
  cat(
    "Pathway ", which(PATHWAYS_FINAL == pathway_name), "/", n_total, ": ",
    pathway_name, "\n", sep = ""
  )
  cat("==============================\n")

  t_pathway_start <- Sys.time()
  row <- run_permtest_for_pathway_external(pathway_name, B = B_FINAL)
  t_pathway_elapsed <- as.numeric(difftime(Sys.time(), t_pathway_start, units = "secs"))

  all_results <- dplyr::bind_rows(all_results, row)

  # Checkpoint: persist intermediate result immediately (CSV + Appendix C)
  write.csv2(all_results, checkpoint_file, row.names = FALSE)
  log_external_permutation_to_md(methods_log_path, row)

  if (identical(row$Status, "OK")) {
    cat(
      "  OK - p_value:", round(row$p_value, 4),
      "| observed_stat:", round(row$observed_stat, 3),
      "| null_mean:", round(row$null_mean, 3),
      "| mean_n_usable_per_draw:", round(row$mean_n_usable_per_draw, 1),
      "| duration:", round(t_pathway_elapsed, 1), "s\n"
    )
  } else {
    cat("  ", row$Status, "\n")
  }

  elapsed_total_min <- as.numeric(difftime(Sys.time(), run_start_time, units = "mins"))
  n_done_this_run <- i
  avg_min_per_pathway <- elapsed_total_min / n_done_this_run
  n_remaining <- length(pathways_to_run) - i
  cat(
    "  Progress: ", n_done_this_run, "/", length(pathways_to_run),
    " (dieser Lauf) - verstrichen: ", round(elapsed_total_min, 1),
    " min - estimated remaining time: ", round(avg_min_per_pathway * n_remaining, 1),
    " min\n", sep = ""
  )
}

log_external_permutation_sessioninfo_to_md(methods_log_path)

cat("\n====================================================\n")
cat("EXTERNAL PERMUTATION TEST -- FULL RUN COMPLETED\n")
print(all_results[, c("pathway", "Status", "n_snps_pathway_drawn", "mean_n_usable_per_draw",
                       "observed_stat", "null_mean", "null_sd", "p_value")])
cat("====================================================\n")
cat(
  "Result file:", checkpoint_file, "\n",
  "Appendix C (Methods-Log) aktualisiert:", methods_log_path, "\n",
  sep = ""
)
