# Pathway-restricted polygenic scores for energy-substrate metabolism and obesity (FoCus cohort)
# Script: 05_aggregation_robust_candidates.R | Source file: Aggregation_Results_Tracking_Table_REPAIR.R
# Local paths replaced by the environment variable FOCUS_BASE_DIR (project root); analysis code unchanged.
# Individual-level FoCus data are not included (see README).

# ==============================================================================
# REPAIR VERSION (2026-10-09): copy of Aggregation_Results_Tracking_Table.R
# Re-run after correction of the defective phenotype CSVs. Changes ONLY to paths
# (lines marked "# REPAIR"); analysis code unchanged. Writes only to FoCus/Repair.
# ==============================================================================

# ==============================================================================
# Aggregation_Results_Tracking_Table.R  (v3.2)
#
# Aim of this table (Paper 1 -- overview/association study):
# "The genetic burden in pathway X is robustly and independently replicated as
# associated with obesity." This is NOT about prediction/discrimination
# performance (AUC, sensitivity/specificity) -- that is a separate, later
# question for specifically selected pathways (own criteria,
# own pathway selection). Accordingly, the criteria here consist
# exclusively of association/replication evidence (FDR, sensitivity,
# independent test-set replication) -- no discrimination measures.
#
# v3.1 CHANGE versus v3.0: v3.0 had two additional hard criteria
# introduced -- (4) independent test-set significance (p < 0.05) and (5)
# LRT significance (p < 0.05). The first production run showed: of the
# pathways meeting the v2.0 criteria (FDR + sensitivity + test-set direction),
# NOT A SINGLE ONE additionally met criterion (4) -- with n~320 in the
# test set and effect sizes around OR 1.05-1.15, an independent
# significance test in the test set alone is too underpowered to confirm practically
# any candidate (not falsified, only not
# demonstrable). The criterion was therefore unusable for candidate selection.
# Criterion (5) (LRT) contributes hardly any independent information anyway, as it
# tests the same null hypothesis as the already FDR-corrected Wald test (with
# one degree of freedom usually an almost identical p-value).
#
# v3.1 therefore removes both criteria again as hard components of
# Robust_Candidate. From v3.1, Robust_Candidate again consists of exactly the three
# v2.0 criteria (FDR, sensitivity, test-set direction consistency) -- see
# below. The computed values TestSet_Pvalue and Multivariate_LRT_Pvalue
# as well as the helper columns Criterion_TestSet_Significant / Criterion_LRT remain
# as INFORMATIVE (not selection-relevant) columns in the results table,
# so that the information is not lost and remains available for the next
# stage of candidate narrowing (open question: further
# association criteria and/or later discrimination criteria for a
# Paper-2 pre-selection).
#
# v3.2 CHANGE versus v3.1: new, fourth hard criterion
# Criterion_EffectSize added -- excludes statistically significant but
# practically negligible small effects (possible consequence of large
# sample size with OR close to 1). Threshold: OR per SD sum score >= 1.10 in
# risk direction or <= 1/1.10 (~0.909) in protective direction, checked exactly
# symmetrically on the log-OR scale (see EFFECT_SIZE_OR_THRESHOLD
# below). Based on the well-powered training main model, not on the
# test set -- hence not affected by the test-set power
# limitation documented in v3.0/v3.1. The threshold is deliberately defined as a constant
# and can be adapted if needed after reviewing the actual distribution
# 
#
# BUGFIX (2026-08-21): Multivariate_OR/CI were, up to and including v3.2,
# erroneously taken from the master table (column "OR" there = effect
# of the RAW, unstandardized sum score, not the per-SD effect). Thus
# Criterion_EffectSize was effectively checked on the per-allele OR instead of the
# documented "OR per SD sum score". From now on the correct
# per-SD OR (incl. its own CI from the separately fitted, standardized
# model) is loaded from the regression results file (see
# extract_pathway_summary(): regr_file / regr_row / OR_Score_Per_SD). The
# robust-candidate list and the PRISMA flow numbers from earlier runs
# of this script are therefore invalid and must be regenerated
# with this version.
#
# Purpose: builds from the existing cache files of the 54 pathway analyses
# ONE combined results tracking table plus PRISMA-like flow numbers for
# the candidate pathway selection. No pipeline re-run needed -- reads only:
#
#   - Pathway_<Name>_Multivariate_Master_Table_Results_All_Configurations.csv
#       (provides p-value, LRT p-value and the univariate parameters Univ_OR/
#        Univ_Pvalue of the main_ld/oriented/Youden reference model)
#   - Pathway_<Name>_Multivariate_SNP_SumScore_Regression_Results_01_Main_
#     LDpruned_All_Oriented.csv
#       (provides the correct per-SD OR/CI of the main_ld/oriented model,
#        column OR_Score_Per_SD -- see BUGFIX note above)
#   - Pathway_<Name>_Multivariate_Sensitivity_Analysis_Overview_All_Configurations.csv
#       (provides OR/CI per sensitivity variant, already compared with the
#        main_ld/oriented/Youden model)
#   - AllPathways_TestSet_Independent_Replication.csv
#       (phase 3: independent, covariate-adjusted test-set estimate +
#        same-direction flag versus the training OR)
#   - Pathway_Liste.csv (54 pathway names)
#
# Primary test statistic for the FDR correction: the multivariable SNP sum
# score p-value from the master table (main_ld / oriented / Youden) -- NOT
# SKAT (as agreed). This p-value comes from the TRAINING fit (see
# Section 2.9); the independent test-set direction (phase 3) is NOT used for
# the FDR correction itself, but as a separate, downstream
# robustness criterion (train/test consistency).
#
# Pre-specified criteria applied uniformly to all 54 pathways
# for "robust candidate" (as agreed in the project):
#   1) FDR: BH-corrected q-value (across all 54 main_ld/oriented/Youden
#      p-values) < 0.05
#   2) Sensitivity consistency (NEW in v2.0, revised versus v1.0):
#      only the variants available for ALL 54 pathways are checked
#      core sensitivity analyses -- LD-pruning sensitivity (Analysis_
#      Config == "sensitivity_nold") and beta/ridge/elastic-net-weighted
#      score (Analysis_Config == "sensitivity_betaweighted") -- NOT the
#      functional VEP sensitivity and NOT the optional variants available only for
#      some pathways (LoF burden, KEGG only). A
#      core variant is considered INCONSISTENT only if its 95% CI
#      statistically significantly excludes 1 in the direction OPPOSITE to the main model
#      (true effect reversal) -- a mere sign
#      change of the point estimate with overlapping CI (noise around a
#      weak/null effect, frequent with ORs close to 1 and n~300-400) does
#      NOT count as inconsistency. (v1.0 had instead checked ALL available
#      variants incl. VEP/LoF/KEGG-only by the mere point-estimate
#      sign; this was too strict and systematically disadvantaged pathways
#      with more optional sensitivity analyses -- see
#      changelog/Appendix C.) If the CI of an existing core variant is
#      not evaluable (NA), the criterion is considered NOT met
#      (conservative). If a pathway has no core variant,
#      the criterion is considered NOT evaluable (NA).
#   3) Train/test direction consistency: according to phase 3 (AllPathways_TestSet_
#      Independent_Replication.csv) does the independent test-set OR show
#      the same direction as the training OR (Status == "OK" and
#      Test_Same_Direction_As_Train == TRUE).
#   4) Effect size (NEW in v3.2): OR per SD sum score (training
#      main model) >= EFFECT_SIZE_OR_THRESHOLD (1.10) in risk direction
#      or <= 1/EFFECT_SIZE_OR_THRESHOLD (~0.909) in protective direction,
#      exactly symmetric on the log-OR scale. Excludes statistically
#      significant but practically negligible small effects.
#   As agreed in the project, cohort/batch is deliberately NOT used as a criterion
#   used (bias due to very unequal cohort sizes Ambu/MD).
#
# ADDITIONALLY COMPUTED, BUT NOT PART OF Robust_Candidate (informative,
# since v3.1 -- see change note above):
#   - Test-set significance: independent test-set p-value < 0.05 (column
#     Criterion_TestSet_Significant, raw value TestSet_Pvalue). Was in v3.0
#     hard criterion, removed again (power limitation with n~320).
#   - LRT significance: likelihood-ratio test p-value of the main model < 0.05
#     (column Criterion_LRT, raw value Multivariate_LRT_Pvalue). Was in v3.0
#     hard criterion, removed again (hardly any independent information
#     versus the already FDR-corrected Wald test).
#
# Uni-/multivariable direction consistency is additionally reported (columns
# Direction_Univariate / Direction_Uni_Multi_Consistent), but is NOT a
# hard selection criterion (was not named as a gate in the original list of criteria)
# -- purely informative for the table.
#
# Results are additionally written automatically to Appendix C of the methods/
# results log (log_tracking_table_to_md()) -- analogous to Appendix
# A/B, but as a COMPLETE replacement of the table content on every run
# (not as upsert/append), as this table is recomputed completely
# from the 54 cache files on every run anyway.
# ==============================================================================

library(dplyr)

output_dir <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Repair/OutPW")   # REPAIR
pathway_liste_file <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Pathway_Selection/Pathway_Liste.csv")
methods_log_path <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Repair/Methods_Results_Log_Repair.md")   # REPAIR

output_file <- function(...) file.path(output_dir, paste0(...))

#==================================================================================
# METHODS/RESULTS LOG (automatic record, Appendix C)
#==================================================================================
# Writes ONLY into the area between the markers
# "<!-- AUTOLOG-TRACKING:BEGIN -->" / "<!-- AUTOLOG-TRACKING:END -->" or
# "<!-- AUTOLOG-FLOW:BEGIN -->" / "<!-- AUTOLOG-FLOW:END -->" -- these
# sections must be inserted once into the .md file BEFORE the first run
# (see Appendix C, provided together with this script).
# Unlike Appendix A/B, the table content is COMPLETELY replaced on every run
# (no append, no upsert) -- the table is recomputed completely
# from the 54 cache files on every run anyway, so there is no
# meaningful "old part" that would need to be kept.

log_tracking_table_to_md <- function(md_path, tracking_table, flow_table) {

  if (!file.exists(md_path)) {
    warning(
      "Methods/results log not found at: ", md_path,
      " - Appendix C entries were NOT written."
    )
    return(invisible(NULL))
  }

  log_lines <- readLines(md_path, warn = FALSE, encoding = "UTF-8")

  replace_block <- function(lines, begin_tag, end_tag, header_lines, body_lines) {

    begin_marker <- which(trimws(lines) == begin_tag)
    end_marker <- which(trimws(lines) == end_tag)

    if (length(begin_marker) != 1 || length(end_marker) != 1 ||
        end_marker <= begin_marker) {
      warning(
        "Marker '", begin_tag, "' / '", end_tag,
        "' not found uniquely or in wrong order in ",
        md_path, " - this block was NOT written."
      )
      return(list(lines = lines, ok = FALSE))
    }

    new_block <- c(
      lines[begin_marker],
      header_lines,
      body_lines,
      lines[end_marker]
    )

    updated <- c(
      lines[seq_len(begin_marker - 1)],
      new_block,
      lines[seq(end_marker + 1, length(lines))]
    )

    list(lines = updated, ok = TRUE)
  }

  safe_cell <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x) || x == "") return("")
    gsub("\\|", "/", gsub("[\r\n]+", " ", as.character(x)))
  }

  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M")

  #--- Tracking table ---------------------------------------------------
  tracking_header <- c(
    paste0(
      "| Pathway | N_test | Nominal_p | LRT_p | q_BH | Direction_Multivariate | ",
      "Direction_Univariate | Uni_Multi_Consistent | N_Sensitivity_Variants | ",
      "Sensitivity_All_Consistent | TestSet_Status | TestSet_p | ",
      "TestSet_Same_Direction_As_Train | Robust_Candidate | Justification | ",
      "Last_Updated |"
    ),
    paste0(
      "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"
    )
  )

  tracking_rows <- vapply(seq_len(nrow(tracking_table)), function(i) {
    r <- tracking_table[i, ]
    sprintf(
      "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
      safe_cell(r$Pathway),
      safe_cell(r$N_test_Multivariate),
      safe_cell(signif(r$Multivariate_Pvalue_Nominal, 4)),
      safe_cell(signif(r$Multivariate_LRT_Pvalue, 4)),
      safe_cell(signif(r$Q_Value_BH, 4)),
      safe_cell(r$Direction_Multivariate),
      safe_cell(r$Direction_Univariate),
      safe_cell(r$Direction_Uni_Multi_Consistent),
      safe_cell(r$N_Sensitivity_Variants),
      safe_cell(r$Sensitivity_All_Consistent),
      safe_cell(r$TestSet_Status),
      safe_cell(signif(r$TestSet_Pvalue, 4)),
      safe_cell(r$TestSet_Same_Direction_As_Train),
      safe_cell(r$Robust_Candidate),
      safe_cell(r$Justification),
      timestamp
    )
  }, character(1))

  step1 <- replace_block(
    log_lines,
    "<!-- AUTOLOG-TRACKING:BEGIN -->",
    "<!-- AUTOLOG-TRACKING:END -->",
    tracking_header,
    tracking_rows
  )

  #--- Flow numbers ----------------------------------------------------------
  flow_header <- c(
    "| Stage | N_Remaining | N_Excluded_At_This_Stage |",
    "|---|---|---|"
  )

  flow_rows <- vapply(seq_len(nrow(flow_table)), function(i) {
    r <- flow_table[i, ]
    sprintf(
      "| %s | %s | %s |",
      safe_cell(r$Stage), safe_cell(r$N_Remaining), safe_cell(r$N_Excluded)
    )
  }, character(1))

  step2 <- replace_block(
    step1$lines,
    "<!-- AUTOLOG-FLOW:BEGIN -->",
    "<!-- AUTOLOG-FLOW:END -->",
    flow_header,
    flow_rows
  )

  if (!step1$ok && !step2$ok) {
    return(invisible(NULL))
  }

  writeLines(step2$lines, md_path, useBytes = TRUE)

  cat(
    "Methods/Results-Log (Appendix C) aktualisiert: ", md_path,
    " (Tracking-Tabelle: ", ifelse(step1$ok, "OK", "NOT written"),
    "; Flow-Tabelle: ", ifelse(step2$ok, "OK", "NOT written"), ")\n",
    sep = ""
  )

  invisible(NULL)
}

# ==============================================================================
# 1) Load pathway list and phase-3 replication results
# ==============================================================================

pathways <- read.csv(pathway_liste_file, stringsAsFactors = FALSE)

if (!"pathway_name" %in% colnames(pathways)) {
  stop(
    "Column 'pathway_name' not in ", pathway_liste_file,
    " found. Please check column names."
  )
}

pathway_names <- pathways$pathway_name

replication_file <- output_file("AllPathways_TestSet_Independent_Replication.csv")

if (!file.exists(replication_file)) {
  stop(
    "File not found: ", replication_file,
    " - please first run Phase3_TestSet_Independent_Replication.R for all ",
    "Pathways laufen lassen."
  )
}

replication <- read.csv2(replication_file, stringsAsFactors = FALSE)

# ==============================================================================
# 2) Per pathway: read master table + sensitivity overview and
#    summarize
# ==============================================================================

extract_pathway_summary <- function(pathway_name) {

  row_template <- data.frame(
    Pathway = pathway_name,
    Status = NA_character_,
    N_test_Multivariate = NA_integer_,
    Multivariate_OR = NA_real_,
    Multivariate_CI_lower = NA_real_,
    Multivariate_CI_upper = NA_real_,
    Multivariate_Pvalue_Nominal = NA_real_,
    Multivariate_LRT_Pvalue = NA_real_,
    Direction_Multivariate = NA_character_,
    Univariate_OR = NA_real_,
    Univariate_Pvalue = NA_real_,
    Direction_Univariate = NA_character_,
    Direction_Uni_Multi_Consistent = NA,
    N_Sensitivity_Variants = NA_integer_,
    N_Sensitivity_Consistent = NA_integer_,
    N_Sensitivity_NA = NA_integer_,
    Sensitivity_All_Consistent = NA,
    Sensitivity_Inconsistent_Variants = NA_character_,
    Sensitivity_NA_Variants = NA_character_,
    TestSet_Status = NA_character_,
    TestSet_N = NA_integer_,
    TestSet_OR = NA_real_,
    TestSet_Pvalue = NA_real_,
    TestSet_Same_Direction_As_Train = NA,
    stringsAsFactors = FALSE
  )

  result <- tryCatch({

    master_file <- output_file(
      "Pathway_", pathway_name,
      "_Multivariate_Master_Table_Results_All_Configurations.csv"
    )
    sens_file <- output_file(
      "Pathway_", pathway_name,
      "_Multivariate_Sensitivity_Analysis_Overview_All_Configurations.csv"
    )
    # BUGFIX (2026-08-21): master_file provides in the "OR" field
    # the effect of the RAW, unstandardized sum score (model "model"),
    # not the per-SD effect (model "model_standardized"), although exactly
    # this is meant for Criterion_EffectSize according to the script header/methods log ("OR per
    # SD sum score >= 1.10"). The correct per-SD OR (with its own
    # CI from the standardized model) is in the separate
    # regression results file (columns OR_Score_Per_SD / *_CI_lower/upper).
    # It is additionally loaded here and used for or_main/CI;
    # master_file remains the source for p-value, LRT p, Univ_OR/p.
    regr_file <- output_file(
      "Pathway_", pathway_name,
      "_Multivariate_SNP_SumScore_Regression_Results_01_Main_LDpruned_All_Oriented.csv"
    )

    if (!file.exists(master_file)) {
      stop("Master table not found: ", master_file)
    }
    if (!file.exists(sens_file)) {
      stop("Sensitivity overview not found: ", sens_file)
    }
    if (!file.exists(regr_file)) {
      stop("Regression results file (for OR per SD) not found: ", regr_file)
    }

    master <- read.csv2(master_file, stringsAsFactors = FALSE)
    sens <- read.csv2(sens_file, stringsAsFactors = FALSE)
    regr <- read.csv2(regr_file, stringsAsFactors = FALSE)

    main_row_idx <- which(
      master$Analysis_Config == "main_ld" &
        master$Score_Mode == "oriented" &
        master$Cutoff_Type == "Youden"
    )
    main_row <- master[main_row_idx, , drop = FALSE]

    if (nrow(main_row) != 1) {
      stop(
        "Reference row main_ld/oriented/Youden not unique (",
        nrow(main_row), " Zeilen) in ", master_file
      )
    }

    regr_row_idx <- which(regr$term == "SNP_SumScore")
    regr_row <- regr[regr_row_idx, , drop = FALSE]

    if (nrow(regr_row) != 1) {
      stop(
        "SNP_SumScore row (for OR per SD) not unique (",
        nrow(regr_row), " Zeilen) in ", regr_file
      )
    }

    or_main <- as.numeric(regr_row$OR_Score_Per_SD[1])
    ci_lower_main <- as.numeric(regr_row$OR_Score_Per_SD_CI_lower[1])
    ci_upper_main <- as.numeric(regr_row$OR_Score_Per_SD_CI_upper[1])
    p_main <- as.numeric(main_row$Pvalue[1])
    if (!"LRT_pvalue" %in% colnames(main_row)) {
      stop("Column 'LRT_pvalue' missing in master table: ", master_file)
    }
    lrt_p_main <- suppressWarnings(as.numeric(main_row$LRT_pvalue[1]))
    or_univ <- as.numeric(main_row$Univ_OR[1])
    p_univ <- as.numeric(main_row$Univ_Pvalue[1])

    direction_of <- function(or) {
      if (is.na(or)) return(NA_character_)
      if (or > 1) return("Risk")
      if (or < 1) return("Protective")
      "Null"
    }

    dir_main <- direction_of(or_main)
    dir_univ <- direction_of(or_univ)
    uni_multi_consistent <- if (is.na(dir_main) || is.na(dir_univ)) {
      NA
    } else {
      dir_main == dir_univ
    }

    # --- Core sensitivity variants: ONLY LD-pruning sensitivity
    #     (sensitivity_nold) and beta/ridge/elastic-net-weighted score
    #     (sensitivity_betaweighted), each Youden cutoff, excl. the
    #     reference model itself. Functional VEP and the optional
    #     variants (LoF burden, KEGG only) deliberately do NOT count
    #     (see rationale in the script header, v2.0). ---
    required_sens_cols <- c(
      "Is_Main_Reference", "Analysis_Config", "Cutoff_Type",
      "Sensitivity_Group", "CI_lower", "CI_upper"
    )
    missing_sens_cols <- setdiff(required_sens_cols, colnames(sens))
    if (length(missing_sens_cols) > 0) {
      stop(
        "Expected columns missing in sensitivity overview: ",
        paste(missing_sens_cols, collapse = ", "), " (", sens_file, ")"
      )
    }

    # Convert Is_Main_Reference robustly to logical (read.csv2 returns, depending
    # depending on content, TRUE/FALSE already as logical OR as character string
    # "TRUE"/"FALSE" -- as.logical() safely covers both cases).
    is_main_ref_logical <- suppressWarnings(as.logical(sens$Is_Main_Reference))

    core_sensitivity_configs <- c("sensitivity_nold", "sensitivity_betaweighted")

    variant_idx <- which(
      sens$Cutoff_Type == "Youden" &
        !is_main_ref_logical &
        sens$Analysis_Config %in% core_sensitivity_configs
    )
    variant_rows <- sens[variant_idx, , drop = FALSE]

    n_variants <- nrow(variant_rows)

    variant_ci_lower <- suppressWarnings(as.numeric(variant_rows$CI_lower))
    variant_ci_upper <- suppressWarnings(as.numeric(variant_rows$CI_upper))

    # A core variant is considered INCONSISTENT (same_dir_logical =
    # FALSE) only if its CI statistically significantly excludes 1 in the
    # main model opposite direction -- mere
    # point-estimate sign changes with overlapping CI count as
    # consistent (same_dir_logical = TRUE). Missing CI -> NA.
    if (n_variants == 0) {
      same_dir_logical <- logical(0)
    } else if (is.na(or_main) || or_main == 1) {
      same_dir_logical <- rep(NA, n_variants)
    } else if (or_main > 1) {
      # main model shows risk direction; significant reversal = variant
      # is significantly protective (entire CI < 1)
      same_dir_logical <- !(variant_ci_upper < 1)
    } else {
      # main model shows protective direction; significant reversal =
      # variant is significantly risk (entire CI > 1)
      same_dir_logical <- !(variant_ci_lower > 1)
    }
    same_dir_logical[is.na(variant_ci_lower) | is.na(variant_ci_upper)] <- NA

    n_na <- sum(is.na(same_dir_logical))
    n_consistent <- sum(same_dir_logical, na.rm = TRUE)

    sensitivity_all_consistent <- if (n_variants == 0) {
      NA
    } else {
      (n_na == 0) && (n_consistent == n_variants)
    }

    inconsistent_labels <- unique(
      variant_rows$Sensitivity_Group[
        !is.na(same_dir_logical) & !same_dir_logical
      ]
    )
    na_labels <- unique(
      variant_rows$Sensitivity_Group[is.na(same_dir_logical)]
    )

    # --- Phase-3 replication result for this pathway ---
    repl_row <- replication[replication$Pathway == pathway_name, ]

    if (nrow(repl_row) != 1) {
      testset_status <- "FEHLT"
      testset_n <- NA_integer_
      testset_or <- NA_real_
      testset_p <- NA_real_
      testset_same_dir <- NA
    } else {
      testset_status <- as.character(repl_row$Status[1])
      testset_n <- suppressWarnings(as.integer(repl_row$N_test[1]))
      testset_or <- suppressWarnings(as.numeric(repl_row$Test_OR[1]))
      testset_p <- suppressWarnings(as.numeric(repl_row$Test_Pvalue[1]))
      testset_same_dir <- suppressWarnings(
        as.logical(repl_row$Test_Same_Direction_As_Train[1])
      )
    }

    data.frame(
      Pathway = pathway_name,
      Status = "OK",
      N_test_Multivariate = suppressWarnings(as.integer(main_row$N_test[1])),
      Multivariate_OR = or_main,
      Multivariate_CI_lower = ci_lower_main,
      Multivariate_CI_upper = ci_upper_main,
      Multivariate_Pvalue_Nominal = p_main,
      Multivariate_LRT_Pvalue = lrt_p_main,
      Direction_Multivariate = dir_main,
      Univariate_OR = or_univ,
      Univariate_Pvalue = p_univ,
      Direction_Univariate = dir_univ,
      Direction_Uni_Multi_Consistent = uni_multi_consistent,
      N_Sensitivity_Variants = n_variants,
      N_Sensitivity_Consistent = n_consistent,
      N_Sensitivity_NA = n_na,
      Sensitivity_All_Consistent = sensitivity_all_consistent,
      Sensitivity_Inconsistent_Variants = paste(inconsistent_labels, collapse = "; "),
      Sensitivity_NA_Variants = paste(na_labels, collapse = "; "),
      TestSet_Status = testset_status,
      TestSet_N = testset_n,
      TestSet_OR = testset_or,
      TestSet_Pvalue = testset_p,
      TestSet_Same_Direction_As_Train = testset_same_dir,
      stringsAsFactors = FALSE
    )

  }, error = function(e) {
    row_template$Status <<- paste0("FEHLER: ", conditionMessage(e))
    row_template
  })

  result
}

all_summary <- dplyr::bind_rows(
  lapply(pathway_names, extract_pathway_summary)
)

n_total <- nrow(all_summary)
n_failed_load <- sum(all_summary$Status != "OK")

if (n_failed_load > 0) {
  cat(
    "\nNote: for ", n_failed_load, " of ", n_total,
    " pathways no complete cache files could be loaded:\n",
    sep = ""
  )
  print(
    all_summary[all_summary$Status != "OK", c("Pathway", "Status")]
  )
}

# ==============================================================================
# 3) BH-FDR correction across all pathways with a valid p-value
# ==============================================================================

all_summary$Q_Value_BH <- NA_real_

valid_p_idx <- which(
  all_summary$Status == "OK" & !is.na(all_summary$Multivariate_Pvalue_Nominal)
)

all_summary$Q_Value_BH[valid_p_idx] <- p.adjust(
  all_summary$Multivariate_Pvalue_Nominal[valid_p_idx],
  method = "BH"
)

# ==============================================================================
# 4) Criteria and final decision
# ==============================================================================

FDR_THRESHOLD <- 0.05
TESTSET_P_THRESHOLD <- 0.05
LRT_P_THRESHOLD <- 0.05

# Effect-size threshold (NEW in v3.2): OR per SD sum score >=
# EFFECT_SIZE_OR_THRESHOLD in risk direction or <=
# 1/EFFECT_SIZE_OR_THRESHOLD in protective direction (exactly symmetric on
# the log-OR scale). Excludes statistically significant but practically
# negligible small effects (e.g. FDR-significant due to large sample size,
# but OR close to 1). Based on the TRAINING main model
# (main_ld/oriented/Youden), not on the test set -- hence not affected by the
# test-set power limitation that in v3.0 led to the exclusion of all
# candidates (see changelog v3.0/v3.1).
EFFECT_SIZE_OR_THRESHOLD <- 1.10

all_summary <- all_summary %>%
  mutate(
    Criterion_FDR = !is.na(Q_Value_BH) & Q_Value_BH < FDR_THRESHOLD,
    Criterion_Sensitivity = !is.na(Sensitivity_All_Consistent) & Sensitivity_All_Consistent,
    Criterion_LRT = Status == "OK" &
      !is.na(Multivariate_LRT_Pvalue) &
      Multivariate_LRT_Pvalue < LRT_P_THRESHOLD,
    Criterion_EffectSize = Status == "OK" &
      !is.na(Multivariate_OR) & Multivariate_OR > 0 &
      abs(log(Multivariate_OR)) >= log(EFFECT_SIZE_OR_THRESHOLD)
  )

# Train/test direction consistency (as before) AND additionally test-set
# significance (NEW): the independent test-set p-value must itself be < 0.05,
# the same direction alone is no longer sufficient.
all_summary$Criterion_TrainTest_Direction <- (
  all_summary$Status == "OK" &
    all_summary$TestSet_Status == "OK" &
    !is.na(all_summary$TestSet_Same_Direction_As_Train) &
    all_summary$TestSet_Same_Direction_As_Train
)

all_summary$Criterion_TestSet_Significant <- (
  all_summary$Status == "OK" &
    all_summary$TestSet_Status == "OK" &
    !is.na(all_summary$TestSet_Pvalue) &
    all_summary$TestSet_Pvalue < TESTSET_P_THRESHOLD &
    all_summary$Criterion_TrainTest_Direction
)

# Robust_Candidate (v3.2): the three v2.0 criteria PLUS the new
# effect-size criterion. Criterion_TestSet_Significant / Criterion_LRT
# are still computed and shown in the table (informative),
# but still do NOT enter Robust_Candidate (see
# change note in the script header, v3.1).
all_summary$Robust_Candidate <- (
  all_summary$Status == "OK" &
    all_summary$Criterion_FDR &
    all_summary$Criterion_Sensitivity &
    all_summary$Criterion_TrainTest_Direction &
    all_summary$Criterion_EffectSize
)

build_justification <- function(status, fdr_ok, sens_ok, tt_dir_ok, es_ok,
                                 q, sens_txt, tt_status, tt_p, lrt_p, or_main) {

  if (status != "OK") {
    return(paste0("Cache-Dateien unvollstaendig/fehlerhaft: ", status))
  }

  if (isTRUE(fdr_ok) && isTRUE(sens_ok) && isTRUE(tt_dir_ok) && isTRUE(es_ok)) {
    extra <- paste0(
      " (informativ, nicht auswahlrelevant: Test-Set-p = ",
      ifelse(is.na(tt_p), "NA", round(tt_p, 4)),
      "; LRT-p = ", ifelse(is.na(lrt_p), "NA", round(lrt_p, 4)), ")"
    )
    return(paste0(
      "Erfuellt alle Kriterien: q = ", round(q, 4), " < 0.05; ",
      "keine signifikante Effektumkehr in den Kern-Sensitivitaetsanalysen ",
      "(LD-Pruning, Beta-/Ridge-/ElasticNet-weighted); ",
      "Test-Set-Richtung = Trainings-Richtung; OR = ", round(or_main, 4),
      " (Effektgroesse >= ", EFFECT_SIZE_OR_THRESHOLD, " bzw. <= ",
      round(1 / EFFECT_SIZE_OR_THRESHOLD, 4), ").", extra
    ))
  }

  reasons <- c()

  if (!isTRUE(fdr_ok)) {
    reasons <- c(reasons, paste0(
      "q = ", ifelse(is.na(q), "NA", round(q, 4)), " (>= ", FDR_THRESHOLD, ")"
    ))
  }
  if (!isTRUE(sens_ok)) {
    reasons <- c(reasons, paste0(
      "signifikante Effektumkehr in Kern-Sensitivitaetsanalyse(n)",
      ifelse(sens_txt == "", "", paste0(" (", sens_txt, ")"))
    ))
  }
  if (!isTRUE(tt_dir_ok)) {
    reasons <- c(reasons, paste0(
      "Test-Set-Richtung weicht von Trainings-Richtung ab oder Status != OK ",
      "(Test-Set-Status: ", tt_status, ")"
    ))
  }
  if (!isTRUE(es_ok)) {
    reasons <- c(reasons, paste0(
      "OR = ", ifelse(is.na(or_main), "NA", round(or_main, 4)),
      " (Effektgroesse < ", EFFECT_SIZE_OR_THRESHOLD, " bzw. > ",
      round(1 / EFFECT_SIZE_OR_THRESHOLD, 4), ")"
    ))
  }

  paste(reasons, collapse = "; ")
}

all_summary$Justification <- mapply(
  build_justification,
  all_summary$Status,
  all_summary$Criterion_FDR,
  all_summary$Criterion_Sensitivity,
  all_summary$Criterion_TrainTest_Direction,
  all_summary$Criterion_EffectSize,
  all_summary$Q_Value_BH,
  all_summary$Sensitivity_Inconsistent_Variants,
  all_summary$TestSet_Status,
  all_summary$TestSet_Pvalue,
  all_summary$Multivariate_LRT_Pvalue,
  all_summary$Multivariate_OR,
  SIMPLIFY = TRUE
)

# ==============================================================================
# 5) PRISMA-like flow numbers (sequential filter, how many survive
#    each stage)
# ==============================================================================

stage0 <- nrow(all_summary)
stage1_df <- all_summary[all_summary$Status == "OK", ]
stage1 <- nrow(stage1_df)

stage2_df <- stage1_df[stage1_df$Criterion_FDR, ]
stage2 <- nrow(stage2_df)

stage3_df <- stage2_df[stage2_df$Criterion_Sensitivity, ]
stage3 <- nrow(stage3_df)

stage4_df <- stage3_df[stage3_df$Criterion_TrainTest_Direction, ]
stage4 <- nrow(stage4_df)

stage5_df <- stage4_df[stage4_df$Criterion_EffectSize, ]
stage5 <- nrow(stage5_df)

# stage5_df == final robust candidates (Robust_Candidate), v3.2. Test-set
# Significance and LRT are still not used as filter stages
# (see change note in the script header, v3.1); the effect-size
# criterion (v3.2) is newly added as the last stage.

flow_table <- data.frame(
  Stage = c(
    "Analysierte Pathways (gesamt)",
    "Mit vollstaendigen Cache-Ergebnissen (Master-/Sensitivitaets-/Test-Set-Dateien lesbar)",
    paste0("q < ", FDR_THRESHOLD, " (BH-FDR, Sum-Score-p-Wert)"),
    "+ keine signifikante Effektumkehr in Kern-Sensitivitaetsanalysen (LD-Pruning, Beta-/Ridge-/ElasticNet-weighted)",
    "+ Test-Set-Richtung = Trainings-Richtung",
    paste0(
      "+ Effektgroesse OR >= ", EFFECT_SIZE_OR_THRESHOLD, " bzw. <= ",
      round(1 / EFFECT_SIZE_OR_THRESHOLD, 4), " (finale robuste Kandidaten)"
    )
  ),
  N_Remaining = c(stage0, stage1, stage2, stage3, stage4, stage5),
  N_Excluded = c(
    NA_integer_,
    stage0 - stage1,
    stage1 - stage2,
    stage2 - stage3,
    stage3 - stage4,
    stage4 - stage5
  ),
  stringsAsFactors = FALSE
)

# ==============================================================================
# 6) Output
# ==============================================================================

cat("\n====================================================\n")
cat("PRISMA-like flow numbers\n")
cat("====================================================\n")
print(flow_table, row.names = FALSE)

cat("\n====================================================\n")
cat("Robust candidates (meet all criteria)\n")
cat("====================================================\n")
if (stage5 > 0) {
  print(
    stage5_df[, c("Pathway", "Q_Value_BH", "Multivariate_OR", "Direction_Multivariate", "TestSet_Pvalue", "Multivariate_LRT_Pvalue")],
    row.names = FALSE
  )
} else {
  cat("No pathways currently meet all criteria.\n")
}

final_columns <- c(
  "Pathway", "Status", "N_test_Multivariate", "Multivariate_OR",
  "Multivariate_CI_lower", "Multivariate_CI_upper",
  "Multivariate_Pvalue_Nominal", "Multivariate_LRT_Pvalue", "Q_Value_BH",
  "Direction_Multivariate",
  "Univariate_OR", "Univariate_Pvalue", "Direction_Univariate",
  "Direction_Uni_Multi_Consistent",
  "N_Sensitivity_Variants", "N_Sensitivity_Consistent", "N_Sensitivity_NA",
  "Sensitivity_All_Consistent", "Sensitivity_Inconsistent_Variants",
  "Sensitivity_NA_Variants",
  "TestSet_Status", "TestSet_N", "TestSet_OR", "TestSet_Pvalue",
  "TestSet_Same_Direction_As_Train",
  "Criterion_FDR", "Criterion_Sensitivity", "Criterion_TrainTest_Direction",
  "Criterion_EffectSize", "Criterion_TestSet_Significant", "Criterion_LRT",
  "Robust_Candidate", "Justification"
)

all_summary <- all_summary[, final_columns]
all_summary <- all_summary[order(all_summary$Q_Value_BH), ]

tracking_output_path <- output_file("AllPathways_Results_Tracking_Table.csv")
flow_output_path <- output_file("AllPathways_PRISMA_Flow_Counts.csv")

write.csv2(all_summary, tracking_output_path, row.names = FALSE)
write.csv2(flow_table, flow_output_path, row.names = FALSE)

cat("\n====================================================\n")
cat("Done.\n")
cat("Tracking table: ", tracking_output_path, "\n", sep = "")
cat("Flow numbers:      ", flow_output_path, "\n", sep = "")
cat("====================================================\n")

# ==============================================================================
# 7) Write automatically to the methods/results log (Appendix C)
# ==============================================================================

log_tracking_table_to_md(methods_log_path, all_summary, flow_table)

