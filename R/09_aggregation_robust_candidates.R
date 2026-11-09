# ==============================================================================
# 09_aggregation_robust_candidates.R
# Aggregation of the 54 pathway analyses: BH-FDR correction, robust-candidate
# criteria and selection flow (Figure 2; Supplementary Tables S2a, S3).
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================

# ==============================================================================
# Aim: identify pathways whose SNP-sum-score association with obesity is
# robust and independently replicated. Discrimination measures (AUC,
# sensitivity/specificity) are not used as criteria.
#
# Inputs (outputs of 07 and 08; no pipeline re-run needed):
#   - Pathway_<Name>_Multivariate_Master_Table_Results_All_Configurations.csv
#       (p-value, LRT p-value and the univariate parameters Univ_OR/
#        Univ_Pvalue of the main_ld/oriented/Youden reference model)
#   - Pathway_<Name>_Multivariate_SNP_SumScore_Regression_Results_01_Main_
#     LDpruned_All_Oriented.csv
#       (per-SD OR/CI of the main_ld/oriented model, column OR_Score_Per_SD)
#   - Pathway_<Name>_Multivariate_Sensitivity_Analysis_Overview_All_Configurations.csv
#       (OR/CI per sensitivity variant)
#   - AllPathways_TestSet_Independent_Replication.csv
#       (08: independent, covariate-adjusted test-set estimate +
#        same-direction flag versus the training OR)
#   - Pathway_Liste.csv (54 pathway names)
#
# Primary test statistic for the FDR correction: the multivariable SNP-sum-
# score p-value from the master table (main_ld / oriented / Youden; training
# fit), not SKAT.
#
# Pre-specified criteria applied uniformly to all 54 pathways
# ("robust candidate" = all four met):
#   1) FDR: BH-corrected q-value (across all 54 main_ld/oriented/Youden
#      p-values) < 0.05
#   2) Sensitivity consistency: the core sensitivity analyses available for
#      all 54 pathways -- LD-pruning sensitivity (Analysis_Config ==
#      "sensitivity_nold") and beta/ridge/elastic-net-weighted score
#      (Analysis_Config == "sensitivity_betaweighted") -- show no
#      significant effect reversal, i.e. no 95% CI that excludes 1 in the
#      direction opposite to the main model. A sign change of the point
#      estimate with overlapping CI does not count as inconsistency. If the
#      CI of an existing core variant is not evaluable (NA), the criterion is
#      considered not met (conservative). The functional-VEP, LoF-burden and
#      KEGG-only variants are not used as criteria.
#   3) Train/test direction consistency: the independent test-set OR (08)
#      has the same direction as the training OR (Status == "OK" and
#      Test_Same_Direction_As_Train == TRUE).
#   4) Effect size: OR per SD sum score (training main model) >=
#      EFFECT_SIZE_OR_THRESHOLD (1.10) in risk direction or
#      <= 1/EFFECT_SIZE_OR_THRESHOLD (~0.909) in protective direction,
#      symmetric on the log-OR scale.
#   Cohort/batch is not used as a criterion (very unequal cohort sizes
#   Ambu/MD).
#
# Additionally reported, not part of Robust_Candidate (informative):
#   - Test-set significance: independent test-set p-value < 0.05
#     (Criterion_TestSet_Significant, raw value TestSet_Pvalue)
#   - LRT significance: likelihood-ratio test p-value of the main model < 0.05
#     (Criterion_LRT, raw value Multivariate_LRT_Pvalue)
#   - Uni-/multivariable direction consistency (Direction_Univariate /
#     Direction_Uni_Multi_Consistent)
#
# Optional: the tracking table is written to Appendix C of the methods/
# results log (complete replacement on every run), if that file exists.
# ==============================================================================

library(dplyr)

output_dir <- file.path(Sys.getenv("FOCUS_RESULTS_DIR", unset = "path/to/results"), "OutPW")
pathway_liste_file <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Pathway_Selection/Pathway_Liste.csv")
methods_log_path <- file.path(Sys.getenv("FOCUS_RESULTS_DIR", unset = "path/to/results"), "Methods_Results_Log.md")

output_file <- function(...) file.path(output_dir, paste0(...))

#==================================================================================
# METHODS/RESULTS LOG (automatic record, Appendix C)
#==================================================================================
# Writes ONLY into the area between the markers
# "<!-- AUTOLOG-TRACKING:BEGIN -->" / "<!-- AUTOLOG-TRACKING:END -->" or
# "<!-- AUTOLOG-FLOW:BEGIN -->" / "<!-- AUTOLOG-FLOW:END -->" -- these
# sections must be inserted once into the .md file BEFORE the first run.
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
    "Methods/results log (Appendix C) updated: ", md_path,
    " (tracking table: ", ifelse(step1$ok, "OK", "NOT written"),
    "; flow table: ", ifelse(step2$ok, "OK", "NOT written"), ")\n",
    sep = ""
  )

  invisible(NULL)
}

# ==============================================================================
# 1) Load pathway list and test-set replication results (08)
# ==============================================================================

pathways <- read.csv(pathway_liste_file, stringsAsFactors = FALSE)

if (!"pathway_name" %in% colnames(pathways)) {
  stop(
    "Column 'pathway_name' not found in ", pathway_liste_file,
    ". Please check column names."
  )
}

pathway_names <- pathways$pathway_name

replication_file <- output_file("AllPathways_TestSet_Independent_Replication.csv")

if (!file.exists(replication_file)) {
  stop(
    "File not found: ", replication_file,
    " - please first run 08_testset_replication.R for all ",
    "pathways."
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
    # The "OR" field of master_file is the effect of the RAW, unstandardized
    # sum score (model "model"), not the per-SD effect (model
    # "model_standardized") required for Criterion_EffectSize ("OR per SD
    # sum score >= 1.10"). The per-SD OR (with its CI from the standardized
    # model) is taken from the separate regression results file (columns
    # OR_Score_Per_SD / *_CI_lower/upper) and used for or_main/CI;
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
        nrow(main_row), " rows) in ", master_file
      )
    }

    regr_row_idx <- which(regr$term == "SNP_SumScore")
    regr_row <- regr[regr_row_idx, , drop = FALSE]

    if (nrow(regr_row) != 1) {
      stop(
        "SNP_SumScore row (for OR per SD) not unique (",
        nrow(regr_row), " rows) in ", regr_file
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
    #     (see rationale in the script header). ---
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

    # Convert Is_Main_Reference robustly to logical (depending on
    # content, read.csv2 returns TRUE/FALSE either as logical or as character
    # string "TRUE"/"FALSE"; as.logical() covers both cases).
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
    # direction opposite to the main model -- mere
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
      # variant is significantly risk-increasing (entire CI > 1)
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

    # --- 08 replication result for this pathway ---
    repl_row <- replication[replication$Pathway == pathway_name, ]

    if (nrow(repl_row) != 1) {
      testset_status <- "MISSING"
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
    row_template$Status <<- paste0("ERROR: ", conditionMessage(e))
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
    " pathways, the cache files could not be loaded completely:\n",
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

# Effect-size threshold: OR per SD sum score >= EFFECT_SIZE_OR_THRESHOLD in
# risk direction or <= 1/EFFECT_SIZE_OR_THRESHOLD in protective direction
# (symmetric on the log-OR scale). Based on the training main model
# (main_ld/oriented/Youden).
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

# Train/test direction consistency (criterion 3). Test-set significance is
# computed separately below (informative only).
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

# Robust_Candidate: criteria 1-4. Criterion_TestSet_Significant /
# Criterion_LRT are shown in the table (informative) but do not enter
# Robust_Candidate.
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
    return(paste0("Cache files incomplete/defective: ", status))
  }

  if (isTRUE(fdr_ok) && isTRUE(sens_ok) && isTRUE(tt_dir_ok) && isTRUE(es_ok)) {
    extra <- paste0(
      " (informative, not selection-relevant: test-set p = ",
      ifelse(is.na(tt_p), "NA", round(tt_p, 4)),
      "; LRT p = ", ifelse(is.na(lrt_p), "NA", round(lrt_p, 4)), ")"
    )
    return(paste0(
      "Meets all criteria: q = ", round(q, 4), " < 0.05; ",
      "no significant effect reversal in the core sensitivity analyses ",
      "(LD pruning, beta-/ridge-/elastic-net-weighted); ",
      "test-set direction = training direction; OR = ", round(or_main, 4),
      " (effect size >= ", EFFECT_SIZE_OR_THRESHOLD, " or <= ",
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
      "significant effect reversal in core sensitivity analysis/analyses",
      ifelse(sens_txt == "", "", paste0(" (", sens_txt, ")"))
    ))
  }
  if (!isTRUE(tt_dir_ok)) {
    reasons <- c(reasons, paste0(
      "test-set direction differs from training direction or status != OK ",
      "(test-set status: ", tt_status, ")"
    ))
  }
  if (!isTRUE(es_ok)) {
    reasons <- c(reasons, paste0(
      "OR = ", ifelse(is.na(or_main), "NA", round(or_main, 4)),
      " (effect size < ", EFFECT_SIZE_OR_THRESHOLD, " or > ",
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

# stage5_df == final robust candidates (Robust_Candidate). Test-set
# significance and LRT are not used as filter stages.

flow_table <- data.frame(
  Stage = c(
    "Pathways analysed (total)",
    "With complete cached results (master/sensitivity/test-set files readable)",
    paste0("q < ", FDR_THRESHOLD, " (BH-FDR, sum-score p-value)"),
    "+ no significant effect reversal in core sensitivity analyses (LD pruning, beta-/ridge-/elastic-net-weighted)",
    "+ test-set direction = training direction",
    paste0(
      "+ effect size OR >= ", EFFECT_SIZE_OR_THRESHOLD, " or <= ",
      round(1 / EFFECT_SIZE_OR_THRESHOLD, 4), " (final robust candidates)"
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

