# ==============================================================================
# 08_testset_replication.R
# Independent test-set replication of the multivariable main model: test-set
# OR per SD, 95% CI, p-value and delta R^2 for all 54 pathways (Table 2).
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================

# ==============================================================================
# For each pathway, the multivariable, covariate-adjusted SNP-sum-score main
# model (main_ld / oriented / Youden) is refitted on the independent test set.
# Risk/protective allele categories, imputation values and the
# standardization parameters (mean/SD of the score) are taken from the
# training set only.
#
# Reference values for the training set (Train_OR_Reference, 95% CI,
# p-value, Score_Mean_Train/Score_SD_Train) are read from the per-pathway
# Regression_Results.csv (row term == "SNP_SumScore", columns
# "OR_Score_Per_SD" / "OR_Score_Per_SD_CI_lower" / "OR_Score_Per_SD_CI_upper"
# / "p.value" / "Score_Mean_Train" / "Score_SD_Train"). The p-values of the
# SD-standardized and of the raw score coefficient are identical (a linear
# rescaling of a predictor does not change the Wald statistic).
#
# Delta R^2: incremental McFadden pseudo-R^2 (pscl::pR2()) of the SNP sum
# score versus a covariates-only model, for the training set (reconstructed
# from global_split$train_ids) and for the test set.
#
# No re-run of QC/PLINK/VEP/SKAT is needed; inputs (outputs of 07):
#   - Pathway_<Name>_Combined_QC_LDpruned.rds
#       (genotypes + phenotype incl. raw covariates, complete sample)
#   - Pathway_<Name>_Global_TrainTest_Split.rds
#       (train_ids / test_ids, identical split as in 07)
#   - Pathway_<Name>_Univariate_SNP_Logistic_Regression_Results_
#       01_Main_LDpruned_All.csv
#       (column "Category": risk/protective per SNP, training set only)
#   - Pathway_<Name>_Univariate_SNP_SumScore_Imputation_Values_
#       01_Main_LDpruned_All_Oriented.csv
#       (column "Training_Mean_Genotype": imputation value per SNP)
#   - Pathway_<Name>_Multivariate_SNP_SumScore_Regression_Results_
#       01_Main_LDpruned_All_Oriented.csv
#       (Score_Mean_Train/Score_SD_Train and OR_Score_Per_SD)
#   - Pathway_<Name>_Multivariate_Master_Table_Results_All_Configurations.csv
#       (N_test plausibility check only)
#   - Global_Genomewide_PCs.rds (06)
#
# Individuals with >10% missing genotypes in the pathway SNP set are excluded
# before imputation (as in 07).
# ==============================================================================

library(dplyr)
library(pscl)

output_dir <- file.path(Sys.getenv("FOCUS_RESULTS_DIR", unset = "path/to/results"), "OutPW")
pathway_liste_file <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Pathway_Selection/Pathway_Liste.csv")
genomewide_pcs_file <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Global_PCA/Global_Genomewide_PCs.rds")

output_file <- function(...) file.path(output_dir, paste0(...))

if (!file.exists(genomewide_pcs_file)) {
  stop("Genome-wide PCs not found: ", genomewide_pcs_file)
}

#==================================================================================
# METHODS/RESULTS LOG (automatic record, Appendix B)
#==================================================================================
# Optional: writes the test-set replication results of all pathways of this
# run into Appendix B of the methods/results log (between the markers
# "<!-- AUTOLOG-PHASE3:BEGIN -->" / "<!-- AUTOLOG-PHASE3:END -->", upsert by
# pathway) plus a session-info block ("<!-- AUTOLOG-PHASE3-SESSIONINFO:BEGIN/END -->"),
# if that file exists.
methods_log_path <- file.path(Sys.getenv("FOCUS_RESULTS_DIR", unset = "path/to/results"), "Methods_Results_Log.md")

log_testset_replication_to_md <- function(md_path, results_table) {

  if (!file.exists(md_path)) {
    warning(
      "Methods/results log not found at: ", md_path,
      " - Appendix B entries were NOT written."
    )
    return(invisible(NULL))
  }

  log_lines <- readLines(md_path, warn = FALSE, encoding = "UTF-8")
  begin_marker <- which(trimws(log_lines) == "<!-- AUTOLOG-PHASE3:BEGIN -->")
  end_marker <- which(trimws(log_lines) == "<!-- AUTOLOG-PHASE3:END -->")

  if (length(begin_marker) != 1 || length(end_marker) != 1) {
    warning(
      "Marker '<!-- AUTOLOG-PHASE3:BEGIN/END -->' not found uniquely in ",
      md_path,
      " - Appendix B entries were NOT written."
    )
    return(invisible(NULL))
  }

  if (end_marker <= begin_marker) {
    warning(
      "Marker '<!-- AUTOLOG-PHASE3:BEGIN/END -->' in unexpected ",
      "order in ", md_path, " - Appendix B entries were NOT ",
      "written."
    )
    return(invisible(NULL))
  }

  # Expected directly after BEGIN: header row, then separator row ("|---|..."),
  # then the actual data rows up to (excluding) END.
  header_end <- begin_marker + 2

  if (header_end > length(log_lines) ||
      !grepl("^\\|\\s*Pathway\\s*\\|", log_lines[begin_marker + 1])) {
    warning(
      "Expected table header after '<!-- AUTOLOG-PHASE3:BEGIN -->' ",
      "not found in ", md_path,
      " - Appendix B entries were NOT written."
    )
    return(invisible(NULL))
  }

  if (!grepl("R2_Train_Full", log_lines[begin_marker + 1])) {
    warning(
      "Table header in ", md_path, " does not contain the ",
      "delta R^2 columns (R2_Train_Full etc.). The rows are written ",
      "anyway -- please check the column assignment afterwards."
    )
  }

  data_start <- header_end + 1
  data_end <- end_marker - 1

  existing_data_rows <- if (data_start <= data_end) {
    log_lines[data_start:data_end]
  } else {
    character(0)
  }

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

  new_pathways <- trimws(as.character(results_table$Pathway))

  new_rows <- vapply(seq_len(nrow(results_table)), function(i) {
    r <- results_table[i, ]
    if (identical(r$Status, "OK")) {
      sprintf(
        "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
        safe_cell(r$Pathway), timestamp,
        safe_cell(r$N_test),
        safe_cell(round(r$Test_OR, 4)),
        safe_cell(round(r$Test_CI_lower, 4)),
        safe_cell(round(r$Test_CI_upper, 4)),
        safe_cell(signif(r$Test_Pvalue, 4)),
        safe_cell(round(r$Train_OR_Reference, 4)),
        safe_cell(signif(r$Train_Pvalue_Reference, 4)),
        safe_cell(r$Test_Same_Direction_As_Train),
        safe_cell(round(r$R2_Train_Full, 4)),
        safe_cell(round(r$Delta_R2_Train, 4)),
        safe_cell(round(r$R2_Test_Full, 4)),
        safe_cell(round(r$Delta_R2_Test, 4))
      )
    } else {
      sprintf(
        "| %s | %s | ERROR | | | | | | | %s | | | | |",
        safe_cell(r$Pathway), timestamp, safe_cell(r$Status)
      )
    }
  }, character(1))

  # Remove existing rows for pathways from THIS run (upsert) --
  # rows of other pathways (e.g. from an earlier partial/test run)
  # remain untouched.
  if (length(existing_data_rows) > 0) {
    existing_pathways <- vapply(
      existing_data_rows, extract_pathway_from_row, character(1)
    )
    kept_existing_rows <- existing_data_rows[
      is.na(existing_pathways) | !(existing_pathways %in% new_pathways)
    ]
  } else {
    kept_existing_rows <- character(0)
  }

  n_replaced <- length(existing_data_rows) - length(kept_existing_rows)

  final_data_rows <- c(kept_existing_rows, new_rows)

  updated_lines <- c(
    log_lines[1:header_end],
    final_data_rows,
    log_lines[end_marker:length(log_lines)]
  )

  writeLines(updated_lines, md_path, useBytes = TRUE)

  cat(
    "Methods/results log (Appendix B) updated: ", md_path,
    " (", nrow(results_table), " pathway rows from this run, ",
    n_replaced, " existing row(s) replaced, ",
    length(kept_existing_rows), " unchanged row(s) of other pathways ",
    "kept)\n", sep = ""
  )

  invisible(NULL)
}

#==================================================================================
# Log session/software info block (R version, packages)
#==================================================================================
# Replaces the content between "<!-- AUTOLOG-PHASE3-SESSIONINFO:BEGIN -->" and
# "...:END -->" completely on every run (no upsert needed, as this is
# ONE block for the respective last run, no per-pathway rows).

log_phase3_sessioninfo_to_md <- function(md_path) {

  if (!file.exists(md_path)) {
    warning(
      "Methods/results log not found at: ", md_path,
      " - session info was NOT written."
    )
    return(invisible(NULL))
  }

  log_lines <- readLines(md_path, warn = FALSE, encoding = "UTF-8")
  begin_marker <- which(
    trimws(log_lines) == "<!-- AUTOLOG-PHASE3-SESSIONINFO:BEGIN -->"
  )
  end_marker <- which(
    trimws(log_lines) == "<!-- AUTOLOG-PHASE3-SESSIONINFO:END -->"
  )

  if (length(begin_marker) != 1 || length(end_marker) != 1 ||
      end_marker <= begin_marker) {
    warning(
      "Marker '<!-- AUTOLOG-PHASE3-SESSIONINFO:BEGIN/END -->' not ",
      "found uniquely in ", md_path, " - session info was NOT ",
      "written."
    )
    return(invisible(NULL))
  }

  pkg_version <- function(pkg) {
    tryCatch(
      as.character(utils::packageVersion(pkg)),
      error = function(e) "not installed"
    )
  }

  info_lines <- c(
    "",
    paste0("- **Last run:** ", format(Sys.time(), "%Y-%m-%d %H:%M")),
    paste0("- **R version:** ", R.version.string),
    paste0("- **Platform:** ", R.version$platform, " (", Sys.info()[["sysname"]], ")"),
    paste0(
      "- **Package versions:** dplyr ", pkg_version("dplyr"),
      ", pscl ", pkg_version("pscl"),
      ", broom ", pkg_version("broom")
    ),
    paste0("- **Script:** 08_testset_replication.R"),
    ""
  )

  updated_lines <- c(
    log_lines[1:begin_marker],
    info_lines,
    log_lines[end_marker:length(log_lines)]
  )

  writeLines(updated_lines, md_path, useBytes = TRUE)

  cat("Methods/results log session info updated: ", md_path, "\n", sep = "")

  invisible(NULL)
}

pathways <- read.csv(pathway_liste_file, stringsAsFactors = FALSE)

if (!"pathway_name" %in% colnames(pathways)) {
  stop(
    "Column 'pathway_name' not found in ", pathway_liste_file,
    ". Please check column names."
  )
}

# ------------------------------------------------------------------
# For a first test run, restrict to 1-2 pathways here if needed,
# e.g.:
#   pathways_to_run <- c("Glycolysis_Gluconeogenesis", "Citrate_Cycle")
# default: all pathways from the list.
# ------------------------------------------------------------------
pathways_to_run <- pathways$pathway_name

# Fixed covariate set of the multivariable main model (identical to
# 07_pathway_pipeline.R, manual_covariate_decisions$default$keep; batch deliberately
# excluded -- sensitivity analysis only there)
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

# Small helper: extract McFadden pseudo-R^2 robustly
safe_mcfadden_r2 <- function(model) {
  tryCatch(
    as.numeric(pscl::pR2(model)["McFadden"]),
    error = function(e) NA_real_
  )
}

# ==============================================================================
# Function: test-set replication + delta R^2 for a single pathway
# ==============================================================================

run_testset_replication <- function(pathway_name) {

  result_row <- data.frame(
    Pathway = pathway_name,
    Status = NA_character_,
    N_full_after_outlier_removal = NA_integer_,
    N_train = NA_integer_,
    N_test = NA_integer_,
    Test_OR = NA_real_,
    Test_CI_lower = NA_real_,
    Test_CI_upper = NA_real_,
    Test_Pvalue = NA_real_,
    Train_OR_Reference = NA_real_,
    Train_OR_Reference_CI_lower = NA_real_,
    Train_OR_Reference_CI_upper = NA_real_,
    Train_Pvalue_Reference = NA_real_,
    N_test_Reference_MasterTable = NA_integer_,
    Test_Same_Direction_As_Train = NA,
    R2_Train_Full = NA_real_,
    R2_Train_CovOnly = NA_real_,
    Delta_R2_Train = NA_real_,
    R2_Test_Full = NA_real_,
    R2_Test_CovOnly = NA_real_,
    Delta_R2_Test = NA_real_,
    stringsAsFactors = FALSE
  )

  tryCatch({

    #----------------------------------------------------------------
    # 1) Load cached raw data (genotypes + phenotype INCL. raw
    #    covariates, complete sample, LD-pruned)
    #----------------------------------------------------------------
    final_dataset_ld_file <- output_file(
      "Pathway_", pathway_name, "_Combined_QC_LDpruned.rds"
    )

    if (!file.exists(final_dataset_ld_file)) {
      stop("File not found: ", final_dataset_ld_file)
    }

    pathway_data <- readRDS(final_dataset_ld_file)
    geno_full  <- as.matrix(pathway_data$genotypes)
    storage.mode(geno_full) <- "numeric"
    pheno_full <- as.data.frame(pathway_data$pheno, stringsAsFactors = FALSE)

    required_pheno_cols <- c(
      "id", "obesity", "gender", "smoking_habits_3m",
      "school_certificate_A", "graduation", "sample_age",
      "sports_activity_min_week", "everyday_activity_min_week",
      "energy_kcal"
    )
    missing_pheno_cols <- setdiff(required_pheno_cols, colnames(pheno_full))
    if (length(missing_pheno_cols) > 0) {
      stop(
        "The following expected covariate/phenotype columns are missing: ",
        paste(missing_pheno_cols, collapse = ", ")
      )
    }

    if (!all(rownames(geno_full) == pheno_full$id)) {
      stop("Order of genotype and phenotype dataset does not match.")
    }

    #----------------------------------------------------------------
    # 1b) Exclude individuals with > 10 % missing SNPs in the pathway set
    #     -- identical rule as in 07_pathway_pipeline.R,
    #     applied BEFORE the training-based imputation.
    #----------------------------------------------------------------
    indiv_missing <- rowMeans(is.na(geno_full))
    keep_missingness <- indiv_missing <= 0.10

    n_before_missingness_filter <- nrow(geno_full)

    geno_full  <- geno_full[keep_missingness, , drop = FALSE]
    pheno_full <- pheno_full[keep_missingness, , drop = FALSE]

    n_after_missingness_filter <- nrow(geno_full)

    if (n_after_missingness_filter < n_before_missingness_filter) {
      cat(
        "Note: ", n_before_missingness_filter - n_after_missingness_filter,
        " individual(s) with > 10% missing genotypes in the SNP set of ",
        pathway_name, " excluded (identical to 07_pathway_pipeline.R).\n",
        sep = ""
      )
    }

    if (!identical(
      as.character(rownames(geno_full)),
      as.character(pheno_full$id)
    )) {
      stop(
        "Order of genotype and phenotype dataset after ",
        "missingness filter does not match."
      )
    }

    #----------------------------------------------------------------
    # 2) Load train/test split (identical, canonical split as
    #    in 07_pathway_pipeline.R)
    #----------------------------------------------------------------
    split_file <- output_file(
      "Pathway_", pathway_name, "_Global_TrainTest_Split.rds"
    )

    if (!file.exists(split_file)) {
      stop("File not found: ", split_file)
    }

    global_split <- readRDS(split_file)

    if (is.null(global_split$train_ids) || is.null(global_split$test_ids)) {
      stop("global_split does not contain both elements 'train_ids'/'test_ids'.")
    }

    train_ids <- as.character(global_split$train_ids)
    test_ids <- as.character(global_split$test_ids)

    #----------------------------------------------------------------
    # 3) Load SNP direction (risk/protective) from the univariate result
    #    -- determined exclusively on the training set
    #----------------------------------------------------------------
    univar_file <- output_file(
      "Pathway_", pathway_name,
      "_Univariate_SNP_Logistic_Regression_Results_01_Main_LDpruned_All.csv"
    )

    if (!file.exists(univar_file)) {
      stop("File not found: ", univar_file)
    }

    univar_results <- read.csv2(univar_file, stringsAsFactors = FALSE)

    if (!all(c("SNP", "Category") %in% colnames(univar_results))) {
      stop("Expected columns 'SNP'/'Category' missing in ", univar_file)
    }

    missing_category_snps <- setdiff(colnames(geno_full), univar_results$SNP)
    if (length(missing_category_snps) > 0) {
      stop(
        "For the following SNPs in the genotype dataset there is no category ",
        "assignment: ", paste(missing_category_snps, collapse = ", ")
      )
    }

    protective_snps <- univar_results$SNP[
      tolower(univar_results$Category) == "protective"
    ]
    protective_snps <- intersect(protective_snps, colnames(geno_full))

    #----------------------------------------------------------------
    # 4) Load imputation values (training-based) and apply to missing
    #    genotypes in the ENTIRE sample (train+test) -- exactly
    #    as in 07_pathway_pipeline.R (Imputation_Source = "Training_Set")
    #----------------------------------------------------------------
    imputation_file <- output_file(
      "Pathway_", pathway_name,
      "_Univariate_SNP_SumScore_Imputation_Values_01_Main_LDpruned_All_Oriented.csv"
    )

    if (!file.exists(imputation_file)) {
      stop("File not found: ", imputation_file)
    }

    imputation_table <- read.csv2(imputation_file, stringsAsFactors = FALSE)

    if (!all(c("SNP", "Training_Mean_Genotype") %in% colnames(imputation_table))) {
      stop(
        "Expected columns 'SNP'/'Training_Mean_Genotype' missing in ",
        imputation_file
      )
    }

    missing_impute_snps <- setdiff(colnames(geno_full), imputation_table$SNP)
    if (length(missing_impute_snps) > 0) {
      stop(
        "For the following SNPs there is no imputation value: ",
        paste(missing_impute_snps, collapse = ", ")
      )
    }

    impute_lookup <- setNames(
      as.numeric(imputation_table$Training_Mean_Genotype),
      imputation_table$SNP
    )

    for (snp in colnames(geno_full)) {
      missing_idx <- is.na(geno_full[, snp])
      if (any(missing_idx)) {
        geno_full[missing_idx, snp] <- impute_lookup[[snp]]
      }
    }

    if (anyNA(geno_full)) {
      stop("Missing genotype values remain after imputation.")
    }

    #----------------------------------------------------------------
    # 5) Apply score orientation (recode protective SNPs) and
    #    build the sum score -- identical logic as in 07_pathway_pipeline.R
    #    (score_mode == "oriented"), on the ENTIRE sample
    #----------------------------------------------------------------
    if (length(protective_snps) > 0) {
      protective_idx <- match(protective_snps, colnames(geno_full))
      geno_full[, protective_idx] <- 2 - geno_full[, protective_idx, drop = FALSE]
    }

    snp_sumscore_full <- rowSums(geno_full)

    #----------------------------------------------------------------
    # 6) Standardize with the training-set parameters
    #    from the per-pathway Regression_Results.csv (the OR field of the
    #    master table refers to the unstandardized score). The master table
    #    is only loaded for the N_test check.
    #----------------------------------------------------------------
    master_file <- output_file(
      "Pathway_", pathway_name,
      "_Multivariate_Master_Table_Results_All_Configurations.csv"
    )

    if (!file.exists(master_file)) {
      stop("File not found: ", master_file)
    }

    master_table <- read.csv2(master_file, stringsAsFactors = FALSE)

    required_master_cols <- c(
      "Analysis_Config", "Score_Mode", "Cutoff_Type", "N_test"
    )
    missing_master_cols <- setdiff(required_master_cols, colnames(master_table))
    if (length(missing_master_cols) > 0) {
      stop(
        "The following columns are missing in the master table: ",
        paste(missing_master_cols, collapse = ", ")
      )
    }

    main_ref <- master_table[
      master_table$Analysis_Config == "main_ld" &
        master_table$Score_Mode == "oriented" &
        master_table$Cutoff_Type == "Youden",
    ]

    if (nrow(main_ref) != 1) {
      stop(
        "Reference row main_ld/oriented/Youden not unique (",
        nrow(main_ref), " rows) in ", master_file
      )
    }

    regression_results_file <- output_file(
      "Pathway_", pathway_name,
      "_Multivariate_SNP_SumScore_Regression_Results_01_Main_LDpruned_All_Oriented.csv"
    )

    if (!file.exists(regression_results_file)) {
      stop("File not found: ", regression_results_file)
    }

    regression_results <- read.csv2(regression_results_file, stringsAsFactors = FALSE)

    required_regression_cols <- c(
      "term", "p.value", "Score_Mean_Train", "Score_SD_Train",
      "OR_Score_Per_SD", "OR_Score_Per_SD_CI_lower", "OR_Score_Per_SD_CI_upper"
    )
    missing_regression_cols <- setdiff(
      required_regression_cols, colnames(regression_results)
    )
    if (length(missing_regression_cols) > 0) {
      stop(
        "The following columns are missing in Regression_Results.csv: ",
        paste(missing_regression_cols, collapse = ", ")
      )
    }

    reg_ref <- regression_results[regression_results$term == "SNP_SumScore", ]

    if (nrow(reg_ref) != 1) {
      stop(
        "Row term=='SNP_SumScore' in Regression_Results.csv not ",
        "unique (", nrow(reg_ref), " rows) in ", regression_results_file
      )
    }

    score_mean_train <- as.numeric(reg_ref$Score_Mean_Train[1])
    score_sd_train   <- as.numeric(reg_ref$Score_SD_Train[1])

    if (!is.finite(score_mean_train) || !is.finite(score_sd_train) ||
        score_sd_train <= 0) {
      stop("Score_Mean_Train/Score_SD_Train from Regression_Results.csv not plausible.")
    }

    train_or_ref <- as.numeric(reg_ref$OR_Score_Per_SD[1])
    train_or_ref_ci_lower <- as.numeric(reg_ref$OR_Score_Per_SD_CI_lower[1])
    train_or_ref_ci_upper <- as.numeric(reg_ref$OR_Score_Per_SD_CI_upper[1])
    train_pvalue_ref <- as.numeric(reg_ref$p.value[1])

    snp_sumscore_full_sd <- (snp_sumscore_full - score_mean_train) / score_sd_train

    #----------------------------------------------------------------
    # 7) Covariate recoding -- identical logic as 07_pathway_pipeline.R
    #    
    #----------------------------------------------------------------
    df <- pheno_full
    df$SNP_SumScore_SD <- snp_sumscore_full_sd

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

    #----------------------------------------------------------------
    # 8) Merge genome-wide PCs -- identical to 07_pathway_pipeline.R
    #    
    #----------------------------------------------------------------
    df$id <- as.character(df$id)

    global_gPCs <- readRDS(genomewide_pcs_file)
    global_gPCs$id <- as.character(global_gPCs$id)

    df <- merge(
      df,
      global_gPCs,
      by = "id",
      all.x = FALSE
    )

    if (nrow(df) == 0) {
      stop("No individuals with genome-wide PCs after merge.")
    }

    #----------------------------------------------------------------
    # 9) Outlier removal -- identical 1.5x IQR rule as
    #    07_pathway_pipeline.R, on the same covariate set
    #----------------------------------------------------------------
    outlier_rows_all <- c()

    for (v in cov_numeric) {

      x <- df[[v]]

      Q1 <- quantile(x, 0.25, na.rm = TRUE)
      Q3 <- quantile(x, 0.75, na.rm = TRUE)
      IQRv <- IQR(x, na.rm = TRUE)

      lower <- Q1 - 1.5 * IQRv
      upper <- Q3 + 1.5 * IQRv

      outlier_idx <- which(x < lower | x > upper)

      outlier_rows_all <- unique(c(outlier_rows_all, outlier_idx))
    }

    if (length(outlier_rows_all) > 0) {
      df <- df[-outlier_rows_all, , drop = FALSE]
    }

    result_row$N_full_after_outlier_removal <- nrow(df)

    #----------------------------------------------------------------
    # 10) Restrict to test set and fit the multivariable model (full
    #     and covariates-only) ONLY on the test set -- identical
    #     covariate set as 07_pathway_pipeline.R main_ld/oriented
    #     (manual_covariate_decisions$default$keep, without batch)
    #----------------------------------------------------------------
    model_vars <- c(
      "id", "obesity", "SNP_SumScore_SD", cov_categorical, cov_numeric
    )

    df_test <- df[df$id %in% test_ids, ]
    df_test_model <- df_test[, model_vars]
    df_test_model <- df_test_model[complete.cases(df_test_model), ]

    if (nrow(df_test_model) < 10) {
      stop(
        "Too few complete test-set cases after covariate ",
        "filtering (n = ", nrow(df_test_model), ")."
      )
    }

    if (length(unique(df_test_model$obesity)) < 2) {
      stop("Test set contains only one outcome class -- model not estimable.")
    }

    model_formula_standardized <- reformulate(
      termlabels = c("SNP_SumScore_SD", cov_categorical, cov_numeric),
      response = "obesity"
    )

    formula_covonly <- reformulate(
      termlabels = c(cov_categorical, cov_numeric),
      response = "obesity"
    )

    model_test <- glm(
      model_formula_standardized,
      data = df_test_model,
      family = binomial
    )

    coefs <- summary(model_test)$coefficients

    if (!"SNP_SumScore_SD" %in% rownames(coefs)) {
      stop("Coefficient 'SNP_SumScore_SD' not found in the test model.")
    }

    ci <- tryCatch(
      confint.default(model_test)["SNP_SumScore_SD", ],
      error = function(e) c(NA_real_, NA_real_)
    )

    test_or <- exp(coefs["SNP_SumScore_SD", "Estimate"])

    model_test_covonly <- glm(
      formula_covonly,
      data = df_test_model,
      family = binomial
    )

    #----------------------------------------------------------------
    # 10b) Analogous full and covariates-only model on the
    #      (internally reconstructed) TRAINING SET for delta R^2 (train)
    #----------------------------------------------------------------
    df_train <- df[df$id %in% train_ids, ]
    df_train_model <- df_train[, model_vars]
    df_train_model <- df_train_model[complete.cases(df_train_model), ]

    if (nrow(df_train_model) < 10) {
      stop(
        "Too few complete training-set cases after covariate ",
        "filtering (n = ", nrow(df_train_model), ")."
      )
    }

    if (length(unique(df_train_model$obesity)) < 2) {
      stop("Training set contains only one outcome class -- model not estimable.")
    }

    model_train_full <- glm(
      model_formula_standardized,
      data = df_train_model,
      family = binomial
    )

    model_train_covonly <- glm(
      formula_covonly,
      data = df_train_model,
      family = binomial
    )

    #----------------------------------------------------------------
    # 10c) McFadden pseudo-R^2 and delta R^2 (train + test)
    #----------------------------------------------------------------
    r2_train_full    <- safe_mcfadden_r2(model_train_full)
    r2_train_covonly <- safe_mcfadden_r2(model_train_covonly)
    r2_test_full      <- safe_mcfadden_r2(model_test)
    r2_test_covonly   <- safe_mcfadden_r2(model_test_covonly)

    delta_r2_train <- if (is.finite(r2_train_full) && is.finite(r2_train_covonly)) {
      r2_train_full - r2_train_covonly
    } else {
      NA_real_
    }

    delta_r2_test <- if (is.finite(r2_test_full) && is.finite(r2_test_covonly)) {
      r2_test_full - r2_test_covonly
    } else {
      NA_real_
    }

    #----------------------------------------------------------------
    # 11) Fill result row
    #----------------------------------------------------------------
    result_row$Status <- "OK"
    result_row$N_train <- nrow(df_train_model)
    result_row$N_test <- nrow(df_test_model)
    result_row$Test_OR <- test_or
    result_row$Test_CI_lower <- exp(ci[1])
    result_row$Test_CI_upper <- exp(ci[2])
    result_row$Test_Pvalue <- coefs["SNP_SumScore_SD", "Pr(>|z|)"]
    result_row$Train_OR_Reference <- train_or_ref
    result_row$Train_OR_Reference_CI_lower <- train_or_ref_ci_lower
    result_row$Train_OR_Reference_CI_upper <- train_or_ref_ci_upper
    result_row$Train_Pvalue_Reference <- train_pvalue_ref
    result_row$N_test_Reference_MasterTable <- as.integer(main_ref$N_test[1])
    result_row$R2_Train_Full <- r2_train_full
    result_row$R2_Train_CovOnly <- r2_train_covonly
    result_row$Delta_R2_Train <- delta_r2_train
    result_row$R2_Test_Full <- r2_test_full
    result_row$R2_Test_CovOnly <- r2_test_covonly
    result_row$Delta_R2_Test <- delta_r2_test

    result_row$Test_Same_Direction_As_Train <- dplyr::case_when(
      is.na(test_or) | is.na(train_or_ref) ~ NA,
      test_or > 1 & train_or_ref > 1 ~ TRUE,
      test_or < 1 & train_or_ref < 1 ~ TRUE,
      test_or == 1 & train_or_ref == 1 ~ TRUE,
      TRUE ~ FALSE
    )

  }, error = function(e) {
    result_row$Status <<- paste0("ERROR: ", conditionMessage(e))
  })

  result_row
}

# ==============================================================================
# Batch run over all (or restricted) pathways
# ==============================================================================

replication_results <- list()

for (pathway_name in pathways_to_run) {

  cat("\n==============================\n")
  cat("Pathway:", pathway_name, "\n")
  cat("==============================\n")

  row <- run_testset_replication(pathway_name)
  replication_results[[pathway_name]] <- row

  if (row$Status == "OK") {
    cat(
      "  OK - N_test:", row$N_test,
      "(master table reference:", row$N_test_Reference_MasterTable, ")",
      "| Test_OR (per SD):", round(row$Test_OR, 3),
      "| Train_OR_Reference (per SD):", round(row$Train_OR_Reference, 3),
      "| same direction:", row$Test_Same_Direction_As_Train,
      "| Delta_R2_Train:", round(row$Delta_R2_Train, 4),
      "| Delta_R2_Test:", round(row$Delta_R2_Test, 4),
      "\n"
    )
  } else {
    cat("  ", row$Status, "\n")
  }
}

final_table <- dplyr::bind_rows(replication_results)

output_path <- output_file("AllPathways_TestSet_Independent_Replication.csv")

write.csv2(
  final_table,
  output_path,
  row.names = FALSE
)

cat("\n====================================================\n")
cat("Done.\n")
cat("Result file:", output_path, "\n")
cat("Status overview:\n")
print(table(ifelse(final_table$Status == "OK", "OK", "ERROR")))
cat("====================================================\n")

#----------------------------------------------------------------
# Additionally write results automatically to the methods/results log
# (Appendix B), if that file exists (see script header).
#----------------------------------------------------------------
log_testset_replication_to_md(methods_log_path, final_table)
log_phase3_sessioninfo_to_md(methods_log_path)
