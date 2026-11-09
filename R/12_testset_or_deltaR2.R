# ==============================================================================
# 12_testset_or_deltaR2.R
# Test-set OR per SD and incremental McFadden pseudo-R^2 (delta R^2) of the
# multivariable main model, for all 54 pathways (Table 2; Figure 3; Results).
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================

#==================================================================================
# Uses the saved pathway outputs of 07_pathway_pipeline.R (df_final RDS +
# Global_TrainTest_Split RDS per pathway); no re-run of genotype QC or score
# construction is needed.
#
# Model specification identical to the main model (function
# run_multivariate_sumscore_analysis(), section "9) Multivariable model
# Train/Test" in 07_pathway_pipeline.R):
#   - Score standardization (SNP_SumScore_SD) with training-set mean/SD,
#     applied unchanged to the test set
#   - Full model: obesity ~ SNP_SumScore_SD + covariates_present (the same
#     11 covariates as in the main model, manual_covariate_decisions$default$keep)
#   - Covariates-only model (for delta R^2): obesity ~ covariates_present,
#     identical to formula_no_snp/model_no_snp in 07_pathway_pipeline.R
#   - Pseudo-R^2: McFadden (pscl::pR2()), the same metric that
#     07_pathway_pipeline.R stores as "McFadden_R2" for the training set
#   - Full model and covariates-only model are fitted separately on the
#     training set and on the test set.
#==================================================================================

library(dplyr)
library(broom)
library(readr)
library(pscl)   # for pR2() / McFadden pseudo-R^2

#----------------------------------------------------------------------------------
# 0) Paths / helper functions - identical to 07_pathway_pipeline.R
#----------------------------------------------------------------------------------

output_dir <- file.path(Sys.getenv("FOCUS_RESULTS_DIR", unset = "path/to/results"), "OutPW")

output_file <- function(...) {
  file.path(output_dir, paste0(...))
}

make_analysis_suffix <- function(config, score_mode = NULL) {
  suffix <- config$label
  if (!is.null(score_mode)) {
    suffix <- paste(suffix, tools::toTitleCase(score_mode), sep = "_")
  }
  suffix
}

make_output_name <- function(pathway_name, analysis_type, analysis_name,
                              result_name, config, score_mode = NULL,
                              ext = "csv") {
  output_file(
    "Pathway_", pathway_name, "_", analysis_type, "_", analysis_name, "_",
    result_name, "_", make_analysis_suffix(config, score_mode), ".", ext
  )
}

# Main-model configuration - identical to analysis_configs$main_ld in 07_pathway_pipeline.R
config <- list(label = "01_Main_LDpruned_All", type = "all_ld")
score_mode <- "oriented"

# Covariates - identical to manual_covariate_decisions$default$keep in 07_pathway_pipeline.R
covariates_present <- c(
  "gender_cat", "smoking_cat", "school_cert_cat", "graduation_cat",
  "sample_age", "sports_activity_min_week", "everyday_activity_min_week",
  "energy_kcal", "gPC1", "gPC2", "gPC3"
)

model_formula_standardized <- reformulate(
  termlabels = c("SNP_SumScore_SD", covariates_present),
  response = "obesity"
)

#----------------------------------------------------------------------------------
# 1) Pathway list - all 54 analysed pathways (Table 2 later filters on
#    Robust_Candidate == TRUE from AllPathways_Results_Tracking_Table.csv)
#----------------------------------------------------------------------------------
pathway_liste_file <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Pathway_Selection/Pathway_Liste.csv")
pathway_names <- read.csv(pathway_liste_file, stringsAsFactors = FALSE)$pathway_name
stopifnot(length(pathway_names) == 54)

#----------------------------------------------------------------------------------
# 2) Core function: test-set OR + delta R^2 for one pathway
#----------------------------------------------------------------------------------

# Covariates-only formula - identical to formula_no_snp in 07_pathway_pipeline.R
formula_no_snp <- reformulate(termlabels = covariates_present, response = "obesity")

compute_test_set_or <- function(pathway_name, config, score_mode,
                                 covariates_present, model_formula_standardized,
                                 formula_no_snp) {

  cat("\n----------------------------------------\n")
  cat("Test-Set-OR:", pathway_name, "\n")
  cat("----------------------------------------\n")

  # --- Load df_final (main-model analysis dataset, already QC/complete-case) ---
  df_final_path <- make_output_name(
    pathway_name, "Multivariate", "SNP_SumScore", "Individuals",
    config, score_mode, ext = "rds"
  )

  if (!file.exists(df_final_path)) {
    warning("df_final RDS not found for ", pathway_name, ": ", df_final_path)
    return(NULL)
  }

  df_final <- readRDS(df_final_path)
  df_final$id <- as.character(df_final$id)

  # --- Load train/test split (pathway-wide, independent of config/score_mode) ---
  split_path <- output_file("Pathway_", pathway_name, "_Global_TrainTest_Split.rds")

  if (!file.exists(split_path)) {
    warning("Global_TrainTest_Split RDS not found for ", pathway_name, ": ", split_path)
    return(NULL)
  }

  global_split <- readRDS(split_path)
  global_split$train_ids <- as.character(global_split$train_ids)
  global_split$test_ids  <- as.character(global_split$test_ids)

  train_set <- df_final %>% filter(id %in% global_split$train_ids)
  test_set  <- df_final %>% filter(id %in% global_split$test_ids)

  if (nrow(train_set) < 20 || nrow(test_set) < 10) {
    warning("Too few individuals in the train/test dataset for ", pathway_name,
            ". Train: ", nrow(train_set), ", Test: ", nrow(test_set))
    return(NULL)
  }

  if (length(unique(test_set$obesity)) < 2) {
    warning("Test dataset contains only one outcome class for ", pathway_name)
    return(NULL)
  }

  # --- Score standardization with training-set mean/SD (identical to 07_pathway_pipeline.R) ---
  score_mean_train <- mean(train_set$SNP_SumScore, na.rm = TRUE)
  score_sd_train   <- sd(train_set$SNP_SumScore,   na.rm = TRUE)

  if (!is.finite(score_sd_train) || score_sd_train <= 0) {
    warning("SNP sum score has insufficient variation for ", pathway_name)
    return(NULL)
  }

  test_set$SNP_SumScore_SD <- (test_set$SNP_SumScore - score_mean_train) / score_sd_train

  # --- Multivariable main model fitted on the test set ---
  model_test <- tryCatch(
    glm(model_formula_standardized, data = test_set, family = binomial),
    error = function(e) {
      warning("Test-set model fit failed for ", pathway_name, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(model_test)) return(NULL)

  test_coef <- broom::tidy(model_test, conf.int = TRUE, conf.level = 0.95) %>%
    filter(term == "SNP_SumScore_SD")

  if (nrow(test_coef) == 0) {
    warning("No SNP_SumScore_SD coefficient in the test-set model for ", pathway_name)
    return(NULL)
  }

  # --- Read training result from existing Regression_Results CSV ---
  # (full model not refitted - identical to the already reported
  # Table 2 values, incl. the McFadden_R2 already stored there for the
  # full training model)
  train_csv_path <- make_output_name(
    pathway_name, "Multivariate", "SNP_SumScore", "Regression_Results",
    config, score_mode, ext = "csv"
  )

  train_or <- NA_real_
  train_ci_low <- NA_real_
  train_ci_high <- NA_real_
  r2_train_full <- NA_real_

  if (file.exists(train_csv_path)) {
    train_res <- readr::read_csv2(train_csv_path, show_col_types = FALSE)
    train_or      <- train_res$OR_Score_Per_SD[1]
    train_ci_low  <- train_res$OR_Score_Per_SD_CI_lower[1]
    train_ci_high <- train_res$OR_Score_Per_SD_CI_upper[1]
    r2_train_full <- train_res$McFadden_R2[1]
  } else {
    warning("Training Regression_Results CSV not found for ", pathway_name,
            " - training values remain NA (only test OR/R^2 are reported).")
  }

  # --- Delta R^2: covariates-only models for train and test ---
  # (full-model R^2 for test comes from model_test above; for train from the
  # existing CSV, see above - no refit of the full training model)

  model_train_no_snp <- tryCatch(
    glm(formula_no_snp, data = train_set, family = binomial),
    error = function(e) {
      warning("Covariates-only model (train) failed for ", pathway_name,
              ": ", conditionMessage(e))
      NULL
    }
  )

  model_test_no_snp <- tryCatch(
    glm(formula_no_snp, data = test_set, family = binomial),
    error = function(e) {
      warning("Covariates-only model (test) failed for ", pathway_name,
              ": ", conditionMessage(e))
      NULL
    }
  )

  r2_train_covonly <- if (!is.null(model_train_no_snp)) {
    tryCatch(as.numeric(pscl::pR2(model_train_no_snp)["McFadden"]), error = function(e) NA_real_)
  } else NA_real_

  r2_test_full <- tryCatch(as.numeric(pscl::pR2(model_test)["McFadden"]), error = function(e) NA_real_)

  r2_test_covonly <- if (!is.null(model_test_no_snp)) {
    tryCatch(as.numeric(pscl::pR2(model_test_no_snp)["McFadden"]), error = function(e) NA_real_)
  } else NA_real_

  data.frame(
    Pathway            = pathway_name,
    N_train            = nrow(train_set),
    N_test             = nrow(test_set),
    OR_train           = train_or,
    CI_low_train       = train_ci_low,
    CI_high_train      = train_ci_high,
    OR_test            = exp(test_coef$estimate),
    CI_low_test        = exp(test_coef$conf.low),
    CI_high_test       = exp(test_coef$conf.high),
    p_test             = test_coef$p.value,
    Direction_train    = sign(log(train_or)),
    Direction_test     = sign(test_coef$estimate),
    R2_train_full      = r2_train_full,
    R2_train_covonly   = r2_train_covonly,
    Delta_R2_train     = r2_train_full - r2_train_covonly,
    R2_test_full       = r2_test_full,
    R2_test_covonly    = r2_test_covonly,
    Delta_R2_test      = r2_test_full - r2_test_covonly,
    stringsAsFactors   = FALSE
  ) %>%
    mutate(Direction_consistent = Direction_train == Direction_test)
}

#----------------------------------------------------------------------------------
# 3) Loop over all pathways
#----------------------------------------------------------------------------------

test_set_results <- lapply(pathway_names, function(pw) {
  tryCatch(
    compute_test_set_or(pw, config, score_mode, covariates_present,
                         model_formula_standardized, formula_no_snp),
    error = function(e) {
      warning("Error for ", pw, ": ", conditionMessage(e))
      NULL
    }
  )
}) %>%
  bind_rows()

cat("\n========================================\n")
cat("Done:", nrow(test_set_results), "of", length(pathway_names), "pathways successful.\n")
cat("========================================\n")

print(test_set_results)

#----------------------------------------------------------------------------------
# 4) Export - Table 2 (test OR) and Results (delta R^2)
#----------------------------------------------------------------------------------

write.csv2(
  test_set_results,
  make_output_name("ALL", "Multivariate", "SNP_SumScore",
                    "TestSet_Regression_Results_DeltaR2", config, score_mode),
  row.names = FALSE
)
