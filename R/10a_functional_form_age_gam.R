# Pathway-restricted polygenic scores for energy-substrate metabolism and obesity (FoCus cohort)
# Script: 10a_functional_form_age_gam.R | Source file: Sensitivity_Age_Functional_Form_GAM_Check_REPAIR.R
# Local paths replaced by the environment variable FOCUS_BASE_DIR (project root); analysis code unchanged.
# Individual-level FoCus data are not included (see README).

# ==============================================================================
# REPAIR VERSION (2026-10-09): copy of Sensitivity_Age_Functional_Form_GAM_Check.R
# Re-run after correction of the defective phenotype CSVs. Changes ONLY to paths
# (lines marked "# REPAIR"); analysis code unchanged. Writes only to FoCus/Repair.
# ==============================================================================

#================================================================================
# SENSITIVITY CHECK: functional form of sample_age in the obesity covariate model
#================================================================================
# Purpose:
#   Checks whether the linearity of `sample_age` on the logit scale assumed in the main model
#   is justified (GAM with s(sample_age) vs. linear GLM) - WITHOUT
#   recomputing the existing pathway main models.
#
#   The covariate set is identical in all pathway models (only the
#   SNP sum score differs per pathway); the functional form of
#   sample_age is therefore a cohort-wide question that needs to be checked exactly once
#   - independent of the pathway-specific genotype extraction.
#
#   Data preparation (pheno merge, covariate derivation) is taken 1:1 from
#   Protocol_SNP_Analysis_FoCus_Pathway_Template_Split_19082026.R,
#   so that exactly the same sample/coding as in the main model is used.
#
# IMPORTANT:
#   `Batch` is, according to the covariate decision table,
#   ("batch only as sensitivity analysis") NOT included in the main model. Models A/B below
#   therefore replicate the main model exactly (without batch). Model C (section 5,
#   optional) additionally adds batch - this answers in passing whether the
#   age form/direction changes when adjusting for Ambu/MD (not only for
#   the SNP score, as in the existing batch sensitivity analysis).
#
#   This file changes NOTHING in existing scripts/outputs - it is a
#   independent, additional check.
#================================================================================

library(data.table)
library(dplyr)
library(mgcv)
library(ggplot2)

#--------------------------------------------------------------------------------
# 0) Paths (identical to the main script)
#--------------------------------------------------------------------------------
base_dir <- Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus")

pheno1_file <- file.path(base_dir, "FoCus Daten/Daten/CSV_repariert/BL_pheno_names_clean_CSV_UTF8_neu.csv")   # REPAIR
pheno2_file <- file.path(base_dir, "FoCus Daten/Daten/CSV_repariert/FO_PA_alldata_CSV_UTF8_neu.csv")   # REPAIR
pheno3_file <- file.path(base_dir, "FoCus Daten/Daten/CSV_repariert/FO_nutri_mima_CSV_UTF8_neu.csv")   # REPAIR

genomewide_pcs_file <- file.path(base_dir, "Global_PCA/Global_Genomewide_PCs.rds")

fam_ambu_file <- file.path(base_dir, "Genetik_Daten_FoCus/FO_GSA_Ambu.fam")
fam_md_file   <- file.path(base_dir, "Genetik_Daten_FoCus/FO_GSA_MD.fam")

output_dir <- file.path(base_dir, "Repair", "Outputs_FunctionalForm")   # REPAIR
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

#--------------------------------------------------------------------------------
# 1) Load and merge pheno data (identical to the logic in the main script)
#--------------------------------------------------------------------------------
pheno1 <- read.csv2(pheno1_file, stringsAsFactors = FALSE)
pheno2 <- read.csv2(pheno2_file, stringsAsFactors = FALSE)
pheno3 <- read.csv2(pheno3_file, stringsAsFactors = FALSE)

stopifnot("id" %in% colnames(pheno1))
stopifnot("ID_2020" %in% colnames(pheno2))
setnames(pheno2, old = "ID_2020", new = "id")

possible_id_cols <- c("id", "ID", "ID_2020", "Id", "subject_id")
id_col_pheno3 <- intersect(possible_id_cols, colnames(pheno3))
if (length(id_col_pheno3) != 1) {
  stop(
    "No unique ID column found in pheno3: ",
    paste(colnames(pheno3), collapse = ", ")
  )
}
setnames(pheno3, old = id_col_pheno3, new = "id")

pheno1$id <- as.character(pheno1$id)
pheno2$id <- as.character(pheno2$id)
pheno3$id <- as.character(pheno3$id)

pheno <- Reduce(
  function(x, y) merge(x, y, by = "id", all = TRUE),
  list(pheno1, pheno2, pheno3)
)
pheno <- as.data.table(pheno)

cat("Pheno merged, n =", nrow(pheno), "\n")

#--------------------------------------------------------------------------------
# 2) Derive covariates (1:1 as in the main script, section "Check global
#    covariate availability" / gender_cat.../obesity)
#--------------------------------------------------------------------------------
pheno[, gender_num := as.numeric(as.character(gender))]
pheno[gender_num %in% c(0, 7777777, 9999999), gender_num := NA]
pheno[, gender_cat := factor(gender_num, levels = c(1, 2), labels = c("male", "female"))]

pheno[, smoking_habits_3m_num := as.numeric(as.character(smoking_habits_3m))]
pheno[smoking_habits_3m_num %in% c(0, 7777777, 9999999), smoking_habits_3m_num := NA]
pheno[, smoking_cat := factor(
  dplyr::case_when(
    smoking_habits_3m_num %in% c(1, 2) ~ "never_smoker",
    smoking_habits_3m_num == 3 ~ "former_smoker",
    smoking_habits_3m_num == 4 ~ "current_smoker",
    TRUE ~ NA_character_
  ),
  levels = c("never_smoker", "former_smoker", "current_smoker")
)]

pheno[, school_certificate_A_num := as.numeric(as.character(school_certificate_A))]
pheno[!school_certificate_A_num %in% 1:5, school_certificate_A_num := NA]
pheno[, school_cert_cat := factor(
  ifelse(school_certificate_A_num %in% c(1, 2, 5), "risk", "reference"),
  levels = c("reference", "risk")
)]

pheno[, graduation_num := as.numeric(as.character(graduation))]
pheno[graduation_num %in% c(6, 7777777, 9999999), graduation_num := NA]
pheno[, graduation_cat := factor(
  ifelse(graduation_num %in% c(1, 2, 5), "risk", "reference"),
  levels = c("reference", "risk")
)]

pheno[, sample_age := as.numeric(as.character(sample_age))]
pheno[, sports_activity_min_week := as.numeric(as.character(sports_activity_min_week))]
pheno[, everyday_activity_min_week := as.numeric(as.character(everyday_activity_min_week))]

if ("energy_kcal" %in% colnames(pheno)) {
  pheno[, energy_kcal := as.numeric(gsub(",", ".", as.character(energy_kcal)))]
} else if ("energy" %in% colnames(pheno)) {
  pheno[, energy_kJ := as.numeric(gsub(",", ".", as.character(energy)))]
  pheno[, energy_kcal := energy_kJ / 4.184]
} else {
  pheno[, energy_kcal := NA_real_]
}

pheno[, BMI := as.numeric(gsub(",", ".", as.character(BMI)))]
pheno[, obesity := ifelse(is.na(BMI), NA, ifelse(BMI >= 30, 1, 0))]

#--------------------------------------------------------------------------------
# 3) Batch (Ambu/MD) - genome-wide, cohort-wide assignment
#    (deliberately NOT via snpStats::read.plink() - we only need the IDs from
#    the .fam file, not the full genotypes; this is faster and
#    independent of any individual pathway extraction)
#--------------------------------------------------------------------------------
fam_ambu <- read.table(fam_ambu_file, header = FALSE, stringsAsFactors = FALSE)
fam_md   <- read.table(fam_md_file,   header = FALSE, stringsAsFactors = FALSE)

# Standard PLINK .fam columns: FID IID PAT MAT SEX PHENOTYPE -> IID = column 2
ids_ambu <- as.character(fam_ambu$V2)
ids_md   <- as.character(fam_md$V2)

pheno[, Batch := NA_character_]
pheno[id %in% ids_ambu, Batch := "Ambu"]
pheno[id %in% ids_md, Batch := "MD"]
pheno[, Batch := factor(Batch, levels = c("MD", "Ambu"))]

cat("\nBatch distribution (entire cohort, before complete-case filter):\n")
print(table(pheno$Batch, useNA = "ifany"))

#--------------------------------------------------------------------------------
# 4) Merge genome-wide PCs (identical to the main script)
#--------------------------------------------------------------------------------
if (!file.exists(genomewide_pcs_file)) {
  stop(
    "Genome-wide PCs not found (", genomewide_pcs_file, "). ",
    "Please run the PCA section in the gene script first."
  )
}

global_gPCs <- readRDS(genomewide_pcs_file)
global_gPCs$id <- as.character(global_gPCs$id)

pheno$id <- as.character(pheno$id)
n_before_gpc_join <- nrow(pheno)

pheno <- merge(pheno, global_gPCs, by = "id", all.x = FALSE)

if (nrow(pheno) < n_before_gpc_join) {
  warning(
    n_before_gpc_join - nrow(pheno),
    " individuals have no genome-wide PCs and were removed."
  )
}

#--------------------------------------------------------------------------------
# 5) Analytical dataset: complete cases on the main-model covariate set
#--------------------------------------------------------------------------------
covariates_main_model <- c(
  "gender_cat", "smoking_cat", "school_cert_cat", "graduation_cat",
  "sample_age", "sports_activity_min_week", "everyday_activity_min_week",
  "energy_kcal", "gPC1", "gPC2", "gPC3"
)

model_vars <- c("obesity", covariates_main_model, "Batch")

df <- pheno[, ..model_vars]
df <- df[complete.cases(df[, ..covariates_main_model]) & !is.na(df$obesity)]

cat("\nAnalytical dataset (complete cases, main-model covariates + obesity):\n")
cat("n =", nrow(df), " | obesity cases:", sum(df$obesity == 1),
    " (", round(100 * mean(df$obesity == 1), 1), "%)\n", sep = "")
cat("Batch distribution in this dataset:\n")
print(table(df$Batch, useNA = "ifany"))

df$obesity <- as.numeric(as.character(df$obesity))

#================================================================================
# 6) MODEL A: linear GLM - replicates the main model exactly (without batch)
#================================================================================
formula_linear <- reformulate(covariates_main_model, response = "obesity")

model_A_linear <- glm(
  formula_linear,
  data   = df,
  family = binomial
)

cat("\n================ MODEL A: linear main model (reference) ================\n")
print(summary(model_A_linear)$coefficients["sample_age", ])
cat("OR per year:", exp(coef(model_A_linear)["sample_age"]), "\n")
cat("OR per decade:", exp(coef(model_A_linear)["sample_age"] * 10), "\n")

#================================================================================
# 7) MODEL B: GAM with s(sample_age) - the actual sensitivity check
#    method = "ML" on both sides so that anova() provides a valid LRT
#================================================================================
covariates_gam_linear <- setdiff(covariates_main_model, "sample_age")

formula_gam_linear <- as.formula(
  paste("obesity ~ sample_age +", paste(covariates_gam_linear, collapse = " + "))
)
formula_gam_smooth <- as.formula(
  paste("obesity ~ s(sample_age) +", paste(covariates_gam_linear, collapse = " + "))
)

model_B0_gam_linear <- gam(
  formula_gam_linear,
  data   = df,
  family = binomial,
  method = "ML"
)

model_B_gam_smooth <- gam(
  formula_gam_smooth,
  data   = df,
  family = binomial,
  method = "ML"
)

cat("\n================ MODEL B: GAM with s(sample_age) ================\n")
print(summary(model_B_gam_smooth)$s.table)
cat(
  "\nInterpretation aid: edf ~= 1 and a non-significant smooth term\n",
  "indicate linearity; edf clearly > 1 and/or a significant\n",
  "smooth term indicate true curvature.\n"
)

cat("\nFormal comparison linear vs. smooth (LRT, method = 'ML'):\n")
print(anova(model_B0_gam_linear, model_B_gam_smooth, test = "Chisq"))

cat("\nAIC comparison (smaller = better):\n")
cat("  Model A (glm, linear):      ", AIC(model_A_linear), "\n")
cat("  Model B0 (gam, linear):     ", AIC(model_B0_gam_linear), "\n")
cat("  Model B (gam, smooth):      ", AIC(model_B_gam_smooth), "\n")

#================================================================================
# 7b) MODEL D: all four continuous covariates simultaneously as smooth
#     -> answers the general question ("does this apply only to age or also to
#     sports_activity_min_week / everyday_activity_min_week / energy_kcal?")
#     in one go, instead of building four individual checks.
#================================================================================
covariates_continuous  <- c(
  "sample_age", "sports_activity_min_week", "everyday_activity_min_week", "energy_kcal"
)
covariates_categorical <- setdiff(covariates_main_model, covariates_continuous)

formula_gam_full_smooth <- as.formula(
  paste(
    "obesity ~",
    paste(sprintf("s(%s)", covariates_continuous), collapse = " + "),
    "+",
    paste(covariates_categorical, collapse = " + ")
  )
)

model_D_gam_full <- gam(
  formula_gam_full_smooth,
  data   = df,
  family = binomial,
  method = "ML"
)

cat("\n================ MODEL D: GAM with smooth for ALL continuous covariates ================\n")
print(summary(model_D_gam_full)$s.table)
cat(
  "\nInterpretation aid (per row/covariate): edf ~= 1 and p >= .05 -> linearity\n",
  "unremarkable; edf clearly > 1 and/or p < .05 -> curvature, look closer.\n"
)

cat("\nAIC: full linear model A vs. full smooth model D:\n")
cat("  Model A (all linear):        ", AIC(model_A_linear), "\n")
cat("  Model D (all 4 continuous as smooth):", AIC(model_D_gam_full), "\n")

# Diagnostic plot: one panel per continuous covariate (mgcv's own plotting
# is entirely sufficient here - no custom ggplot per variable needed)
png(
  filename = file.path(output_dir, "Sensitivity_AllCovariates_Functional_Form_GAM_Plot.png"),
  width = 2400, height = 2000, res = 300
)
par(mfrow = c(2, 2))
plot(model_D_gam_full, shade = TRUE, seWithMean = TRUE, residuals = FALSE)
dev.off()

cat(
  "\nPlot gespeichert:",
  file.path(output_dir, "Sensitivity_AllCovariates_Functional_Form_GAM_Plot.png"), "\n"
)

#================================================================================
# 8) MODEL C (optional/bonus): as B, but additionally with batch
#    -> answers in passing whether the form/direction of sample_age changes,
#       when adjusting for Ambu/MD (complement to the existing
#       SNP score batch sensitivity)
#================================================================================
df_batch <- df[!is.na(df$Batch)]

if (length(unique(df_batch$Batch)) >= 2) {

  formula_gam_smooth_batch <- as.formula(
    paste("obesity ~ s(sample_age) + Batch +", paste(covariates_gam_linear, collapse = " + "))
  )

  model_C_gam_smooth_batch <- gam(
    formula_gam_smooth_batch,
    data   = df_batch,
    family = binomial,
    method = "ML"
  )

  cat("\n================ MODEL C (bonus): GAM + batch ================\n")
  cat("n =", nrow(df_batch), "\n")
  print(summary(model_C_gam_smooth_batch)$s.table)
  print(summary(model_C_gam_smooth_batch)$p.table["BatchAmbu", ])
  cat(
    "\nComparison with model B: does the form/significance of s(sample_age) remain stable after\n",
    "batch adjustment? If so, this argues against Ambu/MD\n",
    "confounding of the age relationship itself (in addition to the already shown\n",
    "Robustheit des SNP-Scores).\n"
  )
} else {
  cat("\nModel C skipped: not enough batch levels in df_batch.\n")
  model_C_gam_smooth_batch <- NULL
}

#================================================================================
# 9) Diagnostic plot: partial effect of sample_age (GAM) vs. linear fit
#================================================================================
age_grid <- data.frame(
  sample_age = seq(min(df$sample_age), max(df$sample_age), length.out = 200)
)
for (v in covariates_gam_linear) {
  if (is.factor(df[[v]])) {
    age_grid[[v]] <- levels(df[[v]])[1]
    age_grid[[v]] <- factor(age_grid[[v]], levels = levels(df[[v]]))
  } else {
    age_grid[[v]] <- mean(df[[v]], na.rm = TRUE)
  }
}

pred_linear <- predict(model_A_linear, newdata = age_grid, type = "link", se.fit = TRUE)
pred_smooth <- predict(model_B_gam_smooth, newdata = age_grid, type = "link", se.fit = TRUE)

plot_df <- rbind(
  data.frame(
    sample_age = age_grid$sample_age,
    fit = pred_linear$fit,
    lower = pred_linear$fit - 1.96 * pred_linear$se.fit,
    upper = pred_linear$fit + 1.96 * pred_linear$se.fit,
    Modell = "A: linear (Hauptmodell)"
  ),
  data.frame(
    sample_age = age_grid$sample_age,
    fit = pred_smooth$fit,
    lower = pred_smooth$fit - 1.96 * pred_smooth$se.fit,
    upper = pred_smooth$fit + 1.96 * pred_smooth$se.fit,
    Modell = "B: GAM s(sample_age)"
  )
)

p_age_form <- ggplot(plot_df, aes(x = sample_age, y = fit, color = Modell, fill = Modell)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  labs(
    x = "sample_age (Jahre)",
    y = "Logit(P(Adipositas))",
    title = "Functional form of sample_age: linear vs. GAM smooth",
    subtitle = "andere Kovariaten am Mittelwert/Referenzlevel gehalten"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(output_dir, "Sensitivity_Age_Functional_Form_GAM_Plot.png"),
  plot = p_age_form, width = 8, height = 5.5, dpi = 300
)

cat("\nPlot saved:",
    file.path(output_dir, "Sensitivity_Age_Functional_Form_GAM_Plot.png"), "\n")

#================================================================================
# 10) Compact summary for the methods/results log
#================================================================================
s_table <- summary(model_B_gam_smooth)$s.table
lrt <- anova(model_B0_gam_linear, model_B_gam_smooth, test = "Chisq")

summary_out <- data.frame(
  N = nrow(df),
  N_obesity = sum(df$obesity == 1),
  OR_per_decade_linear_model = exp(coef(model_A_linear)["sample_age"] * 10),
  edf_smooth_age = s_table[1, "edf"],
  p_smooth_age = s_table[1, "p-value"],
  LRT_Chisq_age = lrt[2, "Deviance"],
  LRT_df_age = lrt[2, "Df"],
  LRT_p_age = lrt[2, "Pr(>Chi)"],
  AIC_linear_glm = AIC(model_A_linear),
  AIC_gam_linear_age = AIC(model_B0_gam_linear),
  AIC_gam_smooth_age = AIC(model_B_gam_smooth)
)

print(summary_out)

write.csv2(
  summary_out,
  file.path(output_dir, "Sensitivity_Age_Functional_Form_GAM_Summary.csv"),
  row.names = FALSE
)

cat("\nSummary (age only) saved:",
    file.path(output_dir, "Sensitivity_Age_Functional_Form_GAM_Summary.csv"), "\n")

# Summary for ALL four continuous covariates (model D)
s_table_all <- summary(model_D_gam_full)$s.table
summary_out_all <- data.frame(
  Variable = covariates_continuous,
  edf      = s_table_all[, "edf"],
  p_value  = s_table_all[, "p-value"]
)

print(summary_out_all)

write.csv2(
  summary_out_all,
  file.path(output_dir, "Sensitivity_AllCovariates_Functional_Form_GAM_Summary.csv"),
  row.names = FALSE
)

cat("\nSummary (all 4 continuous covariates) saved:",
    file.path(output_dir, "Sensitivity_AllCovariates_Functional_Form_GAM_Summary.csv"), "\n")

cat("\n================================================================\n")
cat("DONE. Existing pathway main models were NOT changed.\n")
cat("================================================================\n")

