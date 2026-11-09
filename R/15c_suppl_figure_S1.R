# ==============================================================================
# 15c_suppl_figure_S1.R
# Supplementary Figure S1: functional form of age (linear vs. GAM), with and
# without adjustment for recruitment arm; age OR before/after adjustment.
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================

#================================================================================
# SUPPLEMENTARY FIGURE: sample_age - linear vs. GAM smooth, without vs. with batch
#================================================================================
# Purpose:
#   Builds the publication-ready 2-panel figure for the supplement
#   (Supplementary Figure S1), which visualizes the result of
#   15a_functional_form_age_gam.R:
#     Panel A (without batch adjustment): linear vs. GAM smooth curve for
#       sample_age.
#     Panel B (with batch adjustment): the same comparison after adjustment
#       for recruitment arm.
#
#   Standalone script; the data preparation is duplicated from 15a.
#================================================================================

library(data.table)
library(dplyr)
library(mgcv)
library(ggplot2)

#--------------------------------------------------------------------------------
# 0) Paths
#--------------------------------------------------------------------------------
base_dir <- Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus")

pheno_dir   <- Sys.getenv("FOCUS_PHENO_DIR", unset = "path/to/phenotypes")
pheno1_file <- file.path(pheno_dir, "BL_pheno_names_clean_CSV_UTF8_neu.csv")
pheno2_file <- file.path(pheno_dir, "FO_PA_alldata_CSV_UTF8_neu.csv")
pheno3_file <- file.path(pheno_dir, "FO_nutri_mima_CSV_UTF8_neu.csv")

genomewide_pcs_file <- file.path(base_dir, "Global_PCA/Global_Genomewide_PCs.rds")

fam_ambu_file <- file.path(base_dir, "Genetik_Daten_FoCus/FO_GSA_Ambu.fam")
fam_md_file   <- file.path(base_dir, "Genetik_Daten_FoCus/FO_GSA_MD.fam")

output_dir <- file.path(Sys.getenv("FOCUS_RESULTS_DIR", unset = "path/to/results"), "Outputs_FunctionalForm")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

#--------------------------------------------------------------------------------
# 1) Load and merge pheno data, derive covariates (as in 15a; see there for
#    detailed comments)
#--------------------------------------------------------------------------------
pheno1 <- read.csv2(pheno1_file, stringsAsFactors = FALSE)
pheno2 <- read.csv2(pheno2_file, stringsAsFactors = FALSE)
pheno3 <- read.csv2(pheno3_file, stringsAsFactors = FALSE)

setnames(pheno2, old = "ID_2020", new = "id")
possible_id_cols <- c("id", "ID", "ID_2020", "Id", "subject_id")
id_col_pheno3 <- intersect(possible_id_cols, colnames(pheno3))
setnames(pheno3, old = id_col_pheno3, new = "id")

pheno1$id <- as.character(pheno1$id)
pheno2$id <- as.character(pheno2$id)
pheno3$id <- as.character(pheno3$id)

pheno <- Reduce(function(x, y) merge(x, y, by = "id", all = TRUE), list(pheno1, pheno2, pheno3))
pheno <- as.data.table(pheno)

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
pheno[, energy_kcal := as.numeric(gsub(",", ".", as.character(energy_kcal)))]
pheno[, BMI := as.numeric(gsub(",", ".", as.character(BMI)))]
pheno[, obesity := ifelse(is.na(BMI), NA, ifelse(BMI >= 30, 1, 0))]

fam_ambu <- read.table(fam_ambu_file, header = FALSE, stringsAsFactors = FALSE)
fam_md   <- read.table(fam_md_file,   header = FALSE, stringsAsFactors = FALSE)
pheno[, Batch := NA_character_]
pheno[id %in% as.character(fam_ambu$V2), Batch := "Ambu"]
pheno[id %in% as.character(fam_md$V2), Batch := "MD"]
pheno[, Batch := factor(Batch, levels = c("MD", "Ambu"))]

global_gPCs <- readRDS(genomewide_pcs_file)
global_gPCs$id <- as.character(global_gPCs$id)
pheno$id <- as.character(pheno$id)
pheno <- merge(pheno, global_gPCs, by = "id", all.x = FALSE)

covariates_main_model <- c(
  "gender_cat", "smoking_cat", "school_cert_cat", "graduation_cat",
  "sample_age", "sports_activity_min_week", "everyday_activity_min_week",
  "energy_kcal", "gPC1", "gPC2", "gPC3"
)
model_vars <- c("obesity", covariates_main_model, "Batch")
df <- pheno[, ..model_vars]
df <- df[complete.cases(df[, ..covariates_main_model]) & !is.na(df$obesity) & !is.na(df$Batch)]
df$obesity <- as.numeric(as.character(df$obesity))

cat("n =", nrow(df), "| Batch:\n")
print(table(df$Batch))

#================================================================================
# 2) Four models: linear/smooth x without/with batch
#================================================================================
covariates_gam_linear <- setdiff(covariates_main_model, "sample_age")

# --- without batch (identical to models A/B in 15a) ---
formula_linear_noBatch <- reformulate(covariates_main_model, response = "obesity")
model_linear_noBatch <- glm(formula_linear_noBatch, data = df, family = binomial)

formula_smooth_noBatch <- as.formula(
  paste("obesity ~ s(sample_age) +", paste(covariates_gam_linear, collapse = " + "))
)
model_smooth_noBatch <- gam(formula_smooth_noBatch, data = df, family = binomial, method = "ML")

# --- with batch (model C family: same specification + batch) ---
formula_linear_withBatch <- reformulate(c(covariates_main_model, "Batch"), response = "obesity")
model_linear_withBatch <- glm(formula_linear_withBatch, data = df, family = binomial)

formula_smooth_withBatch <- as.formula(
  paste("obesity ~ s(sample_age) + Batch +", paste(covariates_gam_linear, collapse = " + "))
)
model_smooth_withBatch <- gam(formula_smooth_withBatch, data = df, family = binomial, method = "ML")

#================================================================================
# 3) Prediction grid (other covariates at mean/reference level, batch = "MD"
#    as reference category for the "with batch" curves, so that both panels are based on
#    the same population)
#================================================================================
age_grid <- data.frame(
  sample_age = seq(min(df$sample_age), max(df$sample_age), length.out = 200)
)
for (v in covariates_gam_linear) {
  if (is.factor(df[[v]])) {
    age_grid[[v]] <- factor(levels(df[[v]])[1], levels = levels(df[[v]]))
  } else {
    age_grid[[v]] <- mean(df[[v]], na.rm = TRUE)
  }
}
age_grid$Batch <- factor("MD", levels = levels(df$Batch))

pred <- function(model, grid) {
  p <- predict(model, newdata = grid, type = "link", se.fit = TRUE)
  data.frame(
    sample_age = grid$sample_age,
    fit = p$fit,
    lower = p$fit - 1.96 * p$se.fit,
    upper = p$fit + 1.96 * p$se.fit
  )
}

plot_df <- rbind(
  cbind(pred(model_linear_noBatch, age_grid), Model = "Linear", Panel = "A: Without adjustment for recruitment arm"),
  cbind(pred(model_smooth_noBatch, age_grid), Model = "GAM smooth", Panel = "A: Without adjustment for recruitment arm"),
  cbind(pred(model_linear_withBatch, age_grid), Model = "Linear", Panel = "B: With adjustment for recruitment arm"),
  cbind(pred(model_smooth_withBatch, age_grid), Model = "GAM smooth", Panel = "B: With adjustment for recruitment arm")
)
plot_df$Model <- factor(plot_df$Model, levels = c("Linear", "GAM smooth"))
plot_df$Panel  <- factor(
  plot_df$Panel,
  levels = c("A: Without adjustment for recruitment arm", "B: With adjustment for recruitment arm")
)

#================================================================================
# 4) Publication-ready 2-panel figure
#================================================================================
p_suppl <- ggplot(plot_df, aes(x = sample_age, y = fit, color = Model, fill = Model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ Panel) +
  scale_color_manual(values = c("Linear" = "#D55E00", "GAM smooth" = "#0072B2")) +
  scale_fill_manual(values = c("Linear" = "#D55E00", "GAM smooth" = "#0072B2")) +
  labs(
    x = "Age (years)",
    y = "Log-odds of obesity",
    color = "Model", fill = "Model"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(output_dir, "SupplFigure_S1_Age_Linearity_Batch.png"),
  plot = p_suppl, width = 9, height = 4.5, dpi = 400
)
ggsave(
  filename = file.path(output_dir, "SupplFigure_S1_Age_Linearity_Batch.pdf"),
  plot = p_suppl, width = 9, height = 4.5
)

cat("\nSaved:\n",
    file.path(output_dir, "SupplFigure_S1_Age_Linearity_Batch.png"), "\n",
    file.path(output_dir, "SupplFigure_S1_Age_Linearity_Batch.pdf"), "\n")

cat("\nedf without batch: ", summary(model_smooth_noBatch)$s.table[1, "edf"], "\n")
cat("edf with batch:    ", summary(model_smooth_withBatch)$s.table[1, "edf"], "\n")

#================================================================================
# 5) Age's OR/year, with vs. without recruitment-arm (Batch) adjustment
#    (linear models only - the parametric comparator to the edf-collapse result)
#================================================================================
or_ci_p <- function(model, term) {
  b   <- coef(summary(model))[term, "Estimate"]
  se  <- coef(summary(model))[term, "Std. Error"]
  p   <- coef(summary(model))[term, "Pr(>|z|)"]
  data.frame(
    OR = exp(b), CI_low = exp(b - 1.96 * se), CI_high = exp(b + 1.96 * se), p = p
  )
}

age_or_noArm   <- or_ci_p(model_linear_noBatch, "sample_age")
age_or_withArm <- or_ci_p(model_linear_withBatch, "sample_age")

age_or_compare <- rbind(
  cbind(Model = "Linear, without recruitment-arm adjustment", age_or_noArm),
  cbind(Model = "Linear, with recruitment-arm adjustment",    age_or_withArm)
)
age_or_compare$OR_per_decade <- age_or_compare$OR^10

cat("\n=== Age OR/year, with vs. without recruitment-arm adjustment ===\n")
print(age_or_compare, digits = 4)

write.csv(age_or_compare,
          file.path(output_dir, "Age_OR_Before_After_Batch.csv"),
          row.names = FALSE)
cat("\nSaved:", file.path(output_dir, "Age_OR_Before_After_Batch.csv"), "\n")
