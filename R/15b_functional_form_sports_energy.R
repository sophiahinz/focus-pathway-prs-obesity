# ==============================================================================
# 15b_functional_form_sports_energy.R
# Functional form of sports activity and energy intake, with and without
# adjustment for recruitment arm (Supplementary Table S4).
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================

#================================================================================
# SENSITIVITY CHECK (EXTENSION): sports_activity_min_week & energy_kcal --
# functional form, with vs. without recruitment-arm ("Batch") adjustment
#================================================================================
# Purpose:
#   Extends 15a_functional_form_age_gam.R / 15c_suppl_figure_S1.R to the two
#   other covariates that showed non-linearity in the unadjusted model D
#   (sports_activity_min_week, energy_kcal). everyday_activity_min_week showed
#   no non-linearity and is not re-checked here.
#
#   For each of the two covariates it is tested whether the curvature -
#   analogous to age (model C) - under adjustment for the recruitment arm
#   (outpatient clinic vs. population-based; coded as "Batch" with levels
#   "Ambu"/"MD") disappears completely or partially.
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
# 1) Load and merge pheno data, derive covariates (as in 15a)
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
# Reference level = population-based arm ("MD"), consistent with the age check.
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

cat("n =", nrow(df), "| obesity cases =", sum(df$obesity), sprintf("(%.1f%%)\n", 100 * mean(df$obesity)))
cat("Recruitment arm distribution:\n")
print(table(df$Batch))

#================================================================================
# Helper: for a given continuous covariate, fit linear/GAM x without/with
# recruitment-arm adjustment, run the LRT in each case, and return a summary
# row plus the four fitted models (for the prediction/plotting step below).
#================================================================================
run_covariate_check <- function(df, covariate, covariates_main_model) {

  covariates_other <- setdiff(covariates_main_model, covariate)

  # --- without recruitment-arm adjustment ---
  f_lin_noArm <- reformulate(covariates_main_model, response = "obesity")
  m_lin_noArm <- glm(f_lin_noArm, data = df, family = binomial)

  f_gam_noArm <- as.formula(
    paste0("obesity ~ s(", covariate, ") + ", paste(covariates_other, collapse = " + "))
  )
  m_gam_noArm <- gam(f_gam_noArm, data = df, family = binomial, method = "ML")

  lrt_noArm <- anova(m_lin_noArm, m_gam_noArm, test = "Chisq")
  edf_noArm <- summary(m_gam_noArm)$s.table[1, "edf"]
  p_noArm   <- summary(m_gam_noArm)$s.table[1, "p-value"]

  # --- with recruitment-arm adjustment ---
  f_lin_withArm <- reformulate(c(covariates_main_model, "Batch"), response = "obesity")
  m_lin_withArm <- glm(f_lin_withArm, data = df, family = binomial)

  f_gam_withArm <- as.formula(
    paste0("obesity ~ s(", covariate, ") + Batch + ", paste(covariates_other, collapse = " + "))
  )
  m_gam_withArm <- gam(f_gam_withArm, data = df, family = binomial, method = "ML")

  lrt_withArm <- anova(m_lin_withArm, m_gam_withArm, test = "Chisq")
  edf_withArm <- summary(m_gam_withArm)$s.table[1, "edf"]
  p_withArm   <- summary(m_gam_withArm)$s.table[1, "p-value"]

  cat("\n===", covariate, "===\n")
  cat(sprintf("  Without recruitment-arm adjustment: edf = %.2f, GAM-term p = %.4g\n", edf_noArm, p_noArm))
  print(lrt_noArm)
  cat(sprintf("  With recruitment-arm adjustment:    edf = %.2f, GAM-term p = %.4g\n", edf_withArm, p_withArm))
  print(lrt_withArm)

  list(
    covariate = covariate,
    summary = data.frame(
      Covariate = covariate,
      edf_no_arm = round(edf_noArm, 2),
      p_no_arm = signif(p_noArm, 3),
      LRT_p_no_arm = signif(lrt_noArm$`Pr(>Chi)`[2], 3),
      edf_with_arm = round(edf_withArm, 2),
      p_with_arm = signif(p_withArm, 3),
      LRT_p_with_arm = signif(lrt_withArm$`Pr(>Chi)`[2], 3)
    ),
    models = list(
      lin_noArm = m_lin_noArm, gam_noArm = m_gam_noArm,
      lin_withArm = m_lin_withArm, gam_withArm = m_gam_withArm
    )
  )
}

res_sports <- run_covariate_check(df, "sports_activity_min_week", covariates_main_model)
res_energy <- run_covariate_check(df, "energy_kcal", covariates_main_model)

summary_table <- rbind(res_sports$summary, res_energy$summary)
cat("\n=== Summary ===\n")
print(summary_table)
write.csv(summary_table,
          file.path(output_dir, "Sensitivity_SportsEnergy_Batch_Summary.csv"),
          row.names = FALSE)

#================================================================================
# Plot helper: 2-panel figure (linear vs. GAM smooth, without vs. with
# recruitment-arm adjustment), analogous to SupplFigure_S1 for age.
#================================================================================
make_two_panel_plot <- function(df, covariate, covariates_main_model, models,
                                 x_label, out_stub, x_range = NULL) {

  covariates_other <- setdiff(covariates_main_model, covariate)

  if (is.null(x_range)) {
    x_range <- range(df[[covariate]], na.rm = TRUE)
  }
  grid <- data.frame(seq(x_range[1], x_range[2], length.out = 200))
  names(grid) <- covariate

  for (v in covariates_other) {
    if (is.factor(df[[v]])) {
      grid[[v]] <- factor(levels(df[[v]])[1], levels = levels(df[[v]]))
    } else {
      grid[[v]] <- mean(df[[v]], na.rm = TRUE)
    }
  }
  grid$Batch <- factor("MD", levels = levels(df$Batch))  # population-based = reference

  pred <- function(model, grid) {
    p <- predict(model, newdata = grid, type = "link", se.fit = TRUE)
    out <- data.frame(x = grid[[covariate]], fit = p$fit,
                       lower = p$fit - 1.96 * p$se.fit, upper = p$fit + 1.96 * p$se.fit)
    out
  }

  plot_df <- rbind(
    cbind(pred(models$lin_noArm, grid),   Model = "linear",            Panel = "A: Without adjustment for recruitment arm"),
    cbind(pred(models$gam_noArm, grid),   Model = "GAM smooth",        Panel = "A: Without adjustment for recruitment arm"),
    cbind(pred(models$lin_withArm, grid), Model = "linear",            Panel = "B: With adjustment for recruitment arm"),
    cbind(pred(models$gam_withArm, grid), Model = "GAM smooth",        Panel = "B: With adjustment for recruitment arm")
  )
  plot_df$Model <- factor(plot_df$Model, levels = c("linear", "GAM smooth"))
  plot_df$Panel  <- factor(plot_df$Panel, levels = c(
    "A: Without adjustment for recruitment arm", "B: With adjustment for recruitment arm"
  ))

  p <- ggplot(plot_df, aes(x = x, y = fit, color = Model, fill = Model)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ Panel) +
    scale_color_manual(values = c("linear" = "#D55E00", "GAM smooth" = "#0072B2")) +
    scale_fill_manual(values = c("linear" = "#D55E00", "GAM smooth" = "#0072B2")) +
    labs(x = x_label, y = "Logit(P(obesity))", color = "Model", fill = "Model") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(output_dir, paste0(out_stub, ".png")), plot = p, width = 9, height = 4.5, dpi = 400)
  ggsave(file.path(output_dir, paste0(out_stub, ".pdf")), plot = p, width = 9, height = 4.5)
  cat("Saved:", file.path(output_dir, paste0(out_stub, ".png")), "\n")
}

# Sports activity: full range shown, but note the sparse high-value tail
# (>6,000 min/week) makes the right edge's confidence band very wide -
# consider restricting x_range to, e.g., the 1st-99th percentile if the
# full-range plot looks dominated by that tail.
make_two_panel_plot(
  df, "sports_activity_min_week", covariates_main_model,
  res_sports$models, "Sports activity (min/week)",
  "SupplFigure_SportsActivity_Linearity_Arm"
)

make_two_panel_plot(
  df, "energy_kcal", covariates_main_model,
  res_energy$models, "Energy intake (kcal/day)",
  "SupplFigure_EnergyIntake_Linearity_Arm"
)

cat("\nDONE.\n")
