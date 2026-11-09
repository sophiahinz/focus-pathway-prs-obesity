# ==============================================================================
# 14b_snp_level_models.R
# SNP-level functional-concentration analysis, step 2: effect-magnitude and
# risk-direction-concentration mixed-effects models (Supplementary Tables S9, S10).
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================

# ============================================================
# Models, pooled across the robust candidate pathways
# ============================================================
#
# Tiers (VEP v116; Supplementary Methods 2), priority for multi-tier SNPs:
#   LoF > Missense > Regulatory/UTR > Reference (intronic/synonymous).
#   SNPs that are "functional" under the binary rule of
#   load_functional_vep_snps() but match none of the named tiers
#   (splice_region_variant / stop_lost / protein_altering_variant) are
#   assigned to Reference; VEP classifies splice_region_variant as LOW
#   impact, i.e. the same impact class as synonymous_variant. The number of
#   reclassified SNPs is reported (console + descriptive summary).
#   SIFT/PolyPhen: any-transcript, worst prediction per predictor, combined
#   via OR (either sufficient) for missense-higher vs. -lower pathogenicity.
#   SNPs occurring in >1 pathway are kept at their first occurrence
#   (iteration order = Robust_Candidate order in the tracking table).
#
# Missense is analysed in two ways, fitted as two separate models:
#   - "Split" : Missense_Higher / Missense_Lower vs. Reference, Regulatory_UTR
#   - "Pooled": Missense (not split) vs. Reference, Regulatory_UTR
#   The p-values of both models form one test family (4 tiers x 2
#   hypotheses = 8 tests) with one joint BH-FDR correction. The Pooled
#   model's Regulatory_UTR estimate is dropped as redundant.
#   LoF is reported descriptively only (too few SNPs for formal models).
#   The risk-direction comparison against a 50% null is reported
#   descriptively (proportion + 95% CI per tier).
#
# Only aggregate per-SNP files are read (VEP annotation, univariate
# regression results); no individual-level data are used.
#
# Output:
#   - Q6_Results_Split.csv       (Split model family)
#   - Q6_Results_Pooled.csv      (Pooled model family)
#   - Q6_Descriptive_Tier_Summary_v2.csv (all tiers incl. LoF and a
#     Missense-pooled row; N, mean|logOR|, risk proportion, 95% CI)
#   - Console: model summaries (lmer/glmer) for both model families.
# ============================================================

library(data.table)
library(dplyr)
library(lme4)
library(lmerTest)   # adds Satterthwaite p-values to lmer() output

# ---- Paths ----

vep_dir         <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Outputs_Pathway")   # VEP annotation files (read only)
output_dir      <- file.path(Sys.getenv("FOCUS_RESULTS_DIR", unset = "path/to/results"), "OutPW")       # tracking table, univariate results, outputs
tracking_table <- file.path(output_dir, "AllPathways_Results_Tracking_Table.csv")

out_split_csv       <- file.path(output_dir, "Q6_Results_Split.csv")
out_pooled_csv       <- file.path(output_dir, "Q6_Results_Pooled.csv")
out_descriptive_csv <- file.path(output_dir, "Q6_Descriptive_Tier_Summary_v2.csv")

# ---- Pathway list: same source/order as the permutation test (11) ----

tracking <- fread(tracking_table, sep = ";", dec = ",", encoding = "UTF-8")
PATHWAYS_FINAL <- tracking[Robust_Candidate == TRUE, ]$Pathway
cat("N robust candidate pathways:", length(PATHWAYS_FINAL), "\n\n")

# ---- Term lists (identical to load_functional_vep_snps() / load_lof_vep_annotation()
#      in 07_pathway_pipeline.R) ----

lof_terms <- c(
  "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
  "stop_gained", "frameshift_variant", "start_lost"
)
missense_term <- "missense_variant"
reg_utr_terms <- c(
  "regulatory_region_variant", "TF_binding_site_variant", "TFBS_ablation",
  "TFBS_amplification", "5_prime_UTR_variant", "3_prime_UTR_variant",
  "mature_miRNA_variant"
)
other_coding_terms <- c(
  "protein_altering_variant", "stop_lost", "splice_region_variant"
)

sift_deleterious_prefixes  <- c("deleterious", "deleterious_low_confidence")
polyphen_damaging_prefixes <- c("probably_damaging", "possibly_damaging")

# ---- Per-pathway VEP classification (Unclassified_Functional SNPs are
#      assigned to Reference, see header) ----

classify_pathway_vep <- function(pathway) {

  vep_file <- file.path(vep_dir, paste0("Pathway_", pathway, "_VEP.txt"))
  vep <- fread(vep_file, encoding = "UTF-8")

  vep[, SNP := trimws(as.character(`#Uploaded_variation`))]
  vep <- vep[!is.na(SNP) & SNP != ""]

  vep[, is_lof         := grepl(paste(lof_terms, collapse = "|"), Consequence, ignore.case = TRUE)]
  vep[, is_missense     := grepl(missense_term, Consequence, ignore.case = TRUE)]
  vep[, is_reg_utr      := grepl(paste(reg_utr_terms, collapse = "|"), Consequence, ignore.case = TRUE)]
  vep[, is_other_coding := grepl(paste(other_coding_terms, collapse = "|"), Consequence, ignore.case = TRUE)]
  vep[, is_impact_hm    := IMPACT %in% c("HIGH", "MODERATE")]
  vep[, is_functional_any := is_lof | is_missense | is_reg_utr | is_other_coding | is_impact_hm]

  vep[, sift_prefix     := sub("\\(.*$", "", as.character(SIFT))]
  vep[, polyphen_prefix := sub("\\(.*$", "", as.character(PolyPhen))]
  vep[, sift_del := sift_prefix %in% sift_deleterious_prefixes]
  vep[, pp_dam   := polyphen_prefix %in% polyphen_damaging_prefixes]

  agg <- vep[, .(
    any_lof          = any(is_lof),
    any_missense      = any(is_missense),
    any_reg_utr      = any(is_reg_utr),
    any_other_coding = any(is_other_coding),
    any_functional   = any(is_functional_any),
    sift_worst_del   = any(sift_del),
    pp_worst_dam     = any(pp_dam)
  ), by = SNP]

  # Base tier: worst wins. Default is Reference; "functional-but-unmapped"
  # SNPs stay Reference. `was_unclassified` is kept as a flag for
  # reporting only (see below).
  agg[, tier := "Reference"]
  agg[any_reg_utr == TRUE, tier := "Regulatory_UTR"]
  agg[any_missense == TRUE, tier := "Missense"]
  agg[any_lof == TRUE, tier := "LoF"]
  agg[, was_unclassified := (tier == "Reference" & (any_other_coding == TRUE | any_functional == TRUE))]

  agg[, missense_higher := NA]
  agg[tier == "Missense", missense_higher := (sift_worst_del | pp_worst_dam)]

  # tier_split: Missense split by pathogenicity (primary model)
  agg[, tier_split := tier]
  agg[tier == "Missense" & missense_higher == TRUE,  tier_split := "Missense_Higher"]
  agg[tier == "Missense" & missense_higher == FALSE, tier_split := "Missense_Lower"]

  # tier_pooled: Missense NOT split by pathogenicity (supplementary model)
  # - identical to `tier` itself, kept as an explicit separate column for clarity.
  agg[, tier_pooled := tier]

  agg[, Pathway := pathway]
  agg[]
}

load_univariate <- function(pathway) {
  f <- file.path(
    output_dir,
    paste0("Pathway_", pathway, "_Univariate_SNP_Logistic_Regression_Results_01_Main_LDpruned_All.csv")
  )
  uni <- fread(f, sep = ";", dec = ",", encoding = "UTF-8")
  uni[, SNP := trimws(as.character(SNP))]
  uni[, .(SNP, Estimate_Add, Category)]
}

# ---- Run across all robust pathways, dedup to first-pathway occurrence ----

all_snps  <- list()
seen_snps <- character(0)

for (p in PATHWAYS_FINAL) {
  agg <- classify_pathway_vep(p)
  uni <- load_univariate(p)
  merged <- merge(agg, uni, by = "SNP", all.x = TRUE)
  new_snps <- merged[!(SNP %in% seen_snps)]
  seen_snps <- union(seen_snps, merged$SNP)
  all_snps[[p]] <- new_snps
}

pooled <- rbindlist(all_snps, fill = TRUE)
pooled[, abs_log_or := abs(Estimate_Add)]
pooled[, risk_binary := fifelse(Category == "Risk", 1L, fifelse(Category == "Protective", 0L, NA_integer_))]

cat("Total deduplicated SNPs pooled across", length(PATHWAYS_FINAL), "pathways:", nrow(pooled), "\n")
cat("Tier distribution (tier_split):\n")
print(pooled[, .N, by = tier_split][order(-N)])
cat("\n")

n_reclassified <- sum(pooled$was_unclassified, na.rm = TRUE)
cat("Of the", sum(pooled$tier == "Reference"), "Reference-tier SNPs,", n_reclassified,
    "were reclassified from the former Unclassified_Functional group",
    "(splice_region_variant / stop_lost / protein_altering_variant, not matching",
    "LoF/Missense/Regulatory-UTR specifically).\n\n")

# ============================================================
# Descriptive summary - Reference, Regulatory_UTR, Missense_Higher,
# Missense_Lower, LoF (5 named tiers; Reference includes the
# functional-but-unmapped SNPs, see header), PLUS an extra Missense
# (pooled) row for direct comparison with the Pooled model below.
# ============================================================

exact_binom_ci <- function(x, n) {
  if (n == 0) return(c(NA, NA))
  bt <- binom.test(x, n)
  as.numeric(bt$conf.int)
}

summarize_tier <- function(dt, tier_col) {
  d <- dt[!is.na(get(tier_col)), .(
    N_total          = .N,
    N_with_estimate  = sum(!is.na(abs_log_or)),
    mean_abs_logOR   = mean(abs_log_or, na.rm = TRUE),
    sd_abs_logOR     = sd(abs_log_or, na.rm = TRUE),
    N_with_direction = sum(!is.na(risk_binary)),
    N_risk           = sum(risk_binary == 1, na.rm = TRUE)
  ), by = tier_col]
  setnames(d, tier_col, "Tier")
  d[, risk_proportion := N_risk / N_with_direction]
  ci_list <- mapply(exact_binom_ci, d$N_risk, d$N_with_direction, SIMPLIFY = FALSE)
  d[, risk_CI_lower := sapply(ci_list, `[`, 1)]
  d[, risk_CI_upper := sapply(ci_list, `[`, 2)]
  d[, p_vs_50pct := mapply(function(x, n) if (n == 0) NA_real_ else binom.test(x, n, p = 0.5)$p.value,
                            N_risk, N_with_direction)]
  d
}

desc_split  <- summarize_tier(pooled, "tier_split")
desc_pooled_missense_only <- summarize_tier(pooled, "tier_pooled")[Tier == "Missense"]
desc_pooled_missense_only[, Tier := "Missense_Pooled"]

desc <- rbindlist(list(desc_split, desc_pooled_missense_only), fill = TRUE)

cat("--- Descriptive tier summary (5 named tiers + Missense_Pooled) ---\n")
print(desc)
fwrite(desc, out_descriptive_csv, sep = ";", dec = ",")

# ============================================================
# Formal model family 1 ("Split", primary) - Reference,
# Missense_Higher, Missense_Lower, Regulatory_UTR (LoF excluded)
# ============================================================

model_data_split <- pooled[tier_split %in% c("Reference", "Missense_Higher", "Missense_Lower", "Regulatory_UTR")]
model_data_split[, tier_split := factor(tier_split, levels = c("Reference", "Missense_Higher", "Missense_Lower", "Regulatory_UTR"))]
model_data_split[, Pathway := factor(Pathway)]

cat("\nN SNPs entering the Split formal model:", nrow(model_data_split), "\n")
cat("(excludes LoF: n=", sum(pooled$tier_split == "LoF", na.rm = TRUE), ")\n\n", sep = "")

fit_h6a <- function(dt, tier_var, levels_no_ref) {
  d <- dt[!is.na(abs_log_or)]
  m <- lmerTest::lmer(as.formula(paste("abs_log_or ~", tier_var, "+ (1 | Pathway)")), data = d)
  print(summary(m))
  ci <- confint(m, method = "Wald", parm = "beta_")
  data.table(
    Hypothesis = "Effect-magnitude",
    Tier       = levels_no_ref,
    Estimate   = summary(m)$coefficients[-1, "Estimate"],
    SE         = summary(m)$coefficients[-1, "Std. Error"],
    CI_lower   = ci[-1, 1],
    CI_upper   = ci[-1, 2],
    p_value    = summary(m)$coefficients[-1, "Pr(>|t|)"]
  )
}

fit_h6b <- function(dt, tier_var, levels_no_ref) {
  d <- dt[!is.na(risk_binary)]
  m <- lme4::glmer(as.formula(paste("risk_binary ~", tier_var, "+ (1 | Pathway)")), data = d, family = binomial)
  print(summary(m))
  ci <- confint(m, method = "Wald", parm = "beta_")
  coef <- summary(m)$coefficients
  r <- data.table(
    Hypothesis = "Risk-direction-concentration",
    Tier       = levels_no_ref,
    Estimate   = coef[-1, "Estimate"],
    SE         = coef[-1, "Std. Error"],
    CI_lower   = ci[-1, 1],
    CI_upper   = ci[-1, 2],
    p_value    = coef[-1, "Pr(>|z|)"]
  )
  r[, OR := exp(Estimate)]
  r[, OR_CI_lower := exp(CI_lower)]
  r[, OR_CI_upper := exp(CI_upper)]
  r
}

cat("============================================================\n")
cat("SPLIT MODEL - Effect-magnitude (abs_log_or ~ tier_split + (1|Pathway))\n")
cat("============================================================\n")
h6a_split <- fit_h6a(model_data_split, "tier_split", c("Missense_Higher", "Missense_Lower", "Regulatory_UTR"))

cat("\n============================================================\n")
cat("SPLIT MODEL - Risk-direction-concentration (risk_binary ~ tier_split + (1|Pathway), binomial)\n")
cat("============================================================\n")
h6b_split <- fit_h6b(model_data_split, "tier_split", c("Missense_Higher", "Missense_Lower", "Regulatory_UTR"))

results_split <- rbindlist(list(h6a_split, h6b_split), fill = TRUE)
setcolorder(results_split, c("Hypothesis", "Tier", "Estimate", "SE", "CI_lower", "CI_upper", "p_value"))

cat("\n--- Split model family results (Missense_Higher, Missense_Lower, Regulatory_UTR; q-values below combine with the Pooled family's Missense_Pooled rows) ---\n")
print(results_split)
fwrite(results_split, out_split_csv, sep = ";", dec = ",")

# ============================================================
# Formal model family 2 ("Pooled", supplementary) - Reference,
# Missense (not split), Regulatory_UTR (LoF excluded)
# ============================================================

model_data_pooled <- pooled[tier_pooled %in% c("Reference", "Missense", "Regulatory_UTR")]
model_data_pooled[, tier_pooled := factor(tier_pooled, levels = c("Reference", "Missense", "Regulatory_UTR"))]
model_data_pooled[, Pathway := factor(Pathway)]

cat("\nN SNPs entering the Pooled formal model:", nrow(model_data_pooled), "\n")
cat("(excludes LoF: n=", sum(pooled$tier_pooled == "LoF", na.rm = TRUE), ")\n\n", sep = "")

cat("============================================================\n")
cat("POOLED MODEL - Effect-magnitude (abs_log_or ~ tier_pooled + (1|Pathway))\n")
cat("============================================================\n")
h6a_pooled_full <- fit_h6a(model_data_pooled, "tier_pooled", c("Missense", "Regulatory_UTR"))

cat("\n============================================================\n")
cat("POOLED MODEL - Risk-direction-concentration (risk_binary ~ tier_pooled + (1|Pathway), binomial)\n")
cat("============================================================\n")
h6b_pooled_full <- fit_h6b(model_data_pooled, "tier_pooled", c("Missense", "Regulatory_UTR"))

# Keep only the Missense_Pooled rows from this second fit - its Regulatory_UTR
# rows are the same comparison as in the Split model above (same Reference,
# same Regulatory_UTR SNPs) and would otherwise double-count that tier in the
# combined FDR family below.
results_pooled_full <- rbindlist(list(h6a_pooled_full, h6b_pooled_full), fill = TRUE)
results_pooled_full[Tier == "Missense", Tier := "Missense_Pooled"]
results_pooled <- results_pooled_full[Tier == "Missense_Pooled"]
setcolorder(results_pooled, c("Hypothesis", "Tier", "Estimate", "SE", "CI_lower", "CI_upper", "p_value"))

cat("\n--- Pooled model family results (Missense_Pooled only kept; Regulatory_UTR from this fit dropped as redundant with the Split model's) ---\n")
print(results_pooled)
fwrite(results_pooled, out_pooled_csv, sep = ";", dec = ",")

# ============================================================
# Combined formal test family: 4 tiers (Missense_Higher, Missense_Lower,
# Missense_Pooled, Regulatory_UTR) x 2 hypotheses = 8 tests, one joint
# BH-FDR correction.
# ============================================================

results_combined <- rbindlist(list(results_split, results_pooled), fill = TRUE)
results_combined[, q_value_BH := p.adjust(p_value, method = "BH")]
setcolorder(results_combined, c("Hypothesis", "Tier", "Estimate", "SE", "CI_lower", "CI_upper", "p_value", "q_value_BH"))

cat("\n============================================================\n")
cat("Combined formal test family (8 tests: 4 tiers x 2 hypotheses, ONE joint BH-FDR correction):\n")
cat("============================================================\n")
print(results_combined)
fwrite(results_combined, file.path(output_dir, "Q6_Results_Combined.csv"), sep = ";", dec = ",")

cat("\nDone. Output written to:\n -", out_split_csv, "(pre-FDR)\n -", out_pooled_csv, "(pre-FDR)\n -",
    file.path(output_dir, "Q6_Results_Combined.csv"), "(final, with joint q-values)\n -", out_descriptive_csv, "\n")
cat("\nNote: the descriptive risk-direction vs-50%-null test (per tier, exact binomial) is in the\n")
cat("descriptive summary file (p_vs_50pct column), NOT in the formal BH-FDR family above.\n")
cat("Note: Missense_Higher, Missense_Lower and Missense_Pooled come from TWO separate model\n")
cat("fits (Missense_Pooled cannot share a single categorical factor with Missense_Higher/\n")
cat("Missense_Lower, since it is their union, not a disjoint category) - only their p-values\n")
cat("are combined into one joint BH-FDR family.\n")
