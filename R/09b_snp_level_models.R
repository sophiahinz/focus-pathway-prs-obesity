# Pathway-restricted polygenic scores for energy-substrate metabolism and obesity (FoCus cohort)
# Script: 09b_snp_level_models.R | Source file: Protocol_SNP_Analysis_FoCus_Pathway_Q6_H6a_H6b_Model_v2_REPAIR.R
# Local paths replaced by the environment variable FOCUS_BASE_DIR (project root); analysis code unchanged.
# Individual-level FoCus data are not included (see README).

# ============================================================
# REPAIR COPY (created 2026-11-03) of Notizen_und_Protokolle/Protocol_SNP_Analysis_FoCus_Pathway_Q6_H6a_H6b_Model_v2.R
# Original in Notizen_und_Protokolle/ remains unchanged (read-only).
#
# Reason: after the phenotype-data correction (Repair log, Section 1)
# the robust-candidate set changed from 26 to 19 pathways. The SNP-level
# analysis must therefore be pooled over the 19 REPAIR candidates.
#
# Changes vs. original (all marked "# REPAIR"):
#   (1) output_dir -> Repair/OutPW: tracking table (Robust_Candidate,
#       19 pathways), univariate SNP results and ALL outputs.
#       Output file names are unchanged, so nothing in the original
#       Outputs_Pathway is overwritten.
#   (2) VEP annotation read read-only from the original Outputs_Pathway
#       via orig_output_dir (not regenerated; Repair log, Section 2.1).
#   (3) Guard: stop if the tracking table does not yield 19 candidates.
# Analysis code is otherwise identical.
#
# Note: the univariate SNP result files are identical in old and REPAIR
# runs (checked for Citrate_Cycle 2026-11-03: byte-identical, N = 1,073;
# univariate models do not use the covariates). Differences to the old
# results therefore reflect only the changed pathway pool (26 -> 19).
# Expected values quoted in the comments below (e.g. tier counts, +0.047,
# p = 0.63) refer to the OLD 26-pathway run and will NOT be reproduced.
# ============================================================

# ============================================================
# Q6: SNP-Level Functional Concentration Analysis
# Full effect-magnitude / risk-direction-concentration mixed-effects
# models, pooled across the 26 robust candidate pathways.
#
# v2 (2026-09-02) - Changes vs. the original Q6_H6a_H6b_Model.R,
# per discussion/decision on 2026-09-02:
#   (1) Unclassified_Functional (splice_region_variant / stop_lost /
#       protein_altering_variant, not matching LoF/Missense/Reg-UTR
#       specifically) is now folded into Reference rather than kept
#       as its own excluded tier. Rationale: VEP's own IMPACT ranking
#       classifies splice_region_variant as LOW (same tier as
#       synonymous_variant, already in Reference); these SNPs also
#       co-occur exclusively with intronic/synonymous/non-coding-
#       transcript annotations; and empirically (old run) this group
#       had the LOWEST mean |log OR| of all six tiers (0.117, lower
#       even than Reference's 0.170) - the opposite of what "functional"
#       would predict. A diagnostic count of how many Reference SNPs
#       were reclassified from Unclassified_Functional is still
#       reported (console + descriptive summary) for transparency.
#   (2) Missense is now analyzed TWO ways, run as two separate models:
#         - "Split"  : Missense_Higher / Missense_Lower (SIFT/PolyPhen
#           OR-rule pathogenicity split) vs. Reference, Regulatory_UTR.
#         - "Pooled" : Missense (not split by pathogenicity) vs.
#           Reference, Regulatory_UTR.
#       Fit as two SEPARATE model objects (Missense_Pooled cannot share a
#       single categorical factor with Missense_Higher/Missense_Lower,
#       since it is their union, not a disjoint category) - but their
#       p-values are combined into ONE joint test family: 4 tiers
#       (Missense_Higher, Missense_Lower, Missense_Pooled, Regulatory_UTR)
#       x 2 hypotheses = 8 tests, ONE joint BH-FDR correction (decision
#       2026-09-02). The Pooled model's own Regulatory_UTR estimate is
#       dropped as redundant (same comparison already in the Split model).
#       Rationale for adding the pooled version: missense-higher (n=121)
#       and missense-lower (n=118) point in the same direction for H6a
#       (+0.047 / +0.032 in the original run) - pooling to n=239 may gain
#       power without diluting the signal (unlike pooling with the much
#       larger, near-null Regulatory_UTR group, which was checked and
#       rejected - see manuscript discussion).
#   LoF (n=10) remains excluded from both formal model families and
#   reported descriptively only (unchanged from v1).
#
# Both models (both versions) use "Reference" (intronic/synonymous,
# now incl. former Unclassified_Functional) as the factor baseline.
#
# Design decisions carried over from v1 (see manuscript / Methods-Log
# Section 2.11 for full detail and rationale):
#   - Tier priority (multi-tier SNPs): LoF > Missense > Reg/UTR > Reference.
#   - SIFT/PolyPhen: any-transcript, worst-prediction-wins per predictor,
#     combined via OR (either sufficient) for missense-higher/-lower.
#   - Dedup: SNPs occurring in >1 pathway kept at first occurrence only
#     (iteration order = Robust_Candidate order in the tracking table).
#   - H6b's "vs. 50% null" comparison is NOT part of either formal family;
#     reported descriptively (proportion + 95% CI per tier) only.
#
# Only already-computed, aggregate per-SNP files are read (VEP
# annotation, univariate regression results). No individual-level
# genotype/phenotype data is touched, and no new score is built or
# reoriented.
#
# Output:
#   - Q6_Results_Split.csv       (6 formal tests: Split model family)
#   - Q6_Results_Pooled.csv      (4 formal tests: Pooled model family)
#   - Q6_Descriptive_Tier_Summary_v2.csv (all tiers incl. LoF and a
#     Missense-pooled row; N, mean|logOR|, risk proportion, 95% CI)
#   - Console: full model summaries (lmer/glmer) for both model
#     families, plus the Unclassified->Reference reclassification count.
#
# ============================================================

library(data.table)
library(dplyr)
library(lme4)
if (!requireNamespace("lmerTest", quietly = TRUE)) {
  install.packages("lmerTest")
}
library(lmerTest)   # adds Satterthwaite p-values to lmer() output

# ---- Paths (same convention as Phase4_Permutation_Test.R) ----

orig_output_dir <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Outputs_Pathway")       # REPAIR: read-only (VEP annotation only)
output_dir      <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Repair/OutPW")           # REPAIR: tracking table (19), univariate results, outputs
tracking_table <- file.path(output_dir, "AllPathways_Results_Tracking_Table.csv")

out_split_csv       <- file.path(output_dir, "Q6_Results_Split.csv")
out_pooled_csv       <- file.path(output_dir, "Q6_Results_Pooled.csv")
out_descriptive_csv <- file.path(output_dir, "Q6_Descriptive_Tier_Summary_v2.csv")

# ---- Pathway list: same source/order as the permutation test (Q1) ----

tracking <- fread(tracking_table, sep = ";", dec = ",", encoding = "UTF-8")
PATHWAYS_FINAL <- tracking[Robust_Candidate == TRUE, ]$Pathway
if (length(PATHWAYS_FINAL) != 19) stop("REPAIR: expected 19 robust candidates, got ", length(PATHWAYS_FINAL))  # REPAIR
cat("N robust candidate pathways:", length(PATHWAYS_FINAL), "\n\n")

# ---- Term lists (identical to load_functional_vep_snps() / load_lof_vep_annotation()
#      in Protocol_SNP_Analysis_FoCus_Pathway_Template_Split_19082026.R) ----

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

# ---- Per-pathway VEP classification (identical logic to v1, except the
#      Unclassified_Functional -> Reference fold, see header) ----

classify_pathway_vep <- function(pathway) {

  vep_file <- file.path(orig_output_dir, paste0("Pathway_", pathway, "_VEP.txt"))  # REPAIR
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

  # Base tier: worst wins. Default is Reference, and - unlike v1 -
  # NOTHING further reassigns "functional-but-unmapped" SNPs away from
  # Reference; they simply stay Reference. `was_unclassified` is kept
  # purely as a flag for transparency reporting (see below).
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

# ---- Run across all 26 pathways, dedup to first-pathway occurrence ----

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
# Missense_Lower, LoF (5 named tiers; Reference now includes the
# former Unclassified_Functional group), PLUS an extra Missense
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
# Missense_Pooled, Regulatory_UTR) x 2 hypotheses = 8 tests, ONE joint
# BH-FDR correction (decision 2026-09-02, supersedes the two-separate-
# families approach from the previous version of this script).
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
cat("\nNote: H6b's descriptive vs-50%-null test (per tier, exact binomial) is in the\n")
cat("descriptive summary file (p_vs_50pct column), NOT in the formal BH-FDR family above.\n")
cat("Note: Missense_Higher, Missense_Lower and Missense_Pooled come from TWO separate model\n")
cat("fits (Missense_Pooled cannot share a single categorical factor with Missense_Higher/\n")
cat("Missense_Lower, since it is their union, not a disjoint category) - only their p-values\n")
cat("are combined into one joint BH-FDR family, per your 2026-09-02 decision.\n")
