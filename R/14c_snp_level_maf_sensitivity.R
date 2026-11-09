# ==============================================================================
# 14c_snp_level_maf_sensitivity.R
# SNP-level functional-concentration analysis, step 3: post hoc MAF adjustment
# of the effect-magnitude model (Supplementary Table S10).
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================

# ============================================================
# Question: missense SNPs are rarer than reference SNPs. The sampling
# variance of a univariate log OR grows as MAF falls, so |log OR| is
# larger for rare SNPs even under the null. Does the effect-magnitude
# difference persist after adjustment for MAF?
#
# Pathway list, VEP tier classification, deduplication rule and
# univariate input are identical to 14b_snp_level_models.R, so the
# unadjusted models (a) reproduce the effect-magnitude results of 14b.
#
# Models (all with pathway random intercept, as in 14b):
#   (a) abs_log_or ~ tier            (replication of 14b)
#   (b) abs_log_or ~ tier + log(MAF) (MAF-adjusted, post hoc)
#   (c) abs_log_or ~ tier + log(SE)  (SE-adjusted, alternative)
#   (d) abs_z      ~ tier            (standardized effect |beta|/SE)
# each for the split and the pooled tier coding.
#
# MAF: from the training-set genotype counts (N_Geno0/1/2) in the
# univariate result files (observed, non-imputed genotypes).
#
# Outputs:
#   Q6_MAF_Sensitivity_Results.csv
#   Q6_MAF_Sensitivity_Tier_MAF_Summary.csv
#   Q6_SNP_Level_Pooled_with_MAF.csv   (per-SNP table, all tiers)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
})

# ---- Paths (identical to 14b) ----
vep_dir         <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Outputs_Pathway")   # VEP annotation files (read only)
output_dir      <- file.path(Sys.getenv("FOCUS_RESULTS_DIR", unset = "path/to/results"), "OutPW")       # tracking table, univariate results, outputs
tracking_table <- file.path(output_dir, "AllPathways_Results_Tracking_Table.csv")

out_results <- file.path(output_dir, "Q6_MAF_Sensitivity_Results.csv")
out_mafsum  <- file.path(output_dir, "Q6_MAF_Sensitivity_Tier_MAF_Summary.csv")
out_snps    <- file.path(output_dir, "Q6_SNP_Level_Pooled_with_MAF.csv")

# ---- 1. Pathway list (identical to 14b) ----
tracking <- fread(tracking_table, sep = ";", dec = ",", encoding = "UTF-8")
PATHWAYS_FINAL <- tracking[Robust_Candidate == TRUE, ]$Pathway
cat("N robust candidate pathways:", length(PATHWAYS_FINAL), "\n\n")

# ---- 2. VEP classification (identical to 14b) ----
lof_terms <- c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
               "stop_gained", "frameshift_variant", "start_lost")
missense_term <- "missense_variant"
reg_utr_terms <- c("regulatory_region_variant", "TF_binding_site_variant", "TFBS_ablation",
                   "TFBS_amplification", "5_prime_UTR_variant", "3_prime_UTR_variant",
                   "mature_miRNA_variant")
other_coding_terms <- c("protein_altering_variant", "stop_lost", "splice_region_variant")
sift_deleterious_prefixes  <- c("deleterious", "deleterious_low_confidence")
polyphen_damaging_prefixes <- c("probably_damaging", "possibly_damaging")

classify_pathway_vep <- function(pathway) {
  vep <- fread(file.path(vep_dir, paste0("Pathway_", pathway, "_VEP.txt")), encoding = "UTF-8")
  vep[, SNP := trimws(as.character(`#Uploaded_variation`))]
  vep <- vep[!is.na(SNP) & SNP != ""]
  vep[, is_lof      := grepl(paste(lof_terms, collapse = "|"), Consequence, ignore.case = TRUE)]
  vep[, is_missense := grepl(missense_term, Consequence, ignore.case = TRUE)]
  vep[, is_reg_utr  := grepl(paste(reg_utr_terms, collapse = "|"), Consequence, ignore.case = TRUE)]
  vep[, sift_del := sub("\\(.*$", "", as.character(SIFT))     %in% sift_deleterious_prefixes]
  vep[, pp_dam   := sub("\\(.*$", "", as.character(PolyPhen)) %in% polyphen_damaging_prefixes]

  agg <- vep[, .(
    any_lof        = any(is_lof),
    any_missense   = any(is_missense),
    any_reg_utr    = any(is_reg_utr),
    sift_worst_del = any(sift_del),
    pp_worst_dam   = any(pp_dam),
    Gene           = paste(setdiff(unique(SYMBOL), "-"), collapse = ","),
    SIFT_all       = paste(setdiff(unique(as.character(SIFT)), "-"), collapse = ";"),
    PolyPhen_all   = paste(setdiff(unique(as.character(PolyPhen)), "-"), collapse = ";")
  ), by = SNP]

  agg[, tier := "Reference"]
  agg[any_reg_utr == TRUE,  tier := "Regulatory_UTR"]
  agg[any_missense == TRUE, tier := "Missense"]
  agg[any_lof == TRUE,      tier := "LoF"]
  agg[, tier_split := tier]
  agg[tier == "Missense" & (sift_worst_del | pp_worst_dam),    tier_split := "Missense_Higher"]
  agg[tier == "Missense" & !(sift_worst_del | pp_worst_dam),   tier_split := "Missense_Lower"]
  agg[, tier_pooled := tier]
  agg[, Pathway := pathway]
  agg[]
}

# Univariate results: same file as 14b, plus SE and genotype counts
load_univariate <- function(pathway) {
  f <- file.path(output_dir, paste0("Pathway_", pathway,
                 "_Univariate_SNP_Logistic_Regression_Results_01_Main_LDpruned_All.csv"))
  uni <- fread(f, sep = ";", dec = ",", encoding = "UTF-8")
  uni[, SNP := trimws(as.character(SNP))]
  uni[, .(SNP, Estimate_Add, StdError_Add, OR_Add, CI_lower_Add, CI_upper_Add, Pvalue_Add,
          Allele1, Allele2, Category, N_Geno0, N_Geno1, N_Geno2)]
}

# ---- 3. Pool + dedup to first-pathway occurrence (identical to 14b) ----
all_snps <- list(); seen_snps <- character(0)
for (p in PATHWAYS_FINAL) {
  merged <- merge(classify_pathway_vep(p), load_univariate(p), by = "SNP", all.x = TRUE)
  all_snps[[p]] <- merged[!(SNP %in% seen_snps)]
  seen_snps <- union(seen_snps, merged$SNP)
}
pooled <- rbindlist(all_snps, fill = TRUE)
pooled[, abs_log_or := abs(Estimate_Add)]
pooled[, abs_z      := abs(Estimate_Add) / StdError_Add]
# MAF from observed genotype counts; Allele2 = counted allele (N_Geno2 = homozygous counted)
pooled[, AF  := (2 * N_Geno2 + N_Geno1) / (2 * (N_Geno0 + N_Geno1 + N_Geno2))]
pooled[, MAF := pmin(AF, 1 - AF)]
pooled[, log_MAF := log(MAF)]
pooled[, log_SE  := log(StdError_Add)]

cat("Tier distribution (expected 4479 / 818 / 121 / 118 / 10):\n")
print(pooled[, .N, by = tier_split][order(-N)])
cat("SNPs with MAF == 0 or missing:", pooled[is.na(MAF) | MAF == 0, .N], "\n\n")

maf_sum <- pooled[, .(N = .N,
                      MAF_median = median(MAF, na.rm = TRUE),
                      MAF_Q1 = quantile(MAF, 0.25, na.rm = TRUE),
                      MAF_Q3 = quantile(MAF, 0.75, na.rm = TRUE),
                      SE_median = median(StdError_Add, na.rm = TRUE),
                      abs_log_or_mean = mean(abs_log_or, na.rm = TRUE),
                      abs_z_mean = mean(abs_z, na.rm = TRUE),
                      N_nominal_p05 = sum(Pvalue_Add < 0.05, na.rm = TRUE)),
                  by = tier_split]
print(maf_sum)
fwrite(maf_sum, out_mafsum, sep = ";", dec = ",")

# ---- 4. Models ----
fit_lmm <- function(dt, formula_str, tier_var, levels_no_ref, label) {
  d <- dt[is.finite(log_MAF) & is.finite(log_SE) & !is.na(abs_log_or)]
  m <- lmerTest::lmer(as.formula(formula_str), data = d)
  co <- summary(m)$coefficients
  rows <- paste0(tier_var, levels_no_ref)
  ci <- confint(m, parm = rows, method = "Wald")
  data.table(Model = label, Formula = formula_str, Tier = levels_no_ref,
             Estimate = co[rows, "Estimate"], SE = co[rows, "Std. Error"],
             CI_lower = ci[rows, 1], CI_upper = ci[rows, 2],
             p_value = co[rows, "Pr(>|t|)"], N = nrow(d), Singular = isSingular(m))
}

d_split <- pooled[tier_split %in% c("Reference", "Missense_Higher", "Missense_Lower", "Regulatory_UTR")]
d_split[, tier_split := factor(tier_split, levels = c("Reference", "Missense_Higher", "Missense_Lower", "Regulatory_UTR"))]
d_pool  <- pooled[tier_pooled %in% c("Reference", "Missense", "Regulatory_UTR")]
d_pool[, tier_pooled := factor(tier_pooled, levels = c("Reference", "Missense", "Regulatory_UTR"))]

lv_s <- c("Missense_Higher", "Missense_Lower", "Regulatory_UTR")
lv_p <- c("Missense", "Regulatory_UTR")

res <- rbindlist(list(
  fit_lmm(d_split, "abs_log_or ~ tier_split + (1 | Pathway)",           "tier_split", lv_s, "Split (a) unadjusted"),
  fit_lmm(d_split, "abs_log_or ~ tier_split + log_MAF + (1 | Pathway)", "tier_split", lv_s, "Split (b) MAF-adjusted"),
  fit_lmm(d_split, "abs_log_or ~ tier_split + log_SE + (1 | Pathway)",  "tier_split", lv_s, "Split (c) SE-adjusted"),
  fit_lmm(d_split, "abs_z ~ tier_split + (1 | Pathway)",                "tier_split", lv_s, "Split (d) |z|"),
  fit_lmm(d_pool,  "abs_log_or ~ tier_pooled + (1 | Pathway)",           "tier_pooled", lv_p, "Pooled (a) unadjusted"),
  fit_lmm(d_pool,  "abs_log_or ~ tier_pooled + log_MAF + (1 | Pathway)", "tier_pooled", lv_p, "Pooled (b) MAF-adjusted"),
  fit_lmm(d_pool,  "abs_log_or ~ tier_pooled + log_SE + (1 | Pathway)",  "tier_pooled", lv_p, "Pooled (c) SE-adjusted"),
  fit_lmm(d_pool,  "abs_z ~ tier_pooled + (1 | Pathway)",                "tier_pooled", lv_p, "Pooled (d) |z|")
))
print(res[, .(Model, Tier, Estimate = round(Estimate, 4), CI_lower = round(CI_lower, 3),
              CI_upper = round(CI_upper, 3), p_value = signif(p_value, 3), Singular)])
fwrite(res, out_results, sep = ";", dec = ",")

# ---- 5. Per-SNP table (for reference; not a manuscript table) ----
setorder(pooled, tier_split, -abs_log_or)
fwrite(pooled[, .(SNP, Pathway, tier_split, Gene, SIFT_all, PolyPhen_all,
                  Counted_allele = Allele2, Other_allele = Allele1, MAF,
                  Estimate_Add, StdError_Add, OR_Add, CI_lower_Add, CI_upper_Add,
                  Pvalue_Add, Category)],
       out_snps, sep = ";", dec = ",")

cat("\nDone. Written:\n ", out_results, "\n ", out_mafsum, "\n ", out_snps, "\n")

