# ==============================================================================
# 14a_snp_level_cell_counts.R
# SNP-level functional-concentration analysis, step 1: pooled, deduplicated
# per-tier SNP cell counts across the robust candidate pathways.
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================

# ============================================================
# Pooled, deduplicated per-tier SNP cell counts
# ============================================================
#
# Purpose: before fitting the mixed-effects models (14b), check whether
# the 4-tier VEP-based classification (loss-of-function / missense-higher /
# missense-lower / regulatory-UTR, vs. an intronic/synonymous reference),
# pooled and deduplicated across the robust candidate pathways, yields
# cell sizes large enough to fit those models.
#
# This script only reads already-computed, aggregate per-SNP files
# (VEP annotation text files; univariate SNP regression result CSVs).
# It does NOT touch individual-level genotype or phenotype data, and
# does NOT fit any statistical model - it only counts.
#
# Design choices (see Supplementary Methods 2):
#   - Tier priority when a SNP qualifies for >1 tier: worst wins,
#     LoF > Missense > Regulatory/UTR > Reference.
#   - SIFT / PolyPhen: any-transcript, worst prediction wins (same
#     "union across transcripts" logic already used for the existing
#     Functional-VEP-filtered / LoF-burden sensitivity scores).
#   - SNPs occurring in >1 pathway are deduplicated to their FIRST
#     pathway occurrence (iteration order = Robust_Candidate order in
#     AllPathways_Results_Tracking_Table.csv, i.e. the same order as
#     Table 2).
#
# Two classification questions are reported explicitly in the console
# output / output files, because they affect the tier cell counts:
#
#   (A) Three VEP consequence terms that count as "functional" under
#       the pipeline's existing binary functional/non-functional rule
#       (load_functional_vep_snps()) - protein_altering_variant,
#       stop_lost, splice_region_variant - do not map cleanly onto
#       any of the four named tiers (they are not loss-of-function
#       per the strict LoF term list, not "missense_variant", and not
#       regulatory/UTR). Also: any SNP that is "functional" purely via
#       IMPACT %in% c("HIGH","MODERATE") without matching any named
#       consequence term at all would fall in the same situation.
#       These are pooled into a separate "Unclassified_Functional"
#       tier here (NOT folded into Reference, NOT force-assigned to
#       Missense or Regulatory_UTR) so their size and raw Consequence
#       values can be reviewed.
#
#   (B) The missense higher-/lower-pathogenicity split uses SIFT
#       (deleterious vs. tolerated) AND PolyPhen (probably/possibly
#       damaging vs. benign). This script reports both an OR-rule
#       (higher pathogenicity if SIFT deleterious OR PolyPhen damaging)
#       and an AND-rule (both must agree).
#
# Output:
#   - Console: tier counts, missense split (both rules), tier x
#     Category (Risk/Protective) cross-tab, Unclassified_Functional
#     raw-Consequence breakdown.
#   - Q6_Diagnostic_Cell_Counts.csv
#   - Q6_Diagnostic_Cell_Counts_by_Category.csv
#   - Q6_Diagnostic_Unclassified_Functional_Consequences.csv
#
# ============================================================


library(data.table)
library(dplyr)

# ---- Paths ----

vep_dir         <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Outputs_Pathway")   # VEP annotation files (read only)
output_dir      <- file.path(Sys.getenv("FOCUS_RESULTS_DIR", unset = "path/to/results"), "OutPW")       # tracking table, univariate results, outputs
tracking_table  <- file.path(output_dir, "AllPathways_Results_Tracking_Table.csv")

out_tier_csv           <- file.path(output_dir, "Q6_Diagnostic_Cell_Counts.csv")
out_tier_category_csv  <- file.path(output_dir, "Q6_Diagnostic_Cell_Counts_by_Category.csv")
out_unclassified_csv   <- file.path(output_dir, "Q6_Diagnostic_Unclassified_Functional_Consequences.csv")

# ---- Pathway list: same source/order as the permutation test (11) ----

tracking <- fread(tracking_table, sep = ";", dec = ",", encoding = "UTF-8")
PATHWAYS_FINAL <- tracking[Robust_Candidate == TRUE, ]$Pathway
cat("N robust candidate pathways:", length(PATHWAYS_FINAL), "\n\n")

# ---- Term lists (identical to load_functional_vep_snps() / load_lof_vep_annotation()
#      in 07_pathway_pipeline.R) ----

lof_terms <- c(
  "transcript_ablation",
  "splice_acceptor_variant",
  "splice_donor_variant",
  "stop_gained",
  "frameshift_variant",
  "start_lost"
)

missense_term <- "missense_variant"

reg_utr_terms <- c(
  "regulatory_region_variant",
  "TF_binding_site_variant",
  "TFBS_ablation",
  "TFBS_amplification",
  "5_prime_UTR_variant",
  "3_prime_UTR_variant",
  "mature_miRNA_variant"
)

# Terms that are "functional" under the existing binary rule but do not
# map onto a named tier - see note (A) above. Kept separate on purpose.
other_coding_terms <- c(
  "protein_altering_variant",
  "stop_lost",
  "splice_region_variant"
)

sift_deleterious_prefixes  <- c("deleterious", "deleterious_low_confidence")
sift_tolerated_prefixes    <- c("tolerated", "tolerated_low_confidence")
polyphen_damaging_prefixes <- c("probably_damaging", "possibly_damaging")
polyphen_benign_prefixes   <- c("benign")
# "unknown" and "-" (missing) count toward neither verdict for either predictor.

# ---- Per-pathway VEP classification, any-transcript / worst-wins ----

classify_pathway_vep <- function(pathway) {

  vep_file <- file.path(vep_dir, paste0("Pathway_", pathway, "_VEP.txt"))
  vep <- fread(vep_file, encoding = "UTF-8")

  vep[, SNP := trimws(as.character(`#Uploaded_variation`))]
  vep <- vep[!is.na(SNP) & SNP != ""]

  vep[, is_lof          := grepl(paste(lof_terms, collapse = "|"), Consequence, ignore.case = TRUE)]
  vep[, is_missense      := grepl(missense_term, Consequence, ignore.case = TRUE)]
  vep[, is_reg_utr       := grepl(paste(reg_utr_terms, collapse = "|"), Consequence, ignore.case = TRUE)]
  vep[, is_other_coding  := grepl(paste(other_coding_terms, collapse = "|"), Consequence, ignore.case = TRUE)]
  vep[, is_impact_hm     := IMPACT %in% c("HIGH", "MODERATE")]
  vep[, is_functional_any:= is_lof | is_missense | is_reg_utr | is_other_coding | is_impact_hm]

  vep[, sift_prefix     := sub("\\(.*$", "", as.character(SIFT))]
  vep[, polyphen_prefix := sub("\\(.*$", "", as.character(PolyPhen))]
  vep[, sift_del := sift_prefix %in% sift_deleterious_prefixes]
  vep[, sift_tol := sift_prefix %in% sift_tolerated_prefixes]
  vep[, pp_dam   := polyphen_prefix %in% polyphen_damaging_prefixes]
  vep[, pp_ben   := polyphen_prefix %in% polyphen_benign_prefixes]

  # raw Consequence strings from rows that are functional ONLY via
  # other_coding terms or impact-only (i.e. not already lof/missense/reg_utr)
  vep[, flag_unclassified_row := is_functional_any & !is_lof & !is_missense & !is_reg_utr]

  agg <- vep[, .(
    any_lof           = any(is_lof),
    any_missense       = any(is_missense),
    any_reg_utr       = any(is_reg_utr),
    any_other_coding  = any(is_other_coding),
    any_functional    = any(is_functional_any),
    sift_worst_del    = any(sift_del),
    pp_worst_dam      = any(pp_dam),
    unclassified_consequences = paste(sort(unique(Consequence[flag_unclassified_row])), collapse = ";")
  ), by = SNP]

  # Tier assignment: worst wins. Order matters - later assignments overwrite earlier ones.
  agg[, tier := "Reference"]
  agg[any_reg_utr == TRUE, tier := "Regulatory_UTR"]
  agg[any_missense == TRUE, tier := "Missense"]
  agg[any_lof == TRUE, tier := "LoF"]
  # Functional-but-unmapped SNPs (note A): only applies where nothing else claimed the SNP.
  agg[tier == "Reference" & (any_other_coding == TRUE | any_functional == TRUE), tier := "Unclassified_Functional"]

  # Missense pathogenicity split - both combination rules (note B)
  agg[, missense_higher_OR  := NA]
  agg[tier == "Missense", missense_higher_OR  := (sift_worst_del | pp_worst_dam)]
  agg[, missense_higher_AND := NA]
  agg[tier == "Missense", missense_higher_AND := (sift_worst_del & pp_worst_dam)]

  agg[, Pathway := pathway]
  agg[]
}

# ---- Univariate per-SNP results (Estimate_Add, Category) ----

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
t0 <- Sys.time()

for (p in PATHWAYS_FINAL) {

  agg <- classify_pathway_vep(p)
  uni <- load_univariate(p)

  merged <- merge(agg, uni, by = "SNP", all.x = TRUE)
  n_missing_uni <- sum(is.na(merged$Category))
  if (n_missing_uni > 0) {
    cat("  NOTE (", p, "): ", n_missing_uni,
        " VEP SNP(s) not found in the univariate results file.\n", sep = "")
  }

  new_snps <- merged[!(SNP %in% seen_snps)]
  seen_snps <- union(seen_snps, merged$SNP)

  cat(sprintf("%-40s VEP entries: %6d | unique SNPs: %5d | new after dedup: %5d\n",
              p, nrow(agg), length(unique(agg$SNP)), nrow(new_snps)))

  all_snps[[p]] <- new_snps
}

pooled <- rbindlist(all_snps, fill = TRUE)

cat("\n============================================================\n")
cat("Total deduplicated SNPs pooled across", length(PATHWAYS_FINAL), "pathways:", nrow(pooled), "\n")
cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
cat("============================================================\n\n")

# ---- 1) Primary tier cell counts ----

tier_counts <- pooled[, .N, by = tier][order(-N)]
cat("--- Tier cell counts (pooled, deduplicated) ---\n")
print(tier_counts)

# ---- 2) Missense pathogenicity split, both combination rules ----

cat("\n--- Missense higher/lower split - OR rule (SIFT deleterious OR PolyPhen damaging) ---\n")
print(pooled[tier == "Missense", .N, by = missense_higher_OR])

cat("\n--- Missense higher/lower split - AND rule (SIFT deleterious AND PolyPhen damaging) ---\n")
print(pooled[tier == "Missense", .N, by = missense_higher_AND])

# ---- 3) Tier x Category (Risk/Protective) cross-tab ----

tier_category <- pooled[, .N, by = .(tier, Category)][order(tier, Category)]
cat("\n--- Tier x Category (Risk/Protective) cross-tab ---\n")
print(tier_category)

# ---- 4) Unclassified_Functional: what is actually landing there? ----

unclassified <- pooled[tier == "Unclassified_Functional"]
cat("\n--- Unclassified_Functional: ", nrow(unclassified), " SNPs. Raw Consequence terms found: ---\n", sep = "")
if (nrow(unclassified) > 0) {
  unclassified_terms <- sort(table(unlist(strsplit(unclassified$unclassified_consequences, ";"))), decreasing = TRUE)
  print(unclassified_terms)
} else {
  unclassified_terms <- table(character(0))
  cat("(none)\n")
}

# ---- Write outputs ----

fwrite(tier_counts, out_tier_csv, sep = ";", dec = ",")
fwrite(tier_category, out_tier_category_csv, sep = ";", dec = ",")
fwrite(as.data.frame(unclassified_terms), out_unclassified_csv, sep = ";", dec = ",")

cat("\nDone. Output written to:\n -", out_tier_csv, "\n -", out_tier_category_csv, "\n -", out_unclassified_csv, "\n")
cat("\nNo model was fit. This script only checks cell sizes.\n")
