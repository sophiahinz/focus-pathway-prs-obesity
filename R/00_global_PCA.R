# Pathway-restricted polygenic scores for energy-substrate metabolism and obesity (FoCus cohort)
# Script: 00_global_PCA.R | Source file: Protocol_SNP_Analysis_FoCus_V23b_06082026_PCs_Multivariat_FIX2.R (lines 1-34, 52-962: genome-wide PCA block only)
# Local paths replaced by the environment variable FOCUS_BASE_DIR (project root); analysis code unchanged.
# Individual-level FoCus data are not included (see README).

#======================================================================================
# INSTALL & LOAD PACKAGES
#====================================================================================
#Clean up R: 
rm(list = ls())     # remove all objects
gc()                # free memory

#(Update R if necessary)
# install.packages("installr")   # not run here; install packages beforehand
# library(installr)   # not run here; install packages beforehand
# updateR()   # not run here; install packages beforehand

#BigSnpr
# install.packages("devtools")   # not run here; install packages beforehand
library(devtools)
# devtools::install_github("privefl/bigsnpr")   # not run here; install packages beforehand
library(bigsnpr)

#snpStats
# install.packages("BiocManager")   # not run here; install packages beforehand
# BiocManager::install("snpStats")   # not run here; install packages beforehand
library(snpStats)

library(broom)
library(ggplot2)
library(car)
library(dplyr)
library(pROC)
library(data.table)

#===============================================================================
# Paths
#===============================================================================
base_dir <- Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus")
# (gene-loop paths of the source script omitted)

#==================================================================================
# GENOME-WIDE PCA (POPULATION STRATIFICATION) - ONE-TIME COMPUTATION
#==================================================================================
# Purpose: genome-wide principal components (PCs) to control for population
# stratification in the multivariable models. Computed from the COMPLETE
# GSA array (Ambu + MD, all SNPs) - INDEPENDENT of individual genes/pathways.
# Computed ONCE here and saved as a separate file so that the
# pathway script can reload the same PCs unchanged (no re-computation,
# no pathway-local PCA).
#
# Workflow:
#   1) Merge Ambu + MD (raw, complete GSA files)
#   2) QC: autosomes only, --maf 0.01, --geno 0.05 (SNP missingness),
#      --mind 0.1 (individual missingness)
#   3) Exclusion of known long-range LD regions (HLA, chr8 inversion etc.)
#   4) LD-Pruning: --indep-pairwise 50 5 0.2
#   5) KING-robust kinship: unrelated core set (--king-cutoff 0.0884,
#      corresponds to third-degree relatedness) - serves here as a sensitivity
#      check of whether sporadic relatedness distorts the PCs
#   6) PCA (3 components) on the FULL pruned QC dataset (all
#      individuals) -> these are the PCs used in the models
#   7) Additional PCA on the unrelated core set -> sensitivity comparison:
#      correlation of PC1-3 between both variants for the shared
#      individuals. If close to 1, sporadic relatedness does not relevantly distort the
#      PCs, and the PCs from step 6 can be used unchanged for
#      ALL individuals (incl. the few possibly related ones).
#
# IMPORTANT: the long-range LD regions below are standard coordinates
# for genome build hg19/GRCh37 (see e.g. Anderson et al. 2010, Nat Protoc,
# or the widely used "high-LD-regions.txt"). If the GSA calling was done on
# hg38, these coordinates must be adapted - please check the
# build of the bim files (FO_GSA_Ambu.bim / FO_GSA_MD.bim) before the production
# run.

genetik_raw_dir <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Genetik_Daten_FoCus")

genomewide_pca_dir <- file.path(base_dir, "Global_PCA")
dir.create(genomewide_pca_dir, showWarnings = FALSE, recursive = TRUE)

genomewide_pcs_file <- file.path(genomewide_pca_dir, "Global_Genomewide_PCs.rds")

# PLINK 1.9 (for merge/QC/pruning/KING cutoff) - as used elsewhere in the
# pathway script.
plink <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Working_Plink/plink.exe")

if (!file.exists(plink)) {
  stop("PLINK not found. Please provide the full path, e.g. C:/plink/plink.exe")
}

# PLINK2 is required for the projection step (SNP loadings + --score),
# because the flags used (allele-wts, --score with header-read /
# variance-standardize / --score-col-nums) are PLINK2-specific. If your
# plink.exe is PLINK 1.9, please obtain plink2.exe in addition
# and adapt the path below.
plink2 <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Working_Plink/plink2.exe")

if (!file.exists(plink2)) {
  stop(
    "PLINK2 not found (", plink2, "). For the PC projection ",
    "PLINK2 is required (not only PLINK 1.9). Please obtain plink2.exe ",
    "and adapt the path in 'plink2'."
  )
}

raw_ambu_prefix <- file.path(genetik_raw_dir, "FO_GSA_Ambu")
raw_md_prefix   <- file.path(genetik_raw_dir, "FO_GSA_MD")

merge_prefix         <- file.path(genomewide_pca_dir, "FO_GSA_Merged")
qc_prefix            <- file.path(genomewide_pca_dir, "FO_GSA_Merged_QC")
pruned_prefix        <- file.path(genomewide_pca_dir, "FO_GSA_Merged_QC_Pruned")
unrelated_prefix     <- file.path(genomewide_pca_dir, "FO_GSA_Merged_QC_Pruned_Unrelated")
pca_full_prefix      <- file.path(genomewide_pca_dir, "FO_GSA_Global_PCA_FullSample")
pca_core_prefix      <- file.path(genomewide_pca_dir, "FO_GSA_Global_PCA_UnrelatedCore")
pca_projected_prefix <- file.path(genomewide_pca_dir, "FO_GSA_Global_PCA_Projected_Sensitivity")

n_pcs <- 3

# Known long-range LD regions (hg19/GRCh37), format: CHR START END LABEL
long_range_ld_file <- file.path(genomewide_pca_dir, "high_LD_regions_hg19.txt")

writeLines(
  c(
    "1\t48000000\t52000000\tregion1",
    "2\t86000000\t100500000\tregion2",
    "2\t134500000\t138000000\tregion3",
    "2\t183000000\t190000000\tregion4",
    "3\t47500000\t50000000\tregion5",
    "3\t83500000\t87000000\tregion6",
    "3\t89000000\t97500000\tregion7",
    "5\t44500000\t50500000\tregion8",
    "5\t98000000\t100500000\tregion9",
    "5\t129000000\t132000000\tregion10",
    "5\t135500000\t138500000\tregion11",
    "6\t25500000\t33500000\tregion12_HLA",
    "6\t57000000\t64000000\tregion13",
    "6\t140000000\t142500000\tregion14",
    "7\t55000000\t66000000\tregion15",
    "8\t8000000\t12000000\tregion16_inversion",
    "8\t43000000\t50000000\tregion17",
    "8\t112000000\t115000000\tregion18",
    "10\t37000000\t43000000\tregion19",
    "11\t46000000\t57000000\tregion20",
    "11\t87500000\t90500000\tregion21",
    "12\t33000000\t40000000\tregion22",
    "12\t109500000\t112000000\tregion23",
    "20\t32000000\t34500000\tregion24"
  ),
  long_range_ld_file
)

#==================================================================================
# HARMONIZATION DIAGNOSTICS (strand/allele conflicts between Ambu and MD)
#==================================================================================
# Purpose: before any merge/PCA/gene loop, systematically check whether and to
# what extent allele/strand inconsistencies exist between the two batches
# (Ambu, MD), and whether these are strand/allele-technical
# (resolvable by flipping) or indicate a true batch effect (SNP is
# called systematically differently in the two batches - flipping does NOT
# resolve this, as flipping only changes the allele LABEL, not the
# underlying genotype call).
#
# Workflow:
#   1a) Count ambiguous SNPs (A/T, C/G) per batch (strand not determinable
#       from alleles alone)
#   1b) Raw merge test run (complete array, NO exclusions) ->
#       Number of conflicts from PLINK .missnp, split into ambiguous/unambiguous
#       (purely informative, baseline)
#   1c) Exclude ambiguous SNPs, retry merge -> remaining
#       (unambiguous) conflicts
#   1d) Flip remaining conflicts, retry merge -> resolved vs.
#       still conflicting (the latter are excluded)
#   1e) SNP-level batch-effect test for the now conflict-free
#       dataset: allele-frequency difference and differential missingness
#       between Ambu and MD (batch as pseudo-phenotype)
#   1f) Overall funnel summary across all steps
#   1g) Final exclusion list (ambiguous + not resolvable by flip + batch effect)
#       save as reusable file
#
# IMPORTANT: this block does NOT alter the actual analysis pipeline
# automatically - it produces exclusively diagnostic tables (console + CSV
# export) and an exclusion list for manual review. Applying
# this list in the genome-wide PCA merge (below) and in the gene loop
# (run_gene_pipeline) happens only in a separate step,
# to be decided.
#==================================================================================

harmonization_dir <- file.path(base_dir, "Harmonization_Diagnostics")
dir.create(harmonization_dir, showWarnings = FALSE, recursive = TRUE)

test_merge_dir <- file.path(harmonization_dir, "test_merges")
dir.create(test_merge_dir, showWarnings = FALSE, recursive = TRUE)

batch_effect_pval_threshold <- 1e-5  # <- adapt if needed

harmon_step_log <- list()  # collects the step-by-step funnel table

log_harmon_step <- function(step, description, n_snps_before, n_snps_removed, n_snps_after) {
  harmon_step_log[[length(harmon_step_log) + 1]] <<- data.table(
    Step = step,
    Description = description,
    SNPs_Before = n_snps_before,
    SNPs_Removed = n_snps_removed,
    SNPs_After = n_snps_after
  )
}

read_bim <- function(prefix) {
  fread(
    paste0(prefix, ".bim"), header = FALSE,
    col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2")
  )
}

is_ambiguous <- function(a1, a2) {
  pair <- paste0(toupper(a1), toupper(a2))
  pair %in% c("AT", "TA", "CG", "GC")
}

cat("\n#####################################################\n")
cat("# START HARMONIZATION DIAGNOSTICS (Ambu vs. MD)\n")
cat("#####################################################\n")

bim_ambu_raw <- read_bim(raw_ambu_prefix)
bim_md_raw   <- read_bim(raw_md_prefix)

n_snps_ambu_raw  <- nrow(bim_ambu_raw)
n_snps_md_raw    <- nrow(bim_md_raw)
n_snps_union_raw <- length(union(bim_ambu_raw$SNP, bim_md_raw$SNP))

cat("\n=====================================\n")
cat("baseline\n")
cat("=====================================\n")
cat("SNPs Ambu (roh):        ", n_snps_ambu_raw, "\n")
cat("SNPs MD (roh):          ", n_snps_md_raw, "\n")
cat("SNPs Union (Ambu u MD): ", n_snps_union_raw, "\n")

#--------------------------------------------------------------------------
# 1a) Ambiguous SNPs (A/T, C/G) per batch
#--------------------------------------------------------------------------
bim_ambu_raw[, ambiguous := is_ambiguous(A1, A2)]
bim_md_raw[,   ambiguous := is_ambiguous(A1, A2)]

n_ambig_ambu <- sum(bim_ambu_raw$ambiguous)
n_ambig_md   <- sum(bim_md_raw$ambiguous)

ambig_snps_ambu <- bim_ambu_raw$SNP[bim_ambu_raw$ambiguous]
ambig_snps_md   <- bim_md_raw$SNP[bim_md_raw$ambiguous]
ambig_snps_all  <- unique(c(ambig_snps_ambu, ambig_snps_md))

tab_1a <- data.table(
  Batch = c("Ambu", "MD"),
  SNPs_Total = c(n_snps_ambu_raw, n_snps_md_raw),
  SNPs_Ambiguous_AT_CG = c(n_ambig_ambu, n_ambig_md),
  Pct_Ambiguous = round(100 * c(n_ambig_ambu / n_snps_ambu_raw, n_ambig_md / n_snps_md_raw), 2)
)

cat("\n--- Step 1a: ambiguous SNPs (A/T, C/G) per batch ---\n")
print(tab_1a)
fwrite(tab_1a, file.path(harmonization_dir, "1a_Ambiguous_SNPs_by_Batch.csv"), sep = ";")

log_harmon_step(
  "1a_ambiguous_identified",
  "Ambige SNPs (A/T, C/G) identifiziert (nur markiert, noch nicht entfernt)",
  n_snps_before = n_snps_union_raw,
  n_snps_removed = 0,
  n_snps_after = n_snps_union_raw
)

#--------------------------------------------------------------------------
# 1b) First raw merge test run (full array, NO exclusions)
#--------------------------------------------------------------------------
merge_test1_prefix <- file.path(test_merge_dir, "Test1_RawMerge")

system2(plink, c(
  "--bfile", raw_ambu_prefix,
  "--bmerge",
  paste0(raw_md_prefix, ".bed"),
  paste0(raw_md_prefix, ".bim"),
  paste0(raw_md_prefix, ".fam"),
  "--make-bed",
  "--out", merge_test1_prefix
))

missnp_test1 <- paste0(merge_test1_prefix, "-merge.missnp")
conflict_snps_1 <- if (file.exists(missnp_test1)) readLines(missnp_test1) else character(0)

n_conflict_1       <- length(conflict_snps_1)
n_conflict_1_ambig <- sum(conflict_snps_1 %in% ambig_snps_all)
n_conflict_1_unamb <- n_conflict_1 - n_conflict_1_ambig

tab_1b <- data.table(
  Metric = c(
    "SNPs gesamt (Union Ambu+MD)",
    "Konflikte gesamt (.missnp)",
    "davon ambig (A/T, C/G)",
    "davon eindeutig (nicht-ambig)",
    "Anteil Konflikte an Gesamt-SNPs (%)"
  ),
  Value = c(
    n_snps_union_raw,
    n_conflict_1,
    n_conflict_1_ambig,
    n_conflict_1_unamb,
    round(100 * n_conflict_1 / n_snps_union_raw, 2)
  )
)

cat("\n--- Step 1b: first raw merge test run (baseline, informative) ---\n")
print(tab_1b)
fwrite(tab_1b, file.path(harmonization_dir, "1b_Initial_Merge_Conflicts.csv"), sep = ";")
writeLines(conflict_snps_1, file.path(harmonization_dir, "1b_Conflict_SNPs_List.txt"))

#--------------------------------------------------------------------------
# 1c) Exclude ambiguous SNPs, retry merge
#--------------------------------------------------------------------------
ambig_exclude_file <- file.path(harmonization_dir, "ambiguous_snps_exclude.txt")
writeLines(ambig_snps_all, ambig_exclude_file)

ambu_noambig_prefix <- file.path(test_merge_dir, "Ambu_NoAmbig")
md_noambig_prefix   <- file.path(test_merge_dir, "MD_NoAmbig")

system2(plink, c(
  "--bfile", raw_ambu_prefix,
  "--exclude", ambig_exclude_file,
  "--make-bed",
  "--out", ambu_noambig_prefix
))
system2(plink, c(
  "--bfile", raw_md_prefix,
  "--exclude", ambig_exclude_file,
  "--make-bed",
  "--out", md_noambig_prefix
))

merge_test2_prefix <- file.path(test_merge_dir, "Test2_MergeNoAmbig")

system2(plink, c(
  "--bfile", ambu_noambig_prefix,
  "--bmerge",
  paste0(md_noambig_prefix, ".bed"),
  paste0(md_noambig_prefix, ".bim"),
  paste0(md_noambig_prefix, ".fam"),
  "--make-bed",
  "--out", merge_test2_prefix
))

missnp_test2 <- paste0(merge_test2_prefix, "-merge.missnp")
conflict_snps_2 <- if (file.exists(missnp_test2)) readLines(missnp_test2) else character(0)
n_conflict_2 <- length(conflict_snps_2)

n_snps_after_ambig_removal <- n_snps_union_raw - length(ambig_snps_all)

tab_1c <- data.table(
  Metric = c(
    "SNPs vor diesem Schritt (Union)",
    "Ambige SNPs entfernt",
    "SNPs nach Ausschluss ambiger SNPs",
    "Verbleibende Konflikte (.missnp) danach",
    "Anteil verbleibender Konflikte (%)"
  ),
  Value = c(
    n_snps_union_raw,
    length(ambig_snps_all),
    n_snps_after_ambig_removal,
    n_conflict_2,
    round(100 * n_conflict_2 / n_snps_after_ambig_removal, 2)
  )
)

cat("\n--- Step 1c: ambiguous SNPs removed, merge retried ---\n")
print(tab_1c)
fwrite(tab_1c, file.path(harmonization_dir, "1c_Merge_After_Ambiguous_Removal.csv"), sep = ";")
writeLines(conflict_snps_2, file.path(harmonization_dir, "1c_Remaining_Conflict_SNPs.txt"))

log_harmon_step(
  "1c_remove_ambiguous",
  "Ambige SNPs (A/T, C/G) entfernt",
  n_snps_before = n_snps_union_raw,
  n_snps_removed = length(ambig_snps_all),
  n_snps_after = n_snps_after_ambig_removal
)

#--------------------------------------------------------------------------
# 1d) Flip remaining conflicts, retry merge
#--------------------------------------------------------------------------
md_flipped_prefix <- file.path(test_merge_dir, "MD_NoAmbig_Flipped")

if (n_conflict_2 > 0) {
  system2(plink, c(
    "--bfile", md_noambig_prefix,
    "--flip", missnp_test2,
    "--make-bed",
    "--out", md_flipped_prefix
  ))
} else {
  file.copy(paste0(md_noambig_prefix, ".bed"), paste0(md_flipped_prefix, ".bed"), overwrite = TRUE)
  file.copy(paste0(md_noambig_prefix, ".bim"), paste0(md_flipped_prefix, ".bim"), overwrite = TRUE)
  file.copy(paste0(md_noambig_prefix, ".fam"), paste0(md_flipped_prefix, ".fam"), overwrite = TRUE)
}

merge_test3_prefix <- file.path(test_merge_dir, "Test3_MergeAfterFlip")

system2(plink, c(
  "--bfile", ambu_noambig_prefix,
  "--bmerge",
  paste0(md_flipped_prefix, ".bed"),
  paste0(md_flipped_prefix, ".bim"),
  paste0(md_flipped_prefix, ".fam"),
  "--make-bed",
  "--out", merge_test3_prefix
))

missnp_test3 <- paste0(merge_test3_prefix, "-merge.missnp")
conflict_snps_3 <- if (file.exists(missnp_test3)) readLines(missnp_test3) else character(0)
n_conflict_3 <- length(conflict_snps_3)
n_resolved_by_flip <- n_conflict_2 - n_conflict_3

# If conflicts remain after flipping, exclude them finally
# and retry the merge one last time so that step 1e can build on a
# truly complete, conflict-free dataset.
if (n_conflict_3 > 0) {
  ambu_final_prefix <- file.path(test_merge_dir, "Ambu_Final_NoConflict")
  md_final_prefix   <- file.path(test_merge_dir, "MD_Final_NoConflict")

  system2(plink, c(
    "--bfile", ambu_noambig_prefix,
    "--exclude", missnp_test3,
    "--make-bed",
    "--out", ambu_final_prefix
  ))
  system2(plink, c(
    "--bfile", md_flipped_prefix,
    "--exclude", missnp_test3,
    "--make-bed",
    "--out", md_final_prefix
  ))

  merge_test4_prefix <- file.path(test_merge_dir, "Test4_MergeFinal")
  system2(plink, c(
    "--bfile", ambu_final_prefix,
    "--bmerge",
    paste0(md_final_prefix, ".bed"),
    paste0(md_final_prefix, ".bim"),
    paste0(md_final_prefix, ".fam"),
    "--make-bed",
    "--out", merge_test4_prefix
  ))
  merge_clean_prefix <- merge_test4_prefix
} else {
  merge_clean_prefix <- merge_test3_prefix
}

n_snps_final_conflictfree <- n_snps_after_ambig_removal - n_conflict_3

tab_1d <- data.table(
  Metric = c(
    "SNPs vor Flip-Versuch",
    "Konflikte vor Flip",
    "Durch Flip geloest (nur Beschriftung geaendert)",
    "Nach Flip weiterhin konfliktbehaftet (werden ausgeschlossen)",
    "SNPs im Datensatz nach diesem Schritt"
  ),
  Value = c(
    n_snps_after_ambig_removal,
    n_conflict_2,
    n_resolved_by_flip,
    n_conflict_3,
    n_snps_final_conflictfree
  )
)

cat("\n--- Step 1d: remaining conflicts flipped ---\n")
print(tab_1d)
fwrite(tab_1d, file.path(harmonization_dir, "1d_Flip_Results.csv"), sep = ";")
writeLines(conflict_snps_3, file.path(harmonization_dir, "1d_Still_Conflicting_After_Flip.txt"))

log_harmon_step(
  "1d_flip_and_exclude_unresolved",
  "Verbleibende eindeutige Konflikte geflippt; danach weiterhin konfliktbehaftete SNPs ausgeschlossen",
  n_snps_before = n_snps_after_ambig_removal,
  n_snps_removed = n_conflict_3,
  n_snps_after = n_snps_final_conflictfree
)

#--------------------------------------------------------------------------
# 1e) SNP-level batch-effect test (allele frequency + missingness by batch)
#--------------------------------------------------------------------------
# Batch (Ambu/MD) is coded here as pseudo-phenotype (1 = Ambu, 2 = MD)
# to test with standard PLINK association tests whether individual SNPs
# show batch-specific allele frequencies or missingness -
# an indication of technical batch artefacts that are NOT fixed by simply
# flipping the allele labels.

fam_ambu_ids <- fread(
  paste0(raw_ambu_prefix, ".fam"), header = FALSE, colClasses = "character"
)$V2
fam_md_ids <- fread(
  paste0(raw_md_prefix, ".fam"), header = FALSE, colClasses = "character"
)$V2

fam_clean <- fread(
  paste0(merge_clean_prefix, ".fam"), header = FALSE, colClasses = "character",
  col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")
)

batch_pheno <- data.table(
  FID = fam_clean$FID,
  IID = fam_clean$IID,
  PHENO = ifelse(
    fam_clean$IID %in% fam_ambu_ids, 1,
    ifelse(fam_clean$IID %in% fam_md_ids, 2, NA)
  )
)

n_match_ambu <- sum(fam_clean$IID %in% fam_ambu_ids)
n_match_md   <- sum(fam_clean$IID %in% fam_md_ids)
n_unmatched  <- sum(is.na(batch_pheno$PHENO))

cat("\n[Diagnostic batch pheno] IID matching before PLINK call:\n")
cat("  Individuals in merge:      ", nrow(fam_clean), "\n")
cat("  Of these recognised as Ambu:   ", n_match_ambu, "\n")
cat("  Of these recognised as MD:     ", n_match_md, "\n")
cat("  Not assigned (NA):    ", n_unmatched, "\n")

if (n_unmatched == nrow(fam_clean)) {
  warning(
    "WARNING: not a single individual could be assigned to a batch via IID ",
    "- probably a type/formatting problem of the IDs (e.g. ",
    "leading zeros). The subsequent batch-effect test (1e) will not deliver ",
    "meaningful results with this assignment - please check IDs ",
    "before table 1e is interpreted."
  )
}

batch_pheno_file <- file.path(harmonization_dir, "batch_pseudo_pheno.txt")

# write.table instead of fwrite: more robust against Windows line breaks
# (CRLF), which PLINK's simple text parser may in the worst case completely
# discards (symptom: "0 phenotype values present after --pheno", although
# FID/IID demonstrably match - this was exactly the previous problem).
write.table(
  batch_pheno,
  file = batch_pheno_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  na = "-9",
  eol = "\n"
)

# Re-read check: interpret the file exactly as PLINK would
# (whitespace/tab-separated) to detect format problems BEFORE the PLINK call
# instead of failing at "0 phenotype values present".
pheno_check <- tryCatch(
  read.table(batch_pheno_file, header = FALSE, sep = "", stringsAsFactors = FALSE),
  error = function(e) NULL
)

cat("\n[Diagnostic batch pheno file] Re-read check (as PLINK would parse):\n")
if (is.null(pheno_check)) {
  warning(
    "Konnte '", batch_pheno_file, "' as a 3-column table - ",
    "file is probably malformatted (e.g. CRLF ",
    "line endings). Please check the file in a text editor before PLINK ",
    "erneut aufgerufen wird."
  )
} else {
  cat("  Rows detected:", nrow(pheno_check), "of", nrow(batch_pheno), "expected\n")
  cat("  Columns detected:", ncol(pheno_check), "(expected: 3)\n")
  cat("  First 3 rows:\n")
  print(head(pheno_check, 3))
  if (nrow(pheno_check) != nrow(batch_pheno) || ncol(pheno_check) != 3) {
    warning(
      "Re-read check shows a deviation from the expected row/",
      "Spaltenzahl - bitte '", batch_pheno_file, "' manually (e.g. ",
      "for CRLF/^M characters or wrong separator) before table 1e ",
      "interpretiert wird."
    )
  } else {
    cat("  -> Format looks plausible (3 columns, expected number of rows).\n")
  }
}

batch_assoc_prefix <- file.path(test_merge_dir, "Batch_Effect_Assoc")

system2(plink, c(
  "--bfile", merge_clean_prefix,
  "--pheno", batch_pheno_file,
  "--allow-no-sex",
  "--assoc",
  "--out", batch_assoc_prefix
))

system2(plink, c(
  "--bfile", merge_clean_prefix,
  "--pheno", batch_pheno_file,
  "--allow-no-sex",
  "--test-missing",
  "--out", batch_assoc_prefix
))

assoc_file   <- paste0(batch_assoc_prefix, ".assoc")
missing_file <- paste0(batch_assoc_prefix, ".missing")

flagged_freq    <- character(0)
flagged_missing <- character(0)

if (file.exists(assoc_file)) {
  assoc_res <- fread(assoc_file)
  if ("P" %in% colnames(assoc_res)) {
    flagged_freq <- assoc_res$SNP[!is.na(assoc_res$P) & assoc_res$P < batch_effect_pval_threshold]
  } else {
    warning("Column 'P' not found in .assoc file - allele-frequency batch test skipped.")
  }
}

if (file.exists(missing_file)) {
  missing_res <- fread(missing_file)
  if ("P" %in% colnames(missing_res)) {
    flagged_missing <- missing_res$SNP[!is.na(missing_res$P) & missing_res$P < batch_effect_pval_threshold]
  } else {
    warning("Column 'P' not found in .missing file - missingness batch test skipped.")
  }
}

batch_effect_flagged_snps <- unique(c(flagged_freq, flagged_missing))

assoc_test_ran    <- file.exists(assoc_file) && exists("assoc_res") && "P" %in% colnames(assoc_res) && sum(!is.na(assoc_res$P)) > 0
missing_test_ran  <- file.exists(missing_file) && exists("missing_res") && "P" %in% colnames(missing_res) && sum(!is.na(missing_res$P)) > 0

tab_1e <- data.table(
  Metric = c(
    "SNPs getestet",
    "Schwellenwert (P <)",
    "Allelfrequenz-Test tatsaechlich durchgefuehrt (nicht 0 gueltige P-Werte)",
    "Missingness-Test tatsaechlich durchgefuehrt (nicht 0 gueltige P-Werte)",
    "Auffaellig: Allelfrequenz-Unterschied Ambu vs. MD",
    "Auffaellig: differenzielle Missingness Ambu vs. MD",
    "Auffaellig gesamt (Vereinigung, werden ausgeschlossen)"
  ),
  Value = c(
    n_snps_final_conflictfree,
    batch_effect_pval_threshold,
    as.character(assoc_test_ran),
    as.character(missing_test_ran),
    length(flagged_freq),
    length(flagged_missing),
    length(batch_effect_flagged_snps)
  )
)

if (!assoc_test_ran || !missing_test_ran) {
  warning(
    "Step 1e: at least one of the batch-effect tests was NOT ",
    "run by PLINK (e.g. due to failed phenotype assignment). ",
    "'0 notable SNPs' in this case does NOT mean 'no batch effect', ",
    "but that the test did not run at all. Please check the console/log output of ",
    "PLINK (--assoc / --test-missing above) and the diagnostic output on ",
    "IID matching before table 1e is interpreted."
  )
}

cat("\n--- Step 1e: batch-effect test (allele frequency + missingness) ---\n")
print(tab_1e)
fwrite(tab_1e, file.path(harmonization_dir, "1e_Batch_Effect_Test.csv"), sep = ";")
writeLines(batch_effect_flagged_snps, file.path(harmonization_dir, "1e_Batch_Effect_Flagged_SNPs.txt"))

log_harmon_step(
  "1e_batch_effect_test",
  paste0("SNPs mit batch-abhaengiger Allelfrequenz oder Missingness (P < ", batch_effect_pval_threshold, ") ausgeschlossen"),
  n_snps_before = n_snps_final_conflictfree,
  n_snps_removed = length(batch_effect_flagged_snps),
  n_snps_after = n_snps_final_conflictfree - length(batch_effect_flagged_snps)
)

#--------------------------------------------------------------------------
# 1f) Overall summary (funnel across all steps)
#--------------------------------------------------------------------------
harmon_funnel_summary <- rbindlist(harmon_step_log)

cat("\n=====================================\n")
cat("Harmonization diagnostics: overall summary (funnel)\n")
cat("=====================================\n")
print(harmon_funnel_summary)

fwrite(
  harmon_funnel_summary,
  file.path(harmonization_dir, "1f_Harmonization_Funnel_Summary.csv"),
  sep = ";"
)

#--------------------------------------------------------------------------
# 1g) Final exclusion list (for later reuse in PCA merge +
#     gene loop - applied only after joint decision)
#--------------------------------------------------------------------------
final_exclude_snps <- unique(c(
  ambig_snps_all,               # A/T, C/G (strand not determinable)
  conflict_snps_3,              # still allele-conflicting after flip
  batch_effect_flagged_snps     # batch-dependent frequency/missingness
))

batch_discordant_snps_file <- file.path(
  genomewide_pca_dir, "Batch_Discordant_SNPs_Exclude.txt"
)
writeLines(final_exclude_snps, batch_discordant_snps_file)

tab_1g <- data.table(
  Kategorie = c(
    "Ambig (A/T, C/G)",
    "Nach Flip weiterhin konfliktbehaftet",
    "Batch-Effekt (Frequenz/Missingness)",
    "Gesamt (eindeutige SNP-IDs, Vereinigung)"
  ),
  N_SNPs = c(
    length(ambig_snps_all),
    length(conflict_snps_3),
    length(batch_effect_flagged_snps),
    length(final_exclude_snps)
  )
)

cat("\n--- Step 1g: final exclusion list ---\n")
print(tab_1g)
cat("Saved to:", batch_discordant_snps_file, "\n")
cat(
  "\nNOTE: this list is so far NOT applied automatically.\n",
  "Use in the genome-wide PCA merge and in the gene loop (run_gene_pipeline)\n",
  "only after joint review of the tables above.\n"
)

fwrite(tab_1g, file.path(harmonization_dir, "1g_Final_Exclusion_List_Summary.csv"), sep = ";")

cat("\n#####################################################\n")
cat("# END HARMONIZATION DIAGNOSTICS\n")
cat("# All tables saved as CSV in:", harmonization_dir, "\n")
cat("#####################################################\n\n")

# --- 1) Merge Ambu + MD --------------------------------------------------
# Replaces the earlier raw --bmerge (which failed due to 321,883 allele conflicts,
# see harmonization diagnostics above, steps 1a-1g). Instead of
# merging raw again, the dataset already created in the diagnostics,
# harmonized dataset taken over (merge_clean_prefix: ambiguous SNPs
# removed, strand conflicts flipped, thereafter still conflicting
# SNPs removed) and additionally cleaned of batch-effect SNPs
# (Batch_Discordant_SNPs_Exclude.txt from step 1g). This list
# also contains the already removed ambiguous/conflict SNPs for safety;
# excluding a SNP twice is harmless (PLINK ignores non-existing
# IDs in --exclude). Result: 690,854 - 64,372 = 626,482 SNPs,
# all 1,552 individuals. The variable merge_prefix remains unchanged so that
# the subsequent steps (QC, pruning, PCA) continue unchanged.
stopifnot(
  "merge_clean_prefix nicht gefunden - Harmonisierungs-Diagnostik (Schritte 1a-1g) muss vorher gelaufen sein" =
    exists("merge_clean_prefix"),
  "batch_discordant_snps_file nicht gefunden - Harmonisierungs-Diagnostik (Schritt 1g) muss vorher gelaufen sein" =
    exists("batch_discordant_snps_file")
)

system2(plink, c(
  "--bfile", merge_clean_prefix,
  "--exclude", batch_discordant_snps_file,
  "--make-bed",
  "--out", merge_prefix
))

n_snps_merge_prefix <- tryCatch(
  length(readLines(paste0(merge_prefix, ".bim"))),
  error = function(e) NA
)
n_inds_merge_prefix <- tryCatch(
  length(readLines(paste0(merge_prefix, ".fam"))),
  error = function(e) NA
)

cat("\n--- Step 1 (harmonized): merge result ---\n")
cat("SNPs in merge_prefix:      ", n_snps_merge_prefix, "\n")
cat("Individuals in merge_prefix:", n_inds_merge_prefix, "\n")

# --- 2) QC: autosomes only, MAF, SNP/individual missingness ----------------
system2(plink, c(
  "--bfile", merge_prefix,
  "--autosome",
  "--maf", "0.01",
  "--geno", "0.05",
  "--mind", "0.1",
  "--make-bed",
  "--out", qc_prefix
))

# --- 3) Exclude long-range LD regions + LD pruning list ----------
system2(plink, c(
  "--bfile", qc_prefix,
  "--exclude", "range", long_range_ld_file,
  "--indep-pairwise", "50", "5", "0.2",
  "--out", pruned_prefix
))

# --- 4) Extract pruned SNP list (full QC dataset, all individuals) --
system2(plink, c(
  "--bfile", qc_prefix,
  "--extract", paste0(pruned_prefix, ".prune.in"),
  "--make-bed",
  "--out", pruned_prefix
))

# --- 5) KING-robust kinship: determine unrelated core set -------------
# Writes <out>.king.cutoff.in.id (unrelated) / .king.cutoff.out.id
# (excluded due to relatedness >= cutoff).
# IMPORTANT: --king-cutoff is a PLINK2 flag, PLINK 1.9 does not know it
# ("Error: Unrecognized flag") - hence plink2 instead of plink here.
system2(plink2, c(
  "--bfile", pruned_prefix,
  "--king-cutoff", "0.0884",
  "--out", unrelated_prefix
))

king_in_id_file  <- paste0(unrelated_prefix, ".king.cutoff.in.id")
king_out_id_file <- paste0(unrelated_prefix, ".king.cutoff.out.id")

if (file.exists(king_out_id_file)) {
  n_related_excluded <- length(readLines(king_out_id_file)) - 1  # subtract header
} else {
  n_related_excluded <- 0
}

cat("\n=====================================\n")
cat("KING kinship filter (cutoff 0.0884, ~3rd-degree relatedness)\n")
cat("=====================================\n")
cat("Excluded due to relatedness (core set):", n_related_excluded, "individuals\n")

# Create bed/bim/fam dataset of the unrelated core set
system2(plink, c(
  "--bfile", pruned_prefix,
  "--keep", king_in_id_file,
  "--make-bed",
  "--out", unrelated_prefix
))

# --- 6) PCA on the FULL pruned dataset (all individuals) -----------
# These are the PCs used in the multivariable models.
system2(plink2, c(
  "--bfile", pruned_prefix,
  "--pca", as.character(n_pcs),
  "--out", pca_full_prefix
))

# --- 7) Sensitivity check: PCA on the unrelated core set + projection ---
# of the FULL sample onto the same axes. Serves exclusively to
# check whether sporadic relatedness distorts the PCs from step 6.
system2(plink2, c(
  "--bfile", unrelated_prefix,
  "--freq",
  "--pca", as.character(n_pcs), "allele-wts",
  "--out", pca_core_prefix
))

# Do NOT hard-code column numbers for ID/A1 and the PC columns; instead
# read them from the actual header of the .eigenvec.allele file. Reason:
# in the previous run (PLINK v2.0.0-a.7.1), the hard-coded
# columns "2" (ID) and "5" (A1) discarded ALL 381,724 entries (= 2 x 190,862
# variants, one row per allele) as "mismatching allele codes"
# -> "No valid variants". This is not a true allele/strand
# problem but a column offset: newer PLINK2 versions write in
# --pca allele-wts additional leading columns (e.g. #CHROM, REF, ALT) before
# A1, so that column 5 no longer corresponds to A1. Fix: column names in the header
# look up instead of guessing positions.
allele_header   <- fread(paste0(pca_core_prefix, ".eigenvec.allele"), nrows = 0)
allele_colnames <- colnames(allele_header)
id_col  <- which(allele_colnames == "ID")
a1_col  <- which(allele_colnames == "A1")
pc_cols <- grep("^PC[0-9]+$", allele_colnames)
stopifnot(
  "ID-Spalte nicht in .eigenvec.allele gefunden" = length(id_col) == 1,
  "A1-Spalte nicht in .eigenvec.allele gefunden" = length(a1_col) == 1,
  "Anzahl PC-Spalten in .eigenvec.allele stimmt nicht mit n_pcs ueberein" =
    length(pc_cols) == n_pcs
)
score_col_range <- paste0(min(pc_cols), "-", max(pc_cols))

system2(plink2, c(
  "--bfile", pruned_prefix,
  "--read-freq", paste0(pca_core_prefix, ".afreq"),
  "--score", paste0(pca_core_prefix, ".eigenvec.allele"),
  as.character(id_col), as.character(a1_col),
  "header-read", "no-mean-imputation", "variance-standardize",
  "--score-col-nums", score_col_range,
  "--out", pca_projected_prefix
))

# --- 8) Read into R, sensitivity comparison, save final PCs -----
# PLINK2 writes a header line with --pca (#FID IID PC1 PC2 PC3) - with
# header = FALSE this would be read as an additional "individual" (see
# previous run: 1538 individuals in the PCA, but 1539 in the saved
# final result). Hence header = TRUE here; setnames() renames the columns
# by position anyway, independent of the original header.
eigenvec_full <- fread(paste0(pca_full_prefix, ".eigenvec"), header = TRUE)
setnames(
  eigenvec_full,
  c("FID", "IID", paste0("gPC", seq_len(n_pcs)))
)

sscore_cols <- fread(paste0(pca_projected_prefix, ".sscore"), header = TRUE)
# The last n_pcs columns of the .sscore file contain the projected PCs
# (column names depend on the PLINK2 version, e.g. "PC1_AVG" or similar).
proj_pc_cols <- tail(colnames(sscore_cols), n_pcs)

eigenvec_projected <- sscore_cols[
  ,
  c("IID", proj_pc_cols),
  with = FALSE
]
setnames(eigenvec_projected, c("IID", paste0("gPC", seq_len(n_pcs), "_projected")))

sensitivity_check <- merge(
  eigenvec_full[, .(IID, gPC1, gPC2, gPC3)],
  eigenvec_projected,
  by = "IID"
)

cat("\n=====================================\n")
cat("Sensitivity check: PCs (full sample) vs. PCs (unrelated core set, projected)\n")
cat("=====================================\n")
for (k in seq_len(n_pcs)) {
  corr_k <- suppressWarnings(cor(
    sensitivity_check[[paste0("gPC", k)]],
    sensitivity_check[[paste0("gPC", k, "_projected")]],
    use = "complete.obs"
  ))
  cat("Correlation gPC", k, ": ", round(corr_k, 4), "\n", sep = "")
}
cat(
  "Note: values close to 1 confirm that sporadic relatedness ",
  "does not relevantly distort the PCs from the full sample AND that the ",
  "PLINK2 projection (step 7) worked technically correctly. In case of ",
  "clearly deviating values, please review before the PCs ",
  "weiterverwendet werden.\n"
)

global_genomewide_pcs <- eigenvec_full[, .(id = IID, gPC1, gPC2, gPC3)]

saveRDS(global_genomewide_pcs, genomewide_pcs_file)

write.csv2(
  global_genomewide_pcs,
  file.path(genomewide_pca_dir, "Global_Genomewide_PCs.csv"),
  row.names = FALSE
)

cat("\nGenome-wide PCs saved:", genomewide_pcs_file, "\n")
cat("Number of individuals with genome-wide PCs:", nrow(global_genomewide_pcs), "\n")
