# Pathway-restricted polygenic scores for energy-substrate metabolism and obesity (FoCus cohort)
# Script: 03_pathway_pipeline.R | Source file: Protocol_SNP_Analysis_FoCus_Pathway_Template_Split_19082026_REPAIR.R
# Local paths replaced by the environment variable FOCUS_BASE_DIR (project root); analysis code unchanged.
# Individual-level FoCus data are not included (see README).

#==================================================================================
# REPAIR VERSION (2026-10-09) - re-run after correction of the defective phenotype CSVs
#==================================================================================
# Copy of Protocol_SNP_Analysis_FoCus_Pathway_Template_Split_19082026.R
# Changes ONLY to paths (all marked '# REPAIR'), analysis code unchanged:
#   - phenotype CSVs from FoCus Daten/Daten/CSV_repariert
#   - all write access to FoCus/Repair/... (OutPW = Outputs_Pathway, Working_Plink,
#     Genetik_Daten_FoCus, Methods_Results_Log_Repair.md)
#   - Global_TrainTest_Split.rds and Pathway_<X>_VEP.txt are ONLY read from the
#     original Outputs_Pathway (never regenerated)
# Original folder and files are not written.
#==================================================================================

#======================================================================================
# R PACKAGES
#====================================================================================
#Clean up R: 
rm(list = ls())     # remove all objects
gc()                # free memory

options(repos = c(CRAN = "https://cloud.r-project.org"))   # CRAN mirror for non-interactive Rscript runs

#BigSnpr
if (!requireNamespace("bigsnpr", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  library(devtools)
  devtools::install_github("privefl/bigsnpr")
}
library(bigsnpr)
#snpStats
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("snpStats", quietly = TRUE)) BiocManager::install("snpStats")
library(snpStats)
if (!requireNamespace("SKAT", quietly = TRUE)) install.packages("SKAT")
library(SKAT)
library(broom)
library(ggplot2)
library(car)
library(dplyr)
library(pROC)
library(data.table)


#==================================================================================
# BATCH-PARAMETER (KOMMANDOZEILE)
#==================================================================================
# For the automated run over several pathways, pathway_name,
# pathway_selection_id and run_phase are passed as command-line arguments;
# z. B.: Rscript Protocol_SNP_Analysis_FoCus_Pathway_Template.R Citrate_Cycle hsa00020 1
#
# run_phase controls the breakpoint before the external VCF creation/VEP
# annotation (see block "VCF + VEP ANNOTATION (EXTERNAL STEP)" further
# below):
#   "1" - script runs up to the breakpoint (incl. LD pruning/PLINK export)
#         and then stops cleanly BEFORE the external VCF/VEP generation is required
#         would be.
#   "2" - script runs completely (default; requires that the
#         VEP annotation file Pathway_[pathway]_VEP.txt has already been
#         generated externally).
#
# If the script is still started interactively (e.g. in RStudio) without arguments,
# the previous default values apply unchanged - nothing changes
# in the usual manual work with a single pathway.
cli_args <- commandArgs(trailingOnly = TRUE)

#==================================================================================
# DEFINE PATHWAY
#==================================================================================
# Name as in "<FOCUS_BASE_DIR>/Pathway_Selection/Pathway_Gene_Selection_final.xlsx"
pathway_name <- if (length(cli_args) >= 1) cli_args[1] else "Citrate_Cycle"

#==================================================================================
# PATHWAY SELECTION: READ DATABASE ORIGIN (KEGG vs. KEGG+Reactome)
#==================================================================================
# TODO: adapt manually per pathway (pathway ID from the pathway-selection
# gene-set table, column "Pathway-ID"; for KEGG pathways WITHOUT "_erweitert"
# suffix, for Reactome pathways the complete R-HSA ID incl. version
# suffix, z. B. "R-HSA-70326.10")
pathway_selection_id <- if (length(cli_args) >= 2) cli_args[2] else "hsa00020"

# Breakpoint control for the batch run (see explanation above)
run_phase <- if (length(cli_args) >= 3) cli_args[3] else "2"

#==================================================================================
# METHODS/RESULTS-LOG (AUTOMATISCHES PROTOKOLL FUER DAS SPAETERE MANUSKRIPT)
#==================================================================================
# At each phase-1 breakpoint and at the end of phase 2, a row with
# R/package/PLINK/VEP version information is automatically inserted into the Appendix A
# table of Methods_Results_Log_FoCus_Pathway_Analysis.md (see
# markers there "<!-- AUTO-LOG:BEGIN -->" / "<!-- AUTO-LOG:END -->").
# No manual entry of versions needed anymore.
# REPAIR: base paths (defined after rm(list = ls()))
repair_base       <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Repair")   # REPAIR
orig_output_dir   <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Outputs_Pathway")   # REPAIR: read only (split, VEP)
pheno_repair_dir  <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "FoCus Daten/Daten/CSV_repariert")   # REPAIR
repair_plink_dir  <- file.path(repair_base, "Working_Plink", "plink_pathway_outputs")   # REPAIR
orig_output_file  <- function(...) file.path(orig_output_dir, paste0(...))   # REPAIR
dir.create(repair_plink_dir, recursive = TRUE, showWarnings = FALSE)   # REPAIR

methods_log_path <- file.path(repair_base, "Methods_Results_Log_Repair.md")   # REPAIR

# Read the PLINK version from the first line of a PLINK-generated .log file
# (PLINK automatically writes a .log file with
# the version identifier in line 1 on every call).
extract_plink_version <- function(log_file) {
  if (!file.exists(log_file)) return(NA_character_)
  trimws(readLines(log_file, n = 1, warn = FALSE))
}

# Read VEP version/cache metadata from the "##" comment lines at the start of the
# Read VEP result file (standard output format of Ensembl VEP; fread()
# skips these lines when reading the data table, but they remain
# in the raw file).
extract_vep_meta <- function(vep_file, max_lines = 50) {
  if (!file.exists(vep_file)) return(NA_character_)
  lines <- readLines(vep_file, n = max_lines, warn = FALSE)
  meta_lines <- lines[startsWith(lines, "##")]
  if (length(meta_lines) == 0) return(NA_character_)
  paste(gsub("\\|", "/", meta_lines), collapse = "; ")
}

# Appends a row with version/context information to the Appendix A
# table of the methods/results log (directly before the marker
# "<!-- AUTO-LOG:END -->", so that subsequent sections such as "Changelog"
# below the table remain untouched).
log_pipeline_run <- function(md_path, pathway_name, phase,
                              plink_version = NA_character_,
                              vep_meta = NA_character_,
                              extra_notes = "") {

  if (!file.exists(md_path)) {
    warning(
      "Methods/results log not found at: ", md_path,
      " - log entry was NOT written."
    )
    return(invisible(NULL))
  }

  log_lines <- readLines(md_path, warn = FALSE, encoding = "UTF-8")
  end_marker <- which(trimws(log_lines) == "<!-- AUTO-LOG:END -->")

  if (length(end_marker) != 1) {
    warning(
      "Marker '<!-- AUTO-LOG:END -->' nicht eindeutig gefunden in ", md_path,
      " - log entry was NOT written. Please check file."
    )
    return(invisible(NULL))
  }

  pkg_names <- c("snpStats", "SKAT", "bigsnpr")
  pkg_versions <- paste(
    sprintf(
      "%s %s", pkg_names,
      sapply(pkg_names, function(p) {
        if (requireNamespace(p, quietly = TRUE)) {
          as.character(packageVersion(p))
        } else {
          "NA"
        }
      })
    ),
    collapse = "; "
  )

  safe_cell <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x) || x == "") return("")
    gsub("\\|", "/", gsub("[\r\n]+", " ", x))
  }

  new_row <- sprintf(
    "| %s | %s | %s | %s | %s | %s | %s | %s |",
    pathway_name, phase, format(Sys.time(), "%Y-%m-%d %H:%M"),
    R.version.string, pkg_versions,
    safe_cell(plink_version), safe_cell(vep_meta), safe_cell(extra_notes)
  )

  updated_lines <- append(log_lines, new_row, after = end_marker - 1)
  writeLines(updated_lines, md_path, useBytes = TRUE)

  cat(
    "Methods/Results-Log aktualisiert: ", md_path,
    " (Pathway: ", pathway_name, ", Phase: ", phase, ")\n", sep = ""
  )

  invisible(NULL)
}

pathway_selection_file <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Pathway_Selection/Pathway_Selection_Supplement_GeneSets.csv")

pathway_selection_data <- read.csv(
  pathway_selection_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Base entry (Datensatz-Art "KEGG" or "Reactome", NOT "KEGG erweitert" -
# that is the separately maintained sensitivity variant, see below)
base_rows <- pathway_selection_data[
  pathway_selection_data$`Pathway-ID` == pathway_selection_id &
    pathway_selection_data$`Datensatz-Art` %in% c("KEGG", "Reactome"),
]

if (nrow(base_rows) == 0) {
  stop(
    "Pathway-ID '", pathway_selection_id,
    "' was not found in the pathway-selection gene-set table. ",
    "Please check pathway_selection_id."
  )
}

pathway_db_type <- unique(base_rows$`Datensatz-Art`)

if (length(pathway_db_type) > 1) {
  stop(
    "Ambiguous Datensatz-Art for pathway ID '", pathway_selection_id,
    "': ", paste(pathway_db_type, collapse = ", ")
  )
}

# ------------------------------------------------------------------------
# pathway_file_id: uniform file-naming convention, identical to the one used
# when creating Pathway_Gene_Selection.xlsx and the range/PLINK files
# sanitization used (pathway name -> file name). Used exclusively
# for matching externally (PowerShell/PLINK) generated file names
# - NOT for this script's own output files, which
# continue to use the short "pathway_name" (e.g. "Glycolysis").
# ------------------------------------------------------------------------
clean_pathway_display_name <- function(name) {
  name <- gsub(" - Homo sapiens \\(human\\)", "", name)
  name <- gsub(" \\(KEGG erweitert\\)", "", name)
  name
}

sanitize_pathway_filename <- function(name) {
  name <- gsub("/", "-", name, fixed = TRUE)
  name <- gsub('[\\\\:*?"<>|]', "-", name)
  name <- gsub("\\s+", "_", trimws(name))
  name <- gsub("_+", "_", name)
  name
}

pathway_display_name <- clean_pathway_display_name(base_rows$`Pathway-Name`[1])
pathway_file_id <- sanitize_pathway_filename(pathway_display_name)

cat(
  "Pathway file ID (for range/PLINK/summary files): '",
  pathway_file_id, "'\n", sep = ""
)

# ------------------------------------------------------------------------
# IMPORTANT - ADDITIONAL MANDATORY USE of pathway_file_id:
# Since the introduction of the SHARED train/test split with the gene script
# (see section "GLOBAL TRAIN/TEST SPLIT" further below),
# pathway_file_id is NO LONGER used only for range/PLINK/summary files,
# but also for the CANONICAL storage location of the split.
# The split is generated exclusively by the gene script; there the
# Variable pathway_name_gene_selection CHARACTER-EXACT (upper/lower case,
# special characters, spaces) MUST match the official pathway name from the
# pathway selection, as derived automatically here from
# base_rows$`Pathway-Name` (Pathway_Selection_Supplement_GeneSets.csv)
# (pathway_display_name). If the names differ,
# the split will NOT be found here (see warning in the split section).
# ------------------------------------------------------------------------


# full_genes = gene set for the MAIN ANALYSIS (KEGG+Reactome if KEGG
# pathway with extension; otherwise KEGG alone or a standalone Reactome set)
kegg_genes <- NULL
full_genes <- unique(base_rows$Gen)
include_db_sensitivity <- FALSE

if (pathway_db_type == "KEGG") {

  kegg_genes <- unique(base_rows$Gen)

  extended_id <- paste0(pathway_selection_id, "_erweitert")

  extended_rows <- pathway_selection_data[
    pathway_selection_data$`Pathway-ID` == extended_id &
      pathway_selection_data$`Datensatz-Art` == "KEGG erweitert",
  ]

  if (nrow(extended_rows) > 0) {

    extended_genes <- unique(extended_rows$Gen)
    added_genes <- setdiff(extended_genes, kegg_genes)

    if (length(added_genes) > 0) {

      full_genes <- extended_genes
      include_db_sensitivity <- TRUE

      cat(
        "Database sensitivity analysis KEGG vs. KEGG+Reactome activated for '",
        pathway_selection_id, "': ", length(added_genes),
        " additional Reactome genes (", paste(added_genes, collapse = ", "),
        ").\n", sep = ""
      )

    } else {

      cat(
        "Reactome supplement for '", pathway_selection_id,
        "' is completely redundant (0 additional genes) - ",
        "database sensitivity analysis omitted.\n", sep = ""
      )
    }

  } else {

    cat(
      "No Reactome supplement for KEGG pathway '", pathway_selection_id,
      "' available - database sensitivity analysis omitted.\n", sep = ""
    )
  }

} else if (pathway_db_type == "Reactome") {

  cat(
    "Pathway '", pathway_selection_id,
    "' is a standalone Reactome pathway (zone 1) - ",
    "no additional database sensitivity analysis.\n", sep = ""
  )

} else {

  stop(
    "Unbekannte Datensatz-Art '", pathway_db_type,
    "' for pathway ID '", pathway_selection_id, "'."
  )
}

#==================================================================================
# SET WORKING DIRECTORY
#==================================================================================
getwd()
output_dir <- file.path(repair_base, "OutPW")   # REPAIR (short name due to Windows MAX_PATH = 260)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_file <- function(...) {
  file.path(output_dir, paste0(...))
}

#==================================================================================
# GENOME-WIDE PCs (POPULATION STRATIFICATION) - CENTRAL FILE
#==================================================================================
# Computed once in the gene script (merge Ambu+MD -> QC -> LD pruning ->
# KING cutoff -> PCA, see the "GENOME-WIDE PCA" section there) and loaded
# here unchanged. Replaces the formerly pathway-locally computed
# PC1-PC3 (see further below: gPC1-gPC3 instead of PC1-PC3).
genomewide_pcs_file <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Global_PCA/Global_Genomewide_PCs.rds")

if (!file.exists(genomewide_pcs_file)) {
  stop(
    "Genome-wide PCs not found (", genomewide_pcs_file, "). ",
    "Please run the PCA section in the gene script first."
  )
}

#==================================================================================
# DEFINITION DER IN- UND OUTPUTDATEIEN
#==================================================================================
# Input from PowerShell / PLINK QC
plink_dir <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Working_Plink/plink_pathway_outputs")

plink_ambu_prefix <- file.path(
  plink_dir,
  paste0("Working_Plink_Pathway_", pathway_file_id, "_Ambu_QC")
)

plink_md_prefix <- file.path(
  plink_dir,
  paste0("Working_Plink_Pathway_", pathway_file_id, "_MD_QC")
)

# Shared summary file from the automated PowerShell script (all
# pathways in one run) - filter by PATHWAY == pathway_file_id instead of
# loading a per-pathway file (no longer generated).
qc_summary_file <- file.path(
  plink_dir,
  "all_pathways_qc_summary.tsv"
)

plink_qc_summary_all <- fread(
  qc_summary_file,
  sep = "\t"
)

plink_qc_summary <- plink_qc_summary_all[PATHWAY == pathway_file_id]

if (nrow(plink_qc_summary) == 0) {
  stop(
    "No rows in '", qc_summary_file, "' for pathway '",
    pathway_file_id, "' found. Please check whether the PLINK extraction ",
    "ran for this pathway and the name matches the ",
    "range file."
  )
}

print(plink_qc_summary)

# Output files
out_prefix <- file.path(
  repair_plink_dir,   # REPAIR (inputs ambu/md/qc_summary still from plink_dir)
  paste0("Working_Plink_Pathway_", pathway_name)
)

combined_prefix <- paste0(out_prefix, "_Combined")
after_r_qc_prefix <- paste0(out_prefix, "_after_R_QC")
ld_prune_prefix <- paste0(out_prefix, "_LDprune_r2_02")
ld_pruned_prefix <- paste0(after_r_qc_prefix, "_LDpruned")

numeric_geno_prefix <- paste0(out_prefix, "_Numeric_Geno")
numeric_geno_ld_prefix <- paste0(out_prefix, "_Numeric_Geno_LDpruned")

final_dataset_ld_prefix <- output_file(
  "Pathway_",
  pathway_name,
  "_Combined_QC_LDpruned"
)

final_dataset_unpruned_prefix <- output_file(
  "Pathway_",
  pathway_name,
  "_Combined_QC_noLDpruned"
)

vep_vcf_prefix <- paste0(ld_pruned_prefix, "_VCF")
vep_results_file <- orig_output_file("Pathway_", pathway_name, "_VEP.txt")   # REPAIR: generated externally, read only
# VEP file has already been generated externally from
# Working_Plink_Pathway_<Pathway>_after_R_QC_LDpruned_VCF.vcf generated

out_extended_prefix <- paste0(out_prefix, "_Extended")
out_core_prefix     <- paste0(out_prefix, "_Core")

#==================================================================================
# CREATE CASE-CONTROL DATASET FOR HWE
#==================================================================================
library(data.table)

ph <- fread(file.path(pheno_repair_dir, "BL_pheno_names_clean_CSV_UTF8_neu.csv"))   # REPAIR

ph[, BMI_numeric := as.numeric(gsub(",", ".", BMI))]
ph[, QC_control := ifelse(!is.na(BMI_numeric) & BMI_numeric < 30, 1, 0)]

ph_plink <- ph[, .(id, BMI_numeric, QC_control)]
ph_plink[, `:=`(FID = id, IID = id)]

old_wd <- getwd()
dir.create(file.path(repair_base, "Genetik_Daten_FoCus"), recursive = TRUE, showWarnings = FALSE)   # REPAIR
setwd(file.path(repair_base, "Genetik_Daten_FoCus"))   # REPAIR

fwrite(ph_plink[QC_control == 1, .(FID, IID)],
       "plink_controls_BMI.txt", col.names = FALSE, sep = "\t")

write.csv2(
  ph_plink,
  "pheno_plink_ready.csv",
  row.names = FALSE
)

setwd(old_wd)

#============================================================
# Check BMI availability
#============================================================
# to decide on further procedure: complete cases vs. imputation (>5 %)
n_total <- nrow(ph_plink)
n_available <- sum(!is.na(ph_plink$BMI_numeric))
n_missing <- sum(is.na(ph_plink$BMI_numeric))

cat("\n========================================\n")
cat("BMI availability\n")
cat("========================================\n")
cat("Total individuals: ", n_total, "\n")
cat("BMI available:         ", n_available, " (",
    round(100 * n_available / n_total, 1), "%)\n", sep = "")
cat("BMI missing:           ", n_missing, " (",
    round(100 * n_missing / n_total, 1), "%)\n", sep = "")

#=================================================================================
# GEMEINSAME SNPs (RAW) – PATHWAY
#=================================================================================
library(snpStats)

# Load PLINK datasets (all pathway genes together)
plink_data_ambu <- read.plink(
  bed = paste0(plink_ambu_prefix, ".bed"),
  bim = paste0(plink_ambu_prefix, ".bim"),
  fam = paste0(plink_ambu_prefix, ".fam")
)

plink_data_MD <- read.plink(
  bed = paste0(plink_md_prefix, ".bed"),
  bim = paste0(plink_md_prefix, ".bim"),
  fam = paste0(plink_md_prefix, ".fam")
)

# determine shared SNPs
common_snps_raw <- intersect(
  colnames(plink_data_ambu$genotypes),
  colnames(plink_data_MD$genotypes)
)

# Reduce datasets to shared SNPs
plink_data_ambu_common <- plink_data_ambu
plink_data_MD_common   <- plink_data_MD

plink_data_ambu_common$genotypes <-
  plink_data_ambu$genotypes[, common_snps_raw, drop = FALSE]

plink_data_MD_common$genotypes <-
  plink_data_MD$genotypes[, common_snps_raw, drop = FALSE]

#=================================================================================
# HWE / QC
#=================================================================================
controls <- ph[QC_control == 1, .(FID = id, IID = id)]
cat("Controls loaded:", nrow(controls), "\n")

hwe_filter <- function(plink_data, controls) {
  
  geno <- plink_data$genotypes
  fam  <- plink_data$fam
  
  control_idx <- which(fam$member %in% controls$IID)
  if (length(control_idx) == 0) stop("No controls found!")
  
  geno_controls <- geno[control_idx, , drop = FALSE]
  n_snps <- ncol(geno)
  
  hwe_thresh <- if (n_snps < 100) 0.01 else if (n_snps <= 1000) 0.001 else 1e-6
  
  pvals <- sapply(colnames(geno_controls), function(snp) {
    x <- as(geno_controls[, snp], "numeric")
    tab <- table(factor(x, levels = 0:2))
    if (sum(tab) == 0) return(NA)
    tryCatch(HWT.chisq(c(tab[1], tab[2], tab[3]), cc = 0)$pval,
             error = function(e) NA)
  })
  
  keep_snps <- names(pvals)[pvals >= hwe_thresh | is.na(pvals)]
  
  list(
    genotypes = geno[, keep_snps, drop = FALSE],
    fam       = fam,
    map       = plink_data$map[plink_data$map$snp.name %in% keep_snps, ]
  )
}

filtered_ambu <- hwe_filter(plink_data_ambu_common, controls)
filtered_MD   <- hwe_filter(plink_data_MD_common, controls)

#=================================================================================
# MERGE AFTER HWE – CUMULATIVE PATHWAY
#=================================================================================
# Filter shared SNPs by HWE
common_snps_hwe <- intersect(
  colnames(filtered_ambu$genotypes),
  colnames(filtered_MD$genotypes)
)

cat("Number of shared SNPs after HWE:", length(common_snps_hwe), "\n")

# Reduce genotypes to shared SNPs
geno_ambu_common <- filtered_ambu$genotypes[, common_snps_hwe, drop = FALSE]
geno_MD_common   <- filtered_MD$genotypes[, common_snps_hwe, drop = FALSE]

# Robustness: sort SNP order identically in both datasets
geno_MD_common <- geno_MD_common[, colnames(geno_ambu_common), drop = FALSE]

# Merge individuals
combined_geno <- rbind(geno_ambu_common, geno_MD_common)
combined_fam  <- rbind(filtered_ambu$fam, filtered_MD$fam)

# Check for duplicates
if (any(duplicated(rownames(combined_geno))))
  stop("ERROR: duplicate IDs in combined_geno!")

if (any(duplicated(combined_fam$member)))
  stop("ERROR: duplicate IDs in combined_fam!")

# Create map object for shared SNPs – align order with combined_geno
combined_map <- filtered_ambu$map[
  match(colnames(combined_geno), filtered_ambu$map$snp.name), 
]

# Save the combined map and genotype matrix
saveRDS(
  combined_map,
  paste0(combined_prefix, "_Map.rds")
)

saveRDS(list(
  genotypes = combined_geno,
  fam       = combined_fam,
  map       = combined_map
), paste0(combined_prefix, ".rds"))

cat("Combined PATHWAY SnpMatrix saved.\n")

#=================================================================================
# REPEATED QUALITY CONTROL AFTER MERGING (GENOTYPES ONLY!)
#=================================================================================
library(snpStats)
library(data.table)

# -------- Helper function: MAF
calc_maf <- function(geno_matrix) {
  x <- as(geno_matrix, "numeric")
  allele_counts <- colSums(x, na.rm = TRUE)
  n_alleles <- 2 * colSums(!is.na(x))
  pmin(allele_counts / n_alleles,
       1 - allele_counts / n_alleles)
}

### -------- 1. Individual missingness
indiv_missing <- rowMeans(is.na(as(combined_geno, "numeric")))
inds_to_keep <- rownames(combined_geno)[indiv_missing <= 0.10]

geno_qc <- combined_geno[inds_to_keep, ]
fam_qc  <- combined_fam[combined_fam$member %in% inds_to_keep, ]
fam_qc  <- fam_qc[match(rownames(geno_qc), fam_qc$member), ]

stopifnot(all(fam_qc$member == rownames(geno_qc)))

### -------- 2. SNP Missingness
snp_missing <- colMeans(is.na(as(geno_qc, "numeric")))
snps_to_keep <- names(snp_missing)[snp_missing <= 0.05]

geno_qc <- geno_qc[, snps_to_keep, drop = FALSE]
map_qc  <- combined_map[combined_map$snp.name %in% snps_to_keep, ]

### -------- 3. MAF Filter
maf <- calc_maf(geno_qc)
snps_keep <- names(maf)[maf >= 0.01]

geno_qc <- geno_qc[, snps_keep, drop = FALSE]
map_qc  <- map_qc[map_qc$snp.name %in% snps_keep, ]

### -------- 4. Save (GENOTYPES ONLY!)
saveRDS(
  list(
    genotypes = geno_qc,
    fam       = fam_qc,
    map       = map_qc
  ),
  paste0(out_prefix, "_GENO_QC.rds")
)

cat("QC genotyped data saved (without phenotypes).\n")

#=================================================================================
# QC OVERVIEW TABLE – GENETIC QC
#=================================================================================
pathway_summary <- as.data.table(plink_qc_summary)

#------------------------------------------------------------
# Original PLINK array size before pathway extraction
#------------------------------------------------------------

array_snps_ambu <- unique(
  pathway_summary[COHORT == "Ambu", SNP_COUNT_RAW]
)

array_inds_ambu <- unique(
  pathway_summary[COHORT == "Ambu", IND_COUNT_RAW]
)

array_snps_md <- unique(
  pathway_summary[COHORT == "MD", SNP_COUNT_RAW]
)

array_inds_md <- unique(
  pathway_summary[COHORT == "MD", IND_COUNT_RAW]
)

# Safety check: exactly one value expected per cohort
if (length(array_snps_ambu) != 1 ||
    length(array_inds_ambu) != 1 ||
    length(array_snps_md) != 1 ||
    length(array_inds_md) != 1) {
  
  stop(
    "The PLINK QC summary does not contain exactly one unique ",
    "SNP or individual value per cohort."
  )
}

array_snps_total <- array_snps_ambu + array_snps_md
array_inds_total <- array_inds_ambu + array_inds_md

#------------------------------------------------------------
# Create QC summary
#------------------------------------------------------------

qc_summary <- data.table(
  Pathway = pathway_name,
  
  Step = c(
    "Raw SNP Dataset (Ambu)",
    "Raw SNP Dataset (MD)",
    "Raw SNP Dataset (Combined)",
    "Initial Gene Dataset (Ambu)",
    "Initial Gene Dataset (MD)",
    "Initial Gene Dataset (Combined)",
    "Common SNPs",
    "QC: HWE",
    "QC: Individual Missingness",
    "QC: SNP Missingness",
    "QC: MAF",
    "Genetic QC Dataset"
  ),
  
  SNPs = c(
    array_snps_ambu,
    array_snps_md,
    array_snps_total,
    
    ncol(plink_data_ambu$genotypes),
    ncol(plink_data_MD$genotypes),
    ncol(plink_data_ambu$genotypes) +
      ncol(plink_data_MD$genotypes),
    
    length(common_snps_raw),
    length(common_snps_hwe),
    ncol(combined_geno),
    length(snps_to_keep),
    length(snps_keep),
    ncol(geno_qc)
  ),
  
  Individuals = c(
    array_inds_ambu,
    array_inds_md,
    array_inds_total,
    
    nrow(plink_data_ambu$genotypes),
    nrow(plink_data_MD$genotypes),
    nrow(plink_data_ambu$genotypes) +
      nrow(plink_data_MD$genotypes),
    
    nrow(plink_data_ambu$genotypes) +
      nrow(plink_data_MD$genotypes),
    
    nrow(combined_geno),
    length(inds_to_keep),
    nrow(geno_qc),
    nrow(geno_qc),
    nrow(geno_qc)
  )
)

cat("\n========================================\n")
cat("QC summary – genetic QC\n")
cat("========================================\n")

print(qc_summary)

#------------------------------------------------------------
# Save
#------------------------------------------------------------

saveRDS(
  qc_summary,
  output_file(
    "Pathway_",
    pathway_name,
    "_Combined_QC_Summary.rds"
  )
)

cat(
  "\nQC summary saved at:\n",
  output_file(
    "Pathway_",
    pathway_name,
    "_Combined_QC_Summary.rds"
  ),
  "\n"
)

#=================================================================================
# WRITE NEW PLINK FILE FROM R QC DATA (for LD pruning)
#=================================================================================
library(snpStats)

# geno_qc = SnpMatrix after HWE/missingness/MAF
# fam_qc  = passende Personeninfos
# map_qc  = passende SNP-Map

# Secure order
fam_qc <- fam_qc[match(rownames(geno_qc), fam_qc$member), ]
map_qc <- map_qc[match(colnames(geno_qc), map_qc$snp.name), ]

stopifnot(all(fam_qc$member == rownames(geno_qc)))
stopifnot(all(map_qc$snp.name == colnames(geno_qc)))

# Subject data for PLINK
subject_data <- data.frame(
  pedigree = as.character(rownames(geno_qc)),
  member   = as.character(rownames(geno_qc)),
  father   = "0",
  mother   = "0",
  sex      = as.numeric(fam_qc$sex),
  affected = -9,
  stringsAsFactors = FALSE
)

subject_data$sex[is.na(subject_data$sex)] <- 0
rownames(subject_data) <- rownames(geno_qc)

stopifnot(!any(is.na(subject_data$pedigree)))
stopifnot(!any(is.na(subject_data$member)))

fam_qc$pedigree <- fam_qc$member
fam_qc$father[is.na(fam_qc$father)] <- "0"
fam_qc$mother[is.na(fam_qc$mother)] <- "0"

snp_data <- data.frame(
  chromosome = as.integer(map_qc$chromosome),
  snp.name   = as.character(map_qc$snp.name),
  cM         = 0,
  position   = as.integer(map_qc$position),
  allele.1   = as.character(map_qc$allele.1),
  allele.2   = as.character(map_qc$allele.2),
  stringsAsFactors = FALSE
)

rownames(snp_data) <- colnames(geno_qc)

stopifnot(!any(is.na(snp_data$chromosome)))
stopifnot(!any(snp_data$chromosome == 0))
stopifnot(!any(is.na(snp_data$position)))
stopifnot(!any(snp_data$position == 0))
stopifnot(!any(is.na(snp_data$allele.1)))
stopifnot(!any(is.na(snp_data$allele.2)))

write.plink(
  file.base = after_r_qc_prefix,
  snps = geno_qc,
  subject.data = subject_data,
  snp.data = snp_data
)

bim_fixed <- data.frame(
  V1 = snp_data$chromosome,
  V2 = snp_data$snp.name,
  V3 = 0,
  V4 = snp_data$position,
  V5 = snp_data$allele.1,
  V6 = snp_data$allele.2
)

write.table(
  bim_fixed,
  file = paste0(after_r_qc_prefix, ".bim"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

head(read.table(paste0(after_r_qc_prefix, ".bim")))

stopifnot(file.exists(paste0(after_r_qc_prefix, ".bed")))
stopifnot(file.exists(paste0(after_r_qc_prefix, ".bim")))
stopifnot(file.exists(paste0(after_r_qc_prefix, ".fam")))

plink_after_qc_summary <- data.table(
  Datensatz = after_r_qc_prefix,
  Individuen = nrow(geno_qc),
  SNPs = ncol(geno_qc),
  BED_existiert = file.exists(paste0(after_r_qc_prefix, ".bed")),
  BIM_existiert = file.exists(paste0(after_r_qc_prefix, ".bim")),
  FAM_existiert = file.exists(paste0(after_r_qc_prefix, ".fam"))
)

print(plink_after_qc_summary)

#=================================================================================
# LD PRUNING MIT PLINK VIA R
#=================================================================================
plink <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Working_Plink/plink.exe")

if (!file.exists(plink)) {
  stop("PLINK not found. Please provide the full path, e.g. C:/plink/plink.exe")
}

system2(plink, c(
  "--bfile", after_r_qc_prefix,
  "--indep-pairwise", "50", "5", "0.2",
  "--out", ld_prune_prefix
))

stopifnot(file.exists(paste0(ld_prune_prefix, ".prune.in")))
  
system2(plink, c(
  "--bfile", after_r_qc_prefix,
  "--extract", paste0(ld_prune_prefix, ".prune.in"),
  "--make-bed",
  "--out", ld_pruned_prefix
))

stopifnot(file.exists(paste0(ld_pruned_prefix, ".bed")))

plink_ld <- read.plink(
  bed = paste0(ld_pruned_prefix, ".bed"),
  bim = paste0(ld_pruned_prefix, ".bim"),
  fam = paste0(ld_pruned_prefix, ".fam")
)

geno_ldpruned <- plink_ld$genotypes
geno_numeric_ldpruned <- as(geno_ldpruned, "numeric")
geno_numeric_ldpruned <- matrix(
  geno_numeric_ldpruned,
  nrow = nrow(geno_ldpruned),
  ncol = ncol(geno_ldpruned),
  dimnames = list(rownames(geno_ldpruned), colnames(geno_ldpruned))
)

ldpruned_snps <- data.table(
  SNP = colnames(geno_ldpruned)
)

removed_snps <- setdiff(colnames(geno_qc), colnames(geno_ldpruned))

ld_pruning_summary <- data.table(
  Schritt = c(
    "Vor LD-Pruning",
    "Nach LD-Pruning",
    "Entfernte SNPs"
  ),
  SNPs = c(
    ncol(geno_qc),
    ncol(geno_ldpruned),
    ncol(geno_qc) - ncol(geno_ldpruned)
  )
)

print(ld_pruning_summary)

cat("SNPs before LD pruning:", ncol(geno_qc), "\n")
cat("SNPs after LD pruning:", ncol(geno_ldpruned), "\n")
cat("SNPs removed by LD pruning:", ncol(geno_qc) - ncol(geno_ldpruned), "\n")

saveRDS(
  list(
    genotypes = geno_ldpruned,
    fam       = plink_ld$fam,
    map       = plink_ld$map
  ),
  paste0(out_prefix, "_GENO_QC_LDpruned.rds")
)

saveRDS(geno_numeric_ldpruned,
        paste0(numeric_geno_ld_prefix, ".rds"))


write.csv2(
  ldpruned_snps,
  output_file("Pathway_", pathway_name, "_LDpruned_SNPs.csv"),
  row.names = FALSE
)

write.csv2(
  data.table(SNP = removed_snps),
  output_file("Pathway_", pathway_name, "_LDpruned_Removed_SNPs.csv"),
  row.names = FALSE
)

write.csv2(
  ld_pruning_summary,
  output_file("Pathway_", pathway_name, "_LDpruning_Summary.csv"),
  row.names = FALSE
)

qc_summary <- readRDS(
  output_file("Pathway_", pathway_name, "_Combined_QC_Summary.rds")
)

qc_summary <- rbind(
  qc_summary,
  data.table(
    Pathway = pathway_name,
    Step = "LD-Pruning",
    SNPs = ncol(geno_ldpruned),
    Individuals = nrow(geno_ldpruned)
  )
)

saveRDS(
  qc_summary,
  output_file("Pathway_", pathway_name, "_Combined_QC_Summary.rds")
)

#=================================================================================
# SAVE GENO DATA AS NUMERIC MATRIX
#=================================================================================
library(snpStats)

# Create robust numeric matrix
geno_numeric <- as(geno_qc, "numeric")         # SnpMatrix → numeric
geno_numeric <- matrix(
  geno_numeric,
  nrow = nrow(geno_qc),
  ncol = ncol(geno_qc),
  dimnames = list(rownames(geno_qc), colnames(geno_qc))
)

cat("Dimension of the numeric genotype matrix:", dim(geno_numeric)[1], "individuals x", dim(geno_numeric)[2], "SNPs\n")

geno_counts <- table(factor(unlist(geno_numeric), levels = c(0,1,2,NA)))
print(geno_counts)

saveRDS(geno_numeric, paste0(numeric_geno_prefix, ".rds"))

cat("Numeric genotype matrix saved.\n")

#================================================================================
# MERGE PHENO + SNP DATA – FINAL DATASET
#================================================================================
library(data.table)
library(snpStats)

# --- 1. Load phenotype data --- 
# Load and merge pheno data 
pheno1 <- read.csv2(file.path(pheno_repair_dir, "BL_pheno_names_clean_CSV_UTF8_neu.csv"), stringsAsFactors = FALSE)   # REPAIR
pheno2 <- read.csv2(file.path(pheno_repair_dir, "FO_PA_alldata_CSV_UTF8_neu.csv"), stringsAsFactors = FALSE)   # REPAIR
pheno3 <- read.csv2(file.path(pheno_repair_dir, "FO_nutri_mima_CSV_UTF8_neu.csv"), stringsAsFactors = FALSE)   # REPAIR

stopifnot("id" %in% colnames(pheno1))
stopifnot("ID_2020" %in% colnames(pheno2))

setnames(pheno2, old = "ID_2020", new = "id")

possible_id_cols <- c("id", "ID", "ID_2020", "Id", "subject_id")
id_col_pheno3 <- intersect(possible_id_cols, colnames(pheno3))

if (length(id_col_pheno3) != 1) {
  stop("No unique ID column found in pheno3: ",
       paste(colnames(pheno3), collapse = ", "))
}

setnames(pheno3, old = id_col_pheno3, new = "id")

pheno1$id <- as.character(pheno1$id)
pheno2$id <- as.character(pheno2$id)
pheno3$id <- as.character(pheno3$id)

pheno <- Reduce(
  function(x, y) merge(x, y, by = "id", all = TRUE),
  list(pheno1, pheno2, pheno3)
)

stopifnot("id" %in% colnames(pheno))

# --- 2. Check global covariate availability
pheno <- as.data.table(pheno)
cov_check <- copy(pheno)

cov_check[, gender_num := as.numeric(as.character(gender))]
cov_check[gender_num %in% c(0, 7777777, 9999999), gender_num := NA]
cov_check[, gender_cat := factor(gender_num, levels = c(1, 2),
                                 labels = c("male", "female"))]

cov_check[, smoking_habits_3m_num := as.numeric(as.character(smoking_habits_3m))]
cov_check[smoking_habits_3m_num %in% c(0, 7777777, 9999999),
          smoking_habits_3m_num := NA]

cov_check[, smoking_cat := factor(
  dplyr::case_when(
    smoking_habits_3m_num %in% c(1, 2) ~ "never_smoker",
    smoking_habits_3m_num == 3 ~ "former_smoker",
    smoking_habits_3m_num == 4 ~ "current_smoker",
    TRUE ~ NA_character_
  ),
  levels = c("never_smoker", "former_smoker", "current_smoker")
)]

cov_check[, school_certificate_A_num := as.numeric(as.character(school_certificate_A))]
cov_check[!school_certificate_A_num %in% 1:5, school_certificate_A_num := NA]
cov_check[, school_cert_cat := factor(
  ifelse(school_certificate_A_num %in% c(1, 2, 5), "risk", "reference"),
  levels = c("reference", "risk")
)]

cov_check[, graduation_num := as.numeric(as.character(graduation))]
cov_check[graduation_num %in% c(6, 7777777, 9999999), graduation_num := NA]
cov_check[, graduation_cat := factor(
  ifelse(graduation_num %in% c(1, 2, 5), "risk", "reference"),
  levels = c("reference", "risk")
)]

cov_check[, sample_age := as.numeric(as.character(sample_age))]
cov_check[, sports_activity_min_week := as.numeric(as.character(sports_activity_min_week))]
cov_check[, everyday_activity_min_week := as.numeric(as.character(everyday_activity_min_week))]

if ("energy_kcal" %in% colnames(cov_check)) {
  cov_check[, energy_kcal := as.numeric(gsub(",", ".", as.character(energy_kcal)))]
} else if ("energy" %in% colnames(cov_check)) {
  cov_check[, energy_kJ := as.numeric(gsub(",", ".", as.character(energy)))]
  cov_check[, energy_kcal := energy_kJ / 4.184]
} else {
  cov_check[, energy_kcal := NA_real_]
}

covariates_keep_global <- c(
  "gender_cat",
  "smoking_cat",
  "school_cert_cat",
  "graduation_cat",
  "sample_age",
  "sports_activity_min_week",
  "everyday_activity_min_week",
  "energy_kcal"
)

n_total <- nrow(cov_check)

cov_missing_summary <- data.frame(
  Variable = covariates_keep_global,
  N_available = sapply(covariates_keep_global, function(v) sum(!is.na(cov_check[[v]]))),
  Percent_available = round(100 * sapply(covariates_keep_global, function(v) mean(!is.na(cov_check[[v]]))), 1),
  N_missing = sapply(covariates_keep_global, function(v) sum(is.na(cov_check[[v]]))),
  Percent_missing = round(100 * sapply(covariates_keep_global, function(v) mean(is.na(cov_check[[v]]))), 1)
)

print(cov_missing_summary)

complete_covariates_global <- complete.cases(
  cov_check[, ..covariates_keep_global]
)

cat("\n========================================\n")
cat("Global complete cases covariates\n")
cat("========================================\n")
cat("Total individuals:          ", n_total, "\n")
cat("All covariates available:      ", sum(complete_covariates_global),
    " (", round(100 * mean(complete_covariates_global), 1), "%)\n", sep = "")
cat("At least one covariate missing:", sum(!complete_covariates_global),
    " (", round(100 * mean(!complete_covariates_global), 1), "%)\n", sep = "")

write.csv2(
  cov_missing_summary,
  output_file("Pathway_", pathway_name, "_Global_Covariates_Missingness_Summary.csv"),
  row.names = FALSE
)

# --- 3. Batch assignment ---
plink_data_ambu <- read.plink(
  bed = paste0(plink_ambu_prefix, ".bed"),
  bim = paste0(plink_ambu_prefix, ".bim"),
  fam = paste0(plink_ambu_prefix, ".fam")
)

plink_data_MD <- read.plink(
  bed = paste0(plink_md_prefix, ".bed"),
  bim = paste0(plink_md_prefix, ".bim"),
  fam = paste0(plink_md_prefix, ".fam")
)

ids_ambu <- plink_data_ambu$fam$member
ids_MD   <- plink_data_MD$fam$member

pheno$Batch <- NA
pheno$Batch[pheno$id %in% ids_ambu] <- "Ambu"
pheno$Batch[pheno$id %in% ids_MD]   <- "MD"

cat("Batch distribution:\n")
print(table(pheno$Batch))

# --- 4. Load genotype data ---
geno_numeric <- readRDS(paste0(numeric_geno_ld_prefix, ".rds"))
geno_ids <- rownames(geno_numeric)

# --- 5. Matching IDs ---
common_ids <- intersect(pheno$id, geno_ids)

cat("Shared IDs:", length(common_ids), "\n")

geno_final  <- geno_numeric[common_ids, , drop = FALSE]
pheno_final <- pheno[match(common_ids, pheno$id), ]

rownames(pheno_final) <- pheno_final$id

# --- 6. Check consistency ---
if (!all(rownames(geno_final) == rownames(pheno_final))) {
  stop("Error: ID order does not match!")
}

#--------------------------------------------------------
# Geno/phenotype ID matching
#--------------------------------------------------------

n_pheno <- nrow(pheno)
n_matched <- length(common_ids)

cat("\n========================================\n")
cat("Genotype/phenotype matching\n")
cat("========================================\n")
cat("Phenotypes total:      ", n_pheno, "\n")
cat("Matched individuals:   ", n_matched,
    " (", round(100 * n_matched / n_pheno, 1), "%)\n", sep = "")
cat("Not matched:         ", n_pheno - n_matched,
    " (", round(100 * (n_pheno - n_matched) / n_pheno, 1), "%)\n", sep = "")

#--------------------------------------------------------
# Definition obesity + Complete Cases
#--------------------------------------------------------

if ("BMI" %in% colnames(pheno_final)) {
  pheno_final$BMI <- as.numeric(gsub(",", ".", pheno_final$BMI))
  
  pheno_final$obesity <- ifelse(
    is.na(pheno_final$BMI),
    NA,
    ifelse(pheno_final$BMI >= 30, 1, 0)
  )
} else {
  stop("BMI missing! 'obesity' cannot be computed.")
}

n_before_obesity <- nrow(pheno_final)

# Keep only complete cases for obesity
pheno_final <- pheno_final %>%
  filter(!is.na(obesity))

geno_final <- geno_final[
  rownames(geno_final) %in% pheno_final$id,
  , drop = FALSE
]

pheno_final <- pheno_final[
  match(rownames(geno_final), pheno_final$id),
]

rownames(pheno_final) <- pheno_final$id

stopifnot(all(rownames(geno_final) == rownames(pheno_final)))

n_after_obesity <- nrow(pheno_final)

cat("\n========================================\n")
cat("Complete Cases obesity\n")
cat("========================================\n")
cat("Matched individuals:   ", n_before_obesity, "\n")
cat("Complete Cases obesity: ", n_after_obesity,
    " (", round(100 * n_after_obesity / n_before_obesity, 1), "%)\n", sep = "")
cat("Excluded:         ", n_before_obesity - n_after_obesity,
    " (", round(100 * (n_before_obesity - n_after_obesity) / n_before_obesity, 1), "%)\n", sep = "")

# --- 7. SAVE FINAL DATASET ---
final_data <- list(
  genotypes = geno_final,
  pheno     = pheno_final
)

saveRDS(final_data, paste0(final_dataset_ld_prefix, ".rds"))

cat(
  "Finaler Datensatz gespeichert:",
  paste0(final_dataset_ld_prefix, ".rds"),
  "\n"
)

# --- 8. SAVE PHENO AS CSV (CONSISTENTLY NAMED) ---
pheno_export <- pheno_final
pheno_export$id <- rownames(pheno_export)

write.csv2(
  pheno_export,
  file = output_file("Pathway_", pathway_name, "_Combined_QC_LDpruned_pheno.csv"),
  row.names = FALSE
)

cat("Pheno dataset saved: Working_Plink_Pathway_*_Combined_QC_pheno.csv\n")

qc_summary <- readRDS(
  output_file("Pathway_", pathway_name, "_Combined_QC_Summary.rds")
)

analysis_snps <- ncol(geno_final)

qc_summary_extended <- rbind(
  qc_summary,
  data.table(
    Pathway = pathway_name,
    Step = "Pheno-Matching",
    SNPs = analysis_snps,
    Individuals = n_matched
  ),
  data.table(
    Pathway = pathway_name,
    Step = "Complete Cases (Obesity)",
    SNPs = analysis_snps,
    Individuals = n_after_obesity
  ),
  data.table(
    Pathway = pathway_name,
    Step = "Final Analysis Dataset",
    SNPs = analysis_snps,
    Individuals = n_after_obesity
  )
)

print(qc_summary_extended)

write.csv2(
  qc_summary_extended,
  output_file("Pathway_", pathway_name, "_Combined_QC_Summary.csv"),
  row.names = FALSE
)

#=================================================================================
# QC FLOWCHARTS
#=================================================================================
library(data.table)
library(dplyr)
library(ggplot2)
library(grid)

plot_flow_summary <- function(summary_df, pathway_name, outfile_prefix, show_plot = TRUE) {
  
  df <- as.data.table(summary_df)
  
  get_row <- function(step_name) {
    
    tmp <- df[Step == step_name]
    
    if (nrow(tmp) == 0) {
      stop("Step not found in QC summary: ", step_name)
    }
    
    tmp[1]
  }
  
  raw_ambu  <- get_row("Raw SNP Dataset (Ambu)")
  raw_md    <- get_row("Raw SNP Dataset (MD)")
  
  init_ambu <- get_row("Initial Gene Dataset (Ambu)")
  init_md   <- get_row("Initial Gene Dataset (MD)")
  init_comb <- get_row("Initial Gene Dataset (Combined)")
  
  rest_steps <- df[
    !Step %in% c(
      "Raw SNP Dataset (Ambu)",
      "Raw SNP Dataset (MD)",
      "Raw SNP Dataset (Combined)",
      "Initial Gene Dataset (Ambu)",
      "Initial Gene Dataset (MD)",
      "Initial Gene Dataset (Combined)"
    )
  ]
  
  df_plot <- rbindlist(
    list(
      
      # Separate Ambu/MD boxes + combined box
      data.table(
        Step_clean = c(
          "Raw SNP Dataset – Ambu",
          "Raw SNP Dataset – Pop",
          "Initial Gene Dataset – Ambu",
          "Initial Gene Dataset – Pop",
          "Initial Gene Dataset – Combined"
        ),
        
        Individuals = c(
          raw_ambu$Individuals,
          raw_md$Individuals,
          init_ambu$Individuals,
          init_md$Individuals,
          init_comb$Individuals
        ),
        
        SNPs = c(
          raw_ambu$SNPs,
          raw_md$SNPs,
          init_ambu$SNPs,
          init_md$SNPs,
          init_comb$SNPs
        ),
        
        x = c(
          0.65, 1.35,
          0.65, 1.35,
          1
        ),
        
        y = c(
          12, 12,
          11, 11,
          10
        ),
        
        Fill = "#DDEBF7"
      ),
      
      # Shared steps from common SNPs
      rest_steps[, .(
        Step_clean = Step,
        Individuals,
        SNPs
      )][
        , `:=`(
          x = 1,
          y = seq(
            from = 9,
            by = -1,
            length.out = .N
          ),
          
          Fill = dplyr::case_when(
            
            Step_clean == "Common SNPs" ~ "#DDEBF7",
            
            Step_clean %in% c(
              "QC: HWE",
              "QC: Individual Missingness",
              "QC: SNP Missingness",
              "QC: MAF",
              "QC: MAC",
              "QC: MAF >= 0.01",
              "QC: MAC >= 3",
              "Genetic QC Dataset",
              "LD-Pruning"
            ) ~ "#E2F0D9",
            
            TRUE ~ "#FFF3E6"
          )
        )
      ]
    ),
    fill = TRUE
  )
  
  df_plot[, Label := paste0(
    Step_clean,
    "\n",
    "Individuals: ", Individuals,
    "    SNPs: ", SNPs
  )]
  
  box_width  <- 0.62
  box_height <- 0.52
  
  p <- ggplot(df_plot, aes(x = x, y = y)) +
    
    # Raw Ambu/MD -> Initial Gene Dataset Ambu/MD
    geom_segment(
      data = data.frame(
        x = c(0.65, 1.35),
        xend = c(0.65, 1.35),
        y = c(11.70, 11.70),
        yend = c(11.31, 11.31)
      ),
      aes(
        x = x,
        xend = xend,
        y = y,
        yend = yend
      ),
      arrow = arrow(
        length = unit(0.14, "cm"),
        type = "closed"
      ),
      linewidth = 0.35,
      color = "grey35",
      inherit.aes = FALSE
    ) +
    
    # Initial Gene Dataset Ambu/MD -> Combined
    geom_segment(
      data = data.frame(
        x = c(0.65, 1.35),
        xend = c(1, 1),
        y = c(10.70, 10.70),
        yend = c(10.31, 10.31)
      ),
      aes(
        x = x,
        xend = xend,
        y = y,
        yend = yend
      ),
      arrow = arrow(
        length = unit(0.14, "cm"),
        type = "closed"
      ),
      linewidth = 0.35,
      color = "grey35",
      inherit.aes = FALSE
    ) +
    
    # Combined and all further shared steps
    geom_segment(
      data = df_plot[x == 1 & y > min(y)],
      aes(
        x = 1,
        xend = 1,
        y = y - box_height / 2 - 0.04,
        yend = y - 1 + box_height / 2 + 0.04
      ),
      arrow = arrow(
        length = unit(0.14, "cm"),
        type = "closed"
      ),
      linewidth = 0.35,
      color = "grey35"
    ) +
    
    geom_rect(
      aes(
        xmin = x - box_width / 2,
        xmax = x + box_width / 2,
        ymin = y - box_height / 2,
        ymax = y + box_height / 2,
        fill = Fill
      ),
      color = "grey25",
      linewidth = 0.4
    ) +
    
    geom_text(
      aes(label = Label),
      size = 3.0,
      lineheight = 0.9
    ) +
    
    scale_fill_identity() +
    
    labs(
      title = paste0(
        pathway_name,
        " – QC Flowchart"
      )
    ) +
    
    coord_cartesian(
      xlim = c(0.2, 1.8),
      ylim = c(
        min(df_plot$y) - 0.6,
        12.6
      ),
      clip = "off"
    ) +
    
    theme_void(base_size = 13) +
    
    theme(
      plot.title = element_text(
        face = "plain",
        size = 13,
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      plot.margin = margin(15, 20, 15, 20)
    )
  
  plot_height <- max(
    7,
    0.72 * nrow(df_plot)
  )
  
  ggsave(
    paste0(outfile_prefix, ".pdf"),
    plot = p,
    width = 7.5,
    height = plot_height,
    device = "pdf"
  )
  
  ggsave(
    paste0(outfile_prefix, ".png"),
    plot = p,
    width = 7.5,
    height = plot_height,
    dpi = 600
  )
  
  if (show_plot) {
    print(p)
  }
  
  return(p)
}

qc_summary <- read.csv2(
  output_file("Pathway_", pathway_name, "_Combined_QC_Summary.csv"),
  stringsAsFactors = FALSE
)

qc_flow_plot <- plot_flow_summary(
  summary_df = qc_summary,
  pathway_name = pathway_name,
  outfile_prefix = output_file("Pathway_", pathway_name, "_QC_Flowchart"),
  show_plot = TRUE
)

#================================================================================
# MERGE PHENO + SNP DATA – FINAL DATASET WITHOUT LD PRUNING
#================================================================================
geno_numeric_noLD <- readRDS(
  paste0(numeric_geno_prefix, ".rds")
)
geno_ids_noLD <- rownames(geno_numeric_noLD)

common_ids_noLD <- intersect(pheno$id, geno_ids_noLD)

geno_final_noLD  <- geno_numeric_noLD[common_ids_noLD, , drop = FALSE]
pheno_final_noLD <- pheno[match(common_ids_noLD, pheno$id), ]

rownames(pheno_final_noLD) <- pheno_final_noLD$id

stopifnot(all(rownames(geno_final_noLD) == rownames(pheno_final_noLD)))

#--------------------------------------------------------
# Definition of obesity + complete cases – without LD pruning
#--------------------------------------------------------

if ("BMI" %in% colnames(pheno_final_noLD)) {
  pheno_final_noLD$BMI <- as.numeric(gsub(",", ".", pheno_final_noLD$BMI))
  
  pheno_final_noLD$obesity <- ifelse(
    is.na(pheno_final_noLD$BMI),
    NA,
    ifelse(pheno_final_noLD$BMI >= 30, 1, 0)
  )
} else {
  stop("BMI missing! 'obesity' cannot be computed.")
}

n_before_obesity_noLD <- nrow(pheno_final_noLD)

pheno_final_noLD <- pheno_final_noLD %>%
  filter(!is.na(obesity))

geno_final_noLD <- geno_final_noLD[
  rownames(geno_final_noLD) %in% pheno_final_noLD$id,
  , drop = FALSE
]

pheno_final_noLD <- pheno_final_noLD[
  match(rownames(geno_final_noLD), pheno_final_noLD$id),
]

rownames(pheno_final_noLD) <- pheno_final_noLD$id

stopifnot(all(rownames(geno_final_noLD) == rownames(pheno_final_noLD)))

n_after_obesity_noLD <- nrow(pheno_final_noLD)

cat("\n========================================\n")
cat("Complete cases obesity without LD pruning\n")
cat("========================================\n")
cat("Matched individuals:   ", n_before_obesity_noLD, "\n")
cat("Complete Cases obesity: ", n_after_obesity_noLD,
    " (", round(100 * n_after_obesity_noLD / n_before_obesity_noLD, 1), "%)\n", sep = "")
cat("Excluded:         ", n_before_obesity_noLD - n_after_obesity_noLD,
    " (", round(100 * (n_before_obesity_noLD - n_after_obesity_noLD) / n_before_obesity_noLD, 1), "%)\n", sep = "")

final_data_noLD <- list(
  genotypes = geno_final_noLD,
  pheno     = pheno_final_noLD
)

saveRDS(
  final_data_noLD,
  paste0(final_dataset_unpruned_prefix, ".rds")
)

write.csv2(
  pheno_final_noLD,
  output_file(
    "Pathway_",
    pathway_name,
    "_Combined_QC_noLDpruned_pheno.csv"
  ),
  row.names = FALSE
)

cat("✅ Final dataset WITHOUT LD pruning saved.\n")

#=================================================================================
# VCF + VEP ANNOTATION (EXTERNAL STEP)
#=================================================================================
#
# Basis:
# Working_Plink_Pathway_[pathway]_after_R_QC_LDpruned.bed/.bim/.fam
#
# These LD-pruned PLINK files are converted outside R using:
#
# plink --bfile Working_Plink_Pathway_[pathway]_after_R_QC_LDpruned \
#       --recode vcf \
#       --out Working_Plink_Pathway_[pathway]_after_R_QC_LDpruned_VCF
#
# into a VCF file.
#
# The VCF file is subsequently annotated in Ensembl VEP
# (GRCh37 / hg19).
#
# Result:
# Pathway_[pathway]_VEP.txt
#
# This file is later used for the sensitivity analysis
# "functional_vep" read:
#
# sensitivity_functional_vep
# -> load_functional_vep_snps()
#
# Important:
# The functional sensitivity analysis is based on the
# LD-pruned SNPs and the VEP file generated from them.
#
#=================================================================================

#==================================================================================
# BREAKPOINT PHASE 1 (BATCH-LAUF)
#==================================================================================
# With run_phase == "1" the script deliberately stops here: the LD-pruned PLINK
# files (.bed/.bim/.fam) for this pathway have now been written; the
# next step (VCF conversion + Ensembl VEP annotation) takes place
# outside R. After this external step for ALL pathways,
# the script is restarted with run_phase == "2" (or without a 3rd argument)
# and then continues normally from here.
if (identical(run_phase, "1")) {
  cat(
    "Phase 1 completed for pathway '", pathway_display_name,
    "' (", pathway_name, "). LD-pruned PLINK files are located at:\n",
    "  ", ld_pruned_prefix, ".bed/.bim/.fam\n",
    "Next step (outside R): VCF conversion + Ensembl VEP ",
    "annotation, result expected at:\n",
    "  ", vep_results_file, "\n",
    "Then restart this script with run_phase = 2 for the same pathway ",
    "starten.\n",
    sep = ""
  )

  log_pipeline_run(
    md_path       = methods_log_path,
    pathway_name  = pathway_name,
    phase         = "1 (bis LD-Pruning, vor externem VCF/VEP-Schritt)",
    plink_version = extract_plink_version(paste0(ld_pruned_prefix, ".log")),
    extra_notes   = paste0(
      "SNPs nach QC: ", ncol(geno_qc),
      "; nach LD-Pruning: ", ncol(geno_ldpruned)
    )
  )

  quit(save = "no", status = 0)
}

#================================================================================
# GLOBAL TRAIN/TEST SPLIT (SHARED WITH THE GENE SCRIPT)
#================================================================================
# IMPORTANT - mandatory naming convention: the split is shared between the gene and
# pathway scripts. The CANONICAL file name is necessarily based on the
# OFFICIAL pathway name from the pathway selection (pathway_file_id, see above
# - derived from pathway_display_name / column "Pathway-Name" in
# Pathway_Selection_Supplement_GeneSets.csv) - NOT on the freely chosen
# short label pathway_name (only for this script's other output files,
# see pathway_name definition further above). In the gene script the
# for this pathway_name_gene_selection CHARACTER-EXACT (upper/lower case,
# special characters, spaces) must match the same official pathway name
# - otherwise the split will not be found here.

library(caret)

global_split_canonical_file <- output_file(
  "Pathway_", pathway_file_id, "_Global_TrainTest_Split.rds"
)

# Legacy/compatibility path: all approx. 9 existing read accesses further
# down in this script use "pathway_name" (not pathway_file_id) in the
# file name - an identical copy is additionally stored under this name here,
# so that nothing had to be changed at those places.
global_split_legacy_file <- output_file(
  "Pathway_", pathway_name, "_Global_TrainTest_Split.rds"
)

# REPAIR: the global split is read exclusively from the ORIGINAL Outputs_Pathway
# and NEVER regenerated (identical split as in the original analysis).
# The original regeneration logic has been deliberately removed here.
global_split_orig_file <- orig_output_file(
  "Pathway_", pathway_file_id, "_Global_TrainTest_Split.rds"
)

if (!file.exists(global_split_orig_file)) {
  stop(
    "REPAIR: global split not found at ", global_split_orig_file,
    " - aborting (split is not regenerated in the REPAIR run)."
  )
}

global_split <- readRDS(global_split_orig_file)

cat("✅ REPAIR: global split (original, unchanged) for pathway '",
    pathway_display_name, "' geladen:\n", sep = "")
cat(" File:", global_split_orig_file, "\n")
cat(" Train:", length(global_split$train_ids), "\n")
cat(" Test: ", length(global_split$test_ids), "\n")

# Store a copy in Repair/Outputs_Pathway (canonical name) so that all
# downstream read accesses (output_file) work unchanged.
saveRDS(global_split, global_split_canonical_file)

# Store an identical copy under the legacy/pathway_name file name (see above) -
# from here on all existing places in the script read unchanged as before.
saveRDS(global_split, global_split_legacy_file)

write.csv2(
  data.frame(
    id = c(global_split$train_ids, global_split$test_ids),
    Set = c(
      rep("train", length(global_split$train_ids)),
      rep("test", length(global_split$test_ids))
    )
  ),
  output_file("Pathway_", pathway_name, "_Global_TrainTest_Split.csv"),
  row.names = FALSE
)

#================================================================================
# TRAIN/TEST SPLIT FLOWCHART
#================================================================================
library(ggplot2)
library(grid)

plot_traintest_flow <- function(global_split, pathway_name, outfile_prefix, show_plot = TRUE) {
  
  n_train <- global_split$n_train
  n_test  <- global_split$n_test
  n_total <- n_train + n_test
  
  box_width  <- 0.85
  box_height <- 0.62
  
  col_train  <- 0.4
  col_test   <- 1.6
  col_center <- 1.0
  
  #---------------------------------------------------------------------------
  # Boxes
  #---------------------------------------------------------------------------
  
  box_df <- data.frame(
    Step = c(
      "Final Analysis Dataset",
      "Global Train/Test Split",
      "Training Set",
      "Test Set",
      "Univariate SNP Regression",
      "Weight/Direction Estimation",
      "SNP Sum Score (Train)",
      "SNP Sum Score (Test)",
      "Model + ROC (Train)",
      "ROC & Performance",
      "Final Performance (Test Set)"
    ),
    Label = c(
      paste0(
        "Final Analysis Dataset\n",
        "(Pathway, LD-pruned, complete cases)\n",
        "N = ", n_total
      ),
      paste0(
        "Global Train/Test Split\n",
        "70 / 30, stratified by obesity, seed = ", global_split$seed
      ),
      paste0("Training Set (~70%)\nn = ", n_train),
      paste0("Test Set (~30%)\nn = ", n_test, "  \u00b7  held out"),
      paste0(
        "Univariate SNP Regression\n",
        "Mean-imputed genotypes\n",
        "Direction (Risk/Protective)"
      ),
      paste0(
        "Weight/Direction Estimation\n",
        "Raw . Oriented . Weighted"
      ),
      "SNP Sum Score",
      "SNP Sum Score\n(Test, using training weights)",
      paste0(
        "Model + ROC\n",
        "Obesity ~ Score (+ covariates)\n",
        "Youden cutoff"
      ),
      paste0(
        "ROC & Performance\n",
        "Youden cutoff"
      ),
      paste0(
        "Final Performance (Test Set)\n",
        "AUC . OR per SD . Sensitivity/Specificity\n",
        "Youden cutoff"
      )
    ),
    x = c(
      col_center, col_center,
      col_train, col_test,
      col_train, col_train,
      col_train, col_test,
      col_train, col_test,
      col_center
    ),
    y = c(
      8, 7,
      6, 6,
      5, 4,
      3, 3,
      2, 2,
      1
    ),
    stringsAsFactors = FALSE
  )
  
  box_df$Fill <- dplyr::case_when(
    box_df$Step %in% c(
      "Final Analysis Dataset",
      "Training Set",
      "Test Set"
    ) ~ "#DDEBF7",
    
    box_df$Step %in% c(
      "Global Train/Test Split",
      "Univariate SNP Regression",
      "Weight/Direction Estimation"
    ) ~ "#FFF3E6",
    
    box_df$Step %in% c(
      "SNP Sum Score (Train)",
      "SNP Sum Score (Test)",
      "Final Performance (Test Set)"
    ) ~ "#FCE4D6",
    
    box_df$Step %in% c(
      "Model + ROC (Train)",
      "ROC & Performance"
    ) ~ "#E2F0D9",
    
    TRUE ~ "#FFF3E6"
  )
  
  #---------------------------------------------------------------------------
  # Held-out box (Test genotypes, spans two grid rows in the test column)
  #---------------------------------------------------------------------------
  
  heldout_df <- data.frame(
    xmin = col_test - box_width / 2,
    xmax = col_test + box_width / 2,
    ymin = 3 + box_height / 2 + 0.10,
    ymax = 6 - box_height / 2 - 0.10
  )
  
  heldout_label <- "Test genotypes\nremain unchanged\nno fitting, no label preview"
  
  #---------------------------------------------------------------------------
  # Straight arrows (explicit, one named row each - safer than parallel vectors)
  #---------------------------------------------------------------------------
  
  seg <- function(x, xend, y, yend, style = "solid") {
    data.frame(x = x, xend = xend, y = y, yend = yend, Style = style)
  }
  
  edge_gap <- 0.02
  
  arrow_df <- rbind(
    # Final Dataset -> Split
    seg(col_center, col_center,
        8 - box_height / 2 - edge_gap, 7 + box_height / 2 + edge_gap),
    
    # Split -> Training Set / Test Set
    seg(col_center, col_train,
        7 - box_height / 2 - edge_gap, 6 + box_height / 2 + edge_gap),
    seg(col_center, col_test,
        7 - box_height / 2 - edge_gap, 6 + box_height / 2 + edge_gap),
    
    # Training branch: Training Set -> Univariate Regression -> Weight Estimation
    # -> SNP Sum Score (Train) -> Model + ROC (Train)
    seg(col_train, col_train,
        6 - box_height / 2 - edge_gap, 5 + box_height / 2 + edge_gap),
    seg(col_train, col_train,
        5 - box_height / 2 - edge_gap, 4 + box_height / 2 + edge_gap),
    seg(col_train, col_train,
        4 - box_height / 2 - edge_gap, 3 + box_height / 2 + edge_gap),
    seg(col_train, col_train,
        3 - box_height / 2 - edge_gap, 2 + box_height / 2 + edge_gap),
    
    # Test branch: Test Set -> held-out genotypes -> SNP Sum Score (Test)
    # -> ROC & Performance
    seg(col_test, col_test,
        6 - box_height / 2 - edge_gap, heldout_df$ymax + edge_gap),
    seg(col_test, col_test,
        heldout_df$ymin - edge_gap, 3 + box_height / 2 + edge_gap),
    seg(col_test, col_test,
        3 - box_height / 2 - edge_gap, 2 + box_height / 2 + edge_gap),
    
    # Convergence into Final Performance
    # Training-side ROC is diagnostic only -> faint, dashed
    seg(col_train, col_center,
        2 - box_height / 2 - edge_gap, 1 + box_height / 2 + 0.08, style = "faint"),
    # Test-side ROC is the primary, reported result -> solid
    seg(col_test, col_center,
        2 - box_height / 2 - edge_gap, 1 + box_height / 2 + edge_gap)
  )
  
  #---------------------------------------------------------------------------
  # Curved accent arrows (train-derived information applied to the test branch)
  #---------------------------------------------------------------------------
  
  curve_df <- data.frame(
    x     = c(col_train + box_width / 2, col_train + box_width / 2),
    y     = c(4, 2),
    xend  = c(col_test - box_width / 2, col_test - box_width / 2),
    yend  = c(3.15, 2),
    label = c("weights applied", "model + cutoff applied"),
    ly    = c(3.55, 2.35)
  )
  
  #---------------------------------------------------------------------------
  # Plot
  #---------------------------------------------------------------------------
  
  p <- ggplot() +
    
    geom_rect(
      data = heldout_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "#F2F2F2",
      color = "#9AA6B2",
      linewidth = 0.4,
      linetype = "dashed"
    ) +
    
    geom_text(
      data = data.frame(
        x = col_test,
        y = mean(c(heldout_df$ymin, heldout_df$ymax)),
        label = heldout_label
      ),
      aes(x = x, y = y, label = label),
      size = 3.0,
      lineheight = 0.9,
      color = "grey35"
    ) +
    
    geom_segment(
      data = arrow_df[arrow_df$Style == "solid", ],
      aes(x = x, xend = xend, y = y, yend = yend),
      arrow = arrow(length = unit(0.14, "cm"), type = "closed"),
      linewidth = 0.35,
      color = "grey35"
    ) +
    
    geom_segment(
      data = arrow_df[arrow_df$Style == "faint", ],
      aes(x = x, xend = xend, y = y, yend = yend),
      arrow = arrow(length = unit(0.14, "cm"), type = "closed"),
      linewidth = 0.3,
      linetype = "dashed",
      color = "#9AA6B2"
    ) +
    
    geom_curve(
      data = curve_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      curvature = -0.3,
      linewidth = 0.45,
      linetype = "dashed",
      color = "#B45309",
      arrow = arrow(length = unit(0.14, "cm"), type = "closed")
    ) +
    
    geom_text(
      data = curve_df,
      aes(x = (x + xend) / 2, y = ly, label = label),
      size = 2.9,
      fontface = "bold",
      color = "#B45309"
    ) +
    
    geom_rect(
      data = box_df,
      aes(
        xmin = x - box_width / 2,
        xmax = x + box_width / 2,
        ymin = y - box_height / 2,
        ymax = y + box_height / 2,
        fill = Fill
      ),
      color = "grey25",
      linewidth = 0.4
    ) +
    
    geom_text(
      data = box_df,
      aes(x = x, y = y, label = Label),
      size = 3.0,
      lineheight = 0.9
    ) +
    
    annotate(
      "text",
      x = col_train, y = 6.65,
      label = "TRAINING BRANCH",
      size = 3.1, fontface = "bold", color = "#B45309"
    ) +
    
    annotate(
      "text",
      x = col_test, y = 6.65,
      label = "TEST BRANCH",
      size = 3.1, fontface = "bold", color = "#B45309"
    ) +
    
    scale_fill_identity() +
    
    labs(
      title = paste0(pathway_name, " - Train/Test Split & Analysis Flowchart")
    ) +
    
    coord_cartesian(
      xlim = c(-0.05, 2.05),
      ylim = c(0.3, 8.6),
      clip = "off"
    ) +
    
    theme_void(base_size = 13) +
    
    theme(
      plot.title = element_text(
        face = "plain",
        size = 13,
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      plot.margin = margin(15, 20, 15, 20)
    )
  
  ggsave(
    paste0(outfile_prefix, ".pdf"),
    plot = p,
    width = 9.5,
    height = 10,
    device = "pdf"
  )
  
  ggsave(
    paste0(outfile_prefix, ".png"),
    plot = p,
    width = 9.5,
    height = 10,
    dpi = 600
  )
  
  if (show_plot) {
    print(p)
  }
  
  return(p)
}

traintest_flow_plot <- plot_traintest_flow(
  global_split = global_split,
  pathway_name = pathway_name,
  outfile_prefix = output_file("Pathway_", pathway_name, "_TrainTest_Flowchart"),
  show_plot = TRUE
)

#==================================================================================
# DEFINE ANALYSIS DATASETS (for sensitivity analysis)
#==================================================================================
snp_gene_map_file <- output_file(
  "Pathway_",
  pathway_name,
  "_Reference_SNP_Gene_Map.csv"
)

analysis_configs <- list(
  
  main_ld = list(
    dataset = paste0(final_dataset_ld_prefix, ".rds"),
    label = "01_Main_LDpruned_All",
    type = "all_ld"
  ),
  
  sensitivity_betaweighted = list(
    dataset = paste0(final_dataset_ld_prefix, ".rds"),
    label = "05_Sensitivity_BetaWeighted_LDpruned",
    type = "beta_weighted"
  ),
  
  sensitivity_nold = list(
    dataset = paste0(final_dataset_unpruned_prefix, ".rds"),
    label = "02_Sensitivity_NoLDpruned_All",
    type = "all_nold"
  ),
  
  sensitivity_functional_vep = list(
    dataset = paste0(final_dataset_ld_prefix, ".rds"),
    vep_file = vep_results_file,
    label = "03_Sensitivity_Functional_VEP_LDpruned",
    type = "functional_vep"
  ),
  
  sensitivity_lof_vep = list(
    dataset = paste0(final_dataset_ld_prefix, ".rds"),
    vep_file = vep_results_file,
    label = "06_Sensitivity_LoF_VEP_LDpruned",
    type = "lof_vep"
  )
)

# Add the database sensitivity analysis (KEGG only vs. KEGG+Reactome) only
# if it was established above (block "PATHWAY SELECTION: DATABASE ORIGIN")
# that additional, non-redundant Reactome genes actually exist
# for this pathway (include_db_sensitivity == TRUE). For pure
# KEGG pathways without extension and for standalone Reactome pathways
# (zone 1) this entry is omitted - the analysis would otherwise be identical to the
# main analysis or undefined.
if (include_db_sensitivity) {

  analysis_configs$sensitivity_database_kegg <- list(
    dataset = paste0(final_dataset_ld_prefix, ".rds"),
    kegg_genes = kegg_genes,
    snp_gene_map = snp_gene_map_file,
    label = "07_Sensitivity_KEGGonly_LDpruned",
    type = "database_kegg_only"
  )
}



make_analysis_suffix <- function(config, score_mode = NULL) {
  
  suffix <- config$label
  
  if (!is.null(score_mode)) {
    suffix <- paste(
      suffix,
      tools::toTitleCase(score_mode),
      sep = "_"
    )
  }
  
  suffix
}

make_output_name <- function(
    pathway_name,
    analysis_type,
    analysis_name,
    result_name,
    config,
    score_mode = NULL,
    ext = "csv") {
  
  output_file(
    "Pathway_",
    pathway_name,
    "_",
    analysis_type,
    "_",
    analysis_name,
    "_",
    result_name,
    "_",
    make_analysis_suffix(config, score_mode),
    ".",
    ext
  )
}

#------------------------------------------------------------
# Broad functional VEP SNP selection
#------------------------------------------------------------

load_functional_vep_snps <- function(vep_file) {
  
  vep <- fread(vep_file)
  
  functional_terms <- c(
    "missense_variant",
    "protein_altering_variant",
    "frameshift_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_region_variant",
    "regulatory_region_variant",
    "TF_binding_site_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "mature_miRNA_variant"
  )
  
  functional <- vep[
    IMPACT %in% c("HIGH", "MODERATE") |
      grepl(
        paste(functional_terms, collapse = "|"),
        Consequence,
        ignore.case = TRUE
      )
  ]
  
  unique(
    functional$`#Uploaded_variation`
  )
}


#------------------------------------------------------------
# LoF VEP annotation including annotated allele
#------------------------------------------------------------

load_lof_vep_annotation <- function(vep_file) {
  
  vep <- fread(
    vep_file,
    data.table = FALSE
  )
  
  lof_terms <- c(
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "start_lost"
  )
  
  lof <- vep %>%
    transmute(
      SNP = trimws(
        as.character(`#Uploaded_variation`)
      ),
      VEP_LoF_Allele = toupper(
        trimws(as.character(Allele))
      ),
      VEP_REF_Allele = toupper(
        trimws(as.character(REF_ALLELE))
      ),
      VEP_Uploaded_Allele = toupper(
        trimws(as.character(UPLOADED_ALLELE))
      ),
      Consequence = as.character(Consequence)
    ) %>%
    filter(
      !is.na(SNP),
      SNP != "",
      !is.na(VEP_LoF_Allele),
      VEP_LoF_Allele != "",
      grepl(
        paste(lof_terms, collapse = "|"),
        Consequence,
        ignore.case = TRUE
      )
    )
  
  ambiguous_snps <- lof %>%
    group_by(SNP) %>%
    summarise(
      N_LoF_Alleles = n_distinct(VEP_LoF_Allele),
      .groups = "drop"
    ) %>%
    filter(
      N_LoF_Alleles > 1
    ) %>%
    pull(SNP)
  
  if (length(ambiguous_snps) > 0) {
    
    warning(
      "The following SNPs have several different ",
      "VEP LoF alleles and are excluded: ",
      paste(ambiguous_snps, collapse = ", ")
    )
    
    lof <- lof %>%
      filter(
        !SNP %in% ambiguous_snps
      )
  }
  
  lof %>%
    group_by(
      SNP,
      VEP_LoF_Allele,
      VEP_REF_Allele,
      VEP_Uploaded_Allele
    ) %>%
    summarise(
      LoF_Consequences = paste(
        sort(unique(Consequence)),
        collapse = ";"
      ),
      .groups = "drop"
    )
}


#------------------------------------------------------------
# Plain list of LoF SNPs
#------------------------------------------------------------

load_lof_vep_snps <- function(vep_file) {
  
  load_lof_vep_annotation(
    vep_file
  ) %>%
    pull(SNP) %>%
    unique()
}


#------------------------------------------------------------
# Allele information of the LD-pruned PLINK dataset
#------------------------------------------------------------

lof_bim <- read.table(
  paste0(ld_pruned_prefix, ".bim"),
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(lof_bim) <- c(
  "CHR",
  "SNP",
  "CM",
  "BP",
  "Allele1",
  "Allele2"
)

lof_bim <- lof_bim %>%
  transmute(
    SNP = trimws(as.character(SNP)),
    Matrix_Allele1 = toupper(
      trimws(as.character(Allele1))
    ),
    Matrix_Allele2 = toupper(
      trimws(as.character(Allele2))
    )
  ) %>%
  distinct(
    SNP,
    .keep_all = TRUE
  )

#------------------------------------------------------------
# Consistency check: read.plink() map versus final BIM file
#------------------------------------------------------------

plink_ld_map_check <- plink_ld$map %>%
  transmute(
    SNP = trimws(as.character(snp.name)),
    ReadPlink_Allele1 = toupper(
      trimws(as.character(allele.1))
    ),
    ReadPlink_Allele2 = toupper(
      trimws(as.character(allele.2))
    )
  ) %>%
  inner_join(
    lof_bim,
    by = "SNP"
  )

if (nrow(plink_ld_map_check) == 0) {
  stop(
    "No shared SNPs between plink_ld$map ",
    "and the final LD-pruned BIM file."
  )
}

allele_order_ok <- all(
  plink_ld_map_check$ReadPlink_Allele1 ==
    plink_ld_map_check$Matrix_Allele1 &
    plink_ld_map_check$ReadPlink_Allele2 ==
    plink_ld_map_check$Matrix_Allele2
)

if (!allele_order_ok) {
  stop(
    "The allele order in plink_ld$map does not match ",
    "the final LD-pruned BIM file."
  )
}

cat(
  "Allele order between read.plink() map and BIM confirmed:",
  nrow(plink_ld_map_check),
  "SNPs checked.\n"
)


#------------------------------------------------------------
# Map VEP LoF allele to matrix alleles
#------------------------------------------------------------

create_lof_orientation_table <- function(
    vep_file,
    bim_table,
    matrix_snps) {
  
  lof_annotation <- load_lof_vep_annotation(
    vep_file
  )
  
  lof_annotation %>%
    inner_join(
      bim_table,
      by = "SNP"
    ) %>%
    filter(
      SNP %in% matrix_snps
    ) %>%
    mutate(
      Allele_Match = case_when(
        
        VEP_LoF_Allele == Matrix_Allele2 ~
          "LoF_is_matrix_counted_allele",
        
        VEP_LoF_Allele == Matrix_Allele1 ~
          "LoF_is_matrix_other_allele",
        
        TRUE ~
          "No_allele_match"
      ),
      
      Invert_For_LoF_Score = case_when(
        Allele_Match == "LoF_is_matrix_counted_allele" ~ FALSE,
        Allele_Match == "LoF_is_matrix_other_allele" ~ TRUE,
        TRUE ~ NA
      )
    )
}


#------------------------------------------------------------
# Short check
#------------------------------------------------------------

functional_snps <- load_functional_vep_snps(
  vep_results_file
)

lof_snps <- load_lof_vep_snps(
  vep_results_file
)

cat(
  "Funktionelle VEP-SNPs:",
  length(functional_snps),
  "\n"
)

cat(
  "LoF-VEP-SNPs:",
  length(lof_snps),
  "\n"
)

# ============================================================
# Document LoF availability
# ============================================================

lof_annotation_check <- load_lof_vep_annotation(
  vep_results_file
)

lof_analysis_available <- nrow(lof_annotation_check) > 0

if (!lof_analysis_available) {
  
  lof_status <- data.frame(
    Pathway = pathway_name,
    Analysis = "sensitivity_lof_vep_lof_burden",
    Analysis_Config = "sensitivity_lof_vep",
    Score_Mode = "lof_burden",
    Status = "Not performed",
    Reason = paste0(
      "No classical loss-of-function variants were present ",
      "in the VEP-annotated LD-pruned SNP set."
    ),
    N_LoF_Variants = 0,
    stringsAsFactors = FALSE
  )
  
  write.csv2(
    lof_status,
    output_file(
      "Pathway_",
      pathway_name,
      "_LoF_Sensitivity_Status.csv"
    ),
    row.names = FALSE
  )
  
  message(
    "LoF sensitivity analysis not performed: ",
    "no classical LoF variants available."
  )
}

# ============================================================
# UNIVARIATE LOGISTIC REGRESSION – PATHWAY (main model only)
# ============================================================
library(dplyr)
library(pROC)

results_ld_file <- make_output_name(
  pathway_name,
  "Univariate",
  "SNP_Logistic_Regression",
  "Results",
  analysis_configs$main_ld
)

univar_ld_prefix <- output_file(
  "Pathway_",
  pathway_name,
  "_Univariate_LDpruned"
)

# 1. Load data
pathway_data <- readRDS(paste0(final_dataset_ld_prefix, ".rds"))

geno_filtered  <- as.matrix(pathway_data$genotypes)
pheno_filtered <- pathway_data$pheno

snps <- colnames(geno_filtered)
stopifnot("obesity" %in% colnames(pheno_filtered))

# TRAINING DATA ONLY FOR SNP DIRECTION ESTIMATION
global_split <- readRDS(
  output_file("Pathway_", pathway_name, "_Global_TrainTest_Split.rds")
)

train_ids <- global_split$train_ids

train_idx <- rownames(geno_filtered) %in% train_ids

geno_filtered <- geno_filtered[train_idx, , drop = FALSE]
pheno_filtered <- pheno_filtered[train_idx, ]

# Safety check
stopifnot(all(rownames(geno_filtered) == pheno_filtered$id))

bim_ambu <- read.table(
  paste0(plink_ambu_prefix, ".bim"),
  stringsAsFactors = FALSE
)

bim_md <- read.table(
  paste0(plink_md_prefix, ".bim"),
  stringsAsFactors = FALSE
)

colnames(bim_ambu) <- colnames(bim_md) <- c("CHR","SNP","CM","BP","A1","A2")

bim_all <- unique(rbind(bim_ambu, bim_md))

# 2. Result structure
results <- data.frame(
  SNP = snps,
  
  Estimate_Add = NA, StdError_Add = NA, Zvalue_Add = NA, Pvalue_Add = NA,
  OR_Add = NA, CI_lower_Add = NA, CI_upper_Add = NA, AIC_Add = NA,
  
  Estimate_Geno_1 = NA, OR_Geno_1 = NA,
  Estimate_Geno_2 = NA, OR_Geno_2 = NA,
  
  AIC_Geno = NA, LRT_pvalue = NA,
  
  Allele1 = NA, Allele2 = NA,
  Category = NA,
  Inheritance_Model = NA,
  
  N_Total = NA, N_Geno0 = NA, N_Geno1 = NA, N_Geno2 = NA,
  N_Imputed = NA,
  LowCountsFlag = FALSE,
  
  stringsAsFactors = FALSE
)

# 3. Mode-of-inheritance function
determine_inheritance <- function(OR1, OR2, AIC_add, AIC_geno, LRT_p) {
  
  if (!is.finite(OR1) | !is.finite(OR2)) return("Unclear")
  
  deltaAIC <- AIC_geno - AIC_add
  
  if (deltaAIC < -2 & !is.na(LRT_p) & LRT_p <= 0.05)
    return("Genotypisch")
  
  if (abs(log(OR2) - 2 * log(OR1)) < 0.3)
    return("Additive")
  
  if (abs(log(OR1) - log(OR2)) < 0.25 & (OR1 > 1.1 | OR1 < 0.9))
    return("Dominant")
  
  if (abs(OR1 - 1) < 0.15 & (OR2 > 1.3 | OR2 < 0.7))
    return("Rezessiv")
  
  return("Unclear")
}

# 4. Loop over SNPs
low_count_threshold <- 10

for (i in seq_along(snps)) {
  
  snp_name <- snps[i]
  snp_values <- geno_filtered[, snp_name]
  
  # Only the phenotype must be available; missing genotypes are imputed
  # (mean imputation per SNP, based exclusively on the training data).
  idx_pheno <- !is.na(pheno_filtered$obesity)
  if (sum(idx_pheno) < 10) next
  
  snp_values      <- snp_values[idx_pheno]
  obesity_values  <- pheno_filtered$obesity[idx_pheno]
  
  # -------------------------
  # Mean imputation of missing genotypes
  # -------------------------
  snp_mean <- mean(snp_values, na.rm = TRUE)
  
  if (is.nan(snp_mean)) next  # SNP 100% missing -> not analysable
  
  n_imputed <- sum(is.na(snp_values))
  
  snp_values_observed <- snp_values  # unimputed, for genotype counts
  snp_values_imputed  <- snp_values
  snp_values_imputed[is.na(snp_values_imputed)] <- snp_mean
  
  df <- data.frame(
    obesity = obesity_values,
    snp = snp_values_imputed
  )
  
  results$N_Imputed[i] <- n_imputed
  
  # -------------------------
  # Allele
  # -------------------------
  bim_row <- bim_all[bim_all$SNP == snp_name, ]
  
  if (nrow(bim_row) >= 1) {
    results$Allele1[i] <- bim_row$A1[1]
    results$Allele2[i] <- bim_row$A2[1]
  }
  
  
  
  # -------------------------
  # Genotype counts (based on observed, non-imputed genotypes)
  # -------------------------
  results$N_Total[i] <- nrow(df)
  
  results$N_Geno0[i] <- sum(snp_values_observed == 0, na.rm = TRUE)
  results$N_Geno1[i] <- sum(snp_values_observed == 1, na.rm = TRUE)
  results$N_Geno2[i] <- sum(snp_values_observed == 2, na.rm = TRUE)
  
  results$LowCountsFlag[i] <-
    any(c(results$N_Geno0[i],
          results$N_Geno1[i],
          results$N_Geno2[i]) < low_count_threshold)
  
  # -------------------------
  # Additive model
  # -------------------------
  model_add <- tryCatch(
    glm(obesity ~ snp, data = df, family = binomial),
    error = function(e) NULL
  )
  
  if (!is.null(model_add)) {
    
    co <- summary(model_add)$coefficients
    
    results$Estimate_Add[i] <- co["snp","Estimate"]
    results$StdError_Add[i] <- co["snp","Std. Error"]
    results$Zvalue_Add[i]   <- co["snp","z value"]
    results$Pvalue_Add[i]   <- co["snp","Pr(>|z|)"]
    
    results$OR_Add[i] <- exp(co["snp","Estimate"])
    
    ci <- tryCatch(
      exp(confint.default(model_add)["snp", ]),
      error = function(e) c(NA, NA)
    )
    
    results$CI_lower_Add[i] <- ci[1]
    results$CI_upper_Add[i] <- ci[2]
    
    results$AIC_Add[i] <- AIC(model_add)
  }
  
  # -------------------------
  # Genotypic model
  # -------------------------
  # Fractional, imputed dosage values cannot be represented as 0/1/2 factor levels
  # -> this model uses exclusively the observed genotypes.
  df_geno <- data.frame(
    obesity    = obesity_values[!is.na(snp_values_observed)],
    snp_factor = factor(
      snp_values_observed[!is.na(snp_values_observed)],
      levels = c(0,1,2)
    )
  )
  
  model_geno <- tryCatch(
    glm(obesity ~ snp_factor, data = df_geno, family = binomial),
    error = function(e) NULL
  )
  
  if (!is.null(model_geno)) {
    
    co <- summary(model_geno)$coefficients
    
    if ("snp_factor1" %in% rownames(co)) {
      results$OR_Geno_1[i] <- exp(co["snp_factor1","Estimate"])
    }
    
    if ("snp_factor2" %in% rownames(co)) {
      results$OR_Geno_2[i] <- exp(co["snp_factor2","Estimate"])
    }
    
    results$AIC_Geno[i] <- AIC(model_geno)
    
    # LRT: the additive model must be compared on the same (observed) data as the
    # genotypic model. model_add may run on more
    # (imputed) observations and thus have a different N -> here
    # an additive model on df_geno (observed values
    # only) is carried along for the comparison.
    model_add_obs <- tryCatch(
      glm(
        obesity ~ snp,
        data = data.frame(
          obesity = df_geno$obesity,
          snp = as.numeric(as.character(df_geno$snp_factor))
        ),
        family = binomial
      ),
      error = function(e) NULL
    )
    
    if (!is.null(model_add_obs)) {
      LRT <- anova(model_add_obs, model_geno, test = "LRT")
      results$LRT_pvalue[i] <- LRT$`Pr(>Chi)`[2]
    }
  }
  
  # -------------------------
  # Inheritance mode
  # -------------------------
  results$Inheritance_Model[i] <- determine_inheritance(
    results$OR_Geno_1[i],
    results$OR_Geno_2[i],
    results$AIC_Add[i],
    results$AIC_Geno[i],
    results$LRT_pvalue[i]
  )
  
  # -------------------------
  # Risk / protective
  # -------------------------
  results$Category[i] <-
    ifelse(is.na(results$OR_Add[i]), NA,
           ifelse(results$OR_Add[i] > 1, "Risk", "Protective"))
  
  cat(
    "SNP", i, "/", length(snps), "-", snp_name,
    "| N:", results$N_Total[i],
    "| Genotypes:", results$N_Geno0[i], "/", results$N_Geno1[i], "/", results$N_Geno2[i],
    "| Model:", results$Inheritance_Model[i],
    "| Category:", results$Category[i],
    "| LowCountsFlag:", results$LowCountsFlag[i],
    "\n"
  )
}

# 5. Export
write.csv2(
  results,
  results_ld_file,
  row.names = FALSE
)

cat("\n>>> DONE (no gene mapping)\n")

#===================================================================================
# QQ-PLOT UND HISTOGRAMM DER LOG(OR) – UNIVARIATE REGRESSION (HAUPTMODELL)
#===================================================================================
# Descriptive classification of the univariate effect estimates (log(OR) = Estimate_Add)
# and their p-values (main_ld only, as this section covers only the main model).
qq_hist_data <- results %>%
  filter(!is.na(Pvalue_Add), !is.na(Estimate_Add))
if (nrow(qq_hist_data) > 0) {
  
  # -------------------------
  # QQ plot of p-values
  # -------------------------
  n_snps_qq <- nrow(qq_hist_data)
  
  qq_df <- data.frame(
    observed = -log10(sort(qq_hist_data$Pvalue_Add)),
    expected = -log10(ppoints(n_snps_qq))
  )
  
  qq_plot <- ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(color = "#065A82", size = 1.8, alpha = 0.75) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    labs(
      title = paste0("QQ plot of p-values - univariate regression (", pathway_name, ")"),
      x = "-log10(erwarteter p-Wert)",
      y = "-log10(beobachteter p-Wert)"
    ) +
    theme_minimal(base_size = 12)
  
  ggsave(
    filename = output_file("Pathway_", pathway_name, "_Univariate_QQPlot_Pvalues.png"),
    plot = qq_plot,
    width = 6, height = 6, dpi = 300
  )
  
  print(qq_plot)
  
  # -------------------------
  # Histogram of log(OR)
  # -------------------------
  logor_hist <- ggplot(qq_hist_data, aes(x = Estimate_Add)) +
    geom_histogram(bins = 30, fill = "#065A82", color = "white", alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    labs(
      title = paste0("Distribution of log(OR) - univariate regression (", pathway_name, ")"),
      x = "log(OR) (Estimate_Add)",
      y = "Anzahl SNPs"
    ) +
    theme_minimal(base_size = 12)
  
  ggsave(
    filename = output_file("Pathway_", pathway_name, "_Univariate_Histogram_LogOR.png"),
    plot = logor_hist,
    width = 6, height = 5, dpi = 300
  )
  
  print(logor_hist)
  
  cat("\nQQ plot and histogram of log(OR) saved:\n")
  cat(" -", output_file("Pathway_", pathway_name, "_Univariate_QQPlot_Pvalues.png"), "\n")
  cat(" -", output_file("Pathway_", pathway_name, "_Univariate_Histogram_LogOR.png"), "\n")
  
} else {
  message("No QQ plot/histogram created: no valid Pvalue_Add/Estimate_Add values available.")
}

#===================================================================================
# SNP-GENE-MAPPING
#===================================================================================

range_files_dir <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Working_Plink")

gene_range_file <- file.path(
  range_files_dir,
  paste0("pathway_", pathway_file_id, "_gene_ranges.txt")
)

gene_ranges <- read.table(
  gene_range_file,
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(gene_ranges) <- c(
  "CHR",
  "START",
  "END",
  "GENE"
)

colnames(bim_all) <- c(
  "CHR",
  "SNP",
  "CM",
  "BP",
  "A1",
  "A2"
)

gene_ranges$CHR <- as.character(gene_ranges$CHR)
bim_all$CHR <- as.character(bim_all$CHR)

snp_gene_map <- data.frame()

for(i in seq_len(nrow(gene_ranges))){
  
  chr_i   <- gene_ranges$CHR[i]
  start_i <- gene_ranges$START[i]
  end_i   <- gene_ranges$END[i]
  gene_i  <- gene_ranges$GENE[i]
  
  snps_in_gene <- bim_all %>%
    filter(
      CHR == chr_i,
      BP >= start_i,
      BP <= end_i
    ) %>%
    mutate(GENE = gene_i)
  
  snp_gene_map <- rbind(
    snp_gene_map,
    snps_in_gene
  )
}

write.csv2(
  snp_gene_map,
  snp_gene_map_file,
  row.names = FALSE
)

# ------------------------------------------------------------------------
# Validation: does the actual dataset (extracted externally via PLINK)
# really cover the intended pathway gene set (full_genes, see above - KEGG only
# or KEGG+Reactome, depending on the two-zone model)? Pure
# check, no intervention in the data - full_genes itself does not control the
# extraction (which runs externally via the range file); serves here
# only as a cross-check against silent discrepancies between pathway selection
# and the actually extracted dataset.
# ------------------------------------------------------------------------
mapped_genes <- unique(snp_gene_map$GENE)
missing_from_map <- setdiff(full_genes, mapped_genes)
unexpected_in_map <- setdiff(mapped_genes, full_genes)

if (length(missing_from_map) > 0) {
  warning(
    "The following genes from full_genes (intended pathway gene set) have ",
    "no SNPs in the SNP-gene map / in the PLINK dataset: ",
    paste(missing_from_map, collapse = ", ")
  )
}

if (length(unexpected_in_map) > 0) {
  warning(
    "The following genes are contained in the SNP-gene map although they do not occur in ",
    "full_genes (intended pathway gene set) - please cross-check the range file: ",
    "",
    paste(unexpected_in_map, collapse = ", ")
  )
}

cat(
  length(mapped_genes), "of", length(full_genes),
  "intended genes represented with at least one SNP in the final dataset.\n"
)

#===================================================================================
# UNIVARIATE SNP SUM SCORE – LOOP OVER ANALYSIS DATASETS AND SCORE MODES
#===================================================================================

analysis_plan <- data.frame(
  config_name = c(
    "main_ld",
    "main_ld",
    "sensitivity_nold",
    "sensitivity_nold",
    "sensitivity_functional_vep",
    "sensitivity_functional_vep",
    "sensitivity_lof_vep",
    "sensitivity_betaweighted",
    "sensitivity_betaweighted",
    "sensitivity_betaweighted"
  ),
  score_mode = c(
    "oriented",
    "raw",
    "oriented",
    "raw",
    "oriented",
    "raw",
    "lof_burden",
    "beta_weighted",
    "ridge_weighted",
    "elasticnet_weighted"
  ),
  stringsAsFactors = FALSE
)

# Include the database sensitivity analysis (KEGG only vs. KEGG+Reactome) in the
# execution plan only if it was established above (block "PATHWAY SELECTION:
# DATABASE ORIGIN") that for this pathway
# a sensitivity_database_kegg entry actually exists
# (include_db_sensitivity == TRUE) - analogous to the conditional addition to
# analysis_configs further above. Score modes (oriented/raw) correspond to the
# remaining structural subset sensitivity analyses (nold, functional_vep).
if (include_db_sensitivity) {

  analysis_plan <- rbind(
    analysis_plan,
    data.frame(
      config_name = c("sensitivity_database_kegg", "sensitivity_database_kegg"),
      score_mode = c("oriented", "raw"),
      stringsAsFactors = FALSE
    )
  )
}

if (!lof_analysis_available) {
  
  analysis_plan <- analysis_plan %>%
    filter(
      !(config_name == "sensitivity_lof_vep" &
          score_mode == "lof_burden")
    )
  
  message(
    "LoF analysis removed from analysis_plan ",
    "(no classical LoF variants available)."
  )
}

print(analysis_plan)

valid_score_modes <- c(
  "oriented",
  "raw",
  "beta_weighted",
  "ridge_weighted",
  "elasticnet_weighted",
  "lof_burden"
)

unknown_configs <- setdiff(
  unique(analysis_plan$config_name),
  names(analysis_configs)
)

unknown_score_modes <- setdiff(
  unique(analysis_plan$score_mode),
  valid_score_modes
)

if (length(unknown_configs) > 0) {
  stop(
    "Unknown configurations in analysis_plan: ",
    paste(unknown_configs, collapse = ", ")
  )
}

if (length(unknown_score_modes) > 0) {
  stop(
    "Unknown score modes in analysis_plan: ",
    paste(unknown_score_modes, collapse = ", ")
  )
}

if (anyDuplicated(analysis_plan) > 0) {
  stop(
    "The analysis_plan contains duplicate combinations of ",
    "config_name and score_mode."
  )
}

run_snp_sumscore <- function(config_name, config, score_mode) {
  
  cat("\n====================================================\n")
  cat("Analysis:", config_name, "| Score mode:", score_mode, "\n")
  cat("====================================================\n")
  
  if (!score_mode %in% c(
    "oriented",
    "raw",
    "beta_weighted",
    "ridge_weighted",
    "elasticnet_weighted",
    "lof_burden"
  )) {
    stop(
      "Unbekannter Score-Modus: ",
      score_mode
    )
  }
  
  if (
    score_mode %in% c("ridge_weighted", "elasticnet_weighted") &&
    config$type != "beta_weighted"
  ) {
    stop(
      "score_mode = 'ridge_weighted'/'elasticnet_weighted' may only be used with ",
      "config$type = 'beta_weighted' verwendet werden (gleicher Datensatz ",
      "as the beta-weighted sensitivity analysis)."
    )
  }
  
  if (
    score_mode == "lof_burden" &&
    config$type != "lof_vep"
  ) {
    stop(
      "score_mode = 'lof_burden' may only be used with ",
      "config$type = 'lof_vep' verwendet werden."
    )
  }
  
  #------------------------------------------------------------
  # Output files
  #------------------------------------------------------------
  
  individuals_file <- make_output_name(
    pathway_name,
    "Univariate",
    "SNP_SumScore",
    "Individuals",
    config,
    score_mode
  )
  
  percentiles_file <- make_output_name(
    pathway_name,
    "Univariate",
    "SNP_SumScore",
    "Percentiles",
    config,
    score_mode
  )
  
  weights_file <- make_output_name(
    pathway_name,
    "Univariate",
    "SNP_SumScore",
    "Weights",
    config,
    score_mode
  )
  
  #------------------------------------------------------------
  # Load data
  #------------------------------------------------------------
  pathway_data <- readRDS(
    config$dataset
  )
  
  geno_filtered <- pathway_data$genotypes
  
  pheno_filtered <- as.data.frame(
    pathway_data$pheno,
    stringsAsFactors = FALSE
  )
  
  if (is.null(rownames(geno_filtered))) {
    stop(
      "The genotype matrix has no row names/IDs."
    )
  }
  
  if (!"id" %in% colnames(pheno_filtered)) {
    stop(
      "The phenotype data contain no column 'id'."
    )
  }
  
  rownames(geno_filtered) <- as.character(
    rownames(geno_filtered)
  )
  
  pheno_filtered$id <- as.character(
    pheno_filtered$id
  )
  
  #------------------------------------------------------------
  # Effect sizes from the univariate LD-pruned main analysis
  #------------------------------------------------------------
  
  regression_file <- make_output_name(
    pathway_name = pathway_name,
    analysis_type = "Univariate",
    analysis_name = "SNP_Logistic_Regression",
    result_name = "Results",
    config = analysis_configs$main_ld
  )
  
  if (!file.exists(regression_file)) {
    stop(
      "Univariate regression results not found:\n",
      regression_file
    )
  }
  
  results <- read.csv2(
    regression_file,
    stringsAsFactors = FALSE
  )
  
  required_result_cols <- c(
    "SNP",
    "Category",
    "Estimate_Add",
    "OR_Add",
    "Pvalue_Add"
  )
  
  missing_result_cols <- setdiff(
    required_result_cols,
    colnames(results)
  )
  
  if (length(missing_result_cols) > 0) {
    stop(
      "The following columns are missing in the regression results: ",
      paste(missing_result_cols, collapse = ", ")
    )
  }
  
  results <- results %>%
    mutate(
      SNP = trimws(as.character(SNP)),
      Estimate_Add = as.numeric(Estimate_Add),
      OR_Add = as.numeric(OR_Add),
      Pvalue_Add = as.numeric(Pvalue_Add)
    )
  
  # Table for equally weighted scores
  res_cat <- results %>%
    filter(
      !is.na(SNP),
      SNP != "",
      !is.na(Category),
      Category != ""
    ) %>%
    dplyr::select(
      SNP,
      Category
    ) %>%
    distinct(
      SNP,
      .keep_all = TRUE
    )
  
  # Table for β-weighted score
  res_beta <- results %>%
    filter(
      !is.na(SNP),
      SNP != "",
      is.finite(Estimate_Add)
    ) %>%
    dplyr::select(
      SNP,
      Beta = Estimate_Add,
      OR = OR_Add,
      Pvalue = Pvalue_Add,
      Category
    ) %>%
    distinct(
      SNP,
      .keep_all = TRUE
    )
  
  #------------------------------------------------------------
  # Match IDs robustly and sort identically
  #------------------------------------------------------------
  
  geno_ids <- as.character(
    rownames(geno_filtered)
  )
  
  pheno_ids <- as.character(
    pheno_filtered$id
  )
  
  if (anyDuplicated(geno_ids) > 0) {
    stop(
      "Duplicate IDs in the genotype matrix for ",
      config_name,
      "."
    )
  }
  
  if (anyDuplicated(pheno_ids) > 0) {
    stop(
      "Duplicate IDs in the phenotype data for ",
      config_name,
      "."
    )
  }
  
  # Keep the order of the genotype matrix
  common_ids <- geno_ids[
    geno_ids %in% pheno_ids
  ]
  
  if (length(common_ids) == 0) {
    
    warning(
      "No shared genotype/phenotype IDs for ",
      config_name
    )
    
    return(NULL)
  }
  
  geno_df <- geno_filtered[
    common_ids,
    ,
    drop = FALSE
  ] |>
    as.data.frame(
      stringsAsFactors = FALSE
    )
  
  pheno_match <- match(
    common_ids,
    pheno_filtered$id
  )
  
  if (anyNA(pheno_match)) {
    stop(
      "At least one genotype ID could not be assigned to the ",
      "phenotype data."
    )
  }
  
  pheno_df <- pheno_filtered[
    pheno_match,
    ,
    drop = FALSE
  ]
  
  pheno_df$id <- as.character(
    pheno_df$id
  )
  
  rownames(geno_df) <- common_ids
  rownames(pheno_df) <- common_ids
  
  if (!identical(
    as.character(rownames(geno_df)),
    as.character(pheno_df$id)
  )) {
    
    stop(
      "ID order after genotype/phenotype matching ",
      "does not match."
    )
  }
  
  #------------------------------------------------------------
  # Base set of SNPs depending on score mode
  #------------------------------------------------------------
  
  if (score_mode %in% c("beta_weighted", "ridge_weighted", "elasticnet_weighted")) {
    
    # Only SNPs with a valid, training-based estimated beta
    # (ridge_weighted/elasticnet_weighted use the same SNP base as
    # beta_weighted, so that all three weighting variants are directly comparable)
    valid_snps <- intersect(
      res_beta$SNP,
      colnames(geno_df)
    )
    
    #----------------------------------------------------------
    # Ridge / elastic net: shrunken weights instead of raw log(OR)
    #----------------------------------------------------------
    # Replaces res_beta$Beta with jointly (not SNP-by-SNP univariately)
    # estimated, regularized coefficients - exclusively on the
    # training data, analogous to the Estimate_Add estimation of the beta_weighted mode.
    if (score_mode %in% c("ridge_weighted", "elasticnet_weighted")) {
      
      library(glmnet)
      
      global_split_glmnet <- readRDS(
        output_file("Pathway_", pathway_name, "_Global_TrainTest_Split.rds")
      )
      
      train_ids_glmnet <- intersect(common_ids, global_split_glmnet$train_ids)
      
      if (length(train_ids_glmnet) < 10) {
        stop(
          "Too few training individuals for ridge/elastic-net weighting: ",
          config_name, " / ", score_mode
        )
      }
      
      x_train_raw <- as.matrix(
        geno_df[train_ids_glmnet, valid_snps, drop = FALSE]
      )
      storage.mode(x_train_raw) <- "numeric"
      
      # Mean imputation of missing genotypes per SNP (training data only),
      # consistent with the imputation in the univariate regression.
      x_train <- apply(x_train_raw, 2, function(col) {
        col[is.na(col)] <- mean(col, na.rm = TRUE)
        col
      })
      
      y_train <- as.numeric(
        as.character(pheno_df[train_ids_glmnet, "obesity"])
      )
      
      glmnet_alpha <- if (score_mode == "ridge_weighted") 0 else 0.5
      
      set.seed(123)
      cv_fit <- cv.glmnet(
        x = x_train,
        y = y_train,
        family = "binomial",
        alpha = glmnet_alpha
      )
      
      glmnet_coefs <- as.matrix(coef(cv_fit, s = "lambda.1se"))
      glmnet_beta  <- glmnet_coefs[valid_snps, 1]
      names(glmnet_beta) <- valid_snps
      
      # Replace the univariate beta column (log(OR)) with shrunken coefficients;
      # category/OR/p-value from the univariate analysis are kept
      # in the weights file for classification/documentation.
      res_beta <- res_beta %>%
        filter(SNP %in% valid_snps) %>%
        mutate(Beta = glmnet_beta[SNP])
      
      cat(
        "Ridge/Elastic-Net-Gewichte berechnet (", score_mode,
        ", alpha = ", glmnet_alpha,
        ", lambda.1se = ", signif(cv_fit$lambda.1se, 4),
        ", n_train = ", length(train_ids_glmnet), ")\n",
        sep = ""
      )
    }
    
  } else if (score_mode == "lof_burden") {
    
    # LoF SNPs are defined exclusively biologically via VEP.
    # No dependence on obesity-based risk/protective categories.
    valid_snps <- colnames(geno_df)
    
  } else if (score_mode %in% c("oriented", "raw")) {
    
    # Use the same SNP base for oriented and raw,
    # so that both analyses differ only by orientation.
    valid_snps <- intersect(
      res_cat$SNP,
      colnames(geno_df)
    )
    
  } else {
    
    stop(
      "Unsupported score mode in SNP selection: ",
      score_mode
    )
  }
  
  #------------------------------------------------------------
  # SNP selection depending on analysis configuration
  #------------------------------------------------------------
  lof_orientation_table <- NULL
  
  if (config$type == "functional_vep") {
    
    functional_snps <- load_functional_vep_snps(
      config$vep_file
    )
    
    valid_snps <- intersect(
      valid_snps,
      functional_snps
    )
  }
  
  if (config$type == "lof_vep") {
    
    if (score_mode != "lof_burden") {
      stop(
        "The configuration 'lof_vep' may only be run with ",
        "score_mode = 'lof_burden'."
      )
    }
    
    lof_orientation_table <- create_lof_orientation_table(
      vep_file = config$vep_file,
      bim_table = lof_bim,
      matrix_snps = colnames(geno_df)
    )
    
    if (nrow(lof_orientation_table) == 0) {
      warning(
        "No LoF SNPs could be matched between VEP, BIM and ",
        "genotype matrix."
      )
      return(NULL)
    }
    
    no_match_snps <- lof_orientation_table %>%
      filter(
        Allele_Match == "No_allele_match"
      ) %>%
      pull(SNP) %>%
      unique()
    
    if (length(no_match_snps) > 0) {
      warning(
        "The following VEP LoF alleles match no matrix allele ",
        "and are excluded: ",
        paste(no_match_snps, collapse = ", ")
      )
    }
    
    lof_orientation_table <- lof_orientation_table %>%
      filter(
        Allele_Match != "No_allele_match",
        !is.na(Invert_For_LoF_Score)
      ) %>%
      distinct(
        SNP,
        .keep_all = TRUE
      )
    
    valid_snps <- intersect(
      valid_snps,
      lof_orientation_table$SNP
    )
  }
  
  if (config$type == "database_kegg_only") {
    
    snp_gene_map <- read.csv2(
      config$snp_gene_map,
      stringsAsFactors = FALSE
    )
    
    kegg_only_snps <- snp_gene_map %>%
      filter(
        GENE %in% config$kegg_genes
      ) %>%
      pull(SNP) %>%
      unique()
    
    valid_snps <- intersect(
      valid_snps,
      kegg_only_snps
    )
  }
  
  # Align order with the genotype matrix and remove duplicates
  valid_snps <- colnames(geno_df)[
    colnames(geno_df) %in% unique(valid_snps)
  ]
  
  if (length(valid_snps) == 0) {
    
    warning(
      "No valid SNPs for ",
      config_name,
      " / ",
      score_mode,
      " found."
    )
    
    return(NULL)
  }
  
 
  
  #------------------------------------------------------------
  # Prepare genotype matrix
  #------------------------------------------------------------
  
  geno_score <- geno_df[
    ,
    valid_snps,
    drop = FALSE
  ]
  
  geno_score[] <- lapply(
    geno_score,
    function(x) {
      
      if (is.list(x)) {
        x <- unlist(x)
      }
      
      as.numeric(x)
    }
  )
  
  geno_matrix_original <- as.matrix(geno_score)
  storage.mode(geno_matrix_original) <- "numeric"
  
  qc_range <- range(
    geno_matrix_original,
    na.rm = TRUE
  )
  
  if (
    !all(is.finite(qc_range)) ||
    qc_range[1] < 0 ||
    qc_range[2] > 2
  ) {
    stop(
      "Genotype values outside the expected range 0–2: ",
      paste(qc_range, collapse = " bis ")
    )
  }
  
  #------------------------------------------------------------
  # Exclude individuals with more than 10 % missing SNPs
  #------------------------------------------------------------
  
  indiv_missing <- rowMeans(
    is.na(geno_matrix_original)
  )
  
  keep <- indiv_missing <= 0.10
  
  ids_before_missingness <- as.character(
    rownames(geno_score)
  )
  
  ids_after_missingness <- ids_before_missingness[
    keep
  ]
  
  geno_score <- geno_score[
    ids_after_missingness,
    ,
    drop = FALSE
  ]
  
  pheno_match_missingness <- match(
    ids_after_missingness,
    pheno_df$id
  )
  
  if (anyNA(pheno_match_missingness)) {
    stop(
      "At least one ID could not be assigned to the phenotype data after the missingness filter ",
      ""
    )
  }
  
  pheno_df <- pheno_df[
    pheno_match_missingness,
    ,
    drop = FALSE
  ]
  
  pheno_df$id <- as.character(
    pheno_df$id
  )
  
  rownames(pheno_df) <- pheno_df$id
  
  if (!identical(
    as.character(rownames(geno_score)),
    as.character(pheno_df$id)
  )) {
    
    stop(
      "ID order after the individual missingness filter ",
      "does not match."
    )
  }
  
  #------------------------------------------------------------
  # Training-based genotype imputation for ALL score modes
  #------------------------------------------------------------
  
  global_split <- readRDS(
    output_file(
      "Pathway_",
      pathway_name,
      "_Global_TrainTest_Split.rds"
    )
  )
  
  train_rows <- rownames(geno_score) %in%
    global_split$train_ids
  
  if (!any(train_rows)) {
    stop(
      "No training individuals available for genotype imputation."
    )
  }
  
  # Numeric original matrix after exclusion of individuals
  # with > 10 % missing genotypes
  geno_matrix <- as.matrix(geno_score)
  storage.mode(geno_matrix) <- "numeric"
  
  # Determine SNP means exclusively in the training dataset
  train_snp_means <- colMeans(
    geno_matrix[
      train_rows,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )
  
  # SNPs without a single observed genotype in the training set
  invalid_imputation_snps <- names(
    train_snp_means
  )[
    !is.finite(train_snp_means)
  ]
  
  if (length(invalid_imputation_snps) > 0) {
    
    warning(
      "The following SNPs have no observed genotype values in the training set ",
      "and are excluded: ",
      paste(
        invalid_imputation_snps,
        collapse = ", "
      )
    )
    
    valid_imputation_snps <- setdiff(
      colnames(geno_matrix),
      invalid_imputation_snps
    )
    
    geno_matrix <- geno_matrix[
      ,
      valid_imputation_snps,
      drop = FALSE
    ]
    
    train_snp_means <- train_snp_means[
      valid_imputation_snps
    ]
    
    valid_snps <- intersect(
      valid_snps,
      valid_imputation_snps
    )
    
    if (length(valid_snps) == 0) {
      warning(
        "After excluding non-imputable SNPs, no SNPs remain."
      )
      return(NULL)
    }
  }
  
  # Replace missing values in train AND test exclusively with the
  # SNP means computed in the training set
  n_imputed_total <- 0L
  
  n_imputed_by_snp <- setNames(
    integer(ncol(geno_matrix)),
    colnames(geno_matrix)
  )
  
  for (j in seq_len(ncol(geno_matrix))) {
    
    missing_j <- is.na(
      geno_matrix[, j]
    )
    
    n_imputed_by_snp[j] <- sum(
      missing_j
    )
    
    if (any(missing_j)) {
      
      geno_matrix[
        missing_j,
        j
      ] <- train_snp_means[j]
      
      n_imputed_total <- n_imputed_total +
        sum(missing_j)
    }
  }
  
  if (anyNA(geno_matrix)) {
    stop(
      "After the training-based imputation, missing ",
      "genotype values remain."
    )
  }
  
  # Documentation of the imputation values used
  imputation_table <- data.frame(
    SNP = colnames(geno_matrix),
    Training_Mean_Genotype = as.numeric(
      train_snp_means[
        colnames(geno_matrix)
      ]
    ),
    N_Imputed_All = as.integer(
      n_imputed_by_snp[
        colnames(geno_matrix)
      ]
    ),
    Imputation_Source = "Training_Set",
    stringsAsFactors = FALSE
  )
  
  imputation_file <- make_output_name(
    pathway_name,
    "Univariate",
    "SNP_SumScore",
    "Imputation_Values",
    config,
    score_mode
  )
  
  write.csv2(
    imputation_table,
    imputation_file,
    row.names = FALSE
  )
  
  #------------------------------------------------------------
  # Compute score
  #------------------------------------------------------------
  
  protective <- character(0)
  beta_table_used <- NULL
  prs_raw <- NULL
  snp_sum <- NULL
  
  if (score_mode == "oriented") {
    
    protective <- res_cat %>%
      filter(
        tolower(Category) == "protective"
      ) %>%
      pull(SNP) %>%
      intersect(
        colnames(geno_matrix)
      )
    
    # Re-orient protective SNPs only after imputation
    if (length(protective) > 0) {
      
      protective_idx <- match(
        protective,
        colnames(geno_matrix)
      )
      
      geno_matrix[
        ,
        protective_idx
      ] <- 2 - geno_matrix[
        ,
        protective_idx,
        drop = FALSE
      ]
    }
    
    snp_sum <- rowSums(
      geno_matrix
    )
  }
  
  if (score_mode == "raw") {
    
    snp_sum <- rowSums(
      geno_matrix
    )
  }
  
  if (score_mode == "lof_burden") {
    
    if (is.null(lof_orientation_table)) {
      stop(
        "For the LoF burden no ",
        "allele orientation table was created."
      )
    }
    
    # Create copy:
    # The original geno_matrix remains unchanged.
    geno_matrix_lof <- geno_matrix
    
    # Align the orientation table exactly with the SNP order
    # of the matrix actually used
    lof_orientation_used <- lof_orientation_table[
      match(
        colnames(geno_matrix_lof),
        lof_orientation_table$SNP
      ),
      ,
      drop = FALSE
    ]
    
    if (
      any(is.na(lof_orientation_used$SNP)) ||
      !all(
        lof_orientation_used$SNP ==
        colnames(geno_matrix_lof)
      )
    ) {
      stop(
        "Not for all LoF SNPs used could an ",
        "unambiguous allele orientation be assigned."
      )
    }
    
    # Identify SNPs for which the original matrix
    # counts the counter-allele of the VEP LoF allele
    invert_snps <- lof_orientation_used$SNP[
      lof_orientation_used$Invert_For_LoF_Score %in% TRUE
    ]
    
    # Re-orient only the temporary LoF matrix
    if (length(invert_snps) > 0) {
      
      invert_idx <- match(
        invert_snps,
        colnames(geno_matrix_lof)
      )
      
      geno_matrix_lof[
        ,
        invert_idx
      ] <- 2 - geno_matrix_lof[
        ,
        invert_idx,
        drop = FALSE
      ]
    }
    
    #----------------------------------------------------------
    # Document LoF allele orientation used
    #----------------------------------------------------------
    
    lof_orientation_file <- make_output_name(
      pathway_name,
      "Univariate",
      "SNP_SumScore",
      "LoF_Allele_Orientation",
      config,
      score_mode
    )
    
    lof_orientation_used$Applied_Transformation <- ifelse(
      lof_orientation_used$Invert_For_LoF_Score %in% TRUE,
      "2_minus_original_dosage",
      "original_dosage"
    )
    
    write.csv2(
      lof_orientation_used,
      lof_orientation_file,
      row.names = FALSE
    )
    
    #----------------------------------------------------------
    # Compute LoF burden
    #----------------------------------------------------------
    
    snp_sum <- rowSums(
      geno_matrix_lof
    )
  }
  
  if (score_mode %in% c("beta_weighted", "ridge_weighted", "elasticnet_weighted")) {
    
    #----------------------------------------------------------
    # Align β weights with the order of the imputed matrix
    #----------------------------------------------------------
    
    beta_table_used <- res_beta[
      match(
        colnames(geno_matrix),
        res_beta$SNP
      ),
      ,
      drop = FALSE
    ]
    
    if (
      any(is.na(beta_table_used$SNP)) ||
      any(!is.finite(beta_table_used$Beta))
    ) {
      stop(
        "Not for all selected SNPs could a valid β value ",
        "zugeordnet werden."
      )
    }
    
    stopifnot(
      all(
        beta_table_used$SNP ==
          colnames(geno_matrix)
      )
    )
    
    beta_vector <- beta_table_used$Beta
    names(beta_vector) <- beta_table_used$SNP
    
    # Additionally document unweighted dosage
    prs_raw <- rowSums(
      geno_matrix
    )
    
    # β-weighted score:
    # Sum of genotype dosage × internally estimated log(OR)
    snp_sum <- as.numeric(
      geno_matrix %*% beta_vector
    )
    
    beta_table_used$Training_Mean_Genotype <-
      train_snp_means[
        beta_table_used$SNP
      ]
    
    beta_table_used$N_Imputed_All <-
      n_imputed_by_snp[
        beta_table_used$SNP
      ]
    
    beta_table_used$Weight_Source <- switch(
      score_mode,
      "ridge_weighted" = "Training_Set_Ridge_CV_Lambda1se",
      "elasticnet_weighted" = "Training_Set_ElasticNet_Alpha05_CV_Lambda1se",
      "Training_Set_Univariate_Logistic_Regression"
    )
    
    beta_table_used$Imputation_Source <-
      "Training_Set"
    
    write.csv2(
      beta_table_used,
      weights_file,
      row.names = FALSE
    )
  }
  
  if (is.null(snp_sum)) {
    stop(
      "For score mode ",
      score_mode,
      " no score was computed."
    )
  }
  
  #------------------------------------------------------------
  # LoF-Burden-Information
  #------------------------------------------------------------
  if (score_mode == "lof_burden") {
    
    cat(
      "LoF burden score: sum of biologically annotated ",
      "LoF allele dosages.\n",
      sep = ""
    )
    
    cat(
      "No orientation based on the obesity association.\n"
    )
    
    cat(
      "SNPs inverted for LoF burden: ",
      length(invert_snps),
      "\n"
    )
    
    cat(
      "LoF orientation file: ",
      lof_orientation_file,
      "\n"
    )
  }
  
  if (
    length(snp_sum) != nrow(pheno_df)
  ) {
    stop(
      "The length of the computed score does not match the ",
      "number of individuals."
    )
  }
  
  if (any(!is.finite(snp_sum))) {
    stop(
      "The computed SNP sum score contains invalid values."
    )
  }
  
  #------------------------------------------------------------
  # Result dataset
  #------------------------------------------------------------
  
  df <- pheno_df %>%
    mutate(
      SNP_SumScore = as.numeric(snp_sum),
      Score_Mode = score_mode,
      Analysis_Config = config_name,
      N_Score_SNPs = ncol(geno_matrix),
      N_Imputed_Genotypes = n_imputed_total,
      Imputation_Method = "Training_Mean_Genotype"
    )
  
  if (score_mode %in% c("beta_weighted", "ridge_weighted", "elasticnet_weighted")) {
    
    # Unweighted dosage sum as additional information
    df$PRS_Raw <- as.numeric(prs_raw)
    
  } else {
    
    df$PRS_Raw <- NA_real_
  }
  
  #------------------------------------------------------------
  # Export individuals
  #------------------------------------------------------------
  
  score_table <- data.frame(
    id = as.character(df$id),
    SNP_SumScore = df$SNP_SumScore,
    PRS_Raw = df$PRS_Raw,
    Score_Mode = score_mode,
    Analysis_Config = config_name,
    N_Score_SNPs = ncol(geno_matrix),
    N_Imputed_Genotypes = n_imputed_total,
    Imputation_Method = "Training_Mean_Genotype",
    stringsAsFactors = FALSE
  )
  
  write.csv2(
    score_table,
    individuals_file,
    row.names = FALSE
  )
  
  #------------------------------------------------------------
  # Percentiles based exclusively on the training set
  #------------------------------------------------------------
  
  global_split <- readRDS(
    output_file(
      "Pathway_",
      pathway_name,
      "_Global_TrainTest_Split.rds"
    )
  )
  
  df_train <- df %>%
    filter(id %in% global_split$train_ids)
  
  if (nrow(df_train) == 0) {
    stop("No training individuals available for percentile calculation.")
  }
  
  percentile_table <- data.frame(
    Percentile = seq(0, 100, by = 1),
    
    SNP_SumScore = as.numeric(
      quantile(
        df_train$SNP_SumScore,
        probs = seq(0, 1, by = 0.01),
        na.rm = TRUE
      )
    ),
    
    Score_Mode = score_mode,
    Analysis_Config = config_name,
    Based_On = "Training_Set",
    stringsAsFactors = FALSE
  )
  
  write.csv2(
    percentile_table,
    percentiles_file,
    row.names = FALSE
  )
  
  #------------------------------------------------------------
  # Console output and attributes
  #------------------------------------------------------------
  
  cat(
    "Verwendete SNPs:",
    ncol(geno_matrix),
    "\n"
  )
  
  cat(
    "Genotype imputation for all score modes: ",
    "SNP means exclusively from the training set.\n",
    sep = ""
  )
  
  cat(
    "Imputierte Genotypwerte insgesamt:",
    n_imputed_total,
    "\n"
  )
  
  cat(
    "Imputation file:",
    imputation_file,
    "\n"
  )
  
  if (score_mode == "oriented") {
    
    cat(
      "Invertierte protective SNPs:",
      length(protective),
      "\n"
    )
    
    cat(
      "No β weighting.\n"
    )
  }
  
  if (score_mode == "raw") {
    
    cat(
      "No SNP inversion and no weighting.\n"
    )
  }
  
  if (score_mode == "beta_weighted") {
    
    cat(
      "Weighting: β = Estimate_Add = log(OR)\n"
    )
    
    cat(
      "No SNP inversion.\n"
    )
    
    cat(
      "Gewichtstabelle:",
      weights_file,
      "\n"
    )
  }
  
  if (score_mode %in% c("ridge_weighted", "elasticnet_weighted")) {
    
    cat(
      "Weighting: β = ",
      ifelse(score_mode == "ridge_weighted", "Ridge", "Elastic Net (alpha=0.5)"),
      " coefficient (cv.glmnet, lambda.1se), training data\n",
      sep = ""
    )
    
    cat(
      "No SNP inversion.\n"
    )
    
    cat(
      "Gewichtstabelle:",
      weights_file,
      "\n"
    )
  }
  
  if (score_mode == "lof_burden") {
    
    cat(
      "LoF burden score: sum of biologically annotated ",
      "LoF allele dosages.\n",
      sep = ""
    )
    
    cat(
      "No orientation based on the obesity association.\n"
    )
  }
  
  cat(
    "Score file:",
    individuals_file,
    "\n"
  )
  
  #------------------------------------------------------------
  # Attributes for subsequent model analysis
  #------------------------------------------------------------
  
  attr(
    df,
    "n_score_snps"
  ) <- ncol(geno_matrix)
  
  attr(
    df,
    "score_mode"
  ) <- score_mode
  
  attr(
    df,
    "analysis_config"
  ) <- config_name
  
  attr(
    df,
    "imputation_method"
  ) <- "Training_Mean_Genotype"
  
  attr(
    df,
    "n_imputed_genotypes"
  ) <- n_imputed_total
  
  attr(
    df,
    "imputation_values"
  ) <- imputation_table
  
  if (score_mode %in% c("beta_weighted", "ridge_weighted", "elasticnet_weighted")) {
    
    attr(
      df,
      "beta_weights"
    ) <- beta_table_used
  }
  
  return(df)
}

#===================================================================================
# EXECUTION
#===================================================================================

all_snp_sumscore_results <- list()

for (i in seq_len(nrow(analysis_plan))) {
  
  config_name <- analysis_plan$config_name[i]
  score_mode  <- analysis_plan$score_mode[i]
  
  cat(
    "\n----------------------------------------------------\n",
    "Starte Analyse ",
    i,
    " of ",
    nrow(analysis_plan),
    ": ",
    config_name,
    " / ",
    score_mode,
    "\n",
    "----------------------------------------------------\n",
    sep = ""
  )
  
  #------------------------------------------------------------
  # Check configuration
  #------------------------------------------------------------
  
  if (!config_name %in% names(analysis_configs)) {
    
    warning(
      "Analysis configuration not found: ",
      config_name
    )
    
    next
  }
  
  config <- analysis_configs[[config_name]]
  
  #------------------------------------------------------------
  # Safeguard permitted LoF combination once more
  #------------------------------------------------------------
  
  if (
    score_mode == "lof_burden" &&
    config$type != "lof_vep"
  ) {
    
    warning(
      "Invalid combination skipped: ",
      config_name,
      " / ",
      score_mode,
      ". 'lof_burden' is only permitted for type = 'lof_vep'."
    )
    
    next
  }
  
  if (
    config$type == "lof_vep" &&
    score_mode != "lof_burden"
  ) {
    
    warning(
      "Invalid combination skipped: ",
      config_name,
      " / ",
      score_mode,
      ". A LoF configuration may only be used with ",
      "score_mode = 'lof_burden'."
    )
    
    next
  }
  
  #------------------------------------------------------------
  # Result name
  #------------------------------------------------------------
  
  result_name <- paste(
    config_name,
    score_mode,
    sep = "_"
  )
  
  #------------------------------------------------------------
  # Run analysis
  #------------------------------------------------------------
  
  result <- tryCatch(
    
    run_snp_sumscore(
      config_name = config_name,
      config = config,
      score_mode = score_mode
    ),
    
    error = function(e) {
      
      warning(
        "Error for ",
        result_name,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  #------------------------------------------------------------
  # Save result
  #------------------------------------------------------------
  
  if (is.null(result)) {
    
    warning(
      "For ",
      result_name,
      " no result was produced."
    )
    
    next
  }
  
  all_snp_sumscore_results[[result_name]] <- result
  
  cat(
    "Analyse erfolgreich abgeschlossen: ",
    result_name,
    "\n",
    sep = ""
  )
}

cat(
  "\n====================================================\n",
  "SNP-SumScore-Analysen abgeschlossen\n",
  "Erfolgreiche Analysen: ",
  length(all_snp_sumscore_results),
  " of ",
  nrow(analysis_plan),
  "\n",
  "====================================================\n",
  sep = ""
)

print(
  names(all_snp_sumscore_results)
)

#===================================================================================
# UNIVARIATE SNP SUM SCORE – LOGISTIC REGRESSION / TRAIN-TEST / ROC / CUTOFF
#===================================================================================
library(dplyr)
library(caret)
library(pROC)

#===================================================================================
# HILFSFUNKTIONEN
#===================================================================================

calc_metrics <- function(pred, obs, cutoff) {
  
  pred <- as.numeric(pred)
  obs  <- as.numeric(obs)
  
  valid <- complete.cases(pred, obs)
  
  pred <- pred[valid]
  obs  <- obs[valid]
  
  if (
    length(pred) == 0 ||
    !is.finite(cutoff)
  ) {
    
    return(
      data.frame(
        TN = NA_integer_,
        FP = NA_integer_,
        FN = NA_integer_,
        TP = NA_integer_,
        Sensitivity = NA_real_,
        Specificity = NA_real_,
        PPV = NA_real_,
        NPV = NA_real_,
        Accuracy = NA_real_
      )
    )
  }
  
  pred_class <- ifelse(
    pred >= cutoff,
    1,
    0
  )
  
  cm <- table(
    Predicted = factor(
      pred_class,
      levels = c(0, 1)
    ),
    Observed = factor(
      obs,
      levels = c(0, 1)
    )
  )
  
  TN <- unname(cm[1, 1])
  FP <- unname(cm[2, 1])
  FN <- unname(cm[1, 2])
  TP <- unname(cm[2, 2])
  
  data.frame(
    TN = TN,
    FP = FP,
    FN = FN,
    TP = TP,
    
    Sensitivity = ifelse(
      (TP + FN) > 0,
      TP / (TP + FN),
      NA_real_
    ),
    
    Specificity = ifelse(
      (TN + FP) > 0,
      TN / (TN + FP),
      NA_real_
    ),
    
    PPV = ifelse(
      (TP + FP) > 0,
      TP / (TP + FP),
      NA_real_
    ),
    
    NPV = ifelse(
      (TN + FN) > 0,
      TN / (TN + FN),
      NA_real_
    ),
    
    Accuracy = ifelse(
      sum(cm) > 0,
      (TP + TN) / sum(cm),
      NA_real_
    )
  )
}


calc_score_cutoff <- function(
    prob_cutoff,
    beta0,
    beta1) {
  
  if (
    !is.finite(prob_cutoff) ||
    prob_cutoff <= 0 ||
    prob_cutoff >= 1 ||
    !is.finite(beta0) ||
    !is.finite(beta1) ||
    beta1 == 0
  ) {
    return(NA_real_)
  }
  
  (
    qlogis(prob_cutoff) - beta0
  ) / beta1
}


get_score_cutoff_direction <- function(beta1) {
  
  if (!is.finite(beta1) || beta1 == 0) {
    return(NA_character_)
  }
  
  if (beta1 > 0) {
    return("Obesity predicted when score >= cutoff")
  }
  
  "Obesity predicted when score <= cutoff"
}


safe_roc <- function(obs, pred) {
  
  valid <- complete.cases(obs, pred)
  
  obs  <- obs[valid]
  pred <- pred[valid]
  
  if (
    length(obs) == 0 ||
    length(unique(obs)) < 2 ||
    length(unique(pred)) < 2
  ) {
    return(NULL)
  }
  
  tryCatch(
    pROC::roc(
      response = obs,
      predictor = pred,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    ),
    error = function(e) NULL
  )
}


get_n_score_snps <- function(df) {
  
  n_score_snps <- attr(
    df,
    "n_score_snps"
  )
  
  if (
    is.null(n_score_snps) ||
    length(n_score_snps) == 0
  ) {
    return(NA_integer_)
  }
  
  as.integer(n_score_snps[1])
}

#===================================================================================
# MODELLFUNKTION
#===================================================================================

run_snp_sumscore_model <- function(
    df,
    result_name) {
  
  cat("\n====================================================\n")
  cat("SNP sum score model:", result_name, "\n")
  cat("====================================================\n")
  
  if (is.null(df)) {
    warning("No dataset for ", result_name)
    return(NULL)
  }
  
  required_columns <- c(
    "id",
    "obesity",
    "SNP_SumScore",
    "Score_Mode",
    "Analysis_Config"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(df)
  )
  
  if (length(missing_columns) > 0) {
    
    warning(
      "The following columns are missing for ",
      result_name,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
    
    return(NULL)
  }
  
  analysis_configs_present <- unique(
    as.character(df$Analysis_Config)
  )
  
  score_modes_present <- unique(
    as.character(df$Score_Mode)
  )
  
  analysis_configs_present <- analysis_configs_present[
    !is.na(analysis_configs_present)
  ]
  
  score_modes_present <- score_modes_present[
    !is.na(score_modes_present)
  ]
  
  if (length(analysis_configs_present) != 1) {
    
    warning(
      "Not exactly one Analysis_Config in ",
      result_name
    )
    
    return(NULL)
  }
  
  if (length(score_modes_present) != 1) {
    
    warning(
      "Not exactly one Score_Mode in ",
      result_name
    )
    
    return(NULL)
  }
  
  analysis_config <- analysis_configs_present[1]
  score_mode      <- score_modes_present[1]
  
  if (!analysis_config %in% names(analysis_configs)) {
    
    warning(
      "Unbekannte Analysis_Config: ",
      analysis_config
    )
    
    return(NULL)
  }
  
  config <- analysis_configs[[analysis_config]]
  
  n_score_snps <- get_n_score_snps(df)
  
  #-------------------------------------------------------------------------
  # Consistency check for the LoF analysis
  #-------------------------------------------------------------------------
  
  if (
    score_mode == "lof_burden" &&
    config$type != "lof_vep"
  ) {
    warning(
      "score_mode = 'lof_burden' is only permitted for ",
      "config$type = 'lof_vep': ",
      result_name
    )
    return(NULL)
  }
  
  if (
    config$type == "lof_vep" &&
    score_mode != "lof_burden"
  ) {
    warning(
      "A LoF VEP configuration may only be used with ",
      "score_mode = 'lof_burden' ausgewertet werden: ",
      result_name
    )
    return(NULL)
  }
  
  #-------------------------------------------------------------------------
  # Complete cases
  #-------------------------------------------------------------------------
  
  df_clean <- df %>%
    mutate(
      id = as.character(id),
      obesity = as.numeric(as.character(obesity)),
      SNP_SumScore = as.numeric(SNP_SumScore)
    ) %>%
    filter(
      !is.na(id),
      id != "",
      !is.na(obesity),
      is.finite(SNP_SumScore)
    ) %>%
    distinct(
      id,
      .keep_all = TRUE
    )
  
  if (nrow(df_clean) < 20) {
    warning("Too few complete cases for ", result_name)
    return(NULL)
  }
  
  df_clean$obesity_bin <- ifelse(
    df_clean$obesity == 1,
    1,
    0
  )
  
  if (length(unique(df_clean$obesity_bin)) < 2) {
    warning("Only one outcome class available for ", result_name)
    return(NULL)
  }
  
  #-------------------------------------------------------------------------
  # Apply global train/test split
  #-------------------------------------------------------------------------
  
  global_split <- readRDS(
    output_file(
      "Pathway_",
      pathway_name,
      "_Global_TrainTest_Split.rds"
    )
  )
  
  train_set <- df_clean %>%
    filter(id %in% global_split$train_ids)
  
  test_set <- df_clean %>%
    filter(id %in% global_split$test_ids)
  
  if (
    nrow(train_set) < 10 ||
    nrow(test_set) < 5
  ) {
    
    warning(
      "Train or test dataset is too small for ",
      result_name
    )
    
    return(NULL)
  }
  
  if (
    length(unique(train_set$obesity_bin)) < 2 ||
    length(unique(test_set$obesity_bin)) < 2
  ) {
    
    warning(
      "Train or test contains only one outcome class for ",
      result_name
    )
    
    return(NULL)
  }
  
  if (
    length(unique(train_set$SNP_SumScore)) < 2
  ) {
    
    warning(
      "The SNP sum score has no variation in the training dataset for ",
      result_name
    )
    
    return(NULL)
  }
  
  #-------------------------------------------------------------------------
  # Score distribution in the training dataset
  #-------------------------------------------------------------------------
  
  score_mean_train <- mean(
    train_set$SNP_SumScore,
    na.rm = TRUE
  )
  
  score_sd_train <- sd(
    train_set$SNP_SumScore,
    na.rm = TRUE
  )
  
  if (
    !is.finite(score_sd_train) ||
    score_sd_train <= 0
  ) {
    
    warning(
      "The standard deviation of the score is invalid for ",
      result_name
    )
    
    return(NULL)
  }
  
  # Standardized version:
  # exclusively with mean and SD of the training dataset
  train_set$SNP_SumScore_SD <- (
    train_set$SNP_SumScore - score_mean_train
  ) / score_sd_train
  
  test_set$SNP_SumScore_SD <- (
    test_set$SNP_SumScore - score_mean_train
  ) / score_sd_train
  
  #-------------------------------------------------------------------------
  # Model with score in original units
  #-------------------------------------------------------------------------
  
  model_train <- tryCatch(
    glm(
      obesity_bin ~ SNP_SumScore,
      data = train_set,
      family = binomial
    ),
    error = function(e) NULL
  )
  
  if (is.null(model_train)) {
    warning("Model could not be estimated for ", result_name)
    return(NULL)
  }
  
  coef_table <- summary(model_train)$coefficients
  
  if (!"SNP_SumScore" %in% rownames(coef_table)) {
    
    warning(
      "No SNP sum score coefficient available for ",
      result_name
    )
    
    return(NULL)
  }
  
  beta0 <- unname(
    coef(model_train)["(Intercept)"]
  )
  
  beta1 <- unname(
    coef_table[
      "SNP_SumScore",
      "Estimate"
    ]
  )
  
  se1 <- unname(
    coef_table[
      "SNP_SumScore",
      "Std. Error"
    ]
  )
  
  wald_pvalue <- unname(
    coef_table[
      "SNP_SumScore",
      "Pr(>|z|)"
    ]
  )
  
  OR <- exp(beta1)
  
  ci <- tryCatch(
    exp(
      confint.default(model_train)[
        "SNP_SumScore",
      ]
    ),
    error = function(e) c(NA_real_, NA_real_)
  )
  
  #-------------------------------------------------------------------------
  # Effect per one standard deviation
  #-------------------------------------------------------------------------
  
  model_train_sd <- tryCatch(
    glm(
      obesity_bin ~ SNP_SumScore_SD,
      data = train_set,
      family = binomial
    ),
    error = function(e) NULL
  )
  
  beta_per_sd <- NA_real_
  se_per_sd   <- NA_real_
  OR_per_sd   <- NA_real_
  ci_per_sd   <- c(NA_real_, NA_real_)
  p_per_sd    <- NA_real_
  
  if (!is.null(model_train_sd)) {
    
    coef_sd <- summary(model_train_sd)$coefficients
    
    if ("SNP_SumScore_SD" %in% rownames(coef_sd)) {
      
      beta_per_sd <- coef_sd[
        "SNP_SumScore_SD",
        "Estimate"
      ]
      
      se_per_sd <- coef_sd[
        "SNP_SumScore_SD",
        "Std. Error"
      ]
      
      p_per_sd <- coef_sd[
        "SNP_SumScore_SD",
        "Pr(>|z|)"
      ]
      
      OR_per_sd <- exp(beta_per_sd)
      
      ci_per_sd <- tryCatch(
        exp(
          confint.default(model_train_sd)[
            "SNP_SumScore_SD",
          ]
        ),
        error = function(e) c(
          NA_real_,
          NA_real_
        )
      )
    }
  }
  
  #-------------------------------------------------------------------------
  # Null model and likelihood-ratio test
  #-------------------------------------------------------------------------
  
  model_null <- glm(
    obesity_bin ~ 1,
    data = train_set,
    family = binomial
  )
  
  lrt_pvalue <- tryCatch(
    anova(
      model_null,
      model_train,
      test = "LRT"
    )$`Pr(>Chi)`[2],
    error = function(e) NA_real_
  )
  
  #-------------------------------------------------------------------------
  # Predictions
  #-------------------------------------------------------------------------
  
  train_set$pred_prob <- predict(
    model_train,
    newdata = train_set,
    type = "response"
  )
  
  test_set$pred_prob <- predict(
    model_train,
    newdata = test_set,
    type = "response"
  )
  
  #-------------------------------------------------------------------------
  # ROC
  #-------------------------------------------------------------------------
  
  roc_train <- safe_roc(
    obs = train_set$obesity_bin,
    pred = train_set$pred_prob
  )
  
  roc_obj <- safe_roc(
    obs = test_set$obesity_bin,
    pred = test_set$pred_prob
  )
  
  if (
    is.null(roc_train) ||
    is.null(roc_obj)
  ) {
    
    warning(
      "ROC could not be computed for ",
      result_name
    )
    
    return(NULL)
  }
  
  auc_val <- as.numeric(
    auc(roc_obj)
  )
  
  auc_ci <- tryCatch(
    as.numeric(
      ci.auc(roc_obj)
    ),
    error = function(e) c(
      NA_real_,
      NA_real_,
      NA_real_
    )
  )
  
  #-------------------------------------------------------------------------
  # Cutoffs
  #-------------------------------------------------------------------------
  
  cutoff_0.5 <- 0.5
  
  cutoff_youden <- tryCatch(
    {
      
      cutoff_result <- coords(
        roc_train,
        x = "best",
        ret = "threshold",
        best.method = "youden",
        transpose = FALSE
      )
      
      as.numeric(
        cutoff_result[1]
      )
    },
    error = function(e) NA_real_
  )
  
  score_cutoff_0.5 <- calc_score_cutoff(
    prob_cutoff = cutoff_0.5,
    beta0 = beta0,
    beta1 = beta1
  )
  
  score_cutoff_youden <- calc_score_cutoff(
    prob_cutoff = cutoff_youden,
    beta0 = beta0,
    beta1 = beta1
  )
  
  cutoff_direction <- get_score_cutoff_direction(
    beta1
  )
  
  metrics_0.5 <- calc_metrics(
    pred = test_set$pred_prob,
    obs = test_set$obesity_bin,
    cutoff = cutoff_0.5
  )
  
  metrics_youden <- calc_metrics(
    pred = test_set$pred_prob,
    obs = test_set$obesity_bin,
    cutoff = cutoff_youden
  )
  
  #-------------------------------------------------------------------------
  # Batch sensitivity analysis
  #
  # Main analysis only, not sensitivity_betaweighted.
  #-------------------------------------------------------------------------
  
  batch_sensitivity_results <- NULL
  
  if (
    analysis_config == "main_ld" &&
    "Batch" %in% colnames(train_set)
  ) {
    
    train_set$Batch <- factor(
      train_set$Batch
    )
    
    test_set$Batch <- factor(
      test_set$Batch,
      levels = levels(train_set$Batch)
    )
    
    if (
      nlevels(train_set$Batch) >= 2
    ) {
      
      model_batch <- tryCatch(
        glm(
          obesity_bin ~ SNP_SumScore + Batch,
          data = train_set,
          family = binomial
        ),
        error = function(e) NULL
      )
      
      if (!is.null(model_batch)) {
        
        batch_coef <- summary(
          model_batch
        )$coefficients
        
        if (
          "SNP_SumScore" %in%
          rownames(batch_coef)
        ) {
          
          batch_ci <- tryCatch(
            exp(
              confint.default(model_batch)[
                "SNP_SumScore",
              ]
            ),
            error = function(e) c(
              NA_real_,
              NA_real_
            )
          )
          
          batch_beta <- batch_coef[
            "SNP_SumScore",
            "Estimate"
          ]
          
          batch_sensitivity_results <- data.frame(
            Pathway = pathway_name,
            Analysis = result_name,
            Score_Mode = score_mode,
            Analysis_Config = analysis_config,
            Sensitivity_Type = "Batch_Adjusted",
            N_Train = nobs(model_batch),
            
            Estimate = batch_beta,
            StdError = batch_coef[
              "SNP_SumScore",
              "Std. Error"
            ],
            
            OR = exp(batch_beta),
            CI_lower = batch_ci[1],
            CI_upper = batch_ci[2],
            
            Pvalue = batch_coef[
              "SNP_SumScore",
              "Pr(>|z|)"
            ],
            
            OR_per_SD = exp(
              batch_beta * score_sd_train
            ),
            
            stringsAsFactors = FALSE
          )
          
          batch_file <- make_output_name(
            pathway_name,
            "Univariate",
            "SNP_SumScore",
            "Batch_Sensitivity_Results",
            config,
            score_mode
          )
          
          write.csv2(
            batch_sensitivity_results,
            batch_file,
            row.names = FALSE
          )
          
          cat(
            "Batch sensitivity analysis saved.\n"
          )
        }
      }
    }
  }
  
  #-------------------------------------------------------------------------
  # Main results
  #-------------------------------------------------------------------------
  
  results_out <- data.frame(
    
    Pathway = rep(
      pathway_name,
      2
    ),
    
    Analysis = rep(
      result_name,
      2
    ),
    
    Score_Mode = rep(
      score_mode,
      2
    ),
    
    Analysis_Config = rep(
      analysis_config,
      2
    ),
    
    Score_Type = rep(
      dplyr::case_when(
        score_mode == "beta_weighted" ~
          "Beta-weighted PRS",
        
        score_mode == "ridge_weighted" ~
          "Ridge-weighted PRS",
        
        score_mode == "elasticnet_weighted" ~
          "Elastic-Net-weighted PRS",
        
        score_mode == "lof_burden" ~
          "LoF allele burden",
        
        score_mode == "oriented" ~
          "Obesity-oriented SNP SumScore",
        
        score_mode == "raw" ~
          "Unweighted raw SNP SumScore",
        
        TRUE ~
          "Other SNP Score"
      ),
      2
    ),
    
    Score_Unit = rep(
      dplyr::case_when(
        score_mode == "lof_burden" ~
          "One additional LoF allele copy",
        
        score_mode == "beta_weighted" ~
          "One unit of beta-weighted score",
        
        score_mode == "ridge_weighted" ~
          "One unit of ridge-weighted score",
        
        score_mode == "elasticnet_weighted" ~
          "One unit of elastic-net-weighted score",
        
        score_mode == "oriented" ~
          "One additional obesity-oriented allele",
        
        score_mode == "raw" ~
          "One additional counted allele",
        
        TRUE ~
          "One score unit"
      ),
      2
    ),
    
    N_Score_SNPs = rep(
      n_score_snps,
      2
    ),
    
    N_Total = rep(
      nrow(df_clean),
      2
    ),
    
    N_Train = rep(
      nrow(train_set),
      2
    ),
    
    N_Test = rep(
      nrow(test_set),
      2
    ),
    
    N_Train_NonObese = rep(
      sum(train_set$obesity_bin == 0),
      2
    ),
    
    N_Train_Obese = rep(
      sum(train_set$obesity_bin == 1),
      2
    ),
    
    N_Test_NonObese = rep(
      sum(test_set$obesity_bin == 0),
      2
    ),
    
    N_Test_Obese = rep(
      sum(test_set$obesity_bin == 1),
      2
    ),
    
    Score_Mean_Train = rep(
      score_mean_train,
      2
    ),
    
    Score_SD_Train = rep(
      score_sd_train,
      2
    ),
    
    Estimate_per_Unit = rep(
      beta1,
      2
    ),
    
    StdError_per_Unit = rep(
      se1,
      2
    ),
    
    OR_per_Unit = rep(
      OR,
      2
    ),
    
    CI_lower_per_Unit = rep(
      ci[1],
      2
    ),
    
    CI_upper_per_Unit = rep(
      ci[2],
      2
    ),
    
    Wald_pvalue = rep(
      wald_pvalue,
      2
    ),
    
    Estimate_per_SD = rep(
      beta_per_sd,
      2
    ),
    
    StdError_per_SD = rep(
      se_per_sd,
      2
    ),
    
    OR_per_SD = rep(
      OR_per_sd,
      2
    ),
    
    CI_lower_per_SD = rep(
      ci_per_sd[1],
      2
    ),
    
    CI_upper_per_SD = rep(
      ci_per_sd[2],
      2
    ),
    
    Pvalue_per_SD = rep(
      p_per_sd,
      2
    ),
    
    LRT_pvalue = rep(
      lrt_pvalue,
      2
    ),
    
    AIC = rep(
      AIC(model_train),
      2
    ),
    
    AUC = rep(
      auc_val,
      2
    ),
    
    AUC_CI_lower = rep(
      auc_ci[1],
      2
    ),
    
    AUC_CI_median = rep(
      auc_ci[2],
      2
    ),
    
    AUC_CI_upper = rep(
      auc_ci[3],
      2
    ),
    
    Cutoff_Type = c(
      "Probability_0.5",
      "Youden_Training"
    ),
    
    Cutoff_Prob = c(
      cutoff_0.5,
      cutoff_youden
    ),
    
    SNP_SumScore_Cutoff = c(
      score_cutoff_0.5,
      score_cutoff_youden
    ),
    
    Score_Cutoff_Direction = rep(
      cutoff_direction,
      2
    ),
    
    TN = c(
      metrics_0.5$TN,
      metrics_youden$TN
    ),
    
    FP = c(
      metrics_0.5$FP,
      metrics_youden$FP
    ),
    
    FN = c(
      metrics_0.5$FN,
      metrics_youden$FN
    ),
    
    TP = c(
      metrics_0.5$TP,
      metrics_youden$TP
    ),
    
    Sensitivity = c(
      metrics_0.5$Sensitivity,
      metrics_youden$Sensitivity
    ),
    
    Specificity = c(
      metrics_0.5$Specificity,
      metrics_youden$Specificity
    ),
    
    PPV = c(
      metrics_0.5$PPV,
      metrics_youden$PPV
    ),
    
    NPV = c(
      metrics_0.5$NPV,
      metrics_youden$NPV
    ),
    
    Accuracy = c(
      metrics_0.5$Accuracy,
      metrics_youden$Accuracy
    ),
    
    stringsAsFactors = FALSE
  )
  
  #-------------------------------------------------------------------------
  # Save regression results
  #-------------------------------------------------------------------------
  
  regression_file <- make_output_name(
    pathway_name,
    "Univariate",
    "SNP_SumScore",
    "Regression_Results",
    config,
    score_mode
  )
  
  write.csv2(
    results_out,
    regression_file,
    row.names = FALSE
  )
  
  #-------------------------------------------------------------------------
  # Save test predictions
  #-------------------------------------------------------------------------
  
  prediction_columns <- c(
    "id",
    "obesity",
    "obesity_bin",
    "SNP_SumScore",
    "SNP_SumScore_SD",
    "pred_prob",
    "Score_Mode",
    "Analysis_Config"
  )
  
  # Additionally keep PRS_Raw for the β-weighted score
  if ("PRS_Raw" %in% colnames(test_set)) {
    prediction_columns <- append(
      prediction_columns,
      "PRS_Raw",
      after = 4
    )
  }
  
  prediction_table <- test_set %>%
    dplyr::select(
      all_of(prediction_columns)
    )
  
  prediction_table$Predicted_Class_0.5 <- ifelse(
    prediction_table$pred_prob >= cutoff_0.5,
    1,
    0
  )
  
  prediction_table$Predicted_Class_Youden <- ifelse(
    prediction_table$pred_prob >= cutoff_youden,
    1,
    0
  )
  
  predictions_file <- make_output_name(
    pathway_name,
    "Univariate",
    "SNP_SumScore",
    "Predictions_TestSet",
    config,
    score_mode
  )
  
  write.csv2(
    prediction_table,
    predictions_file,
    row.names = FALSE
  )
  
  #-------------------------------------------------------------------------
  # Save ROC data
  #-------------------------------------------------------------------------
  
  roc_df <- data.frame(
    Threshold = roc_obj$thresholds,
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    Specificity = roc_obj$specificities,
    Sensitivity = roc_obj$sensitivities,
    Analysis = result_name,
    Score_Mode = score_mode,
    Analysis_Config = analysis_config,
    stringsAsFactors = FALSE
  )
  
  roc_file <- make_output_name(
    pathway_name,
    "Univariate",
    "SNP_SumScore",
    "ROC_Data",
    config,
    score_mode
  )
  
  write.csv2(
    roc_df,
    roc_file,
    row.names = FALSE
  )
  
  #-------------------------------------------------------------------------
  # Console output
  #-------------------------------------------------------------------------
  if (score_mode == "lof_burden") {
    
    cat(
      "LoF-Modell: biologisch ausgerichteter Allel-Burden.\n"
    )
    
    cat(
      "Interpretation OR per unit: OR per additional ",
      "LoF-Allelkopie.\n",
      sep = ""
    )
  }
  cat("Score mode:", score_mode, "\n")
  cat("Analysis configuration:", analysis_config, "\n")
  cat("Score-SNPs:", n_score_snps, "\n")
  cat("AUC Test:", round(auc_val, 3), "\n")
  cat("OR per score unit:", round(OR, 3), "\n")
  cat("OR pro SD:", round(OR_per_sd, 3), "\n")
  cat("Saved:", regression_file, "\n")
  
  return(
    list(
      results = results_out,
      model = model_train,
      model_standardized = model_train_sd,
      roc_train = roc_train,
      roc = roc_obj,
      predictions = prediction_table,
      roc_data = roc_df,
      batch_sensitivity = batch_sensitivity_results,
      score_mean_train = score_mean_train,
      score_sd_train = score_sd_train,
      n_score_snps = n_score_snps
    )
  )
}

#===================================================================================
# LOOP OVER ALL SNP SUM SCORE RESULTS
#===================================================================================

all_snp_sumscore_model_results <- list()

for (result_name in names(all_snp_sumscore_results)) {
  
  df_current <- all_snp_sumscore_results[[result_name]]
  
  cat(
    "\n----------------------------------------------------\n",
    "Starte SNP-SumScore-Modell: ",
    result_name,
    "\n",
    "----------------------------------------------------\n",
    sep = ""
  )
  
  if (is.null(df_current)) {
    
    warning(
      "No SNP sum score dataset for ",
      result_name,
      ". Analysis is skipped."
    )
    
    next
  }
  
  model_result <- tryCatch(
    
    run_snp_sumscore_model(
      df = df_current,
      result_name = result_name
    ),
    
    error = function(e) {
      
      warning(
        "Fehler in ",
        result_name,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  if (is.null(model_result)) {
    
    warning(
      "For ",
      result_name,
      " no model result was produced."
    )
    
    next
  }
  
  all_snp_sumscore_model_results[[result_name]] <-
    model_result
  
  cat(
    "Modellanalyse erfolgreich abgeschlossen: ",
    result_name,
    "\n",
    sep = ""
  )
}

#===================================================================================
# EXPORT OVERALL SUMMARY
#===================================================================================

snp_sumscore_model_summary <- bind_rows(
  lapply(
    all_snp_sumscore_model_results,
    function(x) {
      
      if (
        is.null(x) ||
        is.null(x$results)
      ) {
        return(NULL)
      }
      
      x$results
    }
  )
)

if (nrow(snp_sumscore_model_summary) > 0) {
  
  summary_file <- output_file(
    "Pathway_",
    pathway_name,
    "_Univariate_SNP_SumScore_Model_Results_All_Configurations.csv"
  )
  
  write.csv2(
    snp_sumscore_model_summary,
    summary_file,
    row.names = FALSE
  )
  
  cat("\n====================================================\n")
  cat("Overall summary saved:\n")
  cat(summary_file, "\n")
  cat("Successful model analyses:\n")
  
  summary_columns <- c(
    "Analysis_Config",
    "Score_Mode",
    "Score_Type"
  )
  
  # If Score_Unit was added to results_out,
  # it is also shown in the summary.
  if (
    "Score_Unit" %in%
    colnames(snp_sumscore_model_summary)
  ) {
    summary_columns <- c(
      summary_columns,
      "Score_Unit"
    )
  }
  
  print(
    unique(
      snp_sumscore_model_summary[
        ,
        summary_columns,
        drop = FALSE
      ]
    )
  )
  
  #------------------------------------------------------------
  # Additional check of the LoF analysis
  #------------------------------------------------------------
  
  lof_summary <- snp_sumscore_model_summary %>%
    filter(
      Score_Mode == "lof_burden"
    )
  
  if (nrow(lof_summary) > 0) {
    
    cat("\nLoF burden analysis successfully included:\n")
    
    print(
      lof_summary[
        ,
        intersect(
          c(
            "Analysis",
            "Analysis_Config",
            "Score_Mode",
            "Score_Type",
            "Score_Unit",
            "N_Score_SNPs",
            "N_Train",
            "N_Test",
            "OR_per_Unit",
            "CI_lower_per_Unit",
            "CI_upper_per_Unit",
            "Wald_pvalue",
            "OR_per_SD",
            "Pvalue_per_SD",
            "AUC"
          ),
          colnames(lof_summary)
        ),
        drop = FALSE
      ]
    )
    
  } else {
    
    warning(
      "The overall summary contains no ",
      "LoF-Burden-Analyse enthalten."
    )
  }
  
  cat(
    "\nAnzahl erfolgreicher Analysen: ",
    length(all_snp_sumscore_model_results),
    "\n",
    sep = ""
  )
  
} else {
  
  warning(
    "No SNP sum score model results could ",
    "be merged."
  )
}

#===================================================================================
# OPTIONAL SEPARATE SUMMARY: ONE RESULT PER ANALYSIS
# Youden row as primary performance row
#===================================================================================

if (
  exists("snp_sumscore_model_summary") &&
  nrow(snp_sumscore_model_summary) > 0
) {
  
  snp_sumscore_model_primary_summary <-
    snp_sumscore_model_summary %>%
    filter(
      Cutoff_Type == "Youden_Training"
    )
  
  if (
    nrow(
      snp_sumscore_model_primary_summary
    ) > 0
  ) {
    
    write.csv2(
      snp_sumscore_model_primary_summary,
      output_file(
        "Pathway_",
        pathway_name,
        "_Univariate_SNP_SumScore_Model_Primary_Results_All_Configurations.csv"
      ),
      row.names = FALSE
    )
  }
}

#===================================================================================
# UNIVARIATE SNP SUM SCORE – PLOTS
#===================================================================================

library(ggplot2)
library(dplyr)

run_snp_sumscore_plots <- function(
    df,
    model_result,
    result_name) {
  
  cat("\n====================================================\n")
  cat("Plots for:", result_name, "\n")
  cat("====================================================\n")
  
  #------------------------------------------------------------
  # Basic checks
  #------------------------------------------------------------
  
  if (
    is.null(df) ||
    is.null(model_result)
  ) {
    
    warning(
      "No plot data for ",
      result_name
    )
    
    return(NULL)
  }
  
  required_columns <- c(
    "SNP_SumScore",
    "Score_Mode",
    "Analysis_Config"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(df)
  )
  
  if (length(missing_columns) > 0) {
    
    warning(
      "The following columns are missing for the plots of ",
      result_name,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
    
    return(NULL)
  }
  
  analysis_config_values <- unique(
    as.character(df$Analysis_Config)
  )
  
  score_mode_values <- unique(
    as.character(df$Score_Mode)
  )
  
  analysis_config_values <-
    analysis_config_values[
      !is.na(analysis_config_values)
    ]
  
  score_mode_values <-
    score_mode_values[
      !is.na(score_mode_values)
    ]
  
  if (
    length(analysis_config_values) != 1 ||
    length(score_mode_values) != 1
  ) {
    
    warning(
      "Analysis_Config or Score_Mode is not unique for ",
      result_name,
      ""
    )
    
    return(NULL)
  }
  
  analysis_config <- analysis_config_values[1]
  score_mode <- score_mode_values[1]
  
  if (
    !analysis_config %in%
    names(analysis_configs)
  ) {
    
    warning(
      "Unbekannte Analysekonfiguration: ",
      analysis_config
    )
    
    return(NULL)
  }
  
  config <- analysis_configs[[analysis_config]]
  
  #------------------------------------------------------------
  # Consistency check for LoF plots
  #------------------------------------------------------------
  
  if (
    score_mode == "lof_burden" &&
    config$type != "lof_vep"
  ) {
    warning(
      "score_mode = 'lof_burden' is only permitted for ",
      "config$type = 'lof_vep': ",
      result_name
    )
    return(NULL)
  }
  
  if (
    config$type == "lof_vep" &&
    score_mode != "lof_burden"
  ) {
    warning(
      "A LoF VEP configuration may only be used with ",
      "score_mode = 'lof_burden' visualisiert werden: ",
      result_name
    )
    return(NULL)
  }
  
  #------------------------------------------------------------
  # Dynamic score labels
  #------------------------------------------------------------
  
  score_label <- switch(
    score_mode,
    
    oriented =
      "Oriented SNP sum score",
    
    raw =
      "Raw SNP sum score",
    
    beta_weighted =
      expression(
        beta * "-weighted genetic score"
      ),
    
    ridge_weighted =
      "Ridge-weighted genetic score",
    
    elasticnet_weighted =
      "Elastic-Net-weighted genetic score",
    
    lof_burden =
      "LoF allele burden",
    
    "Genetic score"
  )
  
  score_label_text <- switch(
    score_mode,
    
    oriented =
      "Oriented SNP SumScore",
    
    raw =
      "Raw SNP SumScore",
    
    beta_weighted =
      "Beta-weighted Genetic Score",
    
    ridge_weighted =
      "Ridge-weighted Genetic Score",
    
    elasticnet_weighted =
      "Elastic-Net-weighted Genetic Score",
    
    lof_burden =
      "Loss-of-Function Allele Burden",
    
    "Genetic Score"
  )
  
  score_description <- switch(
    score_mode,
    
    oriented =
      "Protective alleles reoriented; unweighted",
    
    raw =
      "Original allele dosage; unweighted",
    
    beta_weighted =
      "Allele dosage weighted by training-set log(OR)",
    
    ridge_weighted =
      "Allele dosage weighted by ridge-regularized (cv.glmnet, lambda.1se) training-set coefficients",
    
    elasticnet_weighted =
      "Allele dosage weighted by elastic-net-regularized (alpha=0.5, cv.glmnet, lambda.1se) training-set coefficients",
    
    lof_burden =
      paste0(
        "VEP-annotated LoF alleles; dosage oriented to the ",
        "LoF allele; unweighted"
      ),
    
    ""
  )
  
  n_score_snps <- attr(
    df,
    "n_score_snps"
  )
  
  if (
    is.null(n_score_snps) &&
    "N_Score_SNPs" %in%
    colnames(df)
  ) {
    
    n_score_snps <- unique(
      df$N_Score_SNPs
    )
    
    n_score_snps <- n_score_snps[
      !is.na(n_score_snps)
    ]
    
    if (length(n_score_snps) > 0) {
      n_score_snps <- n_score_snps[1]
    }
  }
  
  if (
    is.null(n_score_snps) ||
    length(n_score_snps) == 0
  ) {
    n_score_snps <- NA_integer_
  }
  
  plot_subtitle <- paste0(
    score_description,
    ifelse(
      is.finite(n_score_snps),
      paste0(" | SNPs: ", n_score_snps),
      ""
    )
  )
  
  #------------------------------------------------------------
  # Prepare plot data
  #------------------------------------------------------------
  
  df_plot <- df %>%
    mutate(
      SNP_SumScore = as.numeric(
        SNP_SumScore
      )
    ) %>%
    filter(
      is.finite(SNP_SumScore)
    )
  
  if (nrow(df_plot) == 0) {
    
    warning(
      "No valid score values for ",
      result_name
    )
    
    return(NULL)
  }
  
  df_obesity <- df_plot
  
  if (
    "obesity" %in%
    colnames(df_obesity)
  ) {
    
    df_obesity <- df_obesity %>%
      filter(
        !is.na(obesity)
      )
    
    df_obesity$obesity_group <- factor(
      df_obesity$obesity,
      levels = c(0, 1),
      labels = c(
        "Non-obese",
        "Obese"
      )
    )
  }
  
  #------------------------------------------------------------
  # 1) Histogram overall
  #------------------------------------------------------------
  
  p_all <- ggplot(
    df_plot,
    aes(
      x = SNP_SumScore
    )
  ) +
    geom_histogram(
      aes(
        y = after_stat(
          count / sum(count) * 100
        )
      ),
      bins = 30,
      fill = "steelblue",
      color = "white",
      alpha = 0.8
    ) +
    theme_minimal(
      base_size = 13
    ) +
    theme(
      plot.title = element_text(
        face = "plain",
        size = 13,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        face = "plain",
        size = 9,
        hjust = 0.5
      ),
      axis.title.x = element_text(
        size = 10
      ),
      axis.text.x = element_text(
        size = 8
      ),
      axis.title.y = element_text(
        size = 10
      ),
      axis.text.y = element_text(
        size = 8
      )
    ) +
    labs(
      title = paste0(
        pathway_name,
        " – ",
        score_label_text
      ),
      subtitle = plot_subtitle,
      x = score_label,
      y = "Percentage (%)"
    )
  
  ggsave(
    make_output_name(
      pathway_name,
      "Univariate",
      "SNP_SumScore",
      "Hist",
      config,
      score_mode,
      ext = "pdf"
    ),
    plot = p_all,
    width = 8,
    height = 6
  )
  
  ggsave(
    make_output_name(
      pathway_name,
      "Univariate",
      "SNP_SumScore",
      "Hist",
      config,
      score_mode,
      ext = "png"
    ),
    plot = p_all,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  #------------------------------------------------------------
  # 2) Histogram by obesity
  # Percentages within the respective group
  #------------------------------------------------------------
  
  p_obesity <- NULL
  
  if (
    nrow(df_obesity) > 0 &&
    "obesity_group" %in%
    colnames(df_obesity)
  ) {
    
    p_obesity <- ggplot(
      df_obesity,
      aes(
        x = SNP_SumScore,
        fill = obesity_group
      )
    ) +
      geom_histogram(
        aes(
          y = after_stat(
            count /
              ave(
                count,
                group,
                FUN = sum
              ) *
              100
          ),
          group = obesity_group
        ),
        bins = 30,
        alpha = 0.6,
        position = "identity",
        color = "white"
      ) +
      scale_fill_manual(
        values = c(
          "steelblue",
          "firebrick"
        ),
        drop = FALSE
      ) +
      theme_minimal(
        base_size = 13
      ) +
      theme(
        plot.title = element_text(
          face = "plain",
          size = 13,
          hjust = 0.5
        ),
        plot.subtitle = element_text(
          face = "plain",
          size = 9,
          hjust = 0.5
        ),
        axis.title.x = element_text(
          size = 10
        ),
        axis.text.x = element_text(
          size = 8
        ),
        axis.title.y = element_text(
          size = 10
        ),
        axis.text.y = element_text(
          size = 8
        ),
        legend.position = "right"
      ) +
      labs(
        title = paste0(
          pathway_name,
          " – ",
          score_label_text,
          " by Obesity"
        ),
        subtitle = plot_subtitle,
        x = score_label,
        y = "Percentage within group (%)",
        fill = "Group"
      )
    
    ggsave(
      make_output_name(
        pathway_name,
        "Univariate",
        "SNP_SumScore",
        "Hist_Obesity",
        config,
        score_mode,
        ext = "pdf"
      ),
      plot = p_obesity,
      width = 8,
      height = 6
    )
    
    ggsave(
      make_output_name(
        pathway_name,
        "Univariate",
        "SNP_SumScore",
        "Hist_Obesity",
        config,
        score_mode,
        ext = "png"
      ),
      plot = p_obesity,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
  
  #------------------------------------------------------------
  # 3) Histogram by batch
  #------------------------------------------------------------
  
  p_batch <- NULL
  
  if (
    "Batch" %in%
    colnames(df_plot)
  ) {
    
    df_batch <- df_plot %>%
      filter(
        !is.na(Batch)
      )
    
    available_batches <- unique(
      as.character(df_batch$Batch)
    )
    
    if (
      nrow(df_batch) > 0 &&
      length(available_batches) >= 1
    ) {
      
      df_batch$Batch_plot <- factor(
        df_batch$Batch,
        levels = c(
          "Ambu",
          "MD"
        ),
        labels = c(
          "Ambulance",
          "Population"
        )
      )
      
      p_batch <- ggplot(
        df_batch,
        aes(
          x = SNP_SumScore,
          fill = Batch_plot
        )
      ) +
        geom_histogram(
          aes(
            y = after_stat(
              count /
                ave(
                  count,
                  group,
                  FUN = sum
                ) *
                100
            ),
            group = Batch_plot
          ),
          bins = 30,
          alpha = 0.6,
          position = "identity",
          color = "white"
        ) +
        scale_fill_manual(
          values = c(
            "steelblue",
            "firebrick"
          ),
          drop = FALSE
        ) +
        theme_minimal(
          base_size = 13
        ) +
        theme(
          plot.title = element_text(
            face = "plain",
            size = 13,
            hjust = 0.5
          ),
          plot.subtitle = element_text(
            face = "plain",
            size = 9,
            hjust = 0.5
          ),
          axis.title.x = element_text(
            size = 10
          ),
          axis.text.x = element_text(
            size = 8
          ),
          axis.title.y = element_text(
            size = 10
          ),
          axis.text.y = element_text(
            size = 8
          ),
          legend.position = "right"
        ) +
        labs(
          title = paste0(
            pathway_name,
            " – ",
            score_label_text,
            " by Batch"
          ),
          subtitle = plot_subtitle,
          x = score_label,
          y = "Percentage within batch (%)",
          fill = "Batch"
        )
      
      ggsave(
        make_output_name(
          pathway_name,
          "Univariate",
          "SNP_SumScore",
          "Hist_Batch",
          config,
          score_mode,
          ext = "pdf"
        ),
        plot = p_batch,
        width = 8,
        height = 6
      )
      
      ggsave(
        make_output_name(
          pathway_name,
          "Univariate",
          "SNP_SumScore",
          "Hist_Batch",
          config,
          score_mode,
          ext = "png"
        ),
        plot = p_batch,
        width = 8,
        height = 6,
        dpi = 300
      )
    }
  }
  
  #------------------------------------------------------------
  # 4) Boxplot by obesity
  #------------------------------------------------------------
  
  p_box <- NULL
  
  if (
    nrow(df_obesity) > 0 &&
    "obesity_group" %in%
    colnames(df_obesity)
  ) {
    
    p_box <- ggplot(
      df_obesity,
      aes(
        x = obesity_group,
        y = SNP_SumScore,
        fill = obesity_group
      )
    ) +
      geom_boxplot(
        alpha = 0.7,
        outlier.color = "red"
      ) +
      scale_fill_manual(
        values = c(
          "steelblue",
          "firebrick"
        ),
        drop = FALSE
      ) +
      theme_minimal(
        base_size = 13
      ) +
      theme(
        plot.title = element_text(
          face = "plain",
          size = 13,
          hjust = 0.5
        ),
        plot.subtitle = element_text(
          face = "plain",
          size = 9,
          hjust = 0.5
        ),
        axis.title.x = element_text(
          size = 10
        ),
        axis.text.x = element_text(
          size = 8
        ),
        axis.title.y = element_text(
          size = 10
        ),
        axis.text.y = element_text(
          size = 8
        ),
        legend.position = "none"
      ) +
      labs(
        title = paste0(
          pathway_name,
          " – ",
          score_label_text,
          " by Obesity"
        ),
        subtitle = plot_subtitle,
        x = "",
        y = score_label
      )
    
    ggsave(
      make_output_name(
        pathway_name,
        "Univariate",
        "SNP_SumScore",
        "Boxplot",
        config,
        score_mode,
        ext = "pdf"
      ),
      plot = p_box,
      width = 8,
      height = 6
    )
    
    ggsave(
      make_output_name(
        pathway_name,
        "Univariate",
        "SNP_SumScore",
        "Boxplot",
        config,
        score_mode,
        ext = "png"
      ),
      plot = p_box,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
  
  #------------------------------------------------------------
  # 5) Forest plot
  #
  # For all score modes, OR per 1 standard deviation.
  # This allows direct comparison between:
  # oriented, raw, beta_weighted and lof_burden.
  #------------------------------------------------------------
  
  p_wald <- NULL
  
  model_standardized <-
    model_result$model_standardized
  
  if (!is.null(model_standardized)) {
    
    standardized_coef <-
      summary(
        model_standardized
      )$coefficients
    
    if (
      "SNP_SumScore_SD" %in%
      rownames(standardized_coef)
    ) {
      
      est_sd <- standardized_coef[
        "SNP_SumScore_SD",
        "Estimate"
      ]
      
      se_sd <- standardized_coef[
        "SNP_SumScore_SD",
        "Std. Error"
      ]
      
      pvalue_sd <- standardized_coef[
        "SNP_SumScore_SD",
        "Pr(>|z|)"
      ]
      
      wald_df <- data.frame(
        term_friendly = paste0(
          score_label_text,
          "\n(per 1 SD)"
        ),
        OR = exp(est_sd),
        CI_lower = exp(
          est_sd - 1.96 * se_sd
        ),
        CI_upper = exp(
          est_sd + 1.96 * se_sd
        ),
        Pvalue = pvalue_sd,
        stringsAsFactors = FALSE
      )
      
      wald_subtitle <- paste0(
        "OR per 1 SD = ",
        round(
          wald_df$OR,
          3
        ),
        " (95% CI ",
        round(
          wald_df$CI_lower,
          3
        ),
        "–",
        round(
          wald_df$CI_upper,
          3
        ),
        "), p = ",
        format.pval(
          wald_df$Pvalue,
          digits = 3,
          eps = 0.001
        )
      )
      
      p_wald <- ggplot(
        wald_df,
        aes(
          x = term_friendly,
          y = OR
        )
      ) +
        geom_point(
          size = 4,
          color = "steelblue"
        ) +
        geom_errorbar(
          aes(
            ymin = CI_lower,
            ymax = CI_upper
          ),
          width = 0.15,
          color = "steelblue",
          linewidth = 1
        ) +
        geom_hline(
          yintercept = 1,
          linetype = "dashed",
          color = "gray40"
        ) +
        theme_minimal(
          base_size = 13
        ) +
        theme(
          plot.title = element_text(
            face = "plain",
            size = 13,
            hjust = 0.5
          ),
          plot.subtitle = element_text(
            face = "plain",
            size = 9,
            hjust = 0.5
          ),
          axis.title.x = element_text(
            size = 10
          ),
          axis.text.x = element_text(
            size = 8
          ),
          axis.title.y = element_text(
            size = 10
          ),
          axis.text.y = element_text(
            size = 8
          ),
          legend.position = "none"
        ) +
        labs(
          title = paste0(
            pathway_name,
            " – Standardized Score Association"
          ),
          subtitle = wald_subtitle,
          x = "",
          y = "Odds Ratio per 1 SD (95% CI)"
        )
      
      ggsave(
        make_output_name(
          pathway_name,
          "Univariate",
          "SNP_SumScore",
          "Waldplot_PerSD",
          config,
          score_mode,
          ext = "pdf"
        ),
        plot = p_wald,
        width = 7,
        height = 5
      )
      
      ggsave(
        make_output_name(
          pathway_name,
          "Univariate",
          "SNP_SumScore",
          "Waldplot_PerSD",
          config,
          score_mode,
          ext = "png"
        ),
        plot = p_wald,
        width = 7,
        height = 5,
        dpi = 300
      )
    }
  }
  
  if (is.null(p_wald)) {
    
    warning(
      "Standardized model for the forest plot not available: ",
      result_name
    )
  }
  
  #------------------------------------------------------------
  # 6) ROC plot – final test-set performance
  #------------------------------------------------------------
  
  p_roc <- NULL
  
  roc_df <- model_result$roc_data
  
  if (
    !is.null(roc_df) &&
    all(
      c(
        "FPR",
        "TPR"
      ) %in%
      colnames(roc_df)
    )
  ) {
    
    roc_df <- roc_df %>%
      filter(
        is.finite(FPR),
        is.finite(TPR)
      ) %>%
      arrange(FPR, TPR)
    
    auc_val <- unique(
      model_result$results$AUC
    )
    
    auc_val <- auc_val[
      is.finite(auc_val)
    ]
    
    if (length(auc_val) > 0) {
      auc_val <- auc_val[1]
    } else {
      auc_val <- NA_real_
    }
    
    auc_ci_lower <- unique(
      model_result$results$AUC_CI_lower
    )
    
    auc_ci_upper <- unique(
      model_result$results$AUC_CI_upper
    )
    
    auc_ci_lower <- auc_ci_lower[
      is.finite(auc_ci_lower)
    ]
    
    auc_ci_upper <- auc_ci_upper[
      is.finite(auc_ci_upper)
    ]
    
    if (
      length(auc_ci_lower) > 0 &&
      length(auc_ci_upper) > 0
    ) {
      
      roc_subtitle <- paste0(
        "Test AUC = ",
        round(
          auc_val,
          3
        ),
        " (95% CI ",
        round(
          auc_ci_lower[1],
          3
        ),
        "–",
        round(
          auc_ci_upper[1],
          3
        ),
        ")"
      )
      
    } else {
      
      roc_subtitle <- paste0(
        "Test AUC = ",
        round(
          auc_val,
          3
        )
      )
    }
    
    p_roc <- ggplot(
      roc_df,
      aes(
        x = FPR,
        y = TPR
      )
    ) +
      geom_line(
        color = "blue",
        linewidth = 1.2
      ) +
      geom_abline(
        intercept = 0,
        slope = 1,
        linetype = "dashed",
        color = "gray50"
      ) +
      coord_equal(
        xlim = c(0, 1),
        ylim = c(0, 1)
      ) +
      theme_minimal(
        base_size = 13
      ) +
      theme(
        plot.title = element_text(
          face = "plain",
          size = 13,
          hjust = 0.5
        ),
        plot.subtitle = element_text(
          face = "plain",
          size = 10,
          hjust = 0.5
        ),
        axis.title.x = element_text(
          size = 10
        ),
        axis.text.x = element_text(
          size = 8
        ),
        axis.title.y = element_text(
          size = 10
        ),
        axis.text.y = element_text(
          size = 8
        ),
        legend.position = "none"
      ) +
      labs(
        title = paste0(
          pathway_name,
          " – ",
          score_label_text,
          " ROC"
        ),
        subtitle = roc_subtitle,
        x = "1 − Specificity",
        y = "Sensitivity"
      )
    
    ggsave(
      make_output_name(
        pathway_name,
        "Univariate",
        "SNP_SumScore",
        "ROC",
        config,
        score_mode,
        ext = "pdf"
      ),
      plot = p_roc,
      width = 7,
      height = 5
    )
    
    ggsave(
      make_output_name(
        pathway_name,
        "Univariate",
        "SNP_SumScore",
        "ROC",
        config,
        score_mode,
        ext = "png"
      ),
      plot = p_roc,
      width = 7,
      height = 5,
      dpi = 300
    )
  }
  
  if (is.null(p_roc)) {
    
    warning(
      "No valid ROC data for ",
      result_name
    )
  }
  
  #------------------------------------------------------------
  # Completion
  #------------------------------------------------------------
  if (score_mode == "lof_burden") {
    
    cat(
      "LoF visualization: distribution of the biologically ",
      "ausgerichteten LoF-Allel-Burdens.\n",
      sep = ""
    )
    
    cat(
      "Waldplot: OR pro 1 SD des LoF-Burdens.\n"
    )
  }
  
  cat(
    "Plots saved for:",
    result_name,
    "\n"
  )
  
  return(
    list(
      hist = p_all,
      hist_obesity = p_obesity,
      hist_batch = p_batch,
      boxplot = p_box,
      waldplot_per_sd = p_wald,
      roc = p_roc
    )
  )
}

#===================================================================================
# LOOP OVER ALL SNP SUM SCORE PLOTS
#===================================================================================

all_snp_sumscore_plots <- list()

for (
  result_name in names(
    all_snp_sumscore_results
  )
) {
  
  df_current <- all_snp_sumscore_results[[result_name]]
  
  model_result_current <-
    all_snp_sumscore_model_results[[result_name]]
  
  plot_result <- tryCatch(
    
    run_snp_sumscore_plots(
      df = df_current,
      model_result = model_result_current,
      result_name = result_name
    ),
    
    error = function(e) {
      
      warning(
        "Error in the plots for ",
        result_name,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  all_snp_sumscore_plots[[result_name]] <- plot_result
}

#===================================================================================
# UNIVARIATE SKAT ANALYSIS – LOOP OVER ANALYSIS DATASETS
#===================================================================================

library(SKAT)
library(data.table)
library(dplyr)

run_skat_analysis <- function(
    config_name,
    config) {
  
  cat("\n====================================================\n")
  cat("SKAT analysis:", config_name, "\n")
  cat("====================================================\n")
  
  #---------------------------------------------------------------------------
  # The beta_weighted configuration belongs only to the SNP sum score.
  # For SKAT it would be identical to main_ld and would duplicate the result.
  #---------------------------------------------------------------------------
  
  if (config$type == "beta_weighted") {
    
    cat(
      "Skipped: beta_weighted is a score configuration",
      "and not a separate SKAT SNP selection.\n",
      sep = ""
    )
    
    return(NULL)
  }
  
  skat_file <- make_output_name(
    pathway_name,
    "Univariate",
    "SKAT",
    "Results",
    config
  )
  
  imputation_file <- make_output_name(
    pathway_name,
    "Univariate",
    "SKAT",
    "Imputation_Values",
    config
  )
  
  #---------------------------------------------------------------------------
  # Load data
  #---------------------------------------------------------------------------
  
  if (!file.exists(config$dataset)) {
    stop(
      "SKAT dataset not found:\n",
      config$dataset
    )
  }
  
  data_pathway <- readRDS(
    config$dataset
  )
  
  if (
    is.null(data_pathway$genotypes) ||
    is.null(data_pathway$pheno)
  ) {
    stop(
      "The dataset must contain the elements 'genotypes' and 'pheno'."
    )
  }
  
  geno <- data_pathway$genotypes
  pheno <- data_pathway$pheno
  
  if (is.null(rownames(geno))) {
    stop(
      "The genotype matrix has no row names or IDs."
    )
  }
  
  if (!"id" %in% colnames(pheno)) {
    stop(
      "The phenotype data contain no column 'id'."
    )
  }
  
  if (!"obesity" %in% colnames(pheno)) {
    stop(
      "The phenotype data contain no column 'obesity'."
    )
  }
  
  pheno$id <- as.character(
    pheno$id
  )
  
  #---------------------------------------------------------------------------
  # Match IDs
  #---------------------------------------------------------------------------
  
  common_ids <- intersect(
    rownames(geno),
    pheno$id
  )
  
  if (length(common_ids) == 0) {
    
    warning(
      "No shared genotype/phenotype IDs for ",
      config_name
    )
    
    return(NULL)
  }
  
  geno <- geno[
    common_ids,
    ,
    drop = FALSE
  ]
  
  pheno <- pheno[
    match(
      common_ids,
      pheno$id
    ),
    ,
    drop = FALSE
  ]
  
  rownames(pheno) <- pheno$id
  
  stopifnot(
    all(
      rownames(geno) ==
        pheno$id
    )
  )
  
  #---------------------------------------------------------------------------
  # SNP selection by analysis configuration
  #---------------------------------------------------------------------------
  
  valid_snps <- colnames(
    geno
  )
  
  # For the LoF configuration, the allele orientation
  # is applied later to the numeric genotype matrix.
  lof_orientation_table <- NULL
  
  if (config$type == "functional_vep") {
    
    if (
      is.null(config$vep_file) ||
      !file.exists(config$vep_file)
    ) {
      
      warning(
        "VEP file for ",
        config_name,
        " not found."
      )
      
      return(NULL)
    }
    
    functional_snps <- load_functional_vep_snps(
      config$vep_file
    )
    
    valid_snps <- intersect(
      valid_snps,
      functional_snps
    )
  }
  
  if (config$type == "lof_vep") {
    
    if (
      is.null(config$vep_file) ||
      !file.exists(config$vep_file)
    ) {
      warning(
        "VEP file for the LoF SKAT analysis ",
        config_name,
        " not found."
      )
      return(NULL)
    }
    
    lof_orientation_table <- create_lof_orientation_table(
      vep_file = config$vep_file,
      bim_table = lof_bim,
      matrix_snps = colnames(geno)
    )
    
    # Could a shared SNP set be created at all between
    # VEP, BIM and genotype matrix?
    if (nrow(lof_orientation_table) == 0) {
      warning(
        "No LoF SNPs could be matched between VEP, BIM and ",
        "genotype matrix."
      )
      return(NULL)
    }
    
    no_match_snps <- lof_orientation_table %>%
      filter(
        Allele_Match == "No_allele_match"
      ) %>%
      pull(SNP) %>%
      unique()
    
    if (length(no_match_snps) > 0) {
      warning(
        "The following VEP LoF alleles match no matrix allele ",
        "and are excluded for SKAT: ",
        paste(no_match_snps, collapse = ", ")
      )
    }
    
    # Keep only unambiguously orientable LoF SNPs
    lof_orientation_table <- lof_orientation_table %>%
      filter(
        Allele_Match != "No_allele_match",
        !is.na(Invert_For_LoF_Score)
      ) %>%
      distinct(
        SNP,
        .keep_all = TRUE
      )
    
    if (nrow(lof_orientation_table) == 0) {
      warning(
        "After allele matching, no unambiguously ",
        "orientable LoF SNPs remain for SKAT."
      )
      return(NULL)
    }
    
    valid_snps <- intersect(
      valid_snps,
      lof_orientation_table$SNP
    )
  }
  
  if (config$type == "database_kegg_only") {
    
    if (
      is.null(config$snp_gene_map) ||
      !file.exists(config$snp_gene_map)
    ) {
      
      warning(
        "SNP-gene mapping for ",
        config_name,
        " not found."
      )
      
      return(NULL)
    }
    
    snp_gene_map <- read.csv2(
      config$snp_gene_map,
      stringsAsFactors = FALSE
    )
    
    required_map_columns <- c(
      "SNP",
      "GENE"
    )
    
    missing_map_columns <- setdiff(
      required_map_columns,
      colnames(snp_gene_map)
    )
    
    if (length(missing_map_columns) > 0) {
      
      warning(
        "Missing in SNP-gene mapping: ",
        paste(
          missing_map_columns,
          collapse = ", "
        )
      )
      
      return(NULL)
    }
    
    kegg_only_snps <- snp_gene_map %>%
      filter(
        GENE %in% config$kegg_genes
      ) %>%
      pull(SNP) %>%
      as.character() %>%
      unique()
    
    valid_snps <- intersect(
      valid_snps,
      kegg_only_snps
    )
  }
  
  valid_snps <- unique(
    valid_snps
  )
  
  n_snps_selected <- length(
    valid_snps
  )
  
  if (n_snps_selected < 2) {
    
    warning(
      "Too few SNPs for SKAT in ",
      config_name,
      ": ",
      n_snps_selected
    )
    
    return(NULL)
  }
  
  geno <- geno[
    ,
    valid_snps,
    drop = FALSE
  ]
  
  #---------------------------------------------------------------------------
  # Prepare outcome
  #---------------------------------------------------------------------------
  
  obesity_numeric <- suppressWarnings(
    as.numeric(
      as.character(
        pheno$obesity
      )
    )
  )
  
  y <- ifelse(
    is.na(obesity_numeric),
    NA_real_,
    ifelse(
      obesity_numeric == 1,
      1,
      0
    )
  )
  
  outcome_keep <- !is.na(
    y
  )
  
  geno <- geno[
    outcome_keep,
    ,
    drop = FALSE
  ]
  
  pheno <- pheno[
    outcome_keep,
    ,
    drop = FALSE
  ]
  
  y <- y[
    outcome_keep
  ]
  
  if (length(y) < 20) {
    
    warning(
      "Too few complete outcome cases for ",
      config_name
    )
    
    return(NULL)
  }
  
  if (length(unique(y)) < 2) {
    
    warning(
      "Only one outcome class for ",
      config_name
    )
    
    return(NULL)
  }
  
  #---------------------------------------------------------------------------
  # Create robust numeric genotype matrix
  #---------------------------------------------------------------------------
  
  geno_df <- as.data.frame(
    geno,
    stringsAsFactors = FALSE
  )
  
  geno_df[] <- lapply(
    geno_df,
    function(x) {
      
      if (is.list(x)) {
        x <- unlist(x)
      }
      
      suppressWarnings(
        as.numeric(x)
      )
    }
  )
  
  geno_mat <- as.matrix(
    geno_df
  )
  
  storage.mode(
    geno_mat
  ) <- "numeric"
  
  rownames(geno_mat) <- rownames(
    geno
  )
  
  colnames(geno_mat) <- colnames(
    geno
  )
  
  #---------------------------------------------------------------------------
  # LoF allele orientation for SKAT and SKAT-O
  #
  # The original data remain unchanged. Within this local
  # SKAT matrix, however, each LoF SNP is coded such that:
  # 0 = no LoF allele
  # 1 = one LoF allele
  # 2 = two LoF alleles
  #---------------------------------------------------------------------------
  
  lof_orientation_used <- NULL
  lof_orientation_file <- NULL
  invert_snps <- character(0)
  
  if (config$type == "lof_vep") {
    
    if (is.null(lof_orientation_table)) {
      stop(
        "For the LoF SKAT analysis no ",
        "allele orientation table was created."
      )
    }
    
    lof_orientation_used <- lof_orientation_table[
      match(
        colnames(geno_mat),
        lof_orientation_table$SNP
      ),
      ,
      drop = FALSE
    ]
    
    if (
      any(is.na(lof_orientation_used$SNP)) ||
      !all(
        lof_orientation_used$SNP ==
        colnames(geno_mat)
      )
    ) {
      stop(
        "Not for all LoF SNPs of the SKAT matrix could ",
        "an unambiguous allele orientation be assigned."
      )
    }
    
    invert_snps <- lof_orientation_used$SNP[
      lof_orientation_used$Invert_For_LoF_Score %in% TRUE
    ]
    
    if (length(invert_snps) > 0) {
      
      invert_idx <- match(
        invert_snps,
        colnames(geno_mat)
      )
      
      geno_mat[
        ,
        invert_idx
      ] <- 2 - geno_mat[
        ,
        invert_idx,
        drop = FALSE
      ]
    }
    
    lof_orientation_used$Applied_Transformation <- ifelse(
      lof_orientation_used$Invert_For_LoF_Score %in% TRUE,
      "2_minus_original_dosage",
      "original_dosage"
    )
    
    lof_orientation_file <- make_output_name(
      pathway_name,
      "Univariate",
      "SKAT",
      "LoF_Allele_Orientation",
      config
    )
    
    write.csv2(
      lof_orientation_used,
      lof_orientation_file,
      row.names = FALSE
    )
  }
  
  
  #---------------------------------------------------------------------------
  # Plausibility check of genotype values
  #---------------------------------------------------------------------------
  
  observed_values <- geno_mat[
    !is.na(geno_mat)
  ]
  
  if (length(observed_values) == 0) {
    
    warning(
      "No observed genotype values for ",
      config_name
    )
    
    return(NULL)
  }
  
  if (
    min(observed_values) < 0 ||
    max(observed_values) > 2
  ) {
    
    stop(
      "Genotype values outside the expected range 0–2 in ",
      config_name
    )
  }
  
  #---------------------------------------------------------------------------
  # Exclude individuals with more than 10 % missing SNPs
  #---------------------------------------------------------------------------
  
  individual_missingness <- rowMeans(
    is.na(geno_mat)
  )
  
  individual_keep <- individual_missingness <= 0.10
  
  n_individuals_before_missingness <- nrow(
    geno_mat
  )
  
  geno_mat <- geno_mat[
    individual_keep,
    ,
    drop = FALSE
  ]
  
  y <- y[
    individual_keep
  ]
  
  n_individuals_excluded_missingness <-
    sum(
      !individual_keep
    )
  
  if (nrow(geno_mat) < 20) {
    
    warning(
      "Too few individuals after the missingness filter for ",
      config_name
    )
    
    return(NULL)
  }
  
  if (length(unique(y)) < 2) {
    
    warning(
      "Only one outcome class after the missingness filter for ",
      config_name
    )
    
    return(NULL)
  }
  
  #---------------------------------------------------------------------------
  # Remove SNPs without observed values or without variation
  #---------------------------------------------------------------------------
  
  n_observed_per_snp <- colSums(
    !is.na(geno_mat)
  )
  
  snp_variance <- apply(
    geno_mat,
    2,
    var,
    na.rm = TRUE
  )
  
  usable_snp_mask <-
    n_observed_per_snp > 0 &
    is.finite(snp_variance) &
    snp_variance > 0
  
  removed_snps <- colnames(
    geno_mat
  )[
    !usable_snp_mask
  ]
  
  geno_mat <- geno_mat[
    ,
    usable_snp_mask,
    drop = FALSE
  ]
  
  if (
    config$type == "lof_vep" &&
    !is.null(lof_orientation_used)
  ) {
    
    lof_orientation_used <- lof_orientation_used[
      match(
        colnames(geno_mat),
        lof_orientation_used$SNP
      ),
      ,
      drop = FALSE
    ]
    
    if (
      anyNA(lof_orientation_used$SNP) ||
      !all(
        lof_orientation_used$SNP ==
        colnames(geno_mat)
      )
    ) {
      stop(
        "The LoF orientation table could not be updated correctly after the ",
        "SNP variation filter."
      )
    }
    
    write.csv2(
      lof_orientation_used,
      lof_orientation_file,
      row.names = FALSE
    )
  }
  
  if (ncol(geno_mat) < 2) {
    
    warning(
      "After removing unusable or ",
      "monomorphic SNPs, too few SNPs for ",
      config_name
    )
    
    return(NULL)
  }
  
  #---------------------------------------------------------------------------
  # Impute missing genotype values with the SNP mean
  #
  # SKAT here is an inferential analysis of the complete
  # analysis dataset and not a train/test prediction analysis.
  # Hence the SNP mean is computed within the respective SKAT
  # analysis dataset.
  #---------------------------------------------------------------------------
  
  snp_means <- colMeans(
    geno_mat,
    na.rm = TRUE
  )
  
  n_imputed_by_snp <- colSums(
    is.na(geno_mat)
  )
  
  for (j in seq_len(ncol(geno_mat))) {
    
    missing_j <- is.na(
      geno_mat[, j]
    )
    
    if (any(missing_j)) {
      
      geno_mat[
        missing_j,
        j
      ] <- snp_means[j]
    }
  }
  
  if (anyNA(geno_mat)) {
    
    stop(
      "After the SKAT imputation, missing values remain ",
      "genotype values remain."
    )
  }
  
  imputation_table <- data.frame(
    SNP = colnames(geno_mat),
    Mean_Genotype = as.numeric(
      snp_means[
        colnames(geno_mat)
      ]
    ),
    N_Imputed = as.integer(
      n_imputed_by_snp[
        colnames(geno_mat)
      ]
    ),
    Imputation_Method = "Analysis_Set_Mean_Genotype",
    stringsAsFactors = FALSE
  )
  
  write.csv2(
    imputation_table,
    imputation_file,
    row.names = FALSE
  )
  
  #---------------------------------------------------------------------------
  # MAF for documentation
  #---------------------------------------------------------------------------
  
  allele_frequency <- colMeans(
    geno_mat
  ) / 2
  
  maf <- pmin(
    allele_frequency,
    1 - allele_frequency
  )
  
  #---------------------------------------------------------------------------
  # Null model
  #---------------------------------------------------------------------------
  
  null_model <- tryCatch(
    
    SKAT_Null_Model(
      y ~ 1,
      out_type = "D"
    ),
    
    error = function(e) {
      
      warning(
        "SKAT null model failed for ",
        config_name,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  if (is.null(null_model)) {
    return(NULL)
  }
  
  #---------------------------------------------------------------------------
  # SKAT and SKAT-O
  #
  # Default weighting of the SKAT package:
  # Beta-MAF weights with Beta(1,25).
  # No use of the univariate log(OR) weights.
  #---------------------------------------------------------------------------
  
  skat_result <- tryCatch(
    
    SKAT(
      geno_mat,
      null_model,
      method = "SKAT"
    ),
    
    error = function(e) {
      
      warning(
        "SKAT failed for ",
        config_name,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  skato_result <- tryCatch(
    
    SKAT(
      geno_mat,
      null_model,
      method = "SKATO"
    ),
    
    error = function(e) {
      
      warning(
        "SKAT-O failed for ",
        config_name,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  skat_p <- if (
    !is.null(skat_result)
  ) {
    skat_result$p.value
  } else {
    NA_real_
  }
  
  skato_p <- if (
    !is.null(skato_result)
  ) {
    skato_result$p.value
  } else {
    NA_real_
  }
  
  #---------------------------------------------------------------------------
  # Results table
  #---------------------------------------------------------------------------
  
  skat_results <- data.frame(
    Pathway = pathway_name,
    Analysis_Config = config_name,
    Analysis_Label = config$label,
    Analysis_Type = config$type,
    SNP_Set_Description = ifelse(
      config$type == "lof_vep",
      "VEP-annotated LoF variants oriented to LoF allele",
      config$label
    ),
    
    Allele_Orientation = ifelse(
      config$type == "lof_vep",
      "Dosage oriented to VEP LoF allele",
      "Original genotype dosage"
    ),
    
    Num_LoF_SNPs_Inverted = ifelse(
      config$type == "lof_vep",
      length(
        intersect(
          invert_snps,
          colnames(geno_mat)
        )
      ),
      NA_integer_
    ),
    
    Num_Individuals_Before_Missingness =
      n_individuals_before_missingness,
    
    Num_Individuals_Excluded_Missingness =
      n_individuals_excluded_missingness,
    
    Num_Individuals = nrow(
      geno_mat
    ),
    
    Num_Cases = sum(
      y == 1
    ),
    
    Num_Controls = sum(
      y == 0
    ),
    
    Num_SNPs_Selected = n_snps_selected,
    
    Num_SNPs_Removed_NoVariation =
      length(
        removed_snps
      ),
    
    Num_SNPs = ncol(
      geno_mat
    ),
    
    Num_Imputed_Genotypes = sum(
      n_imputed_by_snp
    ),
    
    Mean_MAF = mean(
      maf,
      na.rm = TRUE
    ),
    
    Median_MAF = median(
      maf,
      na.rm = TRUE
    ),
    
    Min_MAF = min(
      maf,
      na.rm = TRUE
    ),
    
    Max_MAF = max(
      maf,
      na.rm = TRUE
    ),
    
    Imputation_Method =
      "Analysis_Set_Mean_Genotype",
    
    SKAT_Weighting =
      "Default Beta(1,25) MAF weights",
    
    SKAT_pvalue = skat_p,
    SKATO_pvalue = skato_p,
    
    Significant_SKAT = ifelse(
      is.finite(skat_p),
      skat_p < 0.05,
      NA
    ),
    
    Significant_SKATO = ifelse(
      is.finite(skato_p),
      skato_p < 0.05,
      NA
    ),
    
    Significant_Any = ifelse(
      is.finite(skat_p) |
        is.finite(skato_p),
      
      (
        is.finite(skat_p) &&
          skat_p < 0.05
      ) |
        (
          is.finite(skato_p) &&
            skato_p < 0.05
        ),
      
      NA
    ),
    
    stringsAsFactors = FALSE
  )
  
  #---------------------------------------------------------------------------
  # Save results
  #---------------------------------------------------------------------------
  
  write.csv2(
    skat_results,
    skat_file,
    row.names = FALSE
  )
  
  cat(
    "Personen:",
    nrow(geno_mat),
    "\n"
  )
  
  cat(
    "Cases / Controls:",
    sum(y == 1),
    "/",
    sum(y == 0),
    "\n"
  )
  
  cat(
    "Verwendete SNPs:",
    ncol(geno_mat),
    "\n"
  )
  
  cat(
    "Imputierte Genotypwerte:",
    sum(n_imputed_by_snp),
    "\n"
  )
  
  cat(
    "SKAT p-value:",
    signif(skat_p, 4),
    "\n"
  )
  
  cat(
    "SKAT-O p-value:",
    signif(skato_p, 4),
    "\n"
  )
  
  if (config$type == "lof_vep") {
    
    cat(
      "LoF SNPs with inverted dosage:",
      length(
        intersect(
          invert_snps,
          colnames(geno_mat)
        )
      ),
      "\n"
    )
    
    cat(
      "LoF orientation file:",
      lof_orientation_file,
      "\n"
    )
  }
  
  cat(
    "Gespeichert:",
    skat_file,
    "\n"
  )
  
  return(
    skat_results
  )
}

#===================================================================================
# LOOP OVER ALL RELEVANT SKAT ANALYSES
#===================================================================================

all_skat_results <- list()

for (
  config_name in names(
    analysis_configs
  )
) {
  
  config <- analysis_configs[[config_name]]
  
  # beta_weighted is exclusively a score sensitivity analysis
  if (config$type == "beta_weighted") {
    
    cat(
      "\nSKAT skipped for ",
      config_name,
      ": no separate SNP selection.\n",
      sep = ""
    )
    
    next
  }
  
  skat_result <- tryCatch(
    
    run_skat_analysis(
      config_name = config_name,
      config = config
    ),
    
    error = function(e) {
      
      warning(
        "Error in the SKAT analysis ",
        config_name,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  all_skat_results[[config_name]] <- skat_result
}

#===================================================================================
# SKAT OVERALL SUMMARY
#===================================================================================

skat_summary <- bind_rows(
  all_skat_results
)

if (nrow(skat_summary) > 0) {
  
  skat_summary_file <- output_file(
    "Pathway_",
    pathway_name,
    "_Univariate_SKAT_Results_All_Configurations.csv"
  )
  
  write.csv2(
    skat_summary,
    skat_summary_file,
    row.names = FALSE
  )
  
  cat(
    "\nSKAT overall summary saved:\n",
    skat_summary_file,
    "\n"
  )
  
} else {
  
  warning(
    "No SKAT results could be merged."
  )
}

#================================================================================
# MULTIVARIABLE ANALYSIS – SNP SUM SCORE LOOP
#================================================================================
library(dplyr)
library(broom)
library(car)
library(vcd)
library(caret)
library(pROC)
library(margins)
library(pscl)
library(ggeffects)
library(ggplot2)

# KOVARIATE PRECHECK 
run_covariate_precheck <- function(df, cov_categorical, cov_numeric, config, score_mode) {
  
  cat("Starting covariate pre-analyses...\n")
  
  #============================================================================
  # Consistency check of the analysis configuration
  #============================================================================
  
  if (
    score_mode == "lof_burden" &&
    config$type != "lof_vep"
  ) {
    stop(
      "score_mode = 'lof_burden' may only be used in the covariate pre-check ",
      "with config$type = 'lof_vep'."
    )
  }
  
  if (
    config$type == "lof_vep" &&
    score_mode != "lof_burden"
  ) {
    stop(
      "A LoF VEP configuration may only be used in the covariate pre-check ",
      "with score_mode = 'lof_burden'."
    )
  }
  
  #============================================================================
  # 1) Remove outliers of numeric covariates directly
  #============================================================================
  
  outlier_log <- list()
  outlier_rows_all <- c()
  
  n_before_outlier_removal <- nrow(df)
  
  for (v in cov_numeric) {
    
    x <- df[[v]]
    
    Q1 <- quantile(x, 0.25, na.rm = TRUE)
    Q3 <- quantile(x, 0.75, na.rm = TRUE)
    IQRv <- IQR(x, na.rm = TRUE)
    
    lower <- Q1 - 1.5 * IQRv
    upper <- Q3 + 1.5 * IQRv
    
    outlier_idx <- which(x < lower | x > upper)
    
    outlier_log[[v]] <- data.frame(
      Variable = v,
      n_outliers = length(outlier_idx),
      Lower_Cutoff = lower,
      Upper_Cutoff = upper
    )
    
    outlier_rows_all <- unique(c(outlier_rows_all, outlier_idx))
  }
  
  if (length(outlier_rows_all) > 0) {
    df <- df[-outlier_rows_all, , drop = FALSE]
  }
  
  n_after_outlier_removal <- nrow(df)
  
  outlier_df <- bind_rows(outlier_log) %>%
    mutate(
      n_before_outlier = n_before_outlier_removal,
      n_after_outlier = n_after_outlier_removal,
      n_removed = n_before_outlier_removal - n_after_outlier_removal
    )
  
  write.csv2(
    outlier_df,
    make_output_name(
      pathway_name,
      "Multivariate",
      "Covariate_Precheck",
      "Outliers",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  cat("\n========================================\n")
  cat("Outlier removal numeric covariates\n")
  cat("========================================\n")
  if (nrow(outlier_df) > 0) {
    
    for (i in seq_len(nrow(outlier_df))) {
      
      cat(
        outlier_df$Variable[i],
        ": ",
        outlier_df$n_outliers[i],
        " Outlier\n",
        sep = ""
      )
    }
  }
  cat("Before outlier removal: ", n_before_outlier_removal, "\n")
  cat("After outlier removal:", n_after_outlier_removal, "\n")
  cat("Removed:               ", n_before_outlier_removal - n_after_outlier_removal, "\n")
  
  #============================================================================
  # 2) Dataset for pre-analyses
  #============================================================================
  
  df_clean <- df %>%
    select(
      id,
      obesity,
      SNP_SumScore,
      all_of(cov_categorical),
      all_of(cov_numeric)
    )
  
  if ("Batch" %in% colnames(df)) {
    
    df_clean$Batch <- factor(
      df$Batch
    )
    
    cov_categorical_all <- unique(
      c(
        cov_categorical,
        "Batch"
      )
    )
    
  } else {
    
    cov_categorical_all <- cov_categorical
  }
  
  #============================================================================
  # 3) Numeric correlation
  #============================================================================
  
  cov_numeric_present <- cov_numeric[
    cov_numeric %in% colnames(df_clean)
  ]
  
  if (length(cov_numeric_present) < 2) {
    
    warning(
      "Fewer than two numeric covariates available. ",
      "Correlation analysis is skipped."
    )
    
    cor_mat <- matrix(
      numeric(0),
      nrow = 0,
      ncol = 0
    )
    
    cor_summary <- data.frame(
      Variable = character(),
      MaxCorrelation = numeric(),
      MaxCor_with = character(),
      HighCorrelation = character(),
      stringsAsFactors = FALSE
    )
    
  } else {
    
    cor_mat <- cor(
      df_clean[
        ,
        cov_numeric_present,
        drop = FALSE
      ],
      use = "pairwise.complete.obs",
      method = "pearson"
    )
    
    write.csv2(
      cor_mat,
      make_output_name(
        pathway_name,
        "Multivariate",
        "Covariate_Precheck",
        "Correlation_Numeric",
        config,
        score_mode
      )
    )
    
    cor_long <- as.data.frame(as.table(cor_mat)) %>%
      filter(Var1 != Var2) %>%
      mutate(abs_cor = abs(Freq))
    
    cor_summary <- cor_long %>%
      group_by(Var1) %>%
      slice_max(abs_cor, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(
        Variable = as.character(Var1),
        MaxCorrelation = Freq,
        MaxCor_with = as.character(Var2),
        HighCorrelation = case_when(
          abs(MaxCorrelation) > 0.90 ~ "CRITICAL",
          abs(MaxCorrelation) > 0.80 ~ "PROBLEMATIC",
          abs(MaxCorrelation) >= 0.50 ~ "OBSERVE",
          TRUE ~ "OK"
        )
      )
  }
  
  #============================================================================
  # 4) Cramer's V
  #============================================================================
  
  cramer_results <- data.frame(
    Var1 = character(),
    Var2 = character(),
    CramersV = numeric(),
    stringsAsFactors = FALSE
  )
  
  if (length(cov_categorical_all) >= 2) {
    for (i in 1:(length(cov_categorical_all) - 1)) {
      for (j in (i + 1):length(cov_categorical_all)) {
        
        tbl <- table(
          df_clean[[cov_categorical_all[i]]],
          df_clean[[cov_categorical_all[j]]]
        )
        
        cv <- tryCatch(
          assocstats(tbl)$cramer,
          error = function(e) NA
        )
        
        cramer_results <- rbind(
          cramer_results,
          data.frame(
            Var1 = cov_categorical_all[i],
            Var2 = cov_categorical_all[j],
            CramersV = cv
          )
        )
      }
    }
  }
  
  write.csv2(
    cramer_results,
    make_output_name(
      pathway_name,
      "Multivariate",
      "Covariate_Precheck",
      "CramersV",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  if (nrow(cramer_results) > 0) {
    
    cramer_summary <- bind_rows(
      
      cramer_results %>%
        transmute(
          Variable = Var1,
          MaxCramersV_with = Var2,
          CramersV = CramersV
        ),
      
      cramer_results %>%
        transmute(
          Variable = Var2,
          MaxCramersV_with = Var1,
          CramersV = CramersV
        )
      
    ) %>%
      filter(
        is.finite(CramersV)
      ) %>%
      group_by(
        Variable
      ) %>%
      slice_max(
        CramersV,
        n = 1,
        with_ties = FALSE
      ) %>%
      ungroup() %>%
      transmute(
        Variable,
        MaxCramersV = CramersV,
        MaxCramersV_with,
        HighCramerV = case_when(
          MaxCramersV > 0.80 ~ "PROBLEMATIC",
          MaxCramersV > 0.50 ~ "CRITICAL_CHECK",
          MaxCramersV >= 0.30 ~ "OBSERVE",
          TRUE ~ "OK"
        )
      )
    
  } else {
    
    cramer_summary <- data.frame(
      Variable = character(),
      MaxCramersV = numeric(),
      MaxCramersV_with = character(),
      HighCramerV = character(),
      stringsAsFactors = FALSE
    )
  }
  
  #============================================================================
  # 5) VIF
  #============================================================================
  
  vif_covariates <- c(cov_categorical, cov_numeric)
  
  if ("Batch" %in% colnames(df_clean)) {
    vif_covariates <- unique(
      c(
        vif_covariates,
        "Batch"
      )
    )
  }
  
  formula_vif <- as.formula(
    paste("obesity ~", paste(vif_covariates, collapse = " + "))
  )
  
  model_vif <- glm(
    formula_vif,
    data = df_clean,
    family = binomial
  )
  
  vif_values <- car::vif(model_vif)
  
  if (is.matrix(vif_values)) {
    vif_df <- data.frame(
      Variable = rownames(vif_values),
      VIF = vif_values[, "GVIF"]^(1 / (2 * vif_values[, "Df"]))
    )
  } else {
    vif_df <- data.frame(
      Variable = names(vif_values),
      VIF = as.numeric(vif_values)
    )
  }
  
  vif_df <- vif_df %>%
    mutate(
      VIF_flag = case_when(
        VIF > 10 ~ "CRITICAL",
        VIF >= 5 ~ "PROBLEMATIC",
        VIF >= 3 ~ "MODERATE",
        TRUE ~ "OK"
      )
    )
  
  write.csv2(
    vif_df,
    make_output_name(
      pathway_name,
      "Multivariate",
      "Covariate_Precheck",
      "VIF",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  #============================================================================
  # 6) Decision Table
  #============================================================================
  
  decision_table <- data.frame(
    Variable = unique(vif_covariates)
  ) %>%
    left_join(vif_df, by = "Variable") %>%
    left_join(cor_summary, by = "Variable") %>%
    left_join(cramer_summary, by = "Variable") %>%
    mutate(
      Recommendation = case_when(
        VIF_flag %in% c("CRITICAL", "PROBLEMATIC") ~ "REVIEW",
        HighCorrelation %in% c("CRITICAL", "PROBLEMATIC") ~ "REVIEW",
        HighCramerV %in% c("CRITICAL_CHECK", "PROBLEMATIC") ~ "REVIEW",
        TRUE ~ "KEEP"
      ),
      Manual_Decision = NA_character_,
      Comment = NA_character_
    )
  
  print(decision_table)
  
  write.csv2(
    decision_table,
    make_output_name(
      pathway_name,
      "Multivariate",
      "Covariate_Precheck",
      "DecisionTable",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  cat("Covariate pre-analyses completed.\n")
  
  return(
    list(
      df = df,
      outlier_df = outlier_df,
      cor_summary = cor_summary,
      cramer_results = cramer_results,
      cramer_summary = cramer_summary,
      vif_df = vif_df,
      decision_table = decision_table,
      cov_categorical_all = cov_categorical_all
    )
  )
}

# KOVARIATEN-PRECHECK LOOP – NUR DECISION TABLES
all_covariate_prechecks <- list()

for (result_name in names(all_snp_sumscore_results)) {
  
  df_score <- all_snp_sumscore_results[[result_name]]
  
  if (is.null(df_score)) next
  
  cat("\n====================================================\n")
  cat("Covariate pre-check:", result_name, "\n")
  cat("====================================================\n")
  
  df <- df_score
  
  # Recoding as later in the analysis
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
  
  #===========================================================================
  # Batch
  #===========================================================================
  
  batch_col <- NULL
  
  # Prefers the spelling "Batch" used in the main dataset.
  # The lower-case variant is recognised for compatibility only.
  if ("Batch" %in% colnames(df)) {
    
    batch_col <- "Batch"
    
  } else if ("batch" %in% colnames(df)) {
    
    batch_col <- "batch"
  }
  
  if (!is.null(batch_col)) {
    
    df$Batch <- factor(
      as.character(df[[batch_col]]),
      levels = c(
        "Ambu",
        "MD"
      )
    )
    
  } else {
    
    df$Batch <- NULL
  }
  
  #===========================================================================
  # Determine score and analysis configuration unambiguously
  #===========================================================================
  
  score_mode <- unique(
    as.character(df$Score_Mode)
  )
  
  score_mode <- score_mode[
    !is.na(score_mode)
  ]
  
  analysis_config <- unique(
    as.character(df$Analysis_Config)
  )
  
  analysis_config <- analysis_config[
    !is.na(analysis_config)
  ]
  
  if (
    length(score_mode) != 1 ||
    length(analysis_config) != 1
  ) {
    
    warning(
      "Score_Mode or Analysis_Config is not unique for ",
      result_name,
      ""
    )
    
    next
  }
  
  score_mode <- score_mode[1]
  analysis_config <- analysis_config[1]
  
  if (!analysis_config %in% names(analysis_configs)) {
    warning(
      "Analysis configuration not found: ",
      analysis_config,
      " for ",
      result_name
    )
    next
  }
  
  config <- analysis_configs[[analysis_config]]
  
  #===========================================================================
  # Consistency check for the LoF configuration
  #===========================================================================
  
  if (
    score_mode == "lof_burden" &&
    config$type != "lof_vep"
  ) {
    warning(
      "Invalid combination in ",
      result_name,
      ": score_mode = 'lof_burden' is only permitted with ",
      "config$type = 'lof_vep'."
    )
    next
  }
  
  if (
    config$type == "lof_vep" &&
    score_mode != "lof_burden"
  ) {
    warning(
      "Invalid combination in ",
      result_name,
      ": a LoF VEP configuration may only be used with ",
      "score_mode = 'lof_burden' verwendet werden."
    )
    next
  }
  
  if (score_mode == "lof_burden") {
    cat(
      "LoF burden analysis: the already biologically oriented ",
      "LoF allele burden is used in the covariate pre-check.\n",
      sep = ""
    )
  }
  
  #===========================================================================
  # Check score for sufficient variation
  #===========================================================================
  
  score_values <- suppressWarnings(
    as.numeric(
      df$SNP_SumScore
    )
  )
  
  score_values <- score_values[
    is.finite(score_values)
  ]
  
  if (
    length(score_values) == 0 ||
    length(unique(score_values)) < 2
  ) {
    warning(
      "The SNP sum score has insufficient variation for ",
      result_name,
      "."
    )
    next
  }
  
  if (score_mode == "lof_burden") {
    
    n_zero_burden <- sum(
      score_values == 0
    )
    
    cat(
      "LoF burden distribution: ",
      n_zero_burden,
      " of ",
      length(score_values),
      " individuals have a burden of 0.\n",
      sep = ""
    )
  }
  
  #===========================================================================
  # Load genome-wide PCs (population stratification)
  #===========================================================================
  # NO LONGER computed pathway-locally from the pathway SNPs, but
  # generated once in the gene script from the complete GSA array (QC + LD pruning +
  # KING-cutoff sensitivity check + PCA) and loaded here unchanged
  # (see genomewide_pcs_file above).
  
  df$id <- as.character(df$id)
  
  global_gPCs <- readRDS(genomewide_pcs_file)
  global_gPCs$id <- as.character(global_gPCs$id)
  
  n_before_gpc_join <- nrow(df)
  
  df <- merge(
    df,
    global_gPCs,
    by = "id",
    all.x = FALSE
  )
  
  if (nrow(df) < n_before_gpc_join) {
    warning(
      n_before_gpc_join - nrow(df),
      " Personen in ", result_name,
      " have no genome-wide PCs and were removed."
    )
  }
  
  if (nrow(df) == 0) {
    warning(
      "No individuals with genome-wide PCs for ",
      result_name
    )
    next
  }
  
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
  
  precheck <- run_covariate_precheck(
    df = df,
    cov_categorical = cov_categorical,
    cov_numeric = cov_numeric,
    config = config,
    score_mode = score_mode
  )
  
  all_covariate_prechecks[[result_name]] <- precheck
  
  cat("\nAUTOMATIC COVARIATE ASSESSMENT FOR:", result_name, "\n")
  print(precheck$decision_table)
}

saveRDS(
  all_covariate_prechecks,
  output_file(
    "Pathway_",
    pathway_name,
    "_Multivariate_Covariate_Precheck_All_Configurations.rds"
  )
)

# MANUAL COVARIATE DECISION AFTER REVIEW OF THE DECISION TABLES
manual_covariate_decisions <- list(
  default = list(
    keep = c(
      "gender_cat",
      "smoking_cat",
      "school_cert_cat",
      "graduation_cat",
      "sample_age",
      "sports_activity_min_week",
      "everyday_activity_min_week",
      "energy_kcal",
      "gPC1",
      "gPC2",
      "gPC3"
    ),
    remove = c(
      "Batch"
    ),
    comment = paste0(
      "Batch nur als Sensitivitätsanalyse; ",
      "school_cert_cat und graduation_cat trotz ",
      "Cramer's V REVIEW inhaltlich beibehalten."
    )
  )
)

#=================================================================================
# ANALYSIS FLOWCHARTS – PATHWAY
#=================================================================================

library(data.table)
library(dplyr)
library(ggplot2)
library(grid)

plot_analysis_flow_summary <- function(
    summary_df,
    pathway_name,
    outfile_prefix,
    score_mode,
    analysis_config,
    show_plot = TRUE) {
  
  #---------------------------------------------------------------------------
  # Check inputs
  #---------------------------------------------------------------------------
  
  df_plot <- as.data.table(
    summary_df
  )
  
  required_columns <- c(
    "Step",
    "SNPs",
    "Individuals"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(df_plot)
  )
  
  if (length(missing_columns) > 0) {
    
    stop(
      "Analysis summary requires the columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  if (nrow(df_plot) == 0) {
    stop(
      "The analysis summary contains no rows."
    )
  }
  
  score_mode <- as.character(
    score_mode
  )[1]
  
  analysis_config <- as.character(
    analysis_config
  )[1]
  
  #---------------------------------------------------------------------------
  # Consistency check for the LoF analysis
  #---------------------------------------------------------------------------
  
  if (!analysis_config %in% names(analysis_configs)) {
    stop(
      "Unknown analysis configuration for the flowchart: ",
      analysis_config
    )
  }
  
  config <- analysis_configs[[analysis_config]]
  
  if (
    score_mode == "lof_burden" &&
    config$type != "lof_vep"
  ) {
    stop(
      "score_mode = 'lof_burden' is only permitted for ",
      "config$type = 'lof_vep'."
    )
  }
  
  if (
    config$type == "lof_vep" &&
    score_mode != "lof_burden"
  ) {
    stop(
      "A LoF VEP configuration may only be used with ",
      "score_mode = 'lof_burden' dargestellt werden."
    )
  }
  
  #---------------------------------------------------------------------------
  # Dynamic label of the score mode
  #---------------------------------------------------------------------------
  
  score_label <- switch(
    score_mode,
    
    oriented =
      "Oriented SNP SumScore",
    
    raw =
      "Raw SNP SumScore",
    
    beta_weighted =
      "Beta-weighted Genetic Score",
    
    ridge_weighted =
      "Ridge-weighted Genetic Score",
    
    elasticnet_weighted =
      "Elastic-Net-weighted Genetic Score",
    
    lof_burden =
      "Loss-of-Function Allele Burden",
    
    "Genetic Score"
  )
  
  score_description <- switch(
    score_mode,
    
    oriented =
      "Protective SNPs reoriented; unweighted score",
    
    raw =
      "Original allele dosage; unweighted score",
    
    beta_weighted =
      "Allele dosage weighted by training-set log(OR)",
    
    ridge_weighted =
      "Allele dosage weighted by ridge-regularized (cv.glmnet, lambda.1se) training-set coefficients",
    
    elasticnet_weighted =
      "Allele dosage weighted by elastic-net-regularized (alpha=0.5, cv.glmnet, lambda.1se) training-set coefficients",
    
    lof_burden =
      paste0(
        "VEP-annotated LoF SNPs; allele dosage oriented ",
        "to the LoF allele; unweighted burden"
      ),
    
    ""
  )
  
  #---------------------------------------------------------------------------
  # Prepare steps
  #---------------------------------------------------------------------------
  
  df_plot[, Step_clean := as.character(Step)]
  
  # Remove only purely technical or redundant steps.
  # The score step is kept so that beta_weighted is visible.
  df_plot <- df_plot[
    !Step_clean %in% c(
      "Covariate Selection",
      "Multivariate Regression"
    )
  ]
  
  if (nrow(df_plot) == 0) {
    stop(
      "After filtering, no flowchart steps remain."
    )
  }
  
  # Rename the first step only if it actually
  # describes the starting dataset
  if (
    df_plot$Step_clean[1] %in% c(
      "Analysis Dataset",
      "Input Dataset",
      "Initial Dataset",
      "Dataset"
    )
  ) {
    
    df_plot[
      1,
      Step_clean := "Final Analysis Dataset"
    ]
  }
  
  # Name the score step mode-specifically
  df_plot[
    Step_clean %in% c(
      "SNP Sum Score",
      "SNP_SumScore",
      "Score Calculation",
      "Genetic Score"
    ),
    Step_clean := score_label
  ]
  
  # If the summary does not yet contain an explicit score step,
  # add one after the starting dataset
  score_step_present <- any(
    df_plot$Step_clean %in% c(
      "Oriented SNP SumScore",
      "Raw SNP SumScore",
      "Beta-weighted Genetic Score",
      "Loss-of-Function Allele Burden",
      "Genetic Score"
    )
  )
  
  if (!score_step_present) {
    
    first_row <- copy(
      df_plot[1]
    )
    
    score_row <- data.table(
      Step = "SNP Sum Score",
      SNPs = first_row$SNPs,
      Individuals = first_row$Individuals,
      Step_clean = score_label
    )
    
    # Fill further columns of the summary
    additional_columns <- setdiff(
      colnames(df_plot),
      colnames(score_row)
    )
    
    for (column_name in additional_columns) {
      score_row[
        ,
        (column_name) := NA
      ]
    }
    
    score_row <- score_row[
      ,
      colnames(df_plot),
      with = FALSE
    ]
    
    df_plot <- rbindlist(
      list(
        df_plot[1],
        score_row,
        df_plot[-1]
      ),
      use.names = TRUE,
      fill = TRUE
    )
  }
  
  #---------------------------------------------------------------------------
  # Positions and colours
  #---------------------------------------------------------------------------
  
  df_plot[, x := 1]
  
  df_plot[, y := rev(
    seq_len(.N)
  )]
  
  df_plot[, Fill := dplyr::case_when(
    
    Step_clean %in% c(
      "PCA"
    ) ~ "#DDEBF7",
    
    Step_clean %in% c(
      "Oriented SNP SumScore",
      "Raw SNP SumScore",
      "Beta-weighted Genetic Score",
      "Loss-of-Function Allele Burden",
      "Genetic Score"
    ) ~ "#FCE4D6",
    
    Step_clean %in% c(
      "Outlier Removed",
      "Complete Cases (Covariates)",
      "Final Analysis Population"
    ) ~ "#E2F0D9",
    
    TRUE ~ "#FFF3E6"
  )]
  
  df_plot[, Label := paste0(
    Step_clean,
    "\n",
    "Individuals: ",
    Individuals,
    "    SNPs: ",
    SNPs
  )]
  
  #---------------------------------------------------------------------------
  # Plot
  #---------------------------------------------------------------------------
  
  box_width <- 0.62
  box_height <- 0.52
  
  arrow_data <- df_plot[
    y > min(y)
  ]
  
  p <- ggplot(
    df_plot,
    aes(
      x = x,
      y = y
    )
  ) +
    geom_segment(
      data = arrow_data,
      aes(
        x = 1,
        xend = 1,
        y = y - box_height / 2 - 0.04,
        yend = y - 1 + box_height / 2 + 0.04
      ),
      arrow = arrow(
        length = unit(
          0.14,
          "cm"
        ),
        type = "closed"
      ),
      linewidth = 0.35,
      color = "grey35"
    ) +
    geom_rect(
      aes(
        xmin = x - box_width / 2,
        xmax = x + box_width / 2,
        ymin = y - box_height / 2,
        ymax = y + box_height / 2,
        fill = Fill
      ),
      color = "grey25",
      linewidth = 0.4
    ) +
    geom_text(
      aes(
        label = Label
      ),
      size = 3,
      lineheight = 0.9
    ) +
    scale_fill_identity() +
    labs(
      title = paste0(
        pathway_name,
        " – Multivariate Analysis Flowchart"
      ),
      subtitle = paste0(
        score_label,
        " | ",
        score_description,
        "\nAnalysis configuration: ",
        analysis_config
      )
    ) +
    coord_cartesian(
      xlim = c(
        0.2,
        1.8
      ),
      ylim = c(
        min(df_plot$y) - 0.6,
        max(df_plot$y) + 0.6
      ),
      clip = "off"
    ) +
    theme_void(
      base_size = 13
    ) +
    theme(
      plot.title = element_text(
        face = "plain",
        size = 13,
        hjust = 0.5,
        margin = margin(
          b = 5
        )
      ),
      plot.subtitle = element_text(
        face = "plain",
        size = 9,
        hjust = 0.5,
        margin = margin(
          b = 10
        )
      ),
      plot.margin = margin(
        15,
        20,
        15,
        20
      )
    )
  
  #---------------------------------------------------------------------------
  # Save
  #---------------------------------------------------------------------------
  
  plot_height <- max(
    6,
    0.75 * nrow(df_plot)
  )
  
  ggsave(
    paste0(
      outfile_prefix,
      ".pdf"
    ),
    plot = p,
    width = 7.5,
    height = plot_height,
    device = "pdf"
  )
  
  ggsave(
    paste0(
      outfile_prefix,
      ".png"
    ),
    plot = p,
    width = 7.5,
    height = plot_height,
    dpi = 600
  )
  
  if (show_plot) {
    print(p)
  }
  
  return(p)
}

# FINAL MULTIVARIABLE ANALYSIS – AFTER MANUAL COVARIATE SELECTION
run_multivariate_sumscore_analysis <- function(df_score, result_name) {
  
  cat("\n====================================================\n")
  cat("Multivariable analysis:", result_name, "\n")
  cat("====================================================\n")
  
  if (is.null(df_score)) {
    warning("No SNP sum score dataset for ", result_name)
    return(NULL)
  }
  
  score_mode <- unique(
    as.character(df_score$Score_Mode)
  )
  
  score_mode <- score_mode[
    !is.na(score_mode)
  ]
  
  analysis_config <- unique(
    as.character(df_score$Analysis_Config)
  )
  
  analysis_config <- analysis_config[
    !is.na(analysis_config)
  ]
  
  if (
    length(score_mode) != 1 ||
    length(analysis_config) != 1
  ) {
    
    warning(
      "Score_Mode or Analysis_Config is not unique for ",
      result_name,
      ""
    )
    
    return(NULL)
  }
  
  score_mode <- score_mode[1]
  analysis_config <- analysis_config[1]
  
  if (!analysis_config %in% names(analysis_configs)) {
    warning(
      "Unbekannte Analysekonfiguration: ",
      analysis_config
    )
    return(NULL)
  }
  
  config <- analysis_configs[[analysis_config]]
  
  #============================================================================
  # Consistency check of analysis and score configuration
  #============================================================================
  
  if (
    score_mode == "lof_burden" &&
    config$type != "lof_vep"
  ) {
    warning(
      "Invalid combination for ",
      result_name,
      ": score_mode = 'lof_burden' is only permitted with ",
      "config$type = 'lof_vep'."
    )
    return(NULL)
  }
  
  if (
    config$type == "lof_vep" &&
    score_mode != "lof_burden"
  ) {
    warning(
      "Invalid combination for ",
      result_name,
      ": a LoF VEP configuration may only be used with ",
      "score_mode = 'lof_burden' analysiert werden."
    )
    return(NULL)
  }
  
  if (score_mode == "lof_burden") {
    cat(
      "LoF burden model: the burden already oriented to the VEP LoF allele ",
      "ausgerichtete Burden wird multivariat analysiert.\n",
      sep = ""
    )
  }
  
  #============================================================================
  # 1) Base data
  #============================================================================
  
  df <- df_score
  
  stopifnot("id" %in% colnames(df))
  stopifnot("obesity" %in% colnames(df))
  stopifnot("SNP_SumScore" %in% colnames(df))
  
  #============================================================================
  # Check score variation
  #============================================================================
  
  score_values <- suppressWarnings(
    as.numeric(
      df$SNP_SumScore
    )
  )
  
  finite_score_values <- score_values[
    is.finite(score_values)
  ]
  
  if (
    length(finite_score_values) == 0 ||
    length(unique(finite_score_values)) < 2
  ) {
    warning(
      "The genetic score has insufficient variation for ",
      result_name,
      "."
    )
    return(NULL)
  }
  
  if (score_mode == "lof_burden") {
    
    cat(
      "LoF burden 0 for ",
      sum(finite_score_values == 0),
      " of ",
      length(finite_score_values),
      " Personen.\n",
      sep = ""
    )
  }
  
  #============================================================================
  # 2) Recoding
  #============================================================================
  
  df$gender_num <- as.numeric(as.character(df$gender))
  df$gender_num[df$gender_num %in% c(0, 7777777, 9999999)] <- NA
  
  df$gender_cat <- factor(
    df$gender_num,
    levels = c(1, 2),
    labels = c("male", "female")
  )
  
  df$smoking_habits_3m_num <- as.numeric(as.character(df$smoking_habits_3m))
  df$smoking_habits_3m_num[
    df$smoking_habits_3m_num %in% c(
      0,
      7777777,
      9999999
    )
  ] <- NA
  
  df$smoking_cat <- factor(
    dplyr::case_when(
      df$smoking_habits_3m_num %in% c(1, 2) ~ "never_smoker",
      df$smoking_habits_3m_num == 3 ~ "former_smoker",
      df$smoking_habits_3m_num == 4 ~ "current_smoker",
      TRUE ~ NA_character_
    ),
    levels = c(
      "never_smoker",
      "former_smoker",
      "current_smoker"
    )
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
  
  #============================================================================
  # 3) Batch
  #============================================================================
  
  batch_col <- NULL
  
  if ("Batch" %in% colnames(df)) {
    
    batch_col <- "Batch"
    
  } else if ("batch" %in% colnames(df)) {
    
    batch_col <- "batch"
  }
  
  if (!is.null(batch_col)) {
    
    df$Batch <- factor(
      as.character(df[[batch_col]]),
      levels = c(
        "Ambu",
        "MD"
      )
    )
    
  } else {
    
    df$Batch <- NULL
  }
  
  #============================================================================
  # 4a) Load genome-wide PCs (population stratification)
  #============================================================================
  # NO LONGER computed pathway-locally from the pathway SNPs, but
  # generated once in the gene script from the complete GSA array (QC + LD pruning +
  # KING-cutoff sensitivity check + PCA) and loaded here unchanged
  # (see genomewide_pcs_file above).
  
  df$id <- as.character(df$id)
  
  global_gPCs <- readRDS(genomewide_pcs_file)
  global_gPCs$id <- as.character(global_gPCs$id)
  
  n_before_gpc_join <- nrow(df)
  
  df <- merge(
    df,
    global_gPCs,
    by = "id",
    all.x = FALSE
  )
  
  if (nrow(df) < n_before_gpc_join) {
    warning(
      n_before_gpc_join - nrow(df),
      " Personen in ", result_name,
      " have no genome-wide PCs and were removed."
    )
  }
  
  if (nrow(df) == 0) {
    stop(
      "No individuals with genome-wide PCs for ",
      result_name,
      "."
    )
  }
  
  #============================================================================
  # 4b) Take over outliers from the covariate pre-check
  #============================================================================
  
  # Number after genotype/phenotype matching and PCA
  n_after_id_matching <- nrow(df)
  
  # Retrieve matching pre-check result
  precheck <- all_covariate_prechecks[[result_name]]
  
  if (is.null(precheck)) {
    stop(
      "No covariate pre-check for ",
      result_name,
      " available."
    )
  }
  
  if (is.null(precheck$df) || !"id" %in% colnames(precheck$df)) {
    stop(
      "The covariate pre-check contains no cleaned dataset ",
      "with ID column for ",
      result_name,
      "."
    )
  }
  
  # Keep IDs that passed the outlier pre-check
  ids_after_outlier_removal <- unique(
    as.character(precheck$df$id)
  )
  
  df <- df %>%
    filter(
      as.character(id) %in% ids_after_outlier_removal
    )
  
  # Number after outlier removal
  n_after_outlier_removal <- nrow(df)
  
  cat("\n========================================\n")
  cat("Outlier removal applied\n")
  cat("========================================\n")
  cat("After ID matching/PCA:   ", n_after_id_matching, "\n")
  cat("After outlier removal:", n_after_outlier_removal, "\n")
  cat("Removed:               ",
      n_after_id_matching - n_after_outlier_removal,
      "\n")
  
  #============================================================================
  # 5) Covariates
  #============================================================================
  
  covariates_keep <- manual_covariate_decisions$default$keep
  covariates_removed <- manual_covariate_decisions$default$remove
  
  # Use only covariates actually present in df
  covariates_present <- intersect(
    covariates_keep,
    colnames(df)
  )
  
  # Document missing covariates
  missing_covariates <- setdiff(
    covariates_keep,
    colnames(df)
  )
  
  if (length(missing_covariates) > 0) {
    warning(
      "The following selected covariates are missing in df and are not used: ",
      paste(missing_covariates, collapse = ", ")
    )
  }
  
  cov_categorical <- intersect(
    covariates_present,
    c(
      "gender_cat",
      "smoking_cat",
      "school_cert_cat",
      "graduation_cat",
      "Batch"
    )
  )
  
  cov_numeric <- setdiff(
    covariates_present,
    cov_categorical
  )
  
  # data.table-conform selection of factor variables
  factor_vars <- c(
    covariates_present[
      vapply(
        df[
          ,
          covariates_present,
          drop = FALSE
        ],
        is.factor,
        logical(1)
      )
    ],
    "obesity"
  )
  
  factor_vars <- unique(
    factor_vars
  )
  
  #============================================================================
  # 6) Covariate pre-analyses
  #============================================================================
  cat("⚠️ Covariate selection based on VIF, correlation, Cramer's V and manual content-based decision.\n")
  cat("\n")
  cat("⚠️ Final covariate selection\n")
  
  cat("KEEP:\n")
  cat(paste(covariates_keep, collapse = ", "), "\n\n")
  
  cat("REMOVE:\n")
  cat(paste(covariates_removed, collapse = ", "), "\n\n")
  
  covariate_selection <- data.frame(
    Variable = c(covariates_keep, covariates_removed),
    Manual_Decision = c(
      rep("KEEP", length(covariates_keep)),
      rep("REMOVE_OR_SENSITIVITY_ONLY", length(covariates_removed))
    ),
    Comment = rep(
      manual_covariate_decisions$default$comment,
      length(c(covariates_keep, covariates_removed))
    )
  )
  
  write.csv2(
    covariate_selection,
    make_output_name(
      pathway_name,
      "Multivariate",
      "Covariate_Precheck",
      "Covariate_Selection",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  #============================================================================
  # 7) Factor-Level-Check
  #============================================================================
  
  # Ensure that only existing variables are checked
  factor_vars <- unique(
    factor_vars[
      factor_vars %in% colnames(df)
    ]
  )
  
  if (length(factor_vars) == 0) {
    
    warning("No factor variables available for the factor-level check.")
    
    factor_level_table <- data.frame(
      Variable = character(),
      Level = character(),
      Count = numeric(),
      stringsAsFactors = FALSE
    )
    
    factor_check_summary <- data.frame(
      Variable = character(),
      IsFactor = logical(),
      N_Levels = integer(),
      Min_Level_Count = numeric(),
      Has_At_Least_2_Levels = logical(),
      Flag = character(),
      stringsAsFactors = FALSE
    )
    
  } else {
    
    factor_level_table <- dplyr::bind_rows(
      lapply(factor_vars, function(v) {
        
        tab <- table(
          df[[v]],
          useNA = "ifany"
        )
        
        data.frame(
          Variable = v,
          Level = names(tab),
          Count = as.numeric(tab),
          stringsAsFactors = FALSE
        )
      })
    )
    
    factor_check_summary <- dplyr::bind_rows(
      lapply(factor_vars, function(v) {
        
        tab <- table(
          df[[v]],
          useNA = "no"
        )
        
        n_levels <- length(tab)
        min_level_count <- if (n_levels > 0) {
          min(as.numeric(tab))
        } else {
          NA_real_
        }
        
        data.frame(
          Variable = v,
          IsFactor = is.factor(df[[v]]),
          N_Levels = n_levels,
          Min_Level_Count = min_level_count,
          Has_At_Least_2_Levels = n_levels >= 2,
          Flag = dplyr::case_when(
            !is.factor(df[[v]]) && v != "obesity" ~ "NOT_FACTOR",
            n_levels < 2 ~ "ONLY_ONE_LEVEL",
            !is.na(min_level_count) && min_level_count < 10 ~ "LOW_LEVEL_COUNT",
            TRUE ~ "OK"
          ),
          stringsAsFactors = FALSE
        )
      })
    )
  }
  
  write.csv2(
    factor_level_table,
    make_output_name(
      pathway_name,
      "Multivariate",
      "Covariate_Precheck",
      "FactorLevels",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  write.csv2(
    factor_check_summary,
    make_output_name(
      pathway_name,
      "Multivariate",
      "Covariate_Precheck",
      "FactorLevel_Check",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  #============================================================================
  # 8) Final analysis dataset
  #============================================================================
  
  required_vars <- c(
    "id",
    "obesity",
    "SNP_SumScore",
    covariates_present
  )
  
  n_before_complete_case <- nrow(df)
  
  df_final <- df %>%
    filter(if_all(all_of(required_vars), ~ !is.na(.)))
  
  n_after_complete_case <- nrow(df_final)
  
  if (nrow(df_final) < 30) {
    warning("Too few complete cases for ", result_name)
    return(NULL)
  }
  
  if (length(unique(df_final$obesity)) < 2) {
    warning("Only one outcome class for ", result_name)
    return(NULL)
  }
  
  write.csv2(
    df_final,
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "Individuals",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  saveRDS(
    df_final,
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "Individuals",
      config,
      score_mode,
      ext = "rds"
    )
  )
  
  #============================================================================
  # Analysis Flow Summary
  #============================================================================
  
  analysis_snps <- attr(df_score, "n_score_snps")
  
  if (is.null(analysis_snps) || is.na(analysis_snps)) {
    analysis_snps <- NA_integer_
  }
  
  analysis_flow <- data.frame(
    Pathway = pathway_name,
    Analysis = result_name,
    Score_Mode = score_mode,
    Analysis_Config = analysis_config,
    Step = c(
      "Final Analysis Dataset",
      "Genomewide PCs Joined",
      "Outlier Removed",
      "Covariate Selection",
      "Complete Cases (Covariates)",
      "Final Analysis Population",
      "Multivariate Regression",
      "SNP Sum Score"
    ),
    SNPs = rep(analysis_snps, 8),
    Individuals = c(
      nrow(df_score),
      n_after_id_matching,
      n_after_outlier_removal,
      n_after_outlier_removal,
      n_after_complete_case,
      nrow(df_final),
      nrow(df_final),
      nrow(df_final)
    ),
    stringsAsFactors = FALSE
  )
  
  print(analysis_flow)
  
  write.csv2(
    analysis_flow,
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "Analysis_Flow_Summary",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  flowchart_prefix <- sub(
    "\\.csv$",
    "",
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "Analysis_Flowchart",
      config,
      score_mode
    )
  )
  
  analysis_flow_plot <- plot_analysis_flow_summary(
    summary_df = analysis_flow,
    pathway_name = pathway_name,
    outfile_prefix = flowchart_prefix,
    score_mode = score_mode,
    analysis_config = analysis_config,
    show_plot = TRUE
  )
  
  numeric_summary <- bind_rows(
    lapply(c("SNP_SumScore", cov_numeric), function(v) {
      
      data.frame(
        Variable = v,
        N = sum(!is.na(df_final[[v]])),
        Mean = mean(df_final[[v]], na.rm = TRUE),
        SD = sd(df_final[[v]], na.rm = TRUE),
        Median = median(df_final[[v]], na.rm = TRUE),
        Q1 = quantile(df_final[[v]], 0.25, na.rm = TRUE),
        Q3 = quantile(df_final[[v]], 0.75, na.rm = TRUE),
        Min = min(df_final[[v]], na.rm = TRUE),
        Max = max(df_final[[v]], na.rm = TRUE),
        Type = "numeric"
      )
    })
  )
  
  categorical_summary <- bind_rows(
    lapply(factor_vars, function(v) {
      
      tab <- table(df_final[[v]])
      
      data.frame(
        Variable = v,
        Level = names(tab),
        N = as.numeric(tab),
        Percent = round(100 * as.numeric(tab) / sum(tab), 2),
        Type = "categorical"
      )
    })
  )
  
  write.csv2(
    numeric_summary,
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "Demographics_Covariates_NumericSummary",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  write.csv2(
    categorical_summary,
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "Demographics_Covariates_CategoricalSummary",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  #============================================================================
  # 9) Multivariable model train/test
  #============================================================================
  
  covariates_formula <- paste(
    covariates_present,
    collapse = " + "
  )
  
  model_formula <- reformulate(
    termlabels = c(
      "SNP_SumScore",
      covariates_present
    ),
    response = "obesity"
  )
  
  model_formula_standardized <- reformulate(
    termlabels = c(
      "SNP_SumScore_SD",
      covariates_present
    ),
    response = "obesity"
  )
  
  global_split <- readRDS(
    output_file(
      "Pathway_",
      pathway_name,
      "_Global_TrainTest_Split.rds"
    )
  )
  
  global_split$train_ids <- as.character(
    global_split$train_ids
  )
  
  global_split$test_ids <- as.character(
    global_split$test_ids
  )
  
  df_final$id <- as.character(
    df_final$id
  )
  
  train_set <- df_final %>%
    filter(id %in% global_split$train_ids)
  
  test_set <- df_final %>%
    filter(id %in% global_split$test_ids)
  
  if (
    nrow(train_set) < 20 ||
    nrow(test_set) < 10
  ) {
    
    warning(
      "Too few individuals in the train or test dataset for ",
      result_name,
      ". Train: ",
      nrow(train_set),
      ", Test: ",
      nrow(test_set)
    )
    
    return(NULL)
  }
  
  if (
    length(unique(train_set$obesity)) < 2 ||
    length(unique(test_set$obesity)) < 2
  ) {
    
    warning(
      "Train or test dataset contains only one outcome class for ",
      result_name
    )
    
    return(NULL)
  }
  
  #------------------------------------------------------------
  # Standardize SNP score using the training set
  #------------------------------------------------------------
  
  score_mean_train <- mean(
    train_set$SNP_SumScore,
    na.rm = TRUE
  )
  
  score_sd_train <- sd(
    train_set$SNP_SumScore,
    na.rm = TRUE
  )
  
  if (
    !is.finite(score_sd_train) ||
    score_sd_train <= 0
  ) {
    
    warning(
      "SNP sum score has insufficient variation for ",
      result_name
    )
    
    return(NULL)
  }
  
  train_set$SNP_SumScore_SD <-
    (
      train_set$SNP_SumScore -
        score_mean_train
    ) /
    score_sd_train
  
  test_set$SNP_SumScore_SD <-
    (
      test_set$SNP_SumScore -
        score_mean_train
    ) /
    score_sd_train
  
  #------------------------------------------------------------
  # Fit models
  #------------------------------------------------------------
  model <- glm(
    model_formula,
    data = train_set,
    family = binomial
  )
  
  model_standardized <- glm(
    model_formula_standardized,
    data = train_set,
    family = binomial
  )
  
  coef_sd <- summary(
    model_standardized
  )$coefficients[
    "SNP_SumScore_SD",
    ,
    drop = FALSE
  ]
  
  ci_sd <- tryCatch(
    exp(
      confint.default(
        model_standardized
      )[
        "SNP_SumScore_SD",
        ,
        drop = TRUE
      ]
    ),
    error = function(e) {
      c(
        NA_real_,
        NA_real_
      )
    }
  )
  
  OR_per_SD <- exp(
    coef_sd[
      "SNP_SumScore_SD",
      "Estimate"
    ]
  )
  
  #------------------------------------------------------------
  # Regression results of the model with original score
  #------------------------------------------------------------
  
  res <- broom::tidy(
    model,
    conf.int = TRUE,
    conf.level = 0.95
  ) %>%
    rename(
      CI_lower = conf.low,
      CI_upper = conf.high
    ) %>%
    mutate(
      OR = exp(estimate),
      OR_CI_lower = exp(CI_lower),
      OR_CI_upper = exp(CI_upper),
      Pathway = pathway_name,
      Analysis = result_name,
      Score_Mode = score_mode,
      Analysis_Config = analysis_config,
      Score_Type = dplyr::case_when(
        score_mode == "oriented" ~
          "Obesity-oriented SNP SumScore",
        
        score_mode == "raw" ~
          "Unweighted raw SNP SumScore",
        
        score_mode == "beta_weighted" ~
          "Beta-weighted genetic score",
        
        score_mode == "ridge_weighted" ~
          "Ridge-weighted genetic score",
        
        score_mode == "elasticnet_weighted" ~
          "Elastic-Net-weighted genetic score",
        
        score_mode == "lof_burden" ~
          "Loss-of-function allele burden",
        
        TRUE ~
          "Genetic score"
      ),
      
      Score_Unit = dplyr::case_when(
        score_mode == "lof_burden" ~
          "One additional LoF allele copy",
        
        score_mode == "oriented" ~
          "One additional obesity-oriented allele",
        
        score_mode == "raw" ~
          "One additional counted allele",
        
        score_mode == "beta_weighted" ~
          "One beta-weighted score unit",
        
        score_mode == "ridge_weighted" ~
          "One ridge-weighted score unit",
        
        score_mode == "elasticnet_weighted" ~
          "One elastic-net-weighted score unit",
        
        TRUE ~
          "One score unit"
      ),
      N_total = nrow(df_final),
      N_train = nrow(train_set),
      N_test = nrow(test_set),
      AIC = AIC(model),
      BIC = BIC(model),
      McFadden_R2 = as.numeric(
        pscl::pR2(model)["McFadden"]
      ),
      Score_Mean_Train = score_mean_train,
      Score_SD_Train = score_sd_train,
      OR_Score_Per_SD = OR_per_SD,
      OR_Score_Per_SD_CI_lower = ci_sd[1],
      OR_Score_Per_SD_CI_upper = ci_sd[2]
    )
  
  write.csv2(
    res,
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "Regression_Results",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  #============================================================================
  # 10) Marginal Effects
  #============================================================================
  
  marg_effects <- tryCatch({
    margins_obj <- margins::margins(model)
    as.data.frame(summary(margins_obj))
  }, error = function(e) NULL)
  
  if (!is.null(marg_effects)) {
    write.csv2(
      marg_effects,
      make_output_name(
        pathway_name,
        "Multivariate",
        "SNP_SumScore",
        "MarginalEffects",
        config,
        score_mode
      ),
      row.names = FALSE
    )
  }
  
  #============================================================================
  # 11) ROC + cutoffs test set
  #============================================================================
  
  train_set$pred_prob <- predict(
    model,
    newdata = train_set,
    type = "response"
  )
  
  test_set$pred_prob <- predict(
    model,
    newdata = test_set,
    type = "response"
  )
  
  roc_train <- roc(
    train_set$obesity,
    train_set$pred_prob,
    quiet = TRUE
  )
  
  roc_obj <- roc(
    test_set$obesity,
    test_set$pred_prob,
    quiet = TRUE
  )
  
  auc_val <- as.numeric(auc(roc_obj))
  
  auc_ci <- tryCatch(
    as.numeric(ci.auc(roc_obj)),
    error = function(e) {
      warning(
        "ci.auc() fehlgeschlagen: ",
        conditionMessage(e),
        " -- AUC confidence interval is set to NA."
      )
      c(NA_real_, auc_val, NA_real_)
    }
  )
  
  youden_cutoff <- tryCatch(
    as.numeric(
      coords(
        roc_train,
        x = "best",
        best.method = "youden",
        ret = "threshold"
      )
    ),
    error = function(e) NA_real_
  )
  
  calc_metrics <- function(obs, pred_prob, cutoff) {
    
    pred_class <- ifelse(pred_prob >= cutoff, 1, 0)
    
    TP <- sum(obs == 1 & pred_class == 1)
    TN <- sum(obs == 0 & pred_class == 0)
    FP <- sum(obs == 0 & pred_class == 1)
    FN <- sum(obs == 1 & pred_class == 0)
    
    data.frame(
      Sensitivity = ifelse((TP + FN) > 0, TP / (TP + FN), NA),
      Specificity = ifelse((TN + FP) > 0, TN / (TN + FP), NA),
      PPV = ifelse((TP + FP) > 0, TP / (TP + FP), NA),
      NPV = ifelse((TN + FN) > 0, TN / (TN + FN), NA),
      Accuracy = (TP + TN) / length(obs)
    )
  }
  
  metrics_05 <- calc_metrics(test_set$obesity, test_set$pred_prob, 0.5)
  metrics_youden <- calc_metrics(test_set$obesity, test_set$pred_prob, youden_cutoff)
  
  beta0 <- coef(model)["(Intercept)"]
  beta_snp <- coef(model)["SNP_SumScore"]
  
  calc_score_cutoff <- function(prob_cutoff, beta0, beta_snp) {
    if (is.na(beta_snp) || beta_snp == 0) return(NA)
    (log(prob_cutoff / (1 - prob_cutoff)) - beta0) / beta_snp
  }
  
  score_cutoff_05 <- NA_real_
  score_cutoff_youden <- NA_real_
  
  formula_no_snp <- reformulate(
    termlabels = covariates_present,
    response = "obesity"
  )
  
  model_no_snp <- glm(
    formula_no_snp,
    data = train_set,
    family = binomial
  )
  
  lrt_pvalue <- anova(model_no_snp, model, test = "LRT")$`Pr(>Chi)`[2]
  
  snp_coef <- summary(model)$coefficients["SNP_SumScore", ]
  snp_ci <- exp(confint.default(model)["SNP_SumScore", ])
  
  base_row <- data.frame(
    Pathway = pathway_name,
    Analysis = result_name,
    Score_Mode = unique(df_final$Score_Mode)[1],
    Analysis_Config = unique(df_final$Analysis_Config)[1],
    Score_Type = dplyr::case_when(
      score_mode == "oriented" ~
        "Obesity-oriented SNP SumScore",
      
      score_mode == "raw" ~
        "Unweighted raw SNP SumScore",
      
      score_mode == "beta_weighted" ~
        "Beta-weighted genetic score",
      
      score_mode == "ridge_weighted" ~
        "Ridge-weighted genetic score",
      
      score_mode == "elasticnet_weighted" ~
        "Elastic-Net-weighted genetic score",
      
      score_mode == "lof_burden" ~
        "Loss-of-function allele burden",
      
      TRUE ~
        "Genetic score"
    ),
    
    Score_Unit = dplyr::case_when(
      score_mode == "lof_burden" ~
        "One additional LoF allele copy",
      
      score_mode == "oriented" ~
        "One additional obesity-oriented allele",
      
      score_mode == "raw" ~
        "One additional counted allele",
      
      score_mode == "beta_weighted" ~
        "One beta-weighted score unit",
      
      score_mode == "ridge_weighted" ~
        "One ridge-weighted score unit",
      
      score_mode == "elasticnet_weighted" ~
        "One elastic-net-weighted score unit",
      
      TRUE ~
        "One score unit"
    ),
    Estimate = snp_coef["Estimate"],
    StdError = snp_coef["Std. Error"],
    OR = exp(snp_coef["Estimate"]),
    CI_lower = snp_ci[1],
    CI_upper = snp_ci[2],
    Pvalue = snp_coef["Pr(>|z|)"],
    LRT_pvalue = lrt_pvalue,
    AUC = auc_val,
    AUC_CI_lower = auc_ci[1],
    AUC_CI_median = auc_ci[2],
    AUC_CI_upper = auc_ci[3],
    N_train = nrow(train_set),
    N_test = nrow(test_set),
    Score_Cutoff_Interpretation =
      "Not uniquely defined in multivariable model"
  )
  
  cutoff_table <- bind_rows(
    cbind(
      base_row,
      Cutoff_Type = "Fixed_0.5",
      Cutoff_Prob = 0.5,
      SNP_SumScore_Cutoff = score_cutoff_05,
      metrics_05
    ),
    cbind(
      base_row,
      Cutoff_Type = "Youden",
      Cutoff_Prob = youden_cutoff,
      SNP_SumScore_Cutoff = score_cutoff_youden,
      metrics_youden
    )
  )
  
  write.csv2(
    cutoff_table,
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "Cutoff_Results",
      config,
      score_mode
    ),
    row.names = FALSE
  )
  
  #============================================================================
  # 12) Batch sensitivity optional – main analysis only
  #============================================================================
  
  batch_results <- NULL
  
  if (
    grepl("^main_ld_", result_name) &&
    "Batch" %in% colnames(train_set) &&
    length(unique(na.omit(train_set$Batch))) > 1
  ) {
    
    
    # Do not include batch twice in the formula
    covariates_without_batch <- setdiff(
      covariates_present,
      "Batch"
    )
    
    formula_no_batch <- reformulate(
      termlabels = c(
        "SNP_SumScore",
        covariates_without_batch
      ),
      response = "obesity"
    )
    
    formula_with_batch <- reformulate(
      termlabels = c(
        "SNP_SumScore",
        "Batch",
        covariates_without_batch
      ),
      response = "obesity"
    )
    
    model_no_batch <- glm(
      formula_no_batch,
      data = train_set,
      family = binomial
    )
    
    model_with_batch <- glm(
      formula_with_batch,
      data = train_set,
      family = binomial
    )
    
    coef_no <- summary(model_no_batch)$coefficients["SNP_SumScore", ]
    coef_ba <- summary(model_with_batch)$coefficients["SNP_SumScore", ]
    
    OR_no <- exp(coef_no["Estimate"])
    OR_ba <- exp(coef_ba["Estimate"])
    
    CI_no <- exp(confint.default(model_no_batch)["SNP_SumScore", ])
    CI_ba <- exp(confint.default(model_with_batch)["SNP_SumScore", ])
    
    batch_results <- data.frame(
      Pathway = pathway_name,
      Analysis = result_name,
      Sensitivity_Type = "Batch_Adjusted_Main_Model_Only",
      OR_no_batch = OR_no,
      CI_low_no_batch = CI_no[1],
      CI_high_no_batch = CI_no[2],
      P_no_batch = coef_no["Pr(>|z|)"],
      OR_with_batch = OR_ba,
      CI_low_with_batch = CI_ba[1],
      CI_high_with_batch = CI_ba[2],
      P_with_batch = coef_ba["Pr(>|z|)"],
      Delta_OR = OR_ba - OR_no,
      Delta_percent = 100 * (OR_ba - OR_no) / OR_no,
      N = nrow(train_set)
    )
    
    write.csv2(
      batch_results,
      make_output_name(
        pathway_name,
        "Multivariate",
        "SNP_SumScore",
        "Batch_Sensitivity",
        config,
        score_mode
      ),
      row.names = FALSE
    )
  }
  
  #============================================================================
  # 13) Plots
  #============================================================================
  
  roc_df <- data.frame(
    TPR = rev(roc_obj$sensitivities),
    FPR = rev(1 - roc_obj$specificities)
  )
  
  roc_plot <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_line(linewidth = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    theme_minimal(base_size = 13) +
    labs(
      title = paste0("Pathway_", pathway_name, " – Multivariate Logistic Regression"),
      subtitle = paste("AUC =", round(auc_val, 3)),
      x = "1 - Specificity",
      y = "Sensitivity"
    )
  
  ggsave(
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "ROC",
      config,
      score_mode,
      ext = "pdf"
    ),
    roc_plot,
    width = 7,
    height = 5
  )
  
  ggsave(
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "ROC",
      config,
      score_mode,
      ext = "png"
    ),
    roc_plot,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  res_standardized <- broom::tidy(
    model_standardized,
    conf.int = TRUE,
    conf.level = 0.95
  ) %>%
    rename(
      CI_lower = conf.low,
      CI_upper = conf.high
    ) %>%
    mutate(
      OR = exp(estimate),
      OR_CI_lower = exp(CI_lower),
      OR_CI_upper = exp(CI_upper)
    )
  
  res_plot <- res_standardized %>%
    filter(term != "(Intercept)") %>%
    filter(
      !is.na(OR),
      !is.na(OR_CI_lower),
      !is.na(OR_CI_upper)
    ) %>%
    mutate(
      term_friendly = case_when(
        term == "SNP_SumScore_SD" ~ case_when(
          score_mode == "oriented" ~
            "Oriented SNP SumScore (per 1 SD)",
          
          score_mode == "raw" ~
            "Raw SNP SumScore (per 1 SD)",
          
          score_mode == "beta_weighted" ~
            "Beta-weighted Genetic Score (per 1 SD)",
          
          score_mode == "ridge_weighted" ~
            "Ridge-weighted Genetic Score (per 1 SD)",
          
          score_mode == "elasticnet_weighted" ~
            "Elastic-Net-weighted Genetic Score (per 1 SD)",
          
          score_mode == "lof_burden" ~
            "Loss-of-Function Allele Burden (per 1 SD)",
          
          TRUE ~
            "Genetic Score (per 1 SD)"
        ),
        term == "gender_catfemale" ~ "Gender female",
        term == "smoking_catformer_smoker" ~ "Former smoker",
        term == "smoking_catcurrent_smoker" ~ "Current smoker",
        term == "school_cert_catrisk" ~ "School certificate risk",
        term == "graduation_catrisk" ~ "Graduation risk",
        term == "sample_age" ~ "Age",
        term == "sports_activity_min_week" ~ "Sports activity",
        term == "everyday_activity_min_week" ~ "Everyday activity",
        term == "energy_kcal" ~ "Energy intake",
        term == "gPC1" ~ "Genomic PC1",
        term == "gPC2" ~ "Genomic PC2",
        term == "gPC3" ~ "Genomic PC3",
        TRUE ~ term
      ),
      is_pc = grepl("^gPC", term),
      direction = case_when(
        OR > 1 ~ "Risk",
        OR < 1 ~ "Protective",
        TRUE ~ "Unclear"
      )
    )

  # Derive sensible axis limits from the non-PC covariates (incl. genetic
  # score) so that the partly extreme PC confidence intervals do not distort the
  # scale and direction/magnitude of the other covariates
  # remain readable. coord_flip(ylim = ...) only zooms the view - no
  # data points are removed; error bars of PCs extending beyond the
  # edge are visibly clipped at the plot margin.
  res_plot_non_pc <- res_plot %>% filter(!is_pc)
  if (nrow(res_plot_non_pc) > 0) {
    or_lo <- min(res_plot_non_pc$OR_CI_lower, na.rm = TRUE)
    or_hi <- max(res_plot_non_pc$OR_CI_upper, na.rm = TRUE)
    or_pad <- (log(or_hi) - log(or_lo)) * 0.2
    if (!is.finite(or_pad) || or_pad <= 0) or_pad <- 0.2
    or_lo <- exp(log(or_lo) - or_pad)
    or_hi <- exp(log(or_hi) + or_pad)
  } else {
    or_lo <- min(res_plot$OR_CI_lower, na.rm = TRUE)
    or_hi <- max(res_plot$OR_CI_upper, na.rm = TRUE)
  }

  forest_plot <- ggplot(res_plot, aes(x = reorder(term_friendly, OR), y = OR, color = direction)) +
    geom_point(aes(shape = is_pc), size = 3) +
    geom_errorbar(
      aes(ymin = OR_CI_lower, ymax = OR_CI_upper),
      width = 0.2
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    scale_y_log10() +
    scale_color_manual(
      values = c("Risk" = "firebrick", "Protective" = "steelblue", "Unclear" = "gray50"),
      guide = "none"
    ) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    coord_flip(ylim = c(or_lo, or_hi)) +
    theme_minimal(base_size = 13) +
    labs(
      title = paste0(pathway_name, " – Multivariate SNP Sum Model"),
      x = "",
      y = "Odds ratio, log scale, 95% CI (axis adapted to non-PC covariates)"
    )

  ggsave(
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "ForestPlot",
      config,
      score_mode,
      ext = "pdf"
    ),
    forest_plot,
    width = 7,
    height = 5
  )
  
  ggsave(
    make_output_name(
      pathway_name,
      "Multivariate",
      "SNP_SumScore",
      "ForestPlot",
      config,
      score_mode,
      ext = "png"
    ),
    forest_plot,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  score_label <- switch(
    score_mode,
    
    oriented =
      "Oriented SNP SumScore",
    
    raw =
      "Raw SNP SumScore",
    
    beta_weighted =
      "Beta-weighted Genetic Score",
    
    ridge_weighted =
      "Ridge-weighted Genetic Score",
    
    elasticnet_weighted =
      "Elastic-Net-weighted Genetic Score",
    
    lof_burden =
      "Loss-of-Function Allele Burden",
    
    "Genetic Score"
  )
  
  partial_plot <- tryCatch({
    
    ggeff <- ggpredict(model, terms = "SNP_SumScore [all]")
    
    ggplot(ggeff, aes(x = x, y = predicted)) +
      geom_line(linewidth = 1) +
      geom_ribbon(
        aes(ymin = conf.low, ymax = conf.high),
        alpha = 0.2
      ) +
      theme_minimal(base_size = 13) +
      labs(
        title = paste0(
          pathway_name,
          " – Adjusted ",
          score_label,
          " Effect"
        ),
        x = score_label,
        y = "Predicted probability of obesity"
      )
    
  }, error = function(e) NULL)
  
  if (!is.null(partial_plot)) {
    
    ggsave(
      make_output_name(
        pathway_name,
        "Multivariate",
        "SNP_SumScore",
        "PartialEffect",
        config,
        score_mode,
        ext = "pdf"
      ),
      partial_plot,
      width = 7,
      height = 5
    )
    
    ggsave(
      make_output_name(
        pathway_name,
        "Multivariate",
        "SNP_SumScore",
        "PartialEffect",
        config,
        score_mode,
        ext = "png"
      ),
      partial_plot,
      width = 7,
      height = 5,
      dpi = 300
    )
  }
  
  if (score_mode == "lof_burden") {
    
    lof_coef <- summary(model)$coefficients[
      "SNP_SumScore",
      ,
      drop = FALSE
    ]
    
    cat(
      "LoF burden OR per additional LoF allele copy:",
      round(
        exp(lof_coef[1, "Estimate"]),
        3
      ),
      "\n"
    )
    
    cat(
      "LoF-Burden OR pro 1 SD:",
      round(
        OR_per_SD,
        3
      ),
      "\n"
    )
  }
  
  cat("✅ Done:", result_name, "\n")
  
  return(
    list(
      regression = res,
      cutoff = cutoff_table,
      batch = batch_results,
      model = model,
      model_standardized = model_standardized,
      roc = roc_obj,
      roc_data = data.frame(
        FPR = 1 - roc_obj$specificities,
        TPR = roc_obj$sensitivities
      ),
      df_final = df_final,
      train_set = train_set,
      test_set = test_set,
      score_mean_train = score_mean_train,
      score_sd_train = score_sd_train,
      covariates_keep = covariates_present
    )
  )
}

options(nwarnings = 10000)

all_multivar_results <- list()

for (result_name in names(all_snp_sumscore_results)) {
  
  df_score <- all_snp_sumscore_results[[result_name]]
  
  if (is.null(df_score)) {
    warning(
      "No score dataset for ",
      result_name,
      "."
    )
    next
  }
  
  analysis_result <- tryCatch(
    
    run_multivariate_sumscore_analysis(
      df_score = df_score,
      result_name = result_name
    ),
    
    error = function(e) {
      
      warning(
        "Error in the multivariable analysis ",
        result_name,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  if (!is.null(analysis_result)) {
    all_multivar_results[[result_name]] <- analysis_result
  }
}

"main_ld_oriented" %in% names(all_multivar_results)

w <- names(warnings())
w[grepl("main_ld_oriented", w)]

names(all_multivar_results)

multivar_summary <- bind_rows(
  lapply(
    all_multivar_results,
    function(x) {
      if (is.null(x)) return(NULL)
      x$cutoff
    }
  )
)

multivar_summary_file <- output_file(
  "Pathway_",
  pathway_name,
  "_Multivariate_SNP_SumScore_Cutoff_Results_All_Configurations.csv"
)

write.csv2(
  multivar_summary,
  multivar_summary_file,
  row.names = FALSE
)

#===================================================================================
# MULTIVARIABLE SKAT ANALYSIS – LOOP OVER ANALYSIS CONFIGURATIONS
#===================================================================================
library(SKAT)
library(dplyr)

run_multivariate_skat <- function(
    df_final,
    covariates_keep,
    result_name) {
  
  cat("\n====================================================\n")
  cat("Multivariable SKAT analysis:", result_name, "\n")
  cat("====================================================\n")
  
  #---------------------------------------------------------------------------
  # Check inputs
  #---------------------------------------------------------------------------
  
  if (is.null(df_final)) {
    warning(
      "No df_final for ",
      result_name
    )
    return(NULL)
  }
  
  if (
    !"Analysis_Config" %in% colnames(df_final) ||
    !"Score_Mode" %in% colnames(df_final)
  ) {
    warning(
      "Analysis_Config or Score_Mode missing for ",
      result_name
    )
    return(NULL)
  }
  
  analysis_config <- unique(
    as.character(df_final$Analysis_Config)
  )
  
  analysis_config <- analysis_config[
    !is.na(analysis_config)
  ]
  
  source_score_mode <- unique(
    as.character(df_final$Score_Mode)
  )
  
  source_score_mode <- source_score_mode[
    !is.na(source_score_mode)
  ]
  
  if (length(analysis_config) != 1) {
    warning(
      "Analysis_Config is not unique for ",
      result_name,
      ""
    )
    return(NULL)
  }
  
  analysis_config <- analysis_config[1]
  
  if (length(source_score_mode) == 1) {
    source_score_mode <- source_score_mode[1]
  } else {
    source_score_mode <- NA_character_
  }
  
  if (!analysis_config %in% names(analysis_configs)) {
    warning(
      "Analysis configuration not found: ",
      analysis_config
    )
    return(NULL)
  }
  
  config <- analysis_configs[[analysis_config]]
  
  # The beta-weighted score is not a separate SNP selection.
  # SKAT is based directly on the genotype matrix and is therefore
  # not recomputed for this pure score sensitivity.
  if (identical(config$type, "beta_weighted")) {
    message(
      "SKAT is ",
      analysis_config,
      " skipped for this configuration, as beta_weighted is only an ",
      "alternative score weighting."
    )
    return(NULL)
  }
  
  #---------------------------------------------------------------------------
  # Output file
  #---------------------------------------------------------------------------
  
  skat_file <- make_output_name(
    pathway_name = pathway_name,
    analysis_type = "Multivariate",
    analysis_name = "SKAT",
    result_name = "Results",
    config = config,
    score_mode = NULL
  )
  
  #---------------------------------------------------------------------------
  # Genotype data and SNP selection
  #---------------------------------------------------------------------------
  
  pathway_data <- readRDS(
    config$dataset
  )
  
  geno <- pathway_data$genotypes
  
  if (is.null(rownames(geno))) {
    stop(
      "The genotype matrix has no row names for ",
      result_name,
      "."
    )
  }
  
  valid_snps <- colnames(geno)
  
  if (config$type == "functional_vep") {
    
    functional_snps <- load_functional_vep_snps(
      config$vep_file
    )
    
    valid_snps <- intersect(
      valid_snps,
      functional_snps
    )
  }
  
  if (config$type == "lof_vep") {
    
    lof_orientation_table <- create_lof_orientation_table(
      vep_file = config$vep_file,
      bim_table = lof_bim,
      matrix_snps = valid_snps
    )
    
    lof_orientation_table <- lof_orientation_table %>%
      filter(
        Allele_Match != "No_allele_match",
        !is.na(Invert_For_LoF_Score)
      ) %>%
      distinct(
        SNP,
        .keep_all = TRUE
      )
    
    valid_snps <- intersect(
      valid_snps,
      lof_orientation_table$SNP
    )
  }
  
  if (config$type == "database_kegg_only") {
    
    snp_gene_map <- read.csv2(
      config$snp_gene_map,
      stringsAsFactors = FALSE
    )
    
    required_map_columns <- c(
      "SNP",
      "GENE"
    )
    
    missing_map_columns <- setdiff(
      required_map_columns,
      colnames(snp_gene_map)
    )
    
    if (length(missing_map_columns) > 0) {
      stop(
        "The following columns are missing in the SNP-gene mapping: ",
        paste(
          missing_map_columns,
          collapse = ", "
        )
      )
    }
    
    kegg_only_snps <- snp_gene_map %>%
      filter(
        GENE %in% config$kegg_genes
      ) %>%
      pull(SNP) %>%
      unique()
    
    valid_snps <- intersect(
      valid_snps,
      kegg_only_snps
    )
  }
  
  valid_snps <- unique(
    valid_snps
  )
  
  if (length(valid_snps) < 2) {
    warning(
      "Too few SNPs for multivariable SKAT in ",
      result_name
    )
    return(NULL)
  }
  
  geno <- geno[
    ,
    valid_snps,
    drop = FALSE
  ]
  
  #---------------------------------------------------------------------------
  # Match IDs
  #---------------------------------------------------------------------------
  
  df_final <- as.data.frame(
    df_final,
    stringsAsFactors = FALSE
  )
  
  if (!"id" %in% colnames(df_final)) {
    stop(
      "ID column missing in the multivariable dataset for ",
      result_name,
      "."
    )
  }
  
  df_final$id <- as.character(
    df_final$id
  )
  
  rownames(geno) <- as.character(
    rownames(geno)
  )
  
  if (anyDuplicated(df_final$id) > 0) {
    stop(
      "Duplicate IDs in the multivariable dataset for ",
      result_name,
      "."
    )
  }
  
  if (anyDuplicated(rownames(geno)) > 0) {
    stop(
      "Duplicate IDs in the genotype matrix for ",
      result_name,
      "."
    )
  }
  
  # Keep order of df_final
  common_ids <- df_final$id[
    df_final$id %in% rownames(geno)
  ]
  
  if (length(common_ids) < 20) {
    warning(
      "Too few shared IDs for multivariable SKAT in ",
      result_name
    )
    return(NULL)
  }
  
  df_clean <- df_final[
    match(
      common_ids,
      df_final$id
    ),
    ,
    drop = FALSE
  ]
  
  geno_mat <- geno[
    common_ids,
    ,
    drop = FALSE
  ]
  
  if (!identical(
    as.character(df_clean$id),
    as.character(rownames(geno_mat))
  )) {
    stop(
      "ID order for multivariable SKAT in ",
      result_name,
      " does not match."
    )
  }
  
  #---------------------------------------------------------------------------
  # Convert genotype matrix to numeric
  #---------------------------------------------------------------------------
  
  geno_mat <- as.matrix(
    geno_mat
  )
  
  storage.mode(
    geno_mat
  ) <- "numeric"
  
  genotype_range <- range(
    geno_mat,
    na.rm = TRUE
  )
  
  if (
    !all(is.finite(genotype_range)) ||
    genotype_range[1] < 0 ||
    genotype_range[2] > 2
  ) {
    stop(
      "Genotype values outside the range 0–2 for ",
      result_name,
      "."
    )
  }
  
  #---------------------------------------------------------------------------
  # Outcome and factor levels
  #---------------------------------------------------------------------------
  
  if (!"obesity" %in% colnames(df_clean)) {
    stop(
      "Outcome obesity missing for ",
      result_name,
      "."
    )
  }
  
  df_clean$obesity <- as.numeric(
    as.character(df_clean$obesity)
  )
  
  if (!all(
    na.omit(df_clean$obesity) %in% c(0, 1)
  )) {
    stop(
      "obesity must contain exclusively 0/1: ",
      result_name
    )
  }
  
  safe_relevel <- function(x, ref) {
    
    if (
      is.factor(x) &&
      ref %in% levels(x)
    ) {
      return(
        relevel(
          x,
          ref = ref
        )
      )
    }
    
    x
  }
  
  if ("gender_cat" %in% colnames(df_clean)) {
    df_clean$gender_cat <- safe_relevel(
      df_clean$gender_cat,
      "male"
    )
  }
  
  if ("smoking_cat" %in% colnames(df_clean)) {
    df_clean$smoking_cat <- safe_relevel(
      df_clean$smoking_cat,
      "never_smoker"
    )
  }
  
  if ("school_cert_cat" %in% colnames(df_clean)) {
    df_clean$school_cert_cat <- safe_relevel(
      df_clean$school_cert_cat,
      "reference"
    )
  }
  
  if ("graduation_cat" %in% colnames(df_clean)) {
    df_clean$graduation_cat <- safe_relevel(
      df_clean$graduation_cat,
      "reference"
    )
  }
  
  #---------------------------------------------------------------------------
  # Covariates and complete cases
  #---------------------------------------------------------------------------
  
  covariates_keep <- intersect(
    covariates_keep,
    colnames(df_clean)
  )
  
  model_vars <- unique(
    c(
      "obesity",
      covariates_keep
    )
  )
  
  keep_complete <- complete.cases(
    df_clean[
      ,
      model_vars,
      drop = FALSE
    ]
  )
  
  df_clean <- df_clean[
    keep_complete,
    ,
    drop = FALSE
  ]
  
  geno_mat <- geno_mat[
    keep_complete,
    ,
    drop = FALSE
  ]
  
  if (nrow(df_clean) < 20) {
    warning(
      "Too few complete cases for multivariable SKAT in ",
      result_name
    )
    return(NULL)
  }
  
  if (length(unique(df_clean$obesity)) < 2) {
    warning(
      "Only one outcome class available for ",
      result_name
    )
    return(NULL)
  }
  
  # Remove factor variables with only one observed level
  single_level_covariates <- covariates_keep[
    vapply(
      covariates_keep,
      function(v) {
        
        x <- df_clean[[v]]
        
        if (is.factor(x)) {
          nlevels(
            droplevels(x)
          ) < 2
        } else {
          length(unique(x)) < 2
        }
      },
      logical(1)
    )
  ]
  
  if (length(single_level_covariates) > 0) {
    
    warning(
      "The following constant or single-level covariates are removed from ",
      "the SKAT null model: ",
      paste(
        single_level_covariates,
        collapse = ", "
      )
    )
    
    covariates_keep <- setdiff(
      covariates_keep,
      single_level_covariates
    )
  }
  
  # Remove unused factor levels
  factor_columns <- names(df_clean)[
    vapply(
      df_clean,
      is.factor,
      logical(1)
    )
  ]
  
  df_clean[
    factor_columns
  ] <- lapply(
    df_clean[
      factor_columns
    ],
    droplevels
  )
  
  #---------------------------------------------------------------------------
  # Remove monomorphic and non-imputable SNPs
  #---------------------------------------------------------------------------
  
  observed_per_snp <- colSums(
    !is.na(geno_mat)
  )
  
  geno_mat <- geno_mat[
    ,
    observed_per_snp > 0,
    drop = FALSE
  ]
  
  if (ncol(geno_mat) < 2) {
    warning(
      "After excluding completely missing SNPs, too few SNPs for ",
      result_name
    )
    return(NULL)
  }
  
  snp_var <- apply(
    geno_mat,
    2,
    var,
    na.rm = TRUE
  )
  
  variable_snp_mask <- is.finite(
    snp_var
  ) &
    snp_var > 0
  
  geno_mat <- geno_mat[
    ,
    variable_snp_mask,
    drop = FALSE
  ]
  
  if (ncol(geno_mat) < 2) {
    warning(
      "After removing monomorphic SNPs, too few SNPs for ",
      result_name
    )
    return(NULL)
  }
  
  # Missing genotypes with means of the current
  # multivariable analysis dataset
  for (j in seq_len(ncol(geno_mat))) {
    
    missing_j <- is.na(
      geno_mat[, j]
    )
    
    if (any(missing_j)) {
      
      snp_mean <- mean(
        geno_mat[, j],
        na.rm = TRUE
      )
      
      if (!is.finite(snp_mean)) {
        stop(
          "SNP ",
          colnames(geno_mat)[j],
          " cannot be imputed for ",
          result_name,
          "."
        )
      }
      
      geno_mat[
        missing_j,
        j
      ] <- snp_mean
    }
  }
  
  if (anyNA(geno_mat)) {
    stop(
      "After genotype imputation, missing values remain for ",
      result_name,
      " available."
    )
  }
  
  #---------------------------------------------------------------------------
  # SKAT null model
  #---------------------------------------------------------------------------
  
  if (length(covariates_keep) > 0) {
    
    formula_skat <- reformulate(
      termlabels = covariates_keep,
      response = "obesity"
    )
    
  } else {
    
    formula_skat <- obesity ~ 1
  }
  
  null_model <- tryCatch(
    SKAT_Null_Model(
      formula_skat,
      data = df_clean,
      out_type = "D"
    ),
    error = function(e) {
      
      warning(
        "SKAT null model could not be created for ",
        result_name,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  if (is.null(null_model)) {
    return(NULL)
  }
  
  #---------------------------------------------------------------------------
  # SKAT and SKAT-O
  #---------------------------------------------------------------------------
  
  skat_p <- tryCatch(
    SKAT(
      geno_mat,
      null_model,
      method = "SKAT"
    )$p.value,
    error = function(e) {
      
      warning(
        "SKAT failed for ",
        result_name,
        ": ",
        conditionMessage(e)
      )
      
      NA_real_
    }
  )
  
  skato_p <- tryCatch(
    SKAT(
      geno_mat,
      null_model,
      method = "SKATO"
    )$p.value,
    error = function(e) {
      
      warning(
        "SKAT-O failed for ",
        result_name,
        ": ",
        conditionMessage(e)
      )
      
      NA_real_
    }
  )
  
  #---------------------------------------------------------------------------
  # Results table
  #---------------------------------------------------------------------------
  
  skat_results <- data.frame(
    Pathway = pathway_name,
    Analysis = result_name,
    Analysis_Config = analysis_config,
    Analysis_Label = config$label,
    Analysis_Type = config$type,
    Source_Score_Mode = source_score_mode,
    Score_Mode = "SKAT",
    Model = "Multivariate_SKAT",
    Covariates = if (
      length(covariates_keep) > 0
    ) {
      paste(
        covariates_keep,
        collapse = ", "
      )
    } else {
      "Intercept only"
    },
    Num_Individuals = nrow(geno_mat),
    Num_SNPs = ncol(geno_mat),
    SKAT_pvalue = skat_p,
    SKATO_pvalue = skato_p,
    Significant_Any = (
      !is.na(skat_p) &&
        skat_p < 0.05
    ) || (
      !is.na(skato_p) &&
        skato_p < 0.05
    ),
    stringsAsFactors = FALSE
  )
  
  write.csv2(
    skat_results,
    skat_file,
    row.names = FALSE
  )
  
  cat("SNPs:", ncol(geno_mat), "\n")
  cat("Individuals:", nrow(geno_mat), "\n")
  cat("SKAT p:", signif(skat_p, 4), "\n")
  cat("SKAT-O p:", signif(skato_p, 4), "\n")
  cat("Saved:", skat_file, "\n")
  
  return(
    skat_results
  )
}

#===================================================================================
# EXECUTION – EXACTLY ONE SKAT ANALYSIS PER SNP SELECTION
#===================================================================================

all_multivar_skat_results <- list()

# For the available SNP selection configurations, the
# raw sum score dataset is used as source for df_final and covariates.
# The sum score itself does not enter SKAT.
skat_result_names <- names(all_multivar_results)[
  grepl("_raw$", names(all_multivar_results)) |
    names(all_multivar_results) ==
    "sensitivity_lof_vep_lof_burden"
]

for (result_name in skat_result_names) {
  
  x <- all_multivar_results[[result_name]]
  
  if (is.null(x)) {
    next
  }
  
  analysis_config <- unique(
    as.character(x$df_final$Analysis_Config)
  )
  
  analysis_config <- analysis_config[
    !is.na(analysis_config)
  ]
  
  if (length(analysis_config) != 1) {
    warning(
      "Analysis_Config is not unique for ",
      result_name,
      ""
    )
    next
  }
  
  config <- analysis_configs[[analysis_config[1]]]
  
  # Pure score weightings do not constitute a separate SKAT analysis
  if (identical(config$type, "beta_weighted")) {
    next
  }
  
  all_multivar_skat_results[[analysis_config[1]]] <- run_multivariate_skat(
    df_final = x$df_final,
    covariates_keep = x$covariates_keep,
    result_name = analysis_config[1]
  )
}

#===================================================================================
# ZUSAMMENFASSUNG
#===================================================================================

multivar_skat_summary <- bind_rows(
  all_multivar_skat_results
)

if (nrow(multivar_skat_summary) > 0) {
  
  write.csv2(
    multivar_skat_summary,
    output_file(
      "Pathway_",
      pathway_name,
      "_Multivariate_SKAT_Results_All_Configurations.csv"
    ),
    row.names = FALSE
  )
  
} else {
  
  warning(
    "No multivariable SKAT results produced."
  )
}

warnings()

#===================================================================================
# MASTER TABLES – CENTRAL RESULTS OVERVIEW BASED ON MULTIVARIABLE ANALYSES
#===================================================================================
library(dplyr)
library(data.table)
library(stringr)

#-----------------------------
# 1) Multivariable SNP sum score results
#-----------------------------

master_multivar <- bind_rows(
  lapply(
    names(all_multivar_results),
    function(result_name) {
      
      x <- all_multivar_results[[result_name]]
      
      if (
        is.null(x) ||
        is.null(x$cutoff)
      ) {
        return(NULL)
      }
      
      x$cutoff %>%
        filter(
          Cutoff_Type == "Youden"
        ) %>%
        mutate(
          Analysis = result_name,
          Model_Type = "Multivariate_SNP_SumScore",
          Multivariable = "Yes",
          Batch_Adjusted = "No"
        )
    }
  )
)

if (nrow(master_multivar) == 0) {
  stop(
    "No multivariable SNP sum score results ",
    "were found for the master table."
  )
}

if (!"Score_Type" %in% colnames(master_multivar)) {
  master_multivar$Score_Type <- NA_character_
}

if (!"Score_Unit" %in% colnames(master_multivar)) {
  master_multivar$Score_Unit <- NA_character_
}

# Add N_Total robustly
if ("N_total" %in% colnames(master_multivar)) {
  
  master_multivar <- master_multivar %>%
    mutate(
      N_Total = N_total
    )
  
} else {
  
  master_multivar <- master_multivar %>%
    mutate(
      N_Total = N_train + N_test
    )
}

#-----------------------------
# 2) Add score and individual information
#-----------------------------

score_info_table <- bind_rows(
  lapply(
    names(all_multivar_results),
    function(result_name) {
      
      x <- all_multivar_results[[result_name]]
      
      if (
        is.null(x) ||
        is.null(x$df_final)
      ) {
        return(NULL)
      }
      
      score_values <- as.numeric(
        x$df_final$SNP_SumScore
      )
      
      score_values_finite <- score_values[
        is.finite(score_values)
      ]
      
      data.frame(
        Analysis = result_name,
        N_Individuals_Final = nrow(x$df_final),
        Score_Min = if (
          length(score_values_finite) > 0
        ) {
          min(score_values_finite)
        } else {
          NA_real_
        },
        Score_Max = if (
          length(score_values_finite) > 0
        ) {
          max(score_values_finite)
        } else {
          NA_real_
        },
        Score_Mean = if (
          length(score_values_finite) > 0
        ) {
          mean(score_values_finite)
        } else {
          NA_real_
        },
        Score_SD = if (
          length(score_values_finite) > 1
        ) {
          sd(score_values_finite)
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }
  )
)

master_multivar <- master_multivar %>%
  left_join(
    score_info_table,
    by = "Analysis"
  )

#-----------------------------
# 3) Add multivariable SKAT/SKAT-O results
#-----------------------------

master_multivar_skat_raw <- bind_rows(
  all_multivar_skat_results
)

if (nrow(master_multivar_skat_raw) > 0) {
  
  master_multivar_skat <- master_multivar_skat_raw %>%
    transmute(
      Analysis_Config,
      SKAT_Source_Analysis = Analysis,
      SKAT_Analysis_Label = Analysis_Label,
      SKAT_Analysis_Type = Analysis_Type,
      SKAT_Source_Score_Mode = Source_Score_Mode,
      SKAT_Covariates = Covariates,
      SKAT_Num_Individuals = Num_Individuals,
      SKAT_Num_SNPs = Num_SNPs,
      SKAT_pvalue,
      SKATO_pvalue,
      SKAT_Significant_Any = Significant_Any
    ) %>%
    distinct(
      Analysis_Config,
      .keep_all = TRUE
    )
  
} else {
  
  master_multivar_skat <- data.frame(
    Analysis_Config = character(),
    SKAT_Source_Analysis = character(),
    SKAT_Analysis_Label = character(),
    SKAT_Analysis_Type = character(),
    SKAT_Source_Score_Mode = character(),
    SKAT_Covariates = character(),
    SKAT_Num_Individuals = integer(),
    SKAT_Num_SNPs = integer(),
    SKAT_pvalue = numeric(),
    SKATO_pvalue = numeric(),
    SKAT_Significant_Any = logical(),
    stringsAsFactors = FALSE
  )
}

#-----------------------------
# 4) Add multivariable batch sensitivity
#-----------------------------

master_batch <- bind_rows(
  lapply(
    names(all_multivar_results),
    function(result_name) {
      
      x <- all_multivar_results[[result_name]]
      
      if (
        is.null(x) ||
        is.null(x$batch)
      ) {
        return(NULL)
      }
      
      x$batch %>%
        mutate(
          Analysis = result_name,
          Batch_Sensitivity_Performed = "Yes"
        )
    }
  )
)

if (nrow(master_batch) > 0) {
  
  master_batch <- master_batch %>%
    transmute(
      Analysis,
      Batch_Sensitivity_Performed,
      OR_no_batch,
      CI_low_no_batch,
      CI_high_no_batch,
      P_no_batch,
      OR_with_batch,
      CI_low_with_batch,
      CI_high_with_batch,
      P_with_batch,
      Delta_OR_Batch = Delta_OR,
      Delta_percent_Batch = Delta_percent
    ) %>%
    distinct(
      Analysis,
      .keep_all = TRUE
    )
  
} else {
  
  master_batch <- data.frame(
    Analysis = character(),
    Batch_Sensitivity_Performed = character(),
    OR_no_batch = numeric(),
    CI_low_no_batch = numeric(),
    CI_high_no_batch = numeric(),
    P_no_batch = numeric(),
    OR_with_batch = numeric(),
    CI_low_with_batch = numeric(),
    CI_high_with_batch = numeric(),
    P_with_batch = numeric(),
    Delta_OR_Batch = numeric(),
    Delta_percent_Batch = numeric(),
    stringsAsFactors = FALSE
  )
}

#-----------------------------
# 5) Add univariate results for comparison
#-----------------------------

sumscore_file <- output_file(
  "Pathway_",
  pathway_name,
  "_Univariate_SNP_SumScore_Model_Results_All_Configurations.csv"
)

if (file.exists(sumscore_file)) {
  
  univar_raw <- read.csv2(
    sumscore_file,
    stringsAsFactors = FALSE
  )
  
  required_univar_columns <- c(
    "Analysis",
    "Cutoff_Type",
    "OR_per_SD",
    "CI_lower_per_SD",
    "CI_upper_per_SD",
    "Pvalue_per_SD",
    "LRT_pvalue",
    "AUC",
    "Sensitivity",
    "Specificity",
    "Accuracy"
  )
  
  missing_univar_columns <- setdiff(
    required_univar_columns,
    colnames(univar_raw)
  )
  
  if (length(missing_univar_columns) > 0) {
    stop(
      "The following columns are missing in the univariate results table: ",
      paste(
        missing_univar_columns,
        collapse = ", "
      )
    )
  }
  
  univar_compare <- univar_raw %>%
    filter(
      Cutoff_Type == "Youden_Training"
    ) %>%
    transmute(
      Analysis,
      Univ_OR = OR_per_SD,
      Univ_CI_lower = CI_lower_per_SD,
      Univ_CI_upper = CI_upper_per_SD,
      Univ_Pvalue = Pvalue_per_SD,
      Univ_LRT_pvalue = LRT_pvalue,
      Univ_AUC = AUC,
      Univ_Sensitivity = Sensitivity,
      Univ_Specificity = Specificity,
      Univ_Accuracy = Accuracy
    ) %>%
    distinct(
      Analysis,
      .keep_all = TRUE
    )
  
} else {
  
  warning(
    "Univariate sum score overall results not found: ",
    sumscore_file
  )
  
  univar_compare <- data.frame(
    Analysis = character(),
    Univ_OR = numeric(),
    Univ_CI_lower = numeric(),
    Univ_CI_upper = numeric(),
    Univ_Pvalue = numeric(),
    Univ_LRT_pvalue = numeric(),
    Univ_AUC = numeric(),
    Univ_Sensitivity = numeric(),
    Univ_Specificity = numeric(),
    Univ_Accuracy = numeric(),
    stringsAsFactors = FALSE
  )
}

#-----------------------------
# 6) Merge master table
#-----------------------------

master_table <- master_multivar %>%
  left_join(
    master_multivar_skat,
    by = "Analysis_Config"
  ) %>%
  left_join(
    master_batch,
    by = "Analysis"
  ) %>%
  left_join(
    univar_compare,
    by = "Analysis"
  ) %>%
  mutate(
    Batch_Sensitivity_Performed = ifelse(
      is.na(Batch_Sensitivity_Performed),
      "No",
      Batch_Sensitivity_Performed
    ),
    
    SKAT_Performed = ifelse(
      is.na(SKAT_pvalue) &
        is.na(SKATO_pvalue),
      "No",
      "Yes"
    ),
    
    Score_Description = case_when(
      Score_Mode == "oriented" ~
        "Protective alleles reoriented; unweighted score",
      
      Score_Mode == "raw" ~
        "Original allele dosage; unweighted score",
      
      Score_Mode == "beta_weighted" ~
        "Allele dosage weighted by training-set log(OR)",
      
      Score_Mode == "ridge_weighted" ~
        "Allele dosage weighted by ridge-regularized (cv.glmnet, lambda.1se) training-set coefficients",
      
      Score_Mode == "elasticnet_weighted" ~
        "Allele dosage weighted by elastic-net-regularized (alpha=0.5, cv.glmnet, lambda.1se) training-set coefficients",
      
      Score_Mode == "lof_burden" ~
        paste0(
          "VEP-annotated loss-of-function variants; ",
          "dosage oriented to the LoF allele; unweighted burden"
        ),
      
      TRUE ~
        "Genetic score"
    )
    
  ) %>%
  dplyr::select(
    Pathway,
    Analysis,
    Analysis_Config,
    Score_Mode,
    Score_Description,
    Score_Type,
    Score_Unit,
    Model_Type,
    Multivariable,
    Batch_Adjusted,
    
    N_Total,
    N_train,
    N_test,
    N_Individuals_Final,
    
    Score_Min,
    Score_Max,
    Score_Mean,
    Score_SD,
    
    Estimate,
    StdError,
    OR,
    CI_lower,
    CI_upper,
    Pvalue,
    LRT_pvalue,
    
    AUC,
    AUC_CI_lower,
    AUC_CI_median,
    AUC_CI_upper,
    
    Sensitivity,
    Specificity,
    PPV,
    NPV,
    Accuracy,
    
    Cutoff_Type,
    Cutoff_Prob,
    SNP_SumScore_Cutoff,
    
    SKAT_Performed,
    SKAT_Source_Analysis,
    SKAT_Source_Score_Mode,
    SKAT_Analysis_Type,
    SKAT_Covariates,
    SKAT_Num_Individuals,
    SKAT_Num_SNPs,
    SKAT_pvalue,
    SKATO_pvalue,
    SKAT_Significant_Any,
    
    Batch_Sensitivity_Performed,
    OR_no_batch,
    CI_low_no_batch,
    CI_high_no_batch,
    P_no_batch,
    OR_with_batch,
    CI_low_with_batch,
    CI_high_with_batch,
    P_with_batch,
    Delta_OR_Batch,
    Delta_percent_Batch,
    
    Univ_OR,
    Univ_CI_lower,
    Univ_CI_upper,
    Univ_Pvalue,
    Univ_LRT_pvalue,
    Univ_AUC,
    Univ_Sensitivity,
    Univ_Specificity,
    Univ_Accuracy
  )

#-----------------------------
# 7) Export master table
#-----------------------------

write.csv2(
  master_table,
  output_file(
    "Pathway_",
    pathway_name,
    "_Multivariate_Master_Table_Results_All_Configurations.csv"
  ),
  row.names = FALSE
)

master_table %>%
  filter(
    Analysis_Config == "main_ld"
  ) %>%
  select(
    Analysis,
    Analysis_Config,
    Score_Mode,
    Cutoff_Type
  )

#===================================================================================
# SENSITIVITY ANALYSIS – CENTRAL OVERVIEW
#===================================================================================

main_reference <- master_table %>%
  filter(
    Analysis_Config == "main_ld",
    Score_Mode == "oriented",
    Cutoff_Type == "Youden"
  )

if (nrow(main_reference) == 0) {
  stop(
    "The reference model main_ld / oriented / Youden was not found."
  )
}

if (nrow(main_reference) > 1) {
  stop(
    paste0(
      "The reference model main_ld / oriented / Youden exists more than once: ",
      nrow(main_reference),
      " Zeilen."
    )
  )
}

# Store reference values as single numeric scalars
main_reference_or <- as.numeric(
  main_reference$OR[1]
)

main_reference_auc <- as.numeric(
  main_reference$AUC[1]
)

sensitivity_overview <- master_table %>%
  mutate(
    
    Is_Main_Reference =
      Analysis_Config == "main_ld" &
      Score_Mode == "oriented" &
      Cutoff_Type == "Youden",
    
    Sensitivity_Group = case_when(
      
      Is_Main_Reference ~
        "Main analysis",
      
      Analysis_Config == "main_ld" &
        Score_Mode == "raw" ~
        "Raw-score sensitivity",
      
      Analysis_Config == "sensitivity_betaweighted" &
        Score_Mode == "beta_weighted" ~
        "Beta-weighted score sensitivity",
      
      Analysis_Config == "sensitivity_betaweighted" &
        Score_Mode == "ridge_weighted" ~
        "Ridge-weighted score sensitivity",
      
      Analysis_Config == "sensitivity_betaweighted" &
        Score_Mode == "elasticnet_weighted" ~
        "Elastic-Net-weighted score sensitivity",
      
      Analysis_Config == "sensitivity_nold" ~
        "LD-pruning sensitivity",
      
      Analysis_Config == "sensitivity_functional_vep" ~
        "Functional SNP sensitivity",
      
      Analysis_Config == "sensitivity_lof_vep" ~
        "Loss-of-function burden sensitivity",
      
      Analysis_Config == "sensitivity_database_kegg" &
        Score_Mode == "oriented" ~
        "Database (KEGG-only) sensitivity",
      
      Analysis_Config == "sensitivity_database_kegg" &
        Score_Mode == "raw" ~
        "Database (KEGG-only) raw-score sensitivity",
      
      TRUE ~
        "Other sensitivity"
    ),
    
    Delta_OR_vs_Main = case_when(
      
      Is_Main_Reference ~
        0,
      
      Score_Mode == "lof_burden" ~
        NA_real_,
      
      is.na(OR) |
        is.na(main_reference_or) ~
        NA_real_,
      
      TRUE ~
        OR - main_reference_or
    ),
    
    Delta_AUC_vs_Main = case_when(
      
      Is_Main_Reference ~
        0,
      
      is.na(AUC) |
        is.na(main_reference_auc) ~
        NA_real_,
      
      TRUE ~
        AUC - main_reference_auc
    ),
    
    Relative_OR_Change_Percent = case_when(
      
      Is_Main_Reference ~
        0,
      
      Score_Mode == "lof_burden" ~
        NA_real_,
      
      is.na(OR) |
        is.na(main_reference_or) |
        main_reference_or == 0 ~
        NA_real_,
      
      TRUE ~
        100 * (OR - main_reference_or) /
        main_reference_or
    ),
    
    Same_Effect_Direction = case_when(
      
      Is_Main_Reference ~
        TRUE,
      
      is.na(OR) |
        is.na(main_reference_or) ~
        NA,
      
      OR > 1 &
        main_reference_or > 1 ~
        TRUE,
      
      OR < 1 &
        main_reference_or < 1 ~
        TRUE,
      
      OR == 1 &
        main_reference_or == 1 ~
        TRUE,
      
      TRUE ~
        FALSE
    ),
    
    Robustness_Interpretation = case_when(
      
      Is_Main_Reference ~
        "Reference model",
      
      Score_Mode == "lof_burden" &
        is.na(Same_Effect_Direction) ~
        "LoF-specific analysis; effect direction unclear",
      
      Score_Mode == "lof_burden" &
        Same_Effect_Direction ~
        paste0(
          "LoF-specific analysis; direction consistent, ",
          "raw OR magnitude not directly comparable"
        ),
      
      Score_Mode == "lof_burden" &
        !Same_Effect_Direction ~
        paste0(
          "LoF-specific analysis; direction differs, ",
          "raw OR magnitude not directly comparable"
        ),
      
      is.na(Same_Effect_Direction) |
        is.na(Delta_AUC_vs_Main) ~
        "Unclear",
      
      Same_Effect_Direction &
        abs(Delta_AUC_vs_Main) < 0.03 ~
        "Robust",
      
      Same_Effect_Direction &
        abs(Delta_AUC_vs_Main) >= 0.03 ~
        "Directionally consistent, AUC differs",
      
      !Same_Effect_Direction ~
        "Direction differs from main model",
      
      TRUE ~
        "Unclear"
    )
  )

write.csv2(
  sensitivity_overview,
  output_file(
    "Pathway_",
    pathway_name,
    "_Multivariate_Sensitivity_Analysis_Overview_All_Configurations.csv"
  ),
  row.names = FALSE
)

#-----------------------------
# Publication table
#-----------------------------

sensitivity_publication_table <- sensitivity_overview %>%
  dplyr::select(
    Analysis,
    Analysis_Config,
    Sensitivity_Group,
    Score_Mode,
    Cutoff_Type,
    Score_Type,
    Score_Unit,
    Score_Description,
    N_Total,
    N_train,
    N_test,
    SKAT_Num_SNPs,
    OR,
    CI_lower,
    CI_upper,
    Pvalue,
    LRT_pvalue,
    AUC,
    AUC_CI_lower,
    AUC_CI_upper,
    Delta_OR_vs_Main,
    Relative_OR_Change_Percent,
    Delta_AUC_vs_Main,
    Same_Effect_Direction,
    SKAT_Performed,
    SKAT_pvalue,
    SKATO_pvalue,
    Batch_Sensitivity_Performed,
    Delta_percent_Batch,
    Univ_OR,
    Univ_AUC,
    Robustness_Interpretation
  )

write.csv2(
  sensitivity_publication_table,
  output_file(
    "Pathway_",
    pathway_name,
    "_Multivariate_Sensitivity_Analysis_Publication_Table_All_Configurations.csv"
  ),
  row.names = FALSE
)

#===================================================================================
# DESKRIPTIVSTATISTIK – TABLE 1
#===================================================================================

library(dplyr)
library(data.table)

cat("\n====================================================\n")
cat("DESKRIPTIVSTATISTIK\n")
cat("====================================================\n")

#--------------------------------------------------
# Load main dataset
#--------------------------------------------------

main_df <- all_snp_sumscore_results[["main_ld_oriented"]]

if (is.null(main_df)) {
  stop(
    "main_ld_oriented wurde in all_snp_sumscore_results ",
    "not found."
  )
}

main_df <- as.data.frame(
  main_df,
  stringsAsFactors = FALSE
)

required_main_columns <- c(
  "id",
  "obesity",
  "SNP_SumScore"
)

missing_main_columns <- setdiff(
  required_main_columns,
  colnames(main_df)
)

if (length(missing_main_columns) > 0) {
  stop(
    "The following variables are missing in the main dataset: ",
    paste(
      missing_main_columns,
      collapse = ", "
    )
  )
}

main_df <- main_df %>%
  filter(
    !is.na(obesity),
    !is.na(SNP_SumScore)
  )

#--------------------------------------------------
# Add train/test and outcome group
#--------------------------------------------------

global_split <- readRDS(
  output_file(
    "Pathway_",
    pathway_name,
    "_Global_TrainTest_Split.rds"
  )
)

global_split$train_ids <- as.character(
  global_split$train_ids
)

global_split$test_ids <- as.character(
  global_split$test_ids
)

main_df$id <- as.character(
  main_df$id
)

main_df$Dataset <- case_when(
  main_df$id %in% global_split$train_ids ~ "Train",
  main_df$id %in% global_split$test_ids ~ "Test",
  TRUE ~ "Unassigned"
)

main_df$Obesity_Group <- case_when(
  main_df$obesity == 1 ~ "Obese",
  main_df$obesity == 0 ~ "Non-obese",
  TRUE ~ NA_character_
)

#--------------------------------------------------
# Helper function: continuous variables
#--------------------------------------------------

summarise_continuous <- function(
    df,
    var,
    group_name) {
  
  if (!var %in% colnames(df)) {
    
    return(
      data.frame(
        Variable = var,
        Group = group_name,
        N = 0,
        Mean = NA_real_,
        SD = NA_real_,
        Median = NA_real_,
        IQR = NA_real_,
        Min = NA_real_,
        Max = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }
  
  x <- suppressWarnings(
    as.numeric(
      gsub(
        ",",
        ".",
        as.character(df[[var]])
      )
    )
  )
  
  x <- x[
    is.finite(x)
  ]
  
  if (length(x) == 0) {
    
    return(
      data.frame(
        Variable = var,
        Group = group_name,
        N = 0,
        Mean = NA_real_,
        SD = NA_real_,
        Median = NA_real_,
        IQR = NA_real_,
        Min = NA_real_,
        Max = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }
  
  data.frame(
    Variable = var,
    Group = group_name,
    N = length(x),
    Mean = round(
      mean(x),
      2
    ),
    SD = if (
      length(x) > 1
    ) {
      round(
        sd(x),
        2
      )
    } else {
      NA_real_
    },
    Median = round(
      median(x),
      2
    ),
    IQR = round(
      IQR(x),
      2
    ),
    Min = round(
      min(x),
      2
    ),
    Max = round(
      max(x),
      2
    ),
    stringsAsFactors = FALSE
  )
}

#--------------------------------------------------
# Continuous variables
#--------------------------------------------------

continuous_vars <- intersect(
  c(
    "BMI",
    "sample_age",
    "SNP_SumScore"
  ),
  colnames(main_df)
)

continuous_groups <- list(
  Overall = main_df,
  Obese = main_df %>%
    filter(
      obesity == 1
    ),
  `Non-obese` = main_df %>%
    filter(
      obesity == 0
    ),
  Train = main_df %>%
    filter(
      Dataset == "Train"
    ),
  Test = main_df %>%
    filter(
      Dataset == "Test"
    )
)

desc_cont <- bind_rows(
  lapply(
    names(continuous_groups),
    function(group_name) {
      
      bind_rows(
        lapply(
          continuous_vars,
          function(v) {
            
            summarise_continuous(
              df = continuous_groups[[group_name]],
              var = v,
              group_name = group_name
            )
          }
        )
      )
    }
  )
)

#--------------------------------------------------
# Categorical variables
#--------------------------------------------------

categorical_vars <- intersect(
  c(
    "gender",
    "gender_cat",
    "Batch",
    "Dataset",
    "Obesity_Group"
  ),
  colnames(main_df)
)

desc_cat <- bind_rows(
  lapply(
    categorical_vars,
    function(v) {
      
      tab <- table(
        main_df[[v]],
        useNA = "ifany"
      )
      
      data.frame(
        Variable = v,
        Category = names(tab),
        N = as.numeric(tab),
        Percent = round(
          100 *
            as.numeric(tab) /
            sum(tab),
          1
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

#--------------------------------------------------
# Export
#--------------------------------------------------

write.csv2(
  desc_cont,
  output_file(
    "Pathway_",
    pathway_name,
    "_Descriptive_Statistics_Continuous.csv"
  ),
  row.names = FALSE
)

write.csv2(
  desc_cat,
  output_file(
    "Pathway_",
    pathway_name,
    "_Descriptive_Statistics_Categorical.csv"
  ),
  row.names = FALSE
)

cat(
  "Deskriptivstatistik gespeichert:\n"
)

cat(
  output_file(
    "Pathway_",
    pathway_name,
    "_Descriptive_Statistics_Continuous.csv"
  ),
  "\n"
)

cat(
  output_file(
    "Pathway_",
    pathway_name,
    "_Descriptive_Statistics_Categorical.csv"
  ),
  "\n"
)

#==================================================================================
# METHODS/RESULTS-LOG: PHASE 2 COMPLETED
#==================================================================================
log_pipeline_run(
  md_path      = methods_log_path,
  pathway_name = pathway_name,
  phase        = "2 (komplett, inkl. VEP-Annotation)",
  vep_meta     = extract_vep_meta(vep_results_file),
  extra_notes  = "Analyse-Pipeline vollstaendig durchgelaufen."
)

