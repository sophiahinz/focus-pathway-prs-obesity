# ==============================================================================
# 13a_external_orientation_prepare.R
# Externally oriented scores, step 1: genome-wide SNP-orientation lookup from
# published BMI GWAS summary statistics (Yengo et al. 2018, primary; Locke et
# al. 2015, secondary cross-check).
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================

#=================================================================================
# Downloads are cached locally; only a small .rds lookup restricted to SNPs
# used in the pipeline (pathway LD-pruned sets + genome-wide permutation pool)
# is kept as input for 13b/13c.
#
# Optional: a run-level summary (software versions, file-access dates, result
# counts) is written to Appendix A of the methods/results log of this
# sensitivity analysis, if that file exists.
#=================================================================================

library(data.table)
library(dplyr)

# ---------------------------------------------------------------------------
# 0) Paths
# ---------------------------------------------------------------------------

base_dir       <- Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus")
external_dir   <- file.path(base_dir, "Sensitivity_External_SNP_Orientation")
raw_dir        <- file.path(external_dir, "Raw_Downloads")
data_dir       <- file.path(external_dir, "External_SNP_Data")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# GIANT consortium data files page:
# https://giant-consortium.web.broadinstitute.org/GIANT_consortium_data_files
# URLs as listed on that page at the time of analysis (accessed 2026-09-04).

yengo_url  <- "https://giant-consortium.web.broadinstitute.org/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz"
yengo_file <- file.path(raw_dir, "Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz")

locke_url  <- "https://giant-consortium.web.broadinstitute.org/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz"
locke_file <- file.path(raw_dir, "SNP_gwas_mc_merge_nogc.tbl.uniq.gz")

# Methods/Results-Log for this sensitivity analysis (auto-logging target, Section 9)
methods_log_path <- file.path(external_dir, "Methods", "Methods_Results_Log_External_SNP_Orientation.md")

# ---------------------------------------------------------------------------
# 0b) Coverage/minimum-count thresholds: 50% coverage AND >=10 SNPs, both required
# ---------------------------------------------------------------------------

MIN_COVERAGE_PCT <- 0.50
MIN_SNP_COUNT    <- 10

pathway_passes_coverage_threshold <- function(n_matched, n_total_main_model) {
  # Returns TRUE only if BOTH the percentage AND absolute-count floor are met.
  # A pathway failing either is marked "n/a" for the externally-oriented
  # sensitivity layer, consistent with the LoF-burden/elastic-net
  # "variant not computable for this pathway" convention.
  (n_matched / n_total_main_model) >= MIN_COVERAGE_PCT & n_matched >= MIN_SNP_COUNT
}

# ---------------------------------------------------------------------------
# 1) Download (once) - both files, primary + secondary cross-check source
# ---------------------------------------------------------------------------

download_if_missing <- function(url, dest) {
  if (!file.exists(dest)) {
    cat("Downloading:", url, "\n")
    download.file(url, destfile = dest, mode = "wb")
  } else {
    cat("Already present, skipping download:", dest, "\n")
  }
  # Log size + mtime for provenance (Appendix-style logging, matches existing
  # pipeline convention of recording what was actually used)
  fi <- file.info(dest)
  cat(sprintf("  size: %.1f MB | mtime: %s\n", fi$size / 1024^2, fi$mtime))
}

download_if_missing(yengo_url, yengo_file)
download_if_missing(locke_url, locke_file)

# ---------------------------------------------------------------------------
# 2) Read raw files - PRINT STRUCTURE FIRST, before assuming column names
#
# Decompression: uses base-R gzfile()/readBin() (bundled zlib) rather than
# fread()'s own gz-handling, which requires the R.utils package.
# ---------------------------------------------------------------------------
decompress_gz_base_r <- function(gz_path) {
  out_path <- sub("\\.gz$", "", gz_path)
  if (!file.exists(out_path)) {
    in_con <- gzfile(gz_path, "rb")
    out_con <- file(out_path, "wb")
    repeat {
      buf <- readBin(in_con, "raw", 1e7)
      if (length(buf) == 0) break
      writeBin(buf, out_con)
    }
    close(in_con); close(out_con)
  }
  out_path
}

yengo_txt <- decompress_gz_base_r(yengo_file)
locke_txt <- decompress_gz_base_r(locke_file)

yengo_raw <- fread(yengo_txt)
cat("\n=== Yengo 2018 file: column names as actually downloaded ===\n")
print(names(yengo_raw))
print(head(yengo_raw, 3))

locke_raw <- fread(locke_txt)
cat("\n=== Locke 2015 file: column names as actually downloaded ===\n")
print(names(locke_raw))
print(head(locke_raw, 3))

yengo_raw %>% filter(SNP %in% c("rs9939609", "rs17782313")) %>% select(SNP, CHR, POS)

# --- Column layout of the downloaded files --------------------------------
#   Yengo 2018:  CHR  POS  SNP  Tested_Allele  Other_Allele
#                Freq_Tested_Allele_in_HRS  BETA  SE  P  N
#   Locke 2015:  SNP  A1  A2  Freq1.Hapmap  b  se  p  N
#
# Genome build: not stated explicitly by the GIANT download page. Both
# meta-analyses were imputed against 1000 Genomes Phase 1 (GRCh37/hg19) per
# their original publications, consistent with the GRCh37/hg19 FoCus
# genotypes. Spot-check:
#   yengo_raw %>% filter(SNP %in% c("rs9939609", "rs17782313")) %>%
#     select(SNP, CHR, POS)
#   Expected under GRCh37: rs9939609 (FTO) chr16:53,820,527;
#   rs17782313 (MC4R region) chr18:57,851,097. (GRCh38 coordinates differ:
#   chr16:53,786,615 / chr18:55,518,281 respectively.)
# -------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 3) Standardize to a common internal schema
# ---------------------------------------------------------------------------

yengo_std <- yengo_raw %>%
  transmute(
    SNP           = SNP,
    Effect_Allele = toupper(Tested_Allele),
    Other_Allele  = toupper(Other_Allele),
    EAF           = Freq_Tested_Allele_in_HRS,
    Beta          = BETA,          # already log-OR / linear-BMI-SD scale per GIANT format
    P              = P,
    N              = N,
    Source        = "Yengo2018"
  )

locke_std <- locke_raw %>%
  transmute(
    SNP           = SNP,
    Effect_Allele = toupper(A1),
    Other_Allele  = toupper(A2),
    EAF           = Freq1.Hapmap,
    Beta          = b,
    P              = p,
    N              = N,
    Source        = "Locke2015"
  )

# ---------------------------------------------------------------------------
# 4) Restrict to GENUINE rsIDs only (explicit, logged step)
#
# Only true dbSNP rsIDs (^rs[0-9]+$) are matched against the external
# files. The FoCus pathway-level .bim files carry a mix of true rsIDs and
# Illumina-internal "GSA-rs..." pseudo-IDs for array probes without a
# confirmed dbSNP rsID. GSA-rs pseudo-IDs
# are Illumina-internal and do not exist in GIANT/Yengo's rsID space, so an
# ID-based join would never match them anyway - this step makes that
# exclusion EXPLICIT and LOGGED per pathway, rather than leaving it as an
# implicit side effect of a failed join (which the QC-flow-count convention
# used elsewhere in this pipeline would otherwise miss).
# ---------------------------------------------------------------------------

is_genuine_rsid <- function(x) grepl("^rs[0-9]+$", x)

# Applied to the FoCus side at usage time (see create_external_orientation_table()
# below) - bim_table is filtered there, not here, since the external files'
# own SNP column is already 100% genuine rsIDs by construction (GIANT/dbSNP).

# ---------------------------------------------------------------------------
# 5) Combine primary + secondary source (Yengo primary; Locke fills gaps /
#    serves as a same-SNP cross-check where both cover a variant)
# ---------------------------------------------------------------------------

external_lookup <- bind_rows(yengo_std, locke_std) %>%
  filter(is_genuine_rsid(SNP)) %>%
  distinct(SNP, Source, .keep_all = TRUE)

# Cross-check table: SNPs covered by BOTH sources, for an internal-consistency
# spot-check (sign concordance) before relying on Yengo alone downstream.
external_crosscheck <- yengo_std %>%
  inner_join(locke_std, by = "SNP", suffix = c("_yengo", "_locke")) %>%
  mutate(Sign_Concordant = sign(Beta_yengo) == sign(
    ifelse(Effect_Allele_yengo == Effect_Allele_locke, Beta_locke, -Beta_locke)
  ))
cat(sprintf(
  "\nCross-check: %d SNPs covered by both Yengo 2018 and Locke 2015; %.1f%% sign-concordant.\n",
  nrow(external_crosscheck), 100 * mean(external_crosscheck$Sign_Concordant, na.rm = TRUE)
))

# ---------------------------------------------------------------------------
# 6) Restrict to the union of SNPs actually used anywhere in the pipeline
#    BEFORE caching (keeps the cached object small - same principle as
#    gene_sets_cache_final.rds).
#
# Why this matters: Yengo 2018 + Locke 2015 together cover on the order of
# ~2-2.5 million genuine-rsID SNPs genome-wide (steps 1-5 keep all of them).
# FoCus's pipeline only ever needs orientation for the SNPs that actually
# appear somewhere in the 54(+) pathway-level LD-pruned sets or the
# genome-wide permutation pool - everything else in external_lookup is dead
# weight that would otherwise be carried into every downstream script's
# loading time and the cached file's disk size for no benefit.
# ---------------------------------------------------------------------------

read_bim_snp_col <- function(f) {
  fread(f, header = FALSE, select = 2, col.names = "SNP")$SNP
}

# Pathway-level, final LD-pruned SNP sets written by 07_pathway_pipeline.R -
# one .bim per screened pathway (plus the KEGG-vs-Reactome database-comparison
# variants). The naming pattern below identifies the FINAL LD-pruned set per
# pathway, as opposed to the intermediate _after_R_QC.bim (pre-pruning) or
# per-batch _Ambu_QC.bim/_MD_QC.bim files.
pathway_bim_dir <- file.path(Sys.getenv("FOCUS_RESULTS_DIR", unset = "path/to/results"), "Working_Plink", "plink_pathway_outputs")
pathway_bim_files <- list.files(
  pathway_bim_dir,
  pattern = "_after_R_QC_LDpruned\\.bim$",
  full.names = TRUE
)
cat(sprintf("\nFound %d pathway-level LD-pruned .bim files in %s\n", length(pathway_bim_files), pathway_bim_dir))
if (length(pathway_bim_files) == 0) {
  warning("No pathway .bim files found - check pathway_bim_dir and the naming pattern before trusting Step 6's filtering.")
}

pathway_snps <- unique(unlist(lapply(pathway_bim_files, read_bim_snp_col)))

# Genome-wide permutation SNP pool - same file used as the random-draw
# source in 11_competitive_permutation_test.R.
permutation_pool_bim  <- file.path(base_dir, "Global_PCA", "FO_GSA_Merged_QC_Pruned.bim")
permutation_pool_snps <- read_bim_snp_col(permutation_pool_bim)

all_used_snps <- unique(c(pathway_snps, permutation_pool_snps))
cat(sprintf(
  "Union of SNPs actually used in the pipeline (pathways + permutation pool): %d\n",
  length(all_used_snps)
))

external_lookup_n_before <- nrow(external_lookup)
external_lookup <- external_lookup %>% filter(SNP %in% all_used_snps)
cat(sprintf(
  "external_lookup restricted from %d to %d rows (%.1f%% retained) before caching.\n",
  external_lookup_n_before, nrow(external_lookup),
  100 * nrow(external_lookup) / external_lookup_n_before
))

# ---------------------------------------------------------------------------
# 7) create_external_orientation_table() - extends create_lof_orientation_table()
#    pattern with palindromic-SNP flagging
# ---------------------------------------------------------------------------

create_external_orientation_table <- function(external_lookup, bim_table, matrix_snps) {

  bim_table_rsid_only <- bim_table %>%
    mutate(Is_Genuine_RsID = is_genuine_rsid(SNP))

  n_total       <- nrow(bim_table_rsid_only)
  n_non_rsid    <- sum(!bim_table_rsid_only$Is_Genuine_RsID)
  cat(sprintf(
    "  SNP-ID QC: %d/%d (%.1f%%) are non-rsID pseudo-IDs (e.g. GSA-rs...) and are excluded from external matching by construction.\n",
    n_non_rsid, n_total, 100 * n_non_rsid / n_total
  ))

  bim_table_rsid_only %>%
    filter(Is_Genuine_RsID) %>%
    left_join(external_lookup, by = "SNP") %>%
    filter(SNP %in% matrix_snps) %>%
    mutate(
      Matched = !is.na(Effect_Allele),

      Is_Palindromic = Matched & (
        (Effect_Allele == "A" & Other_Allele == "T") |
        (Effect_Allele == "T" & Other_Allele == "A") |
        (Effect_Allele == "C" & Other_Allele == "G") |
        (Effect_Allele == "G" & Other_Allele == "C")
      ),

      Allele_Match = case_when(
        !Matched ~ "No_external_match",
        Is_Palindromic ~ "Palindromic_excluded",
        Effect_Allele == Matrix_Allele2 ~ "Effect_is_matrix_counted_allele",
        Effect_Allele == Matrix_Allele1 ~ "Effect_is_matrix_other_allele",
        TRUE ~ "No_allele_match"
      ),

      Invert_For_External_Score = case_when(
        Allele_Match == "Effect_is_matrix_counted_allele" ~ FALSE,
        Allele_Match == "Effect_is_matrix_other_allele" ~ TRUE,
        TRUE ~ NA
      )
    )
}

# ---------------------------------------------------------------------------
# 8) Save compact, indexed lookup
# ---------------------------------------------------------------------------

output_rds <- file.path(data_dir, "External_Orientation_GIANT_Yengo2018_Locke2015.rds")
saveRDS(external_lookup, output_rds)
saveRDS(external_crosscheck, file.path(data_dir, "External_Orientation_CrossCheck_Yengo_vs_Locke.rds"))

cat("\nDone. Cached lookup written to Sensitivity_External_SNP_Orientation/External_SNP_Data/.\n")
cat("Next step: 13b_external_orientation_score.R.\n")

#==================================================================================
# 9) Methods/results log (optional): Appendix A -- upsert + session info
#    (same pattern as in 11_competitive_permutation_test.R)
#==================================================================================

log_external_orientation_to_md <- function(md_path, run_summary) {

  if (!file.exists(md_path)) {
    warning("Methods/results log not found: ", md_path, " - not written.")
    return(invisible(NULL))
  }

  log_lines <- readLines(md_path, warn = FALSE, encoding = "UTF-8")
  begin_marker <- which(trimws(log_lines) == "<!-- AUTOLOG-EXTORIENT-PREP:BEGIN -->")
  end_marker   <- which(trimws(log_lines) == "<!-- AUTOLOG-EXTORIENT-PREP:END -->")

  if (length(begin_marker) != 1 || length(end_marker) != 1 || end_marker <= begin_marker) {
    warning("Marker for Appendix A not found in ", md_path, " - not written.")
    return(invisible(NULL))
  }

  header_end <- begin_marker + 2
  data_start <- header_end + 1
  data_end   <- end_marker - 1
  existing_data_rows <- if (data_start <= data_end) log_lines[data_start:data_end] else character(0)

  safe_cell <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x) || x == "") return("")
    gsub("\\|", "/", gsub("[\r\n]+", " ", as.character(x)))
  }
  extract_rundate_from_row <- function(row) {
    cells <- strsplit(row, "\\|")[[1]]
    if (length(cells) < 2) return(NA_character_)
    trimws(cells[2])
  }

  run_date <- format(Sys.Date(), "%Y-%m-%d")

  new_row <- sprintf(
    "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
    run_date,
    safe_cell(run_summary$yengo_accessed), safe_cell(round(run_summary$yengo_size_mb, 1)),
    safe_cell(run_summary$locke_accessed), safe_cell(round(run_summary$locke_size_mb, 1)),
    safe_cell(run_summary$n_yengo_rsid), safe_cell(run_summary$n_locke_rsid),
    safe_cell(run_summary$n_crosscheck), safe_cell(round(run_summary$crosscheck_concordance_pct, 1)),
    safe_cell(run_summary$output_rds)
  )

  if (length(existing_data_rows) > 0) {
    existing_rundates <- vapply(existing_data_rows, extract_rundate_from_row, character(1))
    kept_existing_rows <- existing_data_rows[is.na(existing_rundates) | existing_rundates != run_date]
  } else {
    kept_existing_rows <- character(0)
  }

  final_data_rows <- c(kept_existing_rows, new_row)
  updated_lines <- c(log_lines[1:header_end], final_data_rows, log_lines[end_marker:length(log_lines)])
  writeLines(updated_lines, md_path, useBytes = TRUE)

  invisible(NULL)
}

log_external_orientation_sessioninfo_to_md <- function(md_path) {

  if (!file.exists(md_path)) return(invisible(NULL))

  log_lines <- readLines(md_path, warn = FALSE, encoding = "UTF-8")
  begin_marker <- which(trimws(log_lines) == "<!-- AUTOLOG-EXTORIENT-PREP-SESSIONINFO:BEGIN -->")
  end_marker   <- which(trimws(log_lines) == "<!-- AUTOLOG-EXTORIENT-PREP-SESSIONINFO:END -->")

  if (length(begin_marker) != 1 || length(end_marker) != 1 || end_marker <= begin_marker) {
    return(invisible(NULL))
  }

  pkg_version <- function(pkg) {
    tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "not installed")
  }

  info_lines <- c(
    "",
    paste0("- **Last run:** ", format(Sys.time(), "%Y-%m-%d %H:%M")),
    paste0("- **R version:** ", R.version.string),
    paste0("- **Platform:** ", R.version$platform, " (", Sys.info()[["sysname"]], ")"),
    paste0(
      "- **Package versions:** data.table ", pkg_version("data.table"),
      ", dplyr ", pkg_version("dplyr")
    ),
    paste0("- **Script:** 13a_external_orientation_prepare.R"),
    paste0("- **Reference files:** Yengo 2018 (", basename(yengo_file), "), Locke 2015 (",
           basename(locke_file), ")"),
    ""
  )

  updated_lines <- c(log_lines[1:begin_marker], info_lines, log_lines[end_marker:length(log_lines)])
  writeLines(updated_lines, md_path, useBytes = TRUE)
  invisible(NULL)
}

run_summary <- list(
  yengo_accessed             = format(file.info(yengo_file)$mtime, "%Y-%m-%d"),
  yengo_size_mb              = file.info(yengo_file)$size / 1024^2,
  locke_accessed             = format(file.info(locke_file)$mtime, "%Y-%m-%d"),
  locke_size_mb              = file.info(locke_file)$size / 1024^2,
  n_yengo_rsid                = sum(external_lookup$Source == "Yengo2018"),
  n_locke_rsid                = sum(external_lookup$Source == "Locke2015"),
  n_crosscheck                = nrow(external_crosscheck),
  crosscheck_concordance_pct  = 100 * mean(external_crosscheck$Sign_Concordant, na.rm = TRUE),
  output_rds                  = basename(output_rds)
)

log_external_orientation_to_md(methods_log_path, run_summary)
log_external_orientation_sessioninfo_to_md(methods_log_path)

cat("\nAppendix A (methods log) updated:", methods_log_path, "\n")
