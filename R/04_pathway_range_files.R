# ==============================================================================
# 04_pathway_range_files.R
# Writes one PLINK range file per pathway (pathway_<Name>_gene_ranges.txt;
# GRCh37 gene coordinates) from Pathway_Gene_Selection_final.xlsx. The range
# files are read by 05_plink_extraction_qc.ps1 (SNP extraction) and by 11/13c.
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================
#
# Format: tab-separated, no header, columns CHR | START | END | GENE; chromosome
# without "chr" prefix (X/MT as text), sorted by gene symbol.
# File name: "pathway_<Pathway with '/' -> '-' and ' ' -> '_'>_gene_ranges.txt".
#
# Genes whose coordinates lie only on alternative contigs/patches of GRCh37
# (no position on the primary assembly, e.g. AQP7B on chr2_gl877871_fix) are
# excluded and logged.
#
# Existing range files are not overwritten unless OVERWRITE_EXISTING <- TRUE.
# ==============================================================================

library(readxl)
library(dplyr)
library(stringr)

base_dir          <- Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus")
working_plink_dir <- file.path(base_dir, "Working_Plink")
excel_file        <- file.path(base_dir, "Pathway_Selection", "Pathway_Gene_Selection_final.xlsx")

OVERWRITE_EXISTING <- FALSE

#======================================================================================
# READ DATA
#======================================================================================
gene_table <- read_excel(excel_file, sheet = 1, col_types = "text")

stopifnot(all(c("Pathway", "Gene", "Chromosome region") %in% names(gene_table)))

cat("Read:", nrow(gene_table), "rows,",
    length(unique(gene_table$Pathway)), "pathways,",
    length(unique(gene_table$Gene)), "unique genes.\n")

#======================================================================================
# STEP 1: split "Chromosome region" (e.g. "chr11:107992243-108018503") into
# CHR / START / END
#======================================================================================
parsed <- gene_table %>%
  mutate(
    Chr   = str_match(`Chromosome region`, "^chr([^:]+):")[, 2],
    Start = str_match(`Chromosome region`, ":(\\d+)-")[, 2],
    End   = str_match(`Chromosome region`, "-(\\d+)$")[, 2]
  )

parse_failed <- parsed %>%
  filter(is.na(Chr) | is.na(Start) | is.na(End) | is.na(`Chromosome region`) | `Chromosome region` == "")

if (nrow(parse_failed) > 0) {
  stop(
    nrow(parse_failed), " row(s) with unparseable 'Chromosome region' - ",
    "no range files written:\n",
    paste(capture.output(print(parse_failed %>% select(Pathway, Gene, `Chromosome region`))),
          collapse = "\n")
  )
}

#======================================================================================
# STEP 2: exclude genes without a position on the primary assembly
#======================================================================================
primary_chromosomes <- c(as.character(1:22), "X", "Y", "MT")
non_primary <- parsed %>% filter(!Chr %in% primary_chromosomes)

if (nrow(non_primary) > 0) {
  cat("Excluded (no position on the GRCh37 primary assembly):\n")
  print(as.data.frame(non_primary %>% select(Pathway, Gene, `Chromosome region`)))
}

parsed <- parsed %>% filter(Chr %in% primary_chromosomes)

#======================================================================================
# STEP 3: check and remove duplicate pathway-gene rows
#======================================================================================
dupl_groups <- parsed %>%
  group_by(Pathway, Gene) %>%
  filter(n() > 1) %>%
  ungroup()

if (nrow(dupl_groups) > 0) {
  inconsistent <- dupl_groups %>%
    group_by(Pathway, Gene) %>%
    filter(n_distinct(Chr, Start, End) > 1) %>%
    ungroup()

  if (nrow(inconsistent) > 0) {
    stop(
      "Duplicate pathway-gene rows with CONFLICTING coordinates in ",
      basename(excel_file), ":\n",
      paste(capture.output(print(inconsistent %>%
              select(Pathway, Gene, `Chromosome region`))), collapse = "\n")
    )
  }

  cat("Note: ", nrow(distinct(dupl_groups, Pathway, Gene)),
      " duplicate pathway-gene combination(s) with identical coordinates removed.\n", sep = "")
}

parsed <- parsed %>% distinct(Pathway, Gene, .keep_all = TRUE)

#======================================================================================
# STEP 4: one range file per pathway
#======================================================================================
pathway_slug <- function(pathway_name) {
  pathway_name %>%
    str_replace_all(fixed("/"), "-") %>%
    str_replace_all(" ", "_")
}

parsed <- parsed %>% mutate(Slug = pathway_slug(Pathway))

log_rows <- list()

for (pw in sort(unique(parsed$Pathway))) {
  df <- parsed %>%
    filter(Pathway == pw) %>%
    arrange(Gene) %>%
    select(Chr, Start, End, Gene)

  slug <- unique(parsed$Slug[parsed$Pathway == pw])
  out_file <- file.path(working_plink_dir, paste0("pathway_", slug, "_gene_ranges.txt"))

  if (file.exists(out_file) && !OVERWRITE_EXISTING) {
    status <- "skipped (exists)"
  } else {
    write.table(df, out_file, sep = "\t", row.names = FALSE, col.names = FALSE,
                quote = FALSE, fileEncoding = "UTF-8")
    status <- "written"
  }

  cat(status, ": ", basename(out_file), " (", nrow(df), " genes)\n", sep = "")
  log_rows[[length(log_rows) + 1]] <- data.frame(
    Pathway = pw, File = basename(out_file), Gene_Count = nrow(df),
    Status = status, stringsAsFactors = FALSE
  )
}

log_df <- do.call(rbind, log_rows)
log_file <- file.path(working_plink_dir, "pathway_range_files_log.csv")
write.csv2(log_df, log_file, row.names = FALSE)

cat("\nDone. ", nrow(log_df), " pathway(s) processed. Log: ", log_file, "\n", sep = "")
