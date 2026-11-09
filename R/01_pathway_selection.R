# ==============================================================================
# 01_pathway_selection.R
# Pathway selection (KEGG-first, two-zone Jaccard overlap model): KEGG base
# list from the BRITE hierarchy, gene-set retrieval and harmonization to HGNC
# symbols, KEGG-Reactome overlap, final gene sets and selection flow counts
# (Figure 1; Supplementary Tables S1, S2b).
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================
#
# The selection combines scripted steps (this file) with manual review steps
# that were documented in a tracking spreadsheet
# (Pathway_Selection_Tracking_Spreadsheet.xlsx, sheet "Tracking"; one row per
# candidate pathway). The spreadsheet holds the inclusion/exclusion decisions
# and is not part of this repository; the resulting gene sets are provided in
# gene_sets/. The workflow is:
#
#   Phase 1 (R)      KEGG base list: all human KEGG pathways in the four
#                    BRITE target categories (carbohydrate, energy, lipid,
#                    amino acid metabolism) -> copied into the spreadsheet.
#   Phase 2 (manual) Gap search across KEGG and review of Reactome candidates
#                    in the four corresponding top-level branches.
#   Phase 3 (manual) Content and hierarchy criteria (inclusion/exclusion
#                    criteria, Supplementary Table S1).
#   Phase 4 (R)      Gene sets per candidate and gene counts (formal and size
#                    criteria, Supplementary Table S1).
#   Phase 5 (R)      Jaccard overlap of each Reactome candidate with its
#                    reference KEGG pathway (two-zone model, J = 0.3).
#   Phase 6 (manual) Review round; decision per candidate in column
#                    "Entscheidung": "in", "in (als KEGG erweitert)" or "out".
#   Phase 7 (R)      Final gene sets (zone 1: independent pathways; zone 2:
#                    KEGG pathways extended by the non-overlapping genes of
#                    the high-overlap Reactome candidates), GMT/CSV export and
#                    flow counts.
#
# Column names and category values of the tracking spreadsheet are German
# (e.g. "Entscheidung", "Referenz-KEGG-Pathway", "Name (final)",
# "Basis-Kandidat") and are used as keys below.
#
# Database versions used: KEGG Release 119.0; Reactome release 97 (via
# ReactomeContentService4R); org.Hs.eg.db 3.22.0 (Bioconductor 3.22).
# ==============================================================================

library(readxl)
library(openxlsx)
library(dplyr)

get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) return(dirname(normalizePath(file_arg)))
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(dirname(ofile))
  "."
}
source(file.path(get_script_dir(), "01_functions_pathway_selection.R"))

export_dir <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Pathway_Selection")
if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

tracking_file  <- file.path(export_dir, "Pathway_Selection_Tracking_Spreadsheet.xlsx")
cache_file     <- file.path(export_dir, "gene_sets_cache_final.rds")
pathway_list   <- file.path(export_dir, "Pathway_Liste.csv")   # short pipeline names (gene_sets/)

read_tracking <- function() read_excel(tracking_file, sheet = "Tracking")

write_tracking <- function(tracking_table) {
  wb <- loadWorkbook(tracking_file)
  writeData(wb, sheet = "Tracking", x = tracking_table)
  saveWorkbook(wb, tracking_file, overwrite = TRUE)
}

if (file.exists(cache_file)) {
  gene_sets_cache <- readRDS(cache_file)
  message("Gene-set cache loaded: ", length(gene_sets_cache), " pathways.")
}

# ============================================================
# PHASE 1 - KEGG base list
# ============================================================
# BRITE hierarchy br08901: B lines = category, C lines = pathway map number
brite_lines <- readLines(url("https://rest.kegg.jp/get/br:br08901"))

target_categories <- c("Carbohydrate metabolism", "Energy metabolism",
                       "Lipid metabolism", "Amino acid metabolism")

current_category <- NA
all_b_headings <- character(0)
hits <- data.frame(map_id = character(0), category = character(0),
                   stringsAsFactors = FALSE)

for (line in brite_lines) {
  if (grepl("^B\\s+\\S", line)) {
    current_category <- trimws(sub("^B\\s+", "", line))
    all_b_headings <- c(all_b_headings, current_category)
  } else if (grepl("^C\\s+\\d{5}", line) && current_category %in% target_categories) {
    hits <- rbind(hits, data.frame(
      map_id   = regmatches(line, regexpr("\\d{5}", line)),
      category = current_category,
      stringsAsFactors = FALSE))
  }
}
hits$hsa_id <- paste0("hsa", hits$map_id)

category_counts <- table(factor(hits$category, levels = target_categories))
print(category_counts)
if (any(category_counts == 0)) {
  warning("0 hits in category/categories: ",
          paste(names(category_counts)[category_counts == 0], collapse = ", "))
}

multiple <- hits[duplicated(hits$hsa_id) | duplicated(hits$hsa_id, fromLast = TRUE), ]
if (nrow(multiple) > 0) {
  warning("Pathway ID(s) in more than one target category: ",
          paste(unique(multiple$hsa_id), collapse = ", "))
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
writeLines(brite_lines, file.path(export_dir, sprintf("br08901_raw_%s.txt", timestamp)))
writeLines(unique(all_b_headings),
           file.path(export_dir, sprintf("br08901_all_categories_%s.txt", timestamp)))
saveRDS(category_counts, file.path(export_dir, "category_counts_phase1.rds"))

all_hsa_pathways <- keggList("pathway", "hsa")
names(all_hsa_pathways) <- sub("^path:", "", names(all_hsa_pathways))

hits <- hits[!duplicated(hits$hsa_id), ]
hits <- hits[hits$hsa_id %in% names(all_hsa_pathways), ]

phase1_export <- data.frame(
  `Pathway-ID`      = hits$hsa_id,
  `Name (Original)` = unname(all_hsa_pathways[hits$hsa_id]),
  `Zielkategorie`   = hits$category,
  check.names = FALSE
)

write.csv2(phase1_export,
           file.path(export_dir, "kegg_base_list_phase1.csv"),
           row.names = FALSE)

# ------------------------------------------------------------
# Phases 2-3 (manual): gap search, Reactome review, content and
# hierarchy criteria -> tracking spreadsheet
# ------------------------------------------------------------

# ============================================================
# PHASE 4 - Gene sets and gene counts per candidate
# ============================================================
tracking_table <- read_tracking()
if (!"Gene-Count" %in% names(tracking_table)) tracking_table$`Gene-Count` <- NA_integer_

rows_to_fill <- which(is.na(tracking_table$`Gene-Count`))
for (i in rows_to_fill) {
  genes <- get_gene_set(tracking_table$`Pathway-ID`[i], tracking_table$DB[i])
  tracking_table$`Gene-Count`[i] <- length(genes)
  if (i %% 20 == 0) saveRDS(gene_sets_cache, file.path(export_dir, "gene_sets_cache_intermediate.rds"))
}

saveRDS(gene_sets_cache, cache_file)
write_tracking(tracking_table)

# ============================================================
# PHASE 5 - Two-zone overlap of Reactome candidates with KEGG
# ============================================================
tracking_table <- read_tracking()

if (!"Overlap zu KEGG-Basis" %in% names(tracking_table)) {
  tracking_table$`Overlap zu KEGG-Basis` <- NA_real_
}

reference_ids <- unique(unlist(lapply(
  na.omit(tracking_table$`Referenz-KEGG-Pathway`), split_kegg_ids)))
for (rid in reference_ids) {
  get_gene_set(rid, "KEGG")
}

reactome_rows <- which(tracking_table$DB == "Reactome" & tracking_table$Entscheidung != "out")

for (i in reactome_rows) {
  reactome_id <- tracking_table$`Pathway-ID`[i]
  kegg_id     <- tracking_table$`Referenz-KEGG-Pathway`[i]

  if (is.na(kegg_id) || kegg_id == "") {
    warning(sprintf("No reference KEGG pathway for %s - overlap cannot be computed.", reactome_id))
    next
  }

  gene_set_reactome <- get_gene_set(reactome_id, "Reactome")
  gene_set_kegg     <- get_gene_set_combined(kegg_id, "KEGG")

  jaccard <- length(intersect(gene_set_kegg, gene_set_reactome)) /
             length(union(gene_set_kegg, gene_set_reactome))

  tracking_table$`Overlap zu KEGG-Basis`[i] <- jaccard
}

saveRDS(gene_sets_cache, cache_file)
write_tracking(tracking_table)

# ------------------------------------------------------------
# Phase 6 (manual): review round, final decisions in "Entscheidung"
# ------------------------------------------------------------

# ============================================================
# PHASE 7 - Final gene sets and flow counts
# ============================================================
tracking_table <- read_tracking()
category_counts <- readRDS(file.path(export_dir, "category_counts_phase1.rds"))

final_tracking_table <- tracking_table %>%
  filter(Entscheidung %in% c("in", "in (als KEGG erweitert)"))

name_lookup <- setNames(
  ifelse(!is.na(tracking_table$`Name (final)`) & tracking_table$`Name (final)` != "",
         tracking_table$`Name (final)`,
         tracking_table$`Name (Original)`),
  tracking_table$`Pathway-ID`
)

# Zone 1: independent pathways (KEGG and low-overlap Reactome)
zone1_ids <- final_tracking_table$`Pathway-ID`[final_tracking_table$Entscheidung == "in"]
gene_sets_zone1 <- gene_sets_cache[zone1_ids]
zone1_names <- name_lookup[zone1_ids]
zone1_db <- setNames(
  final_tracking_table$DB[final_tracking_table$Entscheidung == "in"],
  zone1_ids
)

# Zone 2: KEGG pathways extended by the non-overlapping genes of
# high-overlap Reactome candidates
zone2_rows <- final_tracking_table %>%
  filter(Entscheidung == "in (als KEGG erweitert)")

gene_sets_zone2 <- list()
zone2_names <- character(0)
gene_origin_zone2 <- list()
for (kegg_id in unique(zone2_rows$`Referenz-KEGG-Pathway`)) {
  reactome_ids_here <- zone2_rows$`Pathway-ID`[zone2_rows$`Referenz-KEGG-Pathway` == kegg_id]

  gene_set_kegg_original <- get_gene_set_combined(kegg_id, "KEGG")
  all_reactome_genes     <- unique(unlist(gene_sets_cache[reactome_ids_here]))
  additional_genes       <- setdiff(all_reactome_genes, gene_set_kegg_original)

  list_key <- paste0(gsub(", ?", "_", kegg_id), "_erweitert")
  combined_set <- union(gene_set_kegg_original, additional_genes)
  gene_sets_zone2[[list_key]] <- combined_set

  gene_origin_zone2[[list_key]] <- setNames(
    ifelse(combined_set %in% gene_set_kegg_original, "KEGG", "Reactome"),
    combined_set
  )

  kegg_display_names <- name_lookup[split_kegg_ids(kegg_id)]
  zone2_names[list_key] <- paste0(paste(kegg_display_names, collapse = " / "), " (KEGG erweitert)")
}

gene_sets <- c(gene_sets_zone1, gene_sets_zone2)
pathway_names <- c(zone1_names, zone2_names)

dataset_type <- c(
  setNames(as.character(zone1_db), names(zone1_db)),
  setNames(rep("KEGG erweitert", length(gene_sets_zone2)), names(gene_sets_zone2))
)

saveRDS(gene_sets, file = file.path(export_dir, "gene_sets_final.rds"))

gene_sets_long <- do.call(rbind, lapply(names(gene_sets), function(pid) {
  genes <- gene_sets[[pid]]
  type  <- dataset_type[[pid]]
  origin <- if (type == "KEGG erweitert") gene_origin_zone2[[pid]][genes] else rep(type, length(genes))
  data.frame(`Pathway-ID` = pid,
             `Pathway-Name` = pathway_names[[pid]],
             Gen = genes,
             `Gen-Ursprung` = unname(origin),
             `Datensatz-Art` = type,
             check.names = FALSE,
             stringsAsFactors = FALSE)
}))
write.csv(gene_sets_long, file = file.path(export_dir, "Supplement_GeneSets.csv"), row.names = FALSE)

write_gmt <- function(gene_sets, file, descriptions = NULL) {
  con <- file(file, "w")
  for (pid in names(gene_sets)) {
    desc <- if (!is.null(descriptions)) descriptions[[pid]] else "NA"
    writeLines(paste(c(pid, desc, gene_sets[[pid]]), collapse = "\t"), con)
  }
  close(con)
}

write_gmt(gene_sets, file.path(export_dir, "Supplement_GeneSets.gmt"), descriptions = pathway_names)

# Pipeline input: same long table with the short pathway names used by all
# downstream scripts (gene_sets/Pathway_Liste.csv, matched by pathway ID;
# "_erweitert" sets share the short name of their KEGG pathway).
short_names <- read.csv(pathway_list, stringsAsFactors = FALSE)
short_lookup <- setNames(short_names$pathway_name, short_names$pathway_selection_id)
pipeline_gene_sets <- gene_sets_long
pipeline_gene_sets$`Pathway-Name` <- unname(short_lookup[sub("_erweitert$", "", pipeline_gene_sets$`Pathway-ID`)])
if (any(is.na(pipeline_gene_sets$`Pathway-Name`))) {
  stop("Pathway IDs without short name in Pathway_Liste.csv: ",
       paste(unique(pipeline_gene_sets$`Pathway-ID`[is.na(pipeline_gene_sets$`Pathway-Name`)]), collapse = ", "))
}
write.csv(pipeline_gene_sets,
          file = file.path(export_dir, "Pathway_Selection_Supplement_GeneSets.csv"),
          row.names = FALSE)

# --- Flow counts (Figure 1 / Supplementary Table S2b) ---
flow_counts <- data.frame(
  Stage = c(
    paste("KEGG base list:", names(category_counts)),
    "KEGG starting number (sum of target categories)",
    "KEGG additions via gap search (found)",
    "KEGG additions via gap search (final)",
    "Final KEGG pathways",
    "Reactome candidates screened",
    "Reactome candidates after inclusion/exclusion criteria",
    "Zone 1: independent Reactome pathways",
    "Zone 2: KEGG pathways extended by Reactome genes",
    "Total final pathway-level gene sets"
  ),
  N = c(
    as.integer(category_counts),
    sum(tracking_table$DB == "KEGG" & tracking_table$Status == "Basis-Kandidat"),
    sum(tracking_table$Status == "KEGG-Ergänzung (Lücke)"),
    sum(tracking_table$Status == "KEGG-Ergänzung (Lücke)" & tracking_table$Entscheidung == "in"),
    sum(tracking_table$DB == "KEGG" & tracking_table$Entscheidung == "in"),
    sum(tracking_table$DB == "Reactome"),
    sum(tracking_table$DB == "Reactome" & tracking_table$Entscheidung %in% c("in", "in (als KEGG erweitert)")),
    sum(tracking_table$DB == "Reactome" & tracking_table$Entscheidung == "in"),
    length(unique(unlist(lapply(
      tracking_table$`Referenz-KEGG-Pathway`[tracking_table$Entscheidung == "in (als KEGG erweitert)"],
      split_kegg_ids)))),
    sum(tracking_table$DB == "KEGG" & tracking_table$Entscheidung == "in") +
      sum(tracking_table$DB == "Reactome" & tracking_table$Entscheidung == "in")
  )
)
print(flow_counts, row.names = FALSE)
write.csv(flow_counts, file.path(export_dir, "Pathway_Selection_Flow_Counts.csv"), row.names = FALSE)

# --- Annotation versions ---
cat("Bioconductor:", as.character(BiocManager::version()), "\n")
print(org.Hs.eg_dbInfo())
sessionInfo()
