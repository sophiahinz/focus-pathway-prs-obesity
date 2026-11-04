# Pathway-restricted polygenic scores for energy-substrate metabolism and obesity (FoCus cohort)
# Script: 02_pathway_gene_overlap_jaccard.R | Source file: Pathway_Gene_Overlap_Jaccard_REPAIR.R
# Local paths replaced by the environment variable FOCUS_BASE_DIR (project root); analysis code unchanged.
# Individual-level FoCus data are not included (see README).

# ==============================================================================
# REPAIR VERSION (2026-10-09): copy of Pathway_Gene_Overlap_Jaccard.R
# Re-run after correction of the defective phenotype CSVs. Changes ONLY to paths
# (lines marked "# REPAIR"); analysis code unchanged. Writes only to FoCus/Repair.
# ==============================================================================

# ==============================================================================
# Pathway_Gene_Overlap_Jaccard.R  (v1.0)
#
# Purpose: checks how strongly the gene sets of the robust candidate pathways
# (Robust_Candidate == TRUE in AllPathways_Results_Tracking_Table.csv, see
# Aggregation_Results_Tracking_Table.R) overlap pairwise, and groups
# them, using the same Jaccard threshold (J >= 0.3) already applied in the
# initial pathway selection (two-zone model, section 2.1),
# into "effectively independent" clusters.
#
# Background/rationale: the two-zone check during pathway selection
# only compared THEMATICALLY CORRESPONDING candidate pairs (e.g. a
# Reactome candidate version against the already listed KEGG version
# of the same pathway) to avoid duplicates between the two databases
# -- NOT a comprehensive comparison of all possible pairs of the
# final 54 pathways among each other. Between biologically neighbouring, but
# differently named/categorized pathways (e.g. different
# sub-areas of energy metabolism); gene overlap was never checked.
# This script closes exactly this gap for the robust candidates of
# Paper 1: "N robust candidates" is only a statement about N
# independent biological signals if the underlying gene sets do
# not overlap substantially.
#
# IMPORTANT: this analysis is NOT an additional exclusion criterion. No pathway
# is removed from Robust_Candidate because of overlap -- each one has
# independently met the four pre-specified criteria (FDR, sensitivity, test-set
# direction, effect size). The result serves purely for
# transparent INTERPRETATION ("how many of the N robust candidates are
# genetically effectively independent signals") in the presentation of results and
# discussion of Paper 1.
#
# Input:
#   - Pathway_Gene_Selection_final.xlsx (Sheet "Pathway_Gene_Selection"):
#     Long-format table with one row per gene per pathway (columns
#     Pathway, Gene, ... , Datensatz-Art). Four pathway labels in
#     of this file do NOT directly follow the pipeline naming convention and
#     are resolved explicitly below (see build_pathway_gene_sets()):
#       - "Pyruvate_Metabolism"        -> split by Datensatz-Art into
#                                          Pyruvate_Metabolism_KEGG /
#                                          Pyruvate_Metabolism_Reactome
#       - "Glucagon_Signaling_Pathway" -> Glucagon_Signaling_KEGG
#       - "Glucagon_Signaling"         -> Glucagon_Signaling_Reactome
#       - "Glucose_Metabolism"         -> Glucose_Metabolism_Reactome
#     All remaining 49 names correspond 1:1 to the pipeline naming convention
#     (incl. the 6 pathways whose entire gene set is labelled "KEGG erweitert"
#     -- this is their regular main definition, see
#     Section 2.1/Zone 2, no special handling needed).
#   - AllPathways_Results_Tracking_Table.csv (from
#     Aggregation_Results_Tracking_Table.R): provides the current list of
#     robust candidates automatically -- no manual pathway list needed.
#
# Output:
#   - AllPathways_Gene_Overlap_Jaccard_RobustCandidates.csv (all pairs among
#     the robust candidates, with Jaccard index and cluster assignment)
#   - AllPathways_Gene_Overlap_Jaccard_All54_HighOverlapPairs.csv (supplementary:
#     all pairs among ALL 54 pathways with J >= threshold -- for the
#     discussion of the effective number of independent tests in the
#     BH-FDR correction; deliberately only the notable pairs, not the full
#     54x54 matrix, to keep the file concise)
#   - Automatic entry in Appendix D of the methods/results log
#     (log_jaccard_to_md()) -- like Appendix C, complete replacement on
#     every run, no append/upsert.
# ==============================================================================
library(dplyr)
library(readxl)

output_dir <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Repair/OutPW")   # REPAIR
gene_selection_file <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Pathway_Selection/Pathway_Gene_Selection_final.xlsx")
pathway_liste_file <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Pathway_Selection/Pathway_Liste.csv")
methods_log_path <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Repair/Methods_Results_Log_Repair.md")   # REPAIR

output_file <- function(...) file.path(output_dir, paste0(...))

JACCARD_CLUSTER_THRESHOLD <- 0.3

# ==============================================================================
# 1) Reconstruct gene sets per pathway from the Excel file (54 pathways)
# ==============================================================================

build_pathway_gene_sets <- function(path) {

  raw <- read_excel(path, sheet = "Pathway_Gene_Selection")

  required_cols <- c("Pathway", "Gene", "Datensatz-Art")
  missing_cols <- setdiff(required_cols, colnames(raw))
  if (length(missing_cols) > 0) {
    stop(
      "Erwartete Spalten fehlen in ", path, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  raw <- raw[!is.na(raw$Pathway) & !is.na(raw$Gene), ]

  special_source_names <- c(
    "Pyruvate_Metabolism", "Glucagon_Signaling_Pathway",
    "Glucagon_Signaling", "Glucose_Metabolism"
  )

  gene_sets <- list()

  standard_rows <- raw[!(raw$Pathway %in% special_source_names), ]
  for (p in unique(standard_rows$Pathway)) {
    gene_sets[[p]] <- unique(standard_rows$Gene[standard_rows$Pathway == p])
  }

  # --- Special cases: name resolution as agreed (see script header) ---
  pyr <- raw[raw$Pathway == "Pyruvate_Metabolism", ]
  gene_sets[["Pyruvate_Metabolism_KEGG"]] <- unique(
    pyr$Gene[pyr$`Datensatz-Art` == "KEGG"]
  )
  gene_sets[["Pyruvate_Metabolism_Reactome"]] <- unique(
    pyr$Gene[pyr$`Datensatz-Art` == "Reactome"]
  )

  gsp <- raw[raw$Pathway == "Glucagon_Signaling_Pathway", ]
  gene_sets[["Glucagon_Signaling_KEGG"]] <- unique(gsp$Gene)

  gs <- raw[raw$Pathway == "Glucagon_Signaling", ]
  gene_sets[["Glucagon_Signaling_Reactome"]] <- unique(gs$Gene)

  gm <- raw[raw$Pathway == "Glucose_Metabolism", ]
  gene_sets[["Glucose_Metabolism_Reactome"]] <- unique(gm$Gene)

  # --- Validation against Pathway_Liste.csv (54 names) ---
  pathways_expected <- read.csv(pathway_liste_file, stringsAsFactors = FALSE)$pathway_name
  missing_sets <- setdiff(pathways_expected, names(gene_sets))
  empty_sets <- names(gene_sets)[vapply(gene_sets, length, integer(1)) == 0]

  if (length(missing_sets) > 0) {
    stop(
      "For the following pathways from Pathway_Liste.csv no gene set could be ",
      path, " rekonstruiert werden: ", paste(missing_sets, collapse = ", ")
    )
  }
  if (length(empty_sets) > 0) {
    warning(
      "The following gene sets are empty (0 genes) -- please check: ",
      paste(empty_sets, collapse = ", ")
    )
  }

  gene_sets[pathways_expected]
}

gene_sets <- build_pathway_gene_sets(gene_selection_file)

cat(
  "Gene sets reconstructed for ", length(gene_sets), " Pathways (erwartet: 54).\n",
  sep = ""
)

# ==============================================================================
# 2) Jaccard index function + cluster formation (union-find, no extra package)
# ==============================================================================

jaccard_index <- function(a, b) {
  if (length(a) == 0 || length(b) == 0) return(NA_real_)
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

make_union_find <- function(names_vec) {
  parent <- setNames(names_vec, names_vec)

  find_root <- function(x) {
    while (parent[[x]] != x) {
      parent[[x]] <<- parent[[parent[[x]]]]
      x <- parent[[x]]
    }
    x
  }

  union_sets <- function(x, y) {
    rx <- find_root(x)
    ry <- find_root(y)
    if (rx != ry) parent[[rx]] <<- ry
  }

  list(find = find_root, union = union_sets)
}

compute_pairwise_jaccard <- function(pathway_names, gene_sets, threshold) {

  if (length(pathway_names) < 2) {
    return(list(
      pairs = data.frame(
        Pathway_A = character(0), Pathway_B = character(0),
        N_Genes_A = integer(0), N_Genes_B = integer(0),
        N_Intersection = integer(0), N_Union = integer(0),
        Jaccard = numeric(0), Above_Threshold = logical(0),
        Shared_Genes = character(0), stringsAsFactors = FALSE
      ),
      clusters = data.frame(
        Pathway = pathway_names, Cluster_ID = seq_along(pathway_names),
        stringsAsFactors = FALSE
      )
    ))
  }

  combos <- combn(pathway_names, 2, simplify = FALSE)
  uf <- make_union_find(pathway_names)

  pair_rows <- lapply(combos, function(pair) {
    a <- pair[1]; b <- pair[2]
    ga <- gene_sets[[a]]; gb <- gene_sets[[b]]
    shared <- intersect(ga, gb)
    j <- jaccard_index(ga, gb)
    above <- !is.na(j) && j >= threshold
    if (above) uf$union(a, b)
    data.frame(
      Pathway_A = a, Pathway_B = b,
      N_Genes_A = length(ga), N_Genes_B = length(gb),
      N_Intersection = length(shared), N_Union = length(union(ga, gb)),
      Jaccard = j, Above_Threshold = above,
      Shared_Genes = paste(sort(shared), collapse = "; "),
      stringsAsFactors = FALSE
    )
  })

  pairs <- dplyr::bind_rows(pair_rows)
  pairs <- pairs[order(-pairs$Jaccard), ]

  roots <- vapply(pathway_names, uf$find, character(1))
  cluster_id_map <- setNames(seq_along(unique(roots)), unique(roots))
  clusters <- data.frame(
    Pathway = pathway_names,
    Cluster_ID = unname(cluster_id_map[roots]),
    stringsAsFactors = FALSE
  )
  clusters <- clusters[order(clusters$Cluster_ID, clusters$Pathway), ]

  list(pairs = pairs, clusters = clusters)
}

# ==============================================================================
# 3) Load robust candidates from the aggregation output
# ==============================================================================

tracking_file <- output_file("AllPathways_Results_Tracking_Table.csv")
if (!file.exists(tracking_file)) {
  stop(
    "File not found: ", tracking_file,
    " - please run Aggregation_Results_Tracking_Table.R first."
  )
}

tracking <- read.csv2(tracking_file, stringsAsFactors = FALSE)
tracking$Robust_Candidate <- suppressWarnings(as.logical(tracking$Robust_Candidate))

robust_pathways <- tracking$Pathway[
  which(!is.na(tracking$Robust_Candidate) & tracking$Robust_Candidate)
]

cat(
  "\nRobust candidates (from ", basename(tracking_file), "): ",
  length(robust_pathways), "\n", sep = ""
)

missing_gene_sets <- setdiff(robust_pathways, names(gene_sets))
if (length(missing_gene_sets) > 0) {
  stop(
    "For the following robust candidates a gene set is missing: ",
    paste(missing_gene_sets, collapse = ", ")
  )
}

# ==============================================================================
# 4) Jaccard analysis among the robust candidates
# ==============================================================================

result_candidates <- compute_pairwise_jaccard(
  robust_pathways, gene_sets, JACCARD_CLUSTER_THRESHOLD
)

n_clusters <- length(unique(result_candidates$clusters$Cluster_ID))

cat("\n====================================================\n")
cat("Pairwise Jaccard indices among the robust candidates\n")
cat("====================================================\n")
print(
  result_candidates$pairs[
    , c("Pathway_A", "Pathway_B", "N_Intersection", "Jaccard", "Above_Threshold")
  ],
  row.names = FALSE
)

cat("\n====================================================\n")
cat("Cluster assignment (threshold J >= ", JACCARD_CLUSTER_THRESHOLD, ")\n", sep = "")
cat("====================================================\n")
print(result_candidates$clusters, row.names = FALSE)
cat(
  "\n-> ", length(robust_pathways), " robuste Kandidaten bilden ",
  n_clusters, " effectively independent cluster(s).\n", sep = ""
)

candidates_output_path <- output_file("AllPathways_Gene_Overlap_Jaccard_RobustCandidates.csv")
write.csv2(result_candidates$pairs, candidates_output_path, row.names = FALSE)

clusters_output_path <- output_file("AllPathways_Gene_Overlap_Jaccard_RobustCandidates_Clusters.csv")
write.csv2(result_candidates$clusters, clusters_output_path, row.names = FALSE)

# ==============================================================================
# 5) Additionally: notable pairs (J >= threshold) across ALL 54 pathways
#    (transparency for the BH-FDR independence assumption, not part of the
#    candidate assessment)
# ==============================================================================

all_pathway_names <- names(gene_sets)
result_all <- compute_pairwise_jaccard(
  all_pathway_names, gene_sets, JACCARD_CLUSTER_THRESHOLD
)
high_overlap_all <- result_all$pairs[result_all$pairs$Above_Threshold, ]
high_overlap_all <- high_overlap_all[order(-high_overlap_all$Jaccard), ]

cat("\n====================================================\n")
cat(
  "Pairs with J >= ", JACCARD_CLUSTER_THRESHOLD, " among ALL 54 pathways: ",
  nrow(high_overlap_all), " of ", nrow(result_all$pairs), " possible pairs\n",
  sep = ""
)
cat("====================================================\n")
if (nrow(high_overlap_all) > 0) {
  print(
    high_overlap_all[, c("Pathway_A", "Pathway_B", "N_Intersection", "Jaccard")],
    row.names = FALSE
  )
}

all54_output_path <- output_file("AllPathways_Gene_Overlap_Jaccard_All54_HighOverlapPairs.csv")
write.csv2(high_overlap_all, all54_output_path, row.names = FALSE)

cat("\n====================================================\n")
cat("Done.\n")
cat("Robust candidate pairs:   ", candidates_output_path, "\n", sep = "")
cat("Cluster assignment:          ", clusters_output_path, "\n", sep = "")
cat("All-54 overlap pairs: ", all54_output_path, "\n", sep = "")
cat("====================================================\n")

# ==============================================================================
# 6) Write automatically to the methods/results log (Appendix D)
# ==============================================================================

log_jaccard_to_md <- function(md_path, pairs_df, clusters_df, n_clusters,
                               n_candidates, threshold, high_overlap_all_df,
                               n_all_pairs_checked) {

  if (!file.exists(md_path)) {
    warning(
      "Methods/results log not found at: ", md_path,
      " - Appendix D entries were NOT written."
    )
    return(invisible(NULL))
  }

  log_lines <- readLines(md_path, warn = FALSE, encoding = "UTF-8")

  replace_block <- function(lines, begin_tag, end_tag, header_lines, body_lines) {
    begin_marker <- which(trimws(lines) == begin_tag)
    end_marker <- which(trimws(lines) == end_tag)

    if (length(begin_marker) != 1 || length(end_marker) != 1 ||
        end_marker <= begin_marker) {
      warning(
        "Marker '", begin_tag, "' / '", end_tag,
        "' not found uniquely or in wrong order in ",
        md_path, " - this block was NOT written."
      )
      return(list(lines = lines, ok = FALSE))
    }

    new_block <- c(lines[begin_marker], header_lines, body_lines, lines[end_marker])
    updated <- c(
      lines[seq_len(begin_marker - 1)], new_block,
      lines[seq(end_marker + 1, length(lines))]
    )
    list(lines = updated, ok = TRUE)
  }

  safe_cell <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x) || x == "") return("")
    gsub("\\|", "/", gsub("[\r\n]+", " ", as.character(x)))
  }

  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M")

  #--- Block 1: pairwise Jaccard values among the robust candidates -------
  pairs_header <- c(
    paste0(
      "| Pathway A | Pathway B | Gemeinsame Gene | Jaccard | ",
      "Ueber Schwelle (", threshold, ") | Last_Updated |"
    ),
    "|---|---|---|---|---|---|"
  )

  if (nrow(pairs_df) == 0) {
    pairs_rows <- c(paste0(
      "| *(< 2 robuste Kandidaten -- keine Paare zu vergleichen)* |  |  |  |  | ",
      timestamp, " |"
    ))
  } else {
    pairs_rows <- vapply(seq_len(nrow(pairs_df)), function(i) {
      r <- pairs_df[i, ]
      sprintf(
        "| %s | %s | %s | %s | %s | %s |",
        safe_cell(r$Pathway_A), safe_cell(r$Pathway_B),
        safe_cell(r$N_Intersection), safe_cell(round(r$Jaccard, 4)),
        safe_cell(r$Above_Threshold), timestamp
      )
    }, character(1))
  }

  step1 <- replace_block(
    log_lines, "<!-- AUTOLOG-JACCARD:BEGIN -->", "<!-- AUTOLOG-JACCARD:END -->",
    pairs_header, pairs_rows
  )

  #--- Block 2: cluster assignment --------------------------------------------
  clusters_header <- c(
    "| Pathway | Cluster_ID | Last_Updated |",
    "|---|---|---|"
  )
  clusters_rows <- vapply(seq_len(nrow(clusters_df)), function(i) {
    r <- clusters_df[i, ]
    sprintf(
      "| %s | %s | %s |", safe_cell(r$Pathway), safe_cell(r$Cluster_ID), timestamp
    )
  }, character(1))

  summary_line <- sprintf(
    "*%d robuste Kandidaten bilden %d effektiv unabhaengige(s) Cluster (Schwelle J >= %s). Zuletzt aktualisiert: %s.*",
    n_candidates, n_clusters, threshold, timestamp
  )

  step2 <- replace_block(
    step1$lines,
    "<!-- AUTOLOG-JACCARD-CLUSTERS:BEGIN -->", "<!-- AUTOLOG-JACCARD-CLUSTERS:END -->",
    c(summary_line, "", clusters_header), clusters_rows
  )

  #--- Block 3: notable pairs among ALL 54 (FDR transparency) ----------
  all_header <- c(
    "| Pathway A | Pathway B | Gemeinsame Gene | Jaccard | Last_Updated |",
    "|---|---|---|---|---|"
  )
  if (nrow(high_overlap_all_df) == 0) {
    all_rows <- c(paste0(
      "| *(keine Paare mit J >= ", threshold, " unter allen 54 Pathways)* |  |  |  | ",
      timestamp, " |"
    ))
  } else {
    all_rows <- vapply(seq_len(nrow(high_overlap_all_df)), function(i) {
      r <- high_overlap_all_df[i, ]
      sprintf(
        "| %s | %s | %s | %s | %s |",
        safe_cell(r$Pathway_A), safe_cell(r$Pathway_B),
        safe_cell(r$N_Intersection), safe_cell(round(r$Jaccard, 4)), timestamp
      )
    }, character(1))
  }
  all_summary_line <- sprintf(
    "*%d von %d moeglichen Paaren unter allen 54 Pathways zeigen J >= %s. Zuletzt aktualisiert: %s.*",
    nrow(high_overlap_all_df), n_all_pairs_checked, threshold, timestamp
  )

  step3 <- replace_block(
    step2$lines,
    "<!-- AUTOLOG-JACCARD-ALL54:BEGIN -->", "<!-- AUTOLOG-JACCARD-ALL54:END -->",
    c(all_summary_line, "", all_header), all_rows
  )

  if (!step1$ok && !step2$ok && !step3$ok) {
    return(invisible(NULL))
  }

  writeLines(step3$lines, md_path, useBytes = TRUE)

  cat(
    "\nMethods/Results-Log (Appendix D) aktualisiert: ", md_path,
    " (Kandidaten-Paare: ", ifelse(step1$ok, "OK", "NOT written"),
    "; Cluster: ", ifelse(step2$ok, "OK", "NOT written"),
    "; all-54 transparency: ", ifelse(step3$ok, "OK", "NOT written"), ")\n",
    sep = ""
  )

  invisible(NULL)
}

log_jaccard_to_md(
  methods_log_path,
  result_candidates$pairs, result_candidates$clusters, n_clusters,
  length(robust_pathways), JACCARD_CLUSTER_THRESHOLD,
  high_overlap_all, nrow(result_all$pairs)
)

