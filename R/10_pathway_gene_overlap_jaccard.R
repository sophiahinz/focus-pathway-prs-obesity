# ==============================================================================
# 10_pathway_gene_overlap_jaccard.R
# Pairwise gene overlap (Jaccard index) between the robust candidate pathways
# and among all 54 pathways (Results).
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================

# ==============================================================================
# Purpose: checks how strongly the gene sets of the robust candidate pathways
# (Robust_Candidate == TRUE in AllPathways_Results_Tracking_Table.csv, see
# 09_aggregation_robust_candidates.R) overlap pairwise, and groups
# them, using the same Jaccard threshold (J >= 0.3) as in the pathway
# selection (two-zone model), into "effectively independent" clusters.
#
# Background/rationale: the two-zone check during pathway selection
# only compared THEMATICALLY CORRESPONDING candidate pairs (e.g. a
# Reactome candidate version against the already listed KEGG version
# of the same pathway) to avoid duplicates between the two databases
# -- NOT a comprehensive comparison of all possible pairs of the
# final 54 pathways among each other. Gene overlap between biologically
# neighbouring but differently named/categorized pathways (e.g. different
# sub-areas of energy metabolism) was therefore never checked.
# This script closes this gap for the robust candidates:
# "N robust candidates" is only a statement about N
# independent biological signals if the underlying gene sets do
# not overlap substantially.
#
# IMPORTANT: this analysis is NOT an additional exclusion criterion. No pathway
# is removed from Robust_Candidate because of overlap -- each one has
# independently met the four pre-specified criteria (FDR, sensitivity, test-set
# direction, effect size). The result serves purely for
# transparent INTERPRETATION ("how many of the N robust candidates are
# genetically effectively independent signals") in the Results and
# Discussion.
#
# Input:
#   - gene_sets/Pathway_Selection_Supplement_GeneSets.csv (long format, one
#     row per gene per pathway). The main-analysis gene set of each pathway is
#     used: the "KEGG erweitert" set where it exists (two-zone model, zone 2),
#     otherwise the KEGG or Reactome set -- identical to 07_pathway_pipeline.R.
#   - AllPathways_Results_Tracking_Table.csv (from
#     09_aggregation_robust_candidates.R): list of robust candidates.
#
# Output:
#   - AllPathways_Gene_Overlap_Jaccard_RobustCandidates.csv (all pairs among
#     the robust candidates, with Jaccard index and cluster assignment)
#   - AllPathways_Gene_Overlap_Jaccard_All54_HighOverlapPairs.csv (supplementary:
#     all pairs among ALL 54 pathways with J >= threshold -- for the
#     discussion of the effective number of independent tests in the
#     BH-FDR correction; deliberately only the notable pairs, not the full
#     54x54 matrix, to keep the file concise)
#   - Optional: entry in Appendix D of the methods/results log
#     (log_jaccard_to_md(); complete replacement on every run), if that
#     file exists.
# ==============================================================================
library(dplyr)

output_dir <- file.path(Sys.getenv("FOCUS_RESULTS_DIR", unset = "path/to/results"), "OutPW")
gene_selection_file <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Pathway_Selection/Pathway_Selection_Supplement_GeneSets.csv")
pathway_liste_file <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Pathway_Selection/Pathway_Liste.csv")
methods_log_path <- file.path(Sys.getenv("FOCUS_RESULTS_DIR", unset = "path/to/results"), "Methods_Results_Log.md")

output_file <- function(...) file.path(output_dir, paste0(...))

JACCARD_CLUSTER_THRESHOLD <- 0.3

# ==============================================================================
# 1) Gene sets per pathway from the gene-set CSV (54 pathways)
# ==============================================================================

build_pathway_gene_sets <- function(path) {

  raw <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                  encoding = "UTF-8")

  required_cols <- c("Pathway-Name", "Gen", "Datensatz-Art")
  missing_cols <- setdiff(required_cols, colnames(raw))
  if (length(missing_cols) > 0) {
    stop(
      "Expected columns missing in ", path, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  raw <- raw[!is.na(raw$`Pathway-Name`) & !is.na(raw$Gen), ]

  # Main-analysis gene set per pathway (identical to 07_pathway_pipeline.R):
  # the "KEGG erweitert" set where it exists, otherwise the KEGG/Reactome set
  gene_sets <- list()
  for (p in unique(raw$`Pathway-Name`)) {
    rows <- raw[raw$`Pathway-Name` == p, ]
    if (any(rows$`Datensatz-Art` == "KEGG erweitert")) {
      rows <- rows[rows$`Datensatz-Art` == "KEGG erweitert", ]
    } else {
      rows <- rows[rows$`Datensatz-Art` %in% c("KEGG", "Reactome"), ]
    }
    gene_sets[[p]] <- unique(rows$Gen)
  }

  # --- Validation against Pathway_Liste.csv (54 names) ---
  pathways_expected <- read.csv(pathway_liste_file, stringsAsFactors = FALSE)$pathway_name
  missing_sets <- setdiff(pathways_expected, names(gene_sets))
  empty_sets <- names(gene_sets)[vapply(gene_sets, length, integer(1)) == 0]

  if (length(missing_sets) > 0) {
    stop(
      "For the following pathways from Pathway_Liste.csv no gene set was found in ",
      path, ": ", paste(missing_sets, collapse = ", ")
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
  "Gene sets reconstructed for ", length(gene_sets), " pathways (expected: 54).\n",
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
    " - please run 09_aggregation_robust_candidates.R first."
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
  "\n-> ", length(robust_pathways), " robust candidates form ",
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
      "| Pathway A | Pathway B | Shared genes | Jaccard | ",
      "Above threshold (", threshold, ") | Last_Updated |"
    ),
    "|---|---|---|---|---|---|"
  )

  if (nrow(pairs_df) == 0) {
    pairs_rows <- c(paste0(
      "| *(< 2 robust candidates -- no pairs to compare)* |  |  |  |  | ",
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
    "*%d robust candidates form %d effectively independent cluster(s) (threshold J >= %s). Last updated: %s.*",
    n_candidates, n_clusters, threshold, timestamp
  )

  step2 <- replace_block(
    step1$lines,
    "<!-- AUTOLOG-JACCARD-CLUSTERS:BEGIN -->", "<!-- AUTOLOG-JACCARD-CLUSTERS:END -->",
    c(summary_line, "", clusters_header), clusters_rows
  )

  #--- Block 3: notable pairs among ALL 54 (FDR transparency) ----------
  all_header <- c(
    "| Pathway A | Pathway B | Shared genes | Jaccard | Last_Updated |",
    "|---|---|---|---|---|"
  )
  if (nrow(high_overlap_all_df) == 0) {
    all_rows <- c(paste0(
      "| *(no pairs with J >= ", threshold, " among all 54 pathways)* |  |  |  | ",
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
    "*%d of %d possible pairs among all 54 pathways show J >= %s. Last updated: %s.*",
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
    "\nMethods/results log (Appendix D) updated: ", md_path,
    " (candidate pairs: ", ifelse(step1$ok, "OK", "NOT written"),
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

