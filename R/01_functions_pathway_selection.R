# ==============================================================================
# 01_functions_pathway_selection.R
# Helper functions for 01_pathway_selection.R: gene-set retrieval from KEGG
# (KEGGREST) and Reactome (ReactomeContentService4R), harmonization to HGNC
# symbols (org.Hs.eg.db) and a per-session cache. Sourced by
# 01_pathway_selection.R; not run on its own.
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================

library(KEGGREST)
library(ReactomeContentService4R)
library(org.Hs.eg.db)
library(AnnotationDbi)

# KEGG: Entrez ID (native) -> HGNC symbol
get_kegg_genes <- function(pathway_id) {
  # pathway_id e.g. "hsa00010"
  entrez_raw <- keggLink("hsa", paste0("path:", pathway_id))
  entrez_ids <- sub("hsa:", "", entrez_raw)
  symbols <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = entrez_ids,
                                   keytype = "ENTREZID",
                                   columns = "SYMBOL")$SYMBOL
  unique(na.omit(symbols))
}

# Reactome: UniProt accession (native) -> HGNC symbol
get_reactome_genes <- function(pathway_id) {
  # pathway_id e.g. "R-HSA-70326"
  participants <- getParticipants(pathway_id,
                                  retrieval = "ReferenceEntities")
  uniprot_ids <- unique(participants$identifier)
  symbols <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = uniprot_ids,
                                   keytype = "UNIPROT",
                                   columns = "SYMBOL")$SYMBOL
  unique(na.omit(symbols))
}

# Cache: avoids repeated downloads of the same pathway ID
gene_sets_cache <- list()

get_gene_set <- function(pathway_id, db) {
  if (!is.null(gene_sets_cache[[pathway_id]])) {
    return(gene_sets_cache[[pathway_id]])
  }
  genes <- if (db == "KEGG") get_kegg_genes(pathway_id) else get_reactome_genes(pathway_id)
  gene_sets_cache[[pathway_id]] <<- genes
  genes
}

# Several comma-separated KEGG IDs in one "Referenz-KEGG-Pathway" cell of the
# tracking spreadsheet (e.g. "hsa00061, hsa00062, hsa00071", when KEGG splits
# one Reactome topic over several maps): the gene sets are combined.
split_kegg_ids <- function(id_string) {
  trimws(strsplit(id_string, ",")[[1]])
}

get_gene_set_combined <- function(id_string, db) {
  ids <- if (db == "KEGG") split_kegg_ids(id_string) else id_string
  unique(unlist(lapply(ids, get_gene_set, db = db)))
}
