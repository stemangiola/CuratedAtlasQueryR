#' Checks whether genes in a list of SCE objects overlap
#'
#' @param sce_list A list of SingleCellExperiment object
#' @return A character vector of genes intersection across SCE objects
#' @examples
#' sce1 <-readRDS("~/projects/caq/import_api_pipelines/12eb5fe25994253c1d320ca590a6e999/se.rds")
#' sce2 <- HeOrganAtlasData(tissue = "Liver", ensembl = FALSE, location = TRUE)
#' sce_list <- list(sce1, sce2)
#' sce_list |> check_gene_overlap()
#'
#' @importFrom purrr map reduce map_int
#' @importFrom cli cli_alert_warning

check_gene_overlap <- function(sce_list) {
  # Extract gene names from each SCE object in the list, find intersection, 
  # calculate the number of unique genes in each SCE object
  gene_lists <- map(sce_list, rownames)
  common_genes <- reduce(gene_lists, intersect)
  unique_genes_count <- map_int(gene_lists, ~length(setdiff(.x, common_genes)))
  
  if(any(unique_genes_count > 0)) {
    single_line_str(
      "CuratedAtlasQuery reports: Not all genes completely overlap across the provided SCE objects"
    ) |> cli_alert_warning()
  } 
  
  common_genes
}
