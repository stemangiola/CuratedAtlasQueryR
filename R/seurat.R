# Functions that relate to the Seurat conversion

#' @importFrom assertthat assert_that
#' @importFrom methods as
#' @importFrom SeuratObject as.sparse
#' @exportS3Method
as.sparse.DelayedMatrix <- function(x) {
    # This is glue to ensure the SCE -> Seurat conversion works properly with
    # DelayedArray types
    as(x, "dgCMatrix")
}

#' Given a data frame of HCA metadata, returns a Seurat object corresponding to
#' the samples in that data frame
#'
#' @inheritDotParams get_SingleCellExperiment
#' @importFrom SeuratObject as.Seurat
#' @export
#' @return A Seurat object containing the same data as a call to
#'   get_SingleCellExperiment.
#' @examples
#' meta <- get_metadata() |> head(2)
#' seurat <- get_seurat(meta)
#'
get_seurat <- function(...) {
    get_SingleCellExperiment(...) |> as.Seurat(data = NULL)
}
