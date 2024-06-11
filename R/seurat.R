# Functions that relate to the Seurat conversion

#' @importFrom assertthat assert_that
#' @importFrom methods as
#' @importFrom SeuratObject as.sparse
#' @exportS3Method
as.sparse.DelayedMatrix <- function(x, ...) {
    # This is glue to ensure the SCE -> Seurat conversion works properly with
    # DelayedArray types
    as(x, "dgCMatrix")
}

#' Given a data frame of HCA metadata, returns a Seurat object corresponding to
#' the samples in that data frame
#'
#' @inheritDotParams get_single_cell_experiment
#' @importFrom SeuratObject as.Seurat
#' @export
#' @return A Seurat object containing the same data as a call to
#'   [get_single_cell_experiment()]
#' @examples
#' meta <- get_metadata() |> head(2)
#' seurat <- get_seurat(meta)
#' @references Mangiola, S., M. Milton, N. Ranathunga, C. S. N. Li-Wai-Suen, 
#'   A. Odainic, E. Yang, W. Hutchison et al. "A multi-organ map of the human 
#'   immune system across age, sex and ethnicity." bioRxiv (2023): 2023-06.
#'   doi:10.1101/2023.06.08.542671.
#' @source [Mangiola et al.,2023](https://www.biorxiv.org/content/10.1101/2023.06.08.542671v3)
get_seurat <- function(...) {
    get_single_cell_experiment(...) |> as.Seurat(data = NULL)
}
