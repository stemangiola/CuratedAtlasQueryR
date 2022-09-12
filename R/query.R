#' sccomp_glm main
#'
#' @description The function for linear modelling takes as input a table of cell counts with three columns containing a cell-group identifier, sample identifier, integer count and the covariates (continuous or discrete). The user can define a linear model with an input R formula, where the first covariate is the factor of interest. Alternatively, sccomp accepts single-cell data containers (Seurat, SingleCellExperiment44, cell metadata or group-size). In this case, sccomp derives the count data from cell metadata.
#'
#'
#' @importFrom dplyr pull
#' @importFrom tidySingleCellExperiment inner_join
#' @importFrom purrr reduce
#' @importFrom purrr map
#' @importFrom zellkonverter readH5AD
#' @importFrom BiocGenerics cbind
#' @importFrom glue glue
#'
#' @export
#'
#'

get_SingleCellExperiment = function(.data, repository = "/vast/scratch/users/mangiola.s/human_cell_atlas/splitted_light_data/"){

	.data |>
		pull(.sample) |>
		unique() |>
		as.character() |>
		map(~ readH5AD(glue("{repository}/{.x}"),	use_hdf5 = TRUE	) ) |>

		# Temporary
		map(~ {
			x = .x
			x@int_colData$colPairs = x@int_colData$colPairs[,0]
			x
		}) %>%

		# Combine
		do.call(cbind, .) |>

		# Attach metadata
		inner_join(.data, by=".cell")
}
