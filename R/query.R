#' Given a data frame of HCA metadata, returns a SingleCellExperiment object corresponding to the samples in that data frame
#'
#' @param .data A data frame containing, at minimum, a `.sample` column, which corresponds to a single cell sample ID.
#' This can be obtained from the [get_metadata()] function.
#' @param repository A character vector of length one, which is a file path to where the single cell data is stored
#'
#' @importFrom dplyr pull
#' @importFrom tidySingleCellExperiment inner_join
#' @importFrom purrr reduce
#' @importFrom purrr map
#' @importFrom BiocGenerics cbind
#' @importFrom glue glue
#' @importFrom dplyr as_tibble
#' @importFrom HDF5Array loadHDF5SummarizedExperiment
#' 
#' @export
#'
#'
get_SingleCellExperiment = function(.data, repository = "/vast/projects/RCP/human_cell_atlas/splitted_light_data/"){
  
  # We have to convert to an in-memory table here, or some of the dplyr operations will fail when passed a database connection
  raw_data = as_tibble(.data)

	files_to_read =
	  raw_data |>
		pull(.sample) |>
		unique() |>
		as.character()

	message(glue("Reading {length(files_to_read)} files."))

	sce =
		files_to_read |>
		map(~ {
			cat(".")
			loadHDF5SummarizedExperiment(glue("{repository}/{.x}")	)
			}
		) |>

		# Temporary
		map(~ {
			x = .x
			x@int_colData$colPairs = x@int_colData$colPairs[,0]
			x
		}) |>

		# Combine
		do.call(cbind, args=_) |>

		# Attach metadata
		inner_join(raw_data, by=".cell")

	cat("\n")

	sce

}


#' Returns a data frame of Human Cell Atlas metadata, which should be filtered
#' and ultimately passed into get_SingleCellExperiment.
#'
#' @param sqlite_path Path to the sqlite database where the metadata can be found.
#' Currently this defaults to an internal location within WEHI's milton system.
#'
#' @export
#' 
#' @importFrom DBI dbConnect
#' @importFrom RSQLite SQLite
#' @importFrom dplyr tbl
#'
get_metadata = function(sqlite_path = NULL){
  if (is.null(sqlite_path)){
    sqlite_path = "/vast/projects/RCP/human_cell_atlas/metadata.sqlite"
  }
  
  SQLite() |>
    dbConnect(drv=_, dbname=sqlite_path) |>
    tbl("metadata")
}