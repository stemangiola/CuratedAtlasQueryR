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
#' @importFrom stringr str_remove
#'
#' @export
#'
#'
get_SingleCellExperiment = function(.data, repository = "/vast/projects/RCP/human_cell_atlas/splitted_light_data/"){

  # We have to convert to an in-memory table here, or some of the dplyr operations will fail when passed a database connection
  raw_data = as_tibble(.data)

	files_to_read =
	  raw_data |>
		pull(file_id_db) |>
		unique() |>
		as.character()

	message(glue("Reading {length(files_to_read)} files."))





	sce =
		files_to_read |>
		map(~ {
			cat(".")
			x = loadHDF5SummarizedExperiment(glue("{repository}/{.x}")	) |>
				inner_join(
					raw_data |>
						filter(file_id_db == .x) |>
						mutate(.cell_original = .cell |> str_remove(.sample) |> str_remove("_$")),
					by=c(".cell" = ".cell_original"))

			}
		)

	# Harmonise genes
	all_genes =
		sce |>
		map(rownames) |>
		unlist() |>
		unique()

	sce =
		sce |>
		map(~ {
			missing_genes = all_genes  |> setdiff(rownames(.x))

			missing_matrix =
				HDF5RealizationSink(c(length(missing_genes),ncol(.x)), as.sparse = TRUE) |>
				as("DelayedArray")

			rownames(missing_matrix) = missing_genes
			colnames(missing_matrix) = colnames(.x)

			missing_sce = SingleCellExperiment(list(X=missing_matrix),  colData=colData(.x))
			missing_sce@int_colData = .x@int_colData

			# Make cell name unique
			.x |> rbind(missing_sce	)
		})

	cat("\n")

	sce |>

		# # Temporary
		# map(~ {
		# 	x = .x
		# 	x@int_colData$colPairs = x@int_colData$colPairs[,0]
		# 	x
		# }) |>

		# Combine
		do.call(cbind, args=_)




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
#' @importFrom RSQLite SQLITE_RO
#' @importFrom dplyr tbl
#'
get_metadata = function(sqlite_path = "/vast/projects/RCP/human_cell_atlas/metadata.sqlite"){

  SQLite() |>
    dbConnect(drv=_, dbname=sqlite_path, flags=SQLITE_RO) |>
    tbl("metadata")
}
