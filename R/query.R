#' Given a data frame of HCA metadata, returns a SingleCellExperiment object corresponding to the samples in that data frame
#'
#' @param .data A data frame containing, at minimum, a `.sample` column, which corresponds to a single cell sample ID.
#' This can be obtained from the [get_metadata()] function.
#' @param repository A character vector of length one, which is a file path to where the single cell data is stored
#' @param genes An optional character vector of genes to return the counts for. By default counts for all genes will be returned.
#'
#' @importFrom dplyr pull filter
#' @importFrom tidySingleCellExperiment inner_join 
#' @importFrom purrr reduce map map_int
#' @importFrom BiocGenerics cbind
#' @importFrom glue glue
#' @importFrom dplyr as_tibble
#' @importFrom HDF5Array loadHDF5SummarizedExperiment HDF5RealizationSink loadHDF5SummarizedExperiment
#' @importFrom stringr str_remove
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData assayNames<-
#'
#' @export
#'
#'
get_SingleCellExperiment = function(
  .data,
  repository = "/vast/projects/RCP/human_cell_atlas/splitted_DB2_data",
  genes = NULL
){
  # We have to convert to an in-memory table here, or some of the dplyr operations will fail when passed a database connection
  raw_data = as_tibble(.data)

	files_to_read =
	  raw_data |>
		pull(file_id_db) |>
		unique() |>
		as.character()

	glue("Reading {length(files_to_read)} files.") |>
	  message()

	# Load each file
	sces =
		files_to_read |>
		map(~ {
			cat(".")
		  
		  sce = glue("{repository}/{.x}") |>
			  loadHDF5SummarizedExperiment()
		  
		  if (!is.null(genes)){
		    # Optionally subset the genes
		    sce = sce[
		      intersect(genes, rownames(sce))  
		    ]
		  }
		  
		  sce |>
			  inner_join(
  				# Needed because cell IDs are not unique outside the file_id or file_id_db
  				filter(raw_data, file_id_db == .x),
  				by=".cell"
  			)
		})

	# Drop files with one cell, which causes
	# the DFrame objects to combine must have the same column names
	sces = sces[map_int(sces, ncol)>1]

	cat("\n")

	# Combine
	sce =
		sces |>
		do.call(cbind, args=_)

	# Rename assay
	assayNames(sce) = "counts"

	# Return
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
#' @importFrom RSQLite SQLite SQLITE_RO
#' @importFrom dplyr tbl
#'
get_metadata = function(sqlite_path = "/vast/projects/RCP/human_cell_atlas/metadata.sqlite"){
  SQLite() |>
    dbConnect(drv=_, dbname=sqlite_path, flags=SQLITE_RO) |>
    tbl("metadata")
}
