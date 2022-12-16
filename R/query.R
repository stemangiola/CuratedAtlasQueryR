#' Given a data frame of HCA metadata, returns a SingleCellExperiment object corresponding to the samples in that data frame
#'
#' @param .data A data frame containing, at minimum, a `.sample` column, which corresponds to a single cell sample ID.
#' This can be obtained from the [get_metadata()] function.
#' @param repository A character vector of length one. It should be either a local file path or an HTTP URL pointing to the location where the single cell data is stored.
#' @param cache_dir An optional character vector of length one. If provided, it should indicate a local file path where any remotely accessed files should be copied.
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
#' @importFrom purrr when
#' @importFrom magrittr equals
#' @importFrom httr parse_url
#'
#' @export
#'
#'
get_SingleCellExperiment = function(
  .data,
  repository,
  cache_dir = get_default_cache_dir(),
  genes = NULL
){
  cache_dir |> dir.create(showWarnings = FALSE)
  
  parsed_repo = parse_url(repository)
  
  # We have to convert to an in-memory table here, or some of the dplyr operations will fail when passed a database connection
  raw_data = as_tibble(.data)

	files_to_read =
	  raw_data |>
		pull(file_id_db) |>
		unique() |>
		as.character()
	
	if(!is.null(parsed_repo$scheme)){
	  if (parsed_repo$scheme %in% c("http", "https")){
	    sync_remote_files(parsed_repo, cache_dir, files_to_read)
	    repository = cache_dir
	  }
	  else {
	    glue('Unknown URL scheme "{parsed_repo$scheme}"') |>
	      stop()
	  }
	}

	glue("Reading {length(files_to_read)} files.") |>
	  message()

	# Load each file
	sces =
		files_to_read |>
		map(~ {
			cat(".")

		  sce = file.path(
		    repository,
		    .x
		  ) |>
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

	# Rename assay THIS WILL NOT BE NEEDED EVENTUALLY
	assayNames(sce) = "counts"

	# Return
	sce
}

#' Synchronises one or more remote files with a local copy
#'
#' @param url A character vector of length one. The base HTTP url from which to obtain the files.
#' @param cache_dir A character vector of length one. The local filepath to synchronise files to.
#' @param files A character vector containing one or more file_id_db entries
#'
#' @return Subject to change
#' @importFrom purrr cross2 walk
#' @importFrom httr modify_url GET write_disk
#' @export
#'
sync_remote_files = function(
  url,
  cache_dir,
  files
){
  c("assays.h5", "se.rds") |>
    cross2(files) |>
    walk(\(path_elements){
      # Path to the file of interest from the root path. We use "/"
      # since URLs must use these regardless of OS
      url_path = paste0(url$path, "/", path_elements[[2]], "/", path_elements[[1]])
      
      # Path to save the file on local disk
      # We use file.path since the file separator will differ on other OSs
      output_dir = file.path(
        cache_dir,
        path_elements[[2]]
      )
      output_file = file.path(
        output_dir,
        path_elements[[1]]
      )
      
      # Create the parent dirs
      dir.create(output_dir, recursive=TRUE, showWarnings = FALSE)
      
      # Perform the download
      modify_url(url, path=url_path) |>
        GET(write_disk(output_file))
    })
}

#' Returns the default cache directory
#'
#' @return A length one character vector.
#' @export
#' @importFrom rappdirs user_cache_dir
#'
get_default_cache_dir = function(){
  file.path(
    user_cache_dir(),
    "hca_harmonised"
  )
}

#' @importFrom SeuratObject as.sparse
#' @exportS3Method
as.sparse.DelayedMatrix = function(x){
  # This is glue to ensure the SCE -> Seurat conversion works properly with
  # DelayedArray types
  as(x, "dgCMatrix")
}

#' Given a data frame of HCA metadata, returns a Seurat object corresponding to the samples in that data frame
#'
#' @inheritDotParams get_SingleCellExperiment
#' @importFrom Seurat as.Seurat
#' @export
get_seurat = function(
  ...
){
  get_SingleCellExperiment(...) |> as.Seurat(data=NULL)
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
