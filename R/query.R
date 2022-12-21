# Maps user provided assay names to their corresponding paths in the repository
assay_map = c(
  raw = "original",
  scaled = "cpm"
)

#' Given a data frame of HCA metadata, returns a SingleCellExperiment object corresponding to the samples in that data frame
#'
#' @param data A data frame containing, at minimum, a `.sample` column, which corresponds to a single cell sample ID.
#' This can be obtained from the [get_metadata()] function.
#' @param assay A character vector whose elements must be either "raw" or "scaled", representing the corresponding assay you want to request.
#' @param repository A character vector of length one. If provided, it should be an HTTP URL pointing to the location where the single cell data is stored.
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
#' @importFrom SingleCellExperiment SingleCellExperiment simplifyToSCE
#' @importFrom SummarizedExperiment colData assayNames<-
#' @importFrom purrr when imap
#' @importFrom magrittr equals
#' @importFrom httr parse_url
#' @importFrom assertthat assert_that has_name
#'
#' @export
#'
#'
get_SingleCellExperiment = function(
  data,
  assay = c("raw", "scaled"),
  cache_dir = get_default_cache_dir(),
  repository = NULL,
  genes = NULL
){
  # Parameter validation
  assay %in% names(assay_map) |> all() |> assert_that(msg='assay must be a character vector containing "raw" and/or "scaled"')
  (!anyDuplicated(assay)) |> assert_that()
  inherits(cache_dir, "character") |> assert_that()
  is.null(repository) || is.character(repository) |> assert_that()
  is.null(genes) || is.character(genes) |> assert_that()
  
  # Data parameter validation (last, because it's slower)
  ## Evaluate the promise now so that we get a sensible error message
  data
  ## We have to convert to an in-memory table here, or some of the dplyr operations will fail when passed a database connection
  raw_data = as_tibble(data)
  inherits(raw_data, "tbl") |> assert_that()
  has_name(raw_data, c(".cell", "file_id_db")) |> assert_that()

  cache_dir |> dir.create(showWarnings = FALSE)
  
  files_to_read =
    raw_data |>
    pull(file_id_db) |>
    unique() |>
    as.character()
  
  subdirs = assay_map[assay]
  
  # The repository is optional. If not provided we load only from the cache
  if (!is.null(repository)){
    parsed_repo = parse_url(repository)
    (parsed_repo$scheme %in% c("http", "https")) |> assert_that()
    sync_remote_files(url = parsed_repo, cache_dir = cache_dir, files = files_to_read, subdirs = subdirs)
  }
  
  subdirs |>
    imap(function(current_subdir, current_assay){
    glue("Reading {length(files_to_read)} files.") |>
      message()
    
    # Load each file
    sces =
      files_to_read |>
      map(function(.x){
        cat(".")
        
        sce_path = file.path(
          cache_dir,
          current_subdir,
          .x
        )
        
        file.exists(sce_path) |>
          assert_that(
            msg="Your cache does not contain a file you attempted to 
            query. Please provide the repository parameter so that
            files can be synchronised from the internet"
          )
        
        sce = loadHDF5SummarizedExperiment(sce_path)
        
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
    assayNames(sce) = current_assay
    
    sce
  }) |>
    simplifyToSCE()
  

	# Return
	sce
}

#' Synchronises one or more remote files with a local copy
#'
#' @param url A character vector of length one. The base HTTP url from which to obtain the files.
#' @param cache_dir A character vector of length one. The local filepath to synchronise files to.
#' @param subdirs A character vector of subdirectories within the root URL to sync. These correspond to assays.
#' @param files A character vector containing one or more file_id_db entries
#'
#' @return A character vector of files that have been downloaded
#' @importFrom purrr pmap_chr transpose
#' @importFrom httr modify_url GET write_disk stop_for_status
#' @importFrom dplyr tibble transmute filter full_join
#' @importFrom glue glue
#' @importFrom assertthat assert_that
#' @export
#'
sync_remote_files = function(
  url,
  cache_dir,
  subdirs,
  files
){
  # Find every combination of file name, sample id, and assay, since each
  # will be a separate file we need to download
  expand.grid(
    filename = c("assays.h5", "se.rds"),
    sample_id = files,
    subdir = subdirs,
    stringsAsFactors = FALSE
  ) |>
    transmute(
      # Path to the file of interest from the root path. We use "/"
      # since URLs must use these regardless of OS
      full_url = paste0(url$path, "/",  subdir, "/", sample_id, "/", filename) |> map(~modify_url(url, path=.)),
      
      # Path to save the file on local disk (and its parent directory)
      # We use file.path since the file separator will differ on other OSs
      output_dir = file.path(
        cache_dir,
        subdir,
        sample_id
      ),
      output_file = file.path(
        output_dir,
        filename
      )
    ) |>
    filter(
      # Don't bother downloading files that don't exist
      # TODO: use some kind of hashing to check if the remote file has changed,
      # and proceed with the download if it has. However this is low importance
      # as the repository is not likely to change often
      !file.exists(output_file)
    ) |>
      pmap_chr(function(full_url, output_dir, output_file){
        dir.create(output_dir, recursive=TRUE, showWarnings = FALSE)
        glue("Downloading {full_url} to {output_file}") |> message()
        
        tryCatch(
          GET(full_url, write_disk(output_file)) |> stop_for_status(),
          error = function(e){
            # Clean up if we had an error
            file.remove(output_file)
            glue("File {full_url} could not be downloaded. {e}") |> stop()
          }
        )
        
        output_file
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
