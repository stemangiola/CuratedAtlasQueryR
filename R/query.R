# Maps user provided assay names to their corresponding paths in the repository
assay_map <- c(
  counts = "original",
  cpm = "cpm"
)

# ' Used in a pipeline to run one or more expressions with side effects, but
# ' return the input value as the output value unaffected
aside <- function(x, ...) {
  # Courtesy of Hadley: https://fosstodon.org/@hadleywickham/109558265769090930
  list(...)
  x
}

REMOTE_URL <- "https://harmonised-human-atlas.s3.amazonaws.com/"

#' Given a data frame of HCA metadata, returns a SingleCellExperiment object corresponding to the samples in that data frame
#'
#' @param data A data frame containing, at minimum, a `.sample` column, which corresponds to a single cell sample ID.
#' This can be obtained from the [get_metadata()] function.
#' @param assays A character vector whose elements must be either "raw" or "scaled", representing the corresponding assay you want to request.
#' @param repository A character vector of length one. If provided, it should be an HTTP URL pointing to the location where the single cell data is stored.
#' @param cache_directory An optional character vector of length one. If provided, it should indicate a local file path where any remotely accessed files should be copied.
#' @param features An optional character vector of features (ie genes) to return the counts for. By default counts for all features will be returned.
#' @returns A SingleCellExperiment object, with one assay for each value in the assays argument
#' @examples
#' meta <- get_metadata() |> head(2)
#' sce <- get_SingleCellExperiment(meta)
#'
#' @importFrom dplyr pull filter as_tibble
#' @importFrom tidySingleCellExperiment inner_join
#' @importFrom purrr reduce map map_int imap keep
#' @importFrom BiocGenerics cbind
#' @importFrom glue glue
#' @importFrom HDF5Array loadHDF5SummarizedExperiment HDF5RealizationSink loadHDF5SummarizedExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment simplifyToSCE
#' @importFrom SummarizedExperiment colData assayNames<-
#' @importFrom httr parse_url
#' @importFrom assertthat assert_that has_name
#' @importFrom cli cli_alert_success cli_alert_info
#' @importFrom rlang .data
#' @importFrom stats setNames
#'
#' @export
#'
#'
get_SingleCellExperiment <- function(data,
                                     assays = c("counts", "cpm"),
                                     cache_directory = get_default_cache_dir(),
                                     repository = REMOTE_URL,
                                     features = NULL) {
  # Parameter validation
  assays %in% names(assay_map) |>
    all() |>
    assert_that(msg = 'assays must be a character vector containing "counts" and/or "cpm"')
  (!anyDuplicated(assays)) |> assert_that()
  inherits(cache_directory, "character") |> assert_that()
  is.null(repository) || is.character(repository) |> assert_that()
  is.null(features) || is.character(features) |> assert_that()

  # Data parameter validation (last, because it's slower)
  ## Evaluate the promise now so that we get a sensible error message
  data
  ## We have to convert to an in-memory table here, or some of the dplyr operations will fail when passed a database connection
  cli_alert_info("Realising metadata.")
  raw_data <- as_tibble(data)
  inherits(raw_data, "tbl") |> assert_that()
  has_name(raw_data, c(".cell", "file_id_db")) |> assert_that()

  cache_directory |> dir.create(showWarnings = FALSE)

  files_to_read <-
    raw_data |>
    pull(.data$file_id_db) |>
    unique() |>
    as.character()

  subdirs <- assay_map[assays]

  # The repository is optional. If not provided we load only from the cache
  if (!is.null(repository)) {
    cli_alert_info("Synchronising files")
    parsed_repo <- parse_url(repository)
    (parsed_repo$scheme %in% c("http", "https")) |> assert_that()
    sync_assay_files(url = parsed_repo, cache_dir = cache_directory, files = files_to_read, subdirs = subdirs)
  }
  files_to_read <-
    raw_data |>
    pull(.data$file_id_db) |>
    unique() |>
    as.character()

  subdirs |>
    imap(function(current_subdir, current_assay) {
      # Load each file
      sces <-
        files_to_read |>
        map(function(.x) {
          sce_path <- file.path(
            cache_directory,
            current_subdir,
            .x
          )

          file.exists(sce_path) |>
            assert_that(
              msg = "Your cache does not contain a file you attempted to
            query. Please provide the repository parameter so that
            files can be synchronised from the internet"
            )

          sce <- loadHDF5SummarizedExperiment(sce_path)

          if (!is.null(features)) {
            # Optionally subset the genes
            sce <- sce[
              rownames(sce) |> intersect(features)
            ]
          }

          sce
        }, .progress = list(name = "Reading files")) |>
        # Drop files with one cell, which causes
        # the DFrame objects to combine must have the same column names
        keep(~ ncol(.) > 1) |>
        # Combine each sce by column, since each sce has a different set of cells
        do.call(cbind, args = _) |>
        # We only need the assay, since we ultimately need to combine them
        # We need to use :: here since we already have an assays argument
        SummarizedExperiment::assays() |>
        setNames(current_assay)
    }) |>
    aside(cli_alert_info("Compiling Single Cell Experiment.")) |>
    # Combine the assays into one list
    reduce(c) |>
    SingleCellExperiment(assays = _) |>
    aside(cli_alert_info("Attaching metadata.")) |>
    # Join back to metadata, which will become coldata annotations
    inner_join(
      # Needed because cell IDs are not unique outside the file_id or file_id_db
      filter(raw_data, .data$file_id_db %in% files_to_read),
      by = ".cell"
    )
}

#' Synchronises one or more remote assays with a local copy
#'
#' @param url A character vector of length one. The base HTTP URL from which to obtain the files.
#' @param cache_dir A character vector of length one. The local filepath to synchronise files to.
#' @param subdirs A character vector of subdirectories within the root URL to sync. These correspond to assays.
#' @param files A character vector containing one or more file_id_db entries
#' @returns A character vector consisting of file paths to all the newly
#' downloaded files
#'
#' @return A character vector of files that have been downloaded
#' @importFrom purrr pmap_chr transpose
#' @importFrom httr modify_url GET write_disk stop_for_status
#' @importFrom dplyr tibble transmute filter full_join
#' @importFrom glue glue
#' @importFrom assertthat assert_that
#' @importFrom cli cli_alert_success cli_alert_info cli_abort
#' @noRd
#'
sync_assay_files <- function(url = httr::parse_url(REMOTE_URL),
                             cache_dir,
                             subdirs,
                             files) {
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
      full_url = paste0(url$path, "/", .data$subdir, "/", .data$sample_id, "/", .data$filename) |> map(~ modify_url(url, path = .)),

      # Path to save the file on local disk (and its parent directory)
      # We use file.path since the file separator will differ on other OSs
      output_dir = file.path(
        cache_dir,
        .data$subdir,
        .data$sample_id
      ),
      output_file = file.path(
        .data$output_dir,
        .data$filename
      )
    ) |>
    filter(
      # Don't bother downloading files that don't exist
      # TODO: use some kind of hashing to check if the remote file has changed,
      # and proceed with the download if it has. However this is low importance
      # as the repository is not likely to change often
      !file.exists(.data$output_file)
    ) |>
    pmap_chr(function(full_url, output_dir, output_file) {
      sync_remote_file(full_url, output_file)
      output_file
    }, .progress = list(name = "Downloading files"))
}

#' Synchronises a single remote file with a local path
#' @noRd
sync_remote_file <- function(full_url, output_file, ...) {
  if (!file.exists(output_file)) {
    output_dir <- dirname(output_file)
    dir.create(output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    cli_alert_info("Downloading {full_url} to {output_file}")

    tryCatch(
      GET(full_url, write_disk(output_file), ...) |> stop_for_status(),
      error = function(e) {
        # Clean up if we had an error
        file.remove(output_file)
        cli_abort("File {full_url} could not be downloaded. {e}")
      }
    )
  }
}

#' Returns the default cache directory
#'
#' @return A length one character vector.
#' @export
#' @importFrom rappdirs user_cache_dir
#'
get_default_cache_dir <- function() {
  file.path(
    user_cache_dir(),
    "hca_harmonised"
  )
}

#' @importFrom SeuratObject as.sparse
#' @importFrom assertthat assert_that
#' @importFrom methods as
#' @exportS3Method
as.sparse.DelayedMatrix <- function(x) {
  # This is glue to ensure the SCE -> Seurat conversion works properly with
  # DelayedArray types
  as(x, "dgCMatrix")
}

#' Given a data frame of HCA metadata, returns a Seurat object corresponding to the samples in that data frame
#'
#' @inheritDotParams get_SingleCellExperiment
#' @importFrom Seurat as.Seurat
#' @export
#' @return A Seurat object containing the same data as a call to get_SingleCellExperiment.
#' @examples
#' meta <- get_metadata() |> head(2)
#' seurat <- get_seurat(meta, repository = Sys.getenv("REMOTE_HCA"), cache_directory)
#'
get_seurat <- function(...) {
  get_SingleCellExperiment(...) |> as.Seurat(data = NULL)
}


#' Returns a data frame of Human Cell Atlas metadata, which should be filtered
#' and ultimately passed into get_SingleCellExperiment.
#'
#' @param repository Optional character vector of length 1. An HTTP URL pointing to the
#' location of the sqlite database.
#' @param cache_directory Optional character vector of length 1. A file path on
#' your local system to a directory (not a file) that will be used to store
#' metadata.sqlite
#' @return A lazy data.frame subclass containing the metadata. You can interact
#' with this object using most standard dplyr functions. However, it is recommended
#' that you use the %LIKE% operator for string matching, as most stringr functions
#' will not work.
#' @export
#' @examples
#' filtered_metadata <- get_metadata() |>
#'   filter(
#'     ethnicity == "African" &
#'       assay %LIKE% "%10x%" &
#'       tissue == "lung parenchyma" &
#'       cell_type %LIKE% "%CD4%"
#'   )
#'
#' @importFrom DBI dbConnect
#' @importFrom RSQLite SQLite SQLITE_RO
#' @importFrom dplyr tbl
#' @importFrom httr progress
#'
get_metadata <- function(repository = "https://harmonised-human-atlas.s3.amazonaws.com/metadata.sqlite",
                         cache_directory = get_default_cache_dir()) {
  sqlite_path <- file.path(cache_directory, "metadata.sqlite")
  sync_remote_file(repository, sqlite_path, progress(type = "down", con = stderr()))
  SQLite() |>
    dbConnect(drv = _, dbname = sqlite_path, flags = SQLITE_RO) |>
    tbl("metadata")
}
