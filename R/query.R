# These are hacks to force the above packages to be loaded, and also to
# satisfy R CMD check. We don't need to attach them at all.
#' @import dbplyr 
#' @import Seurat
NULL

# Maps user provided assay names to their corresponding paths in the repository
assay_map <- c(
    counts = "original",
    cpm = "cpm"
)

#' Used in a pipeline to run one or more expressions with side effects, but
#' return the input value as the output value unaffected
aside <- function(x, ...) {
    # Courtesy of Hadley: https://fosstodon.org/@hadleywickham/109558265769090930
    list(...)
    x
}

REMOTE_URL <- "https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas"

#' Given a data frame of HCA metadata, returns a SingleCellExperiment object
#' corresponding to the samples in that data frame
#'
#' @param data A data frame containing, at minimum, a `.sample` column, which
#'   corresponds to a single cell sample ID. This can be obtained from the
#'   [get_metadata()] function.
#' @param assays A character vector whose elements must be either "raw" or
#'   "scaled", representing the corresponding assay you want to request.
#' @param repository A character vector of length one. If provided, it should be
#'   an HTTP URL pointing to the location where the single cell data is stored.
#' @param cache_directory An optional character vector of length one. If
#'   provided, it should indicate a local file path where any remotely accessed
#'   files should be copied.
#' @param features An optional character vector of features (ie genes) to return
#'   the counts for. By default counts for all features will be returned.
#' @returns A SingleCellExperiment object, with one assay for each value in the
#'   assays argument
#' @examples
#' meta <- get_metadata() |> head(2)
#' sce <- get_SingleCellExperiment(meta)
#'
#' @importFrom dplyr pull filter as_tibble inner_join collect
#' @importFrom tibble column_to_rownames
#' @importFrom purrr reduce map map_int imap keep
#' @importFrom BiocGenerics cbind
#' @importFrom glue glue
#' @importFrom HDF5Array loadHDF5SummarizedExperiment HDF5RealizationSink
#'   loadHDF5SummarizedExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment simplifyToSCE
#' @importFrom SummarizedExperiment colData assayNames<-
#' @importFrom httr parse_url
#' @importFrom assertthat assert_that has_name
#' @importFrom cli cli_alert_success cli_alert_info
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @importFrom S4Vectors DataFrame
#'
#' @export
#'
#'
get_SingleCellExperiment <- function(
    data,
    assays = c("counts", "cpm"),
    cache_directory = get_default_cache_dir(),
    repository = REMOTE_URL,
    features = NULL
) {
    # Parameter validation
    assays %in% names(assay_map) |>
        all() |>
        assert_that(
            msg = 'assays must be a character vector containing "counts" and/or
          "cpm"'
        )
    (!anyDuplicated(assays)) |> assert_that()
    inherits(cache_directory, "character") |> assert_that()
    is.null(repository) || is.character(repository) |> assert_that()
    is.null(features) || is.character(features) |> assert_that()

    # Data parameter validation (last, because it's slower)
    ## Evaluate the promise now so that we get a sensible error message
    data
    ## We have to convert to an in-memory table here, or some of the dplyr
    ## operations will fail when passed a database connection
    cli_alert_info("Realising metadata.")
    raw_data <- collect(data)
    inherits(raw_data, "tbl") |> assert_that()
    has_name(raw_data, c(".cell", "file_id_db")) |> assert_that()

    cache_directory |> dir.create(showWarnings = FALSE)

    cells_of_interest <- raw_data |>
        pull(.data$.cell) |>
        unique() |>
        as.character()

    subdirs <- assay_map[assays]

    # The repository is optional. If not provided we load only from the cache
    if (!is.null(repository)) {
        cli_alert_info("Synchronising files")
        files_to_read <-
            raw_data |>
            pull(.data$file_id_db) |>
            unique() |>
            as.character()
        parsed_repo <- parse_url(repository)
        (parsed_repo$scheme %in% c("http", "https")) |> assert_that()
        sync_assay_files(
            url = parsed_repo,
            cache_dir = cache_directory,
            files = files_to_read,
            subdirs = subdirs
        )
    }

    cli_alert_info("Reading files.")
    sces <- subdirs |>
        imap(function(current_subdir, current_assay) {
            # Build up an SCE for each assay
            dir_prefix <- file.path(
                cache_directory,
                current_subdir
            )

            raw_data |>
                dplyr::group_by(.data$file_id_db) |>
                # Load each file and attach metadata
                dplyr::summarise(sces = list(group_to_sce(
                    dplyr::cur_group_id(),
                    dplyr::cur_data_all(),
                    dir_prefix,
                    features
                ))) |>
                dplyr::pull(sces) |>
                # Combine each sce by column, since each sce has a different set
                # of cells
                do.call(cbind, args = _)
        })

    cli_alert_info("Compiling Single Cell Experiment.")
    # Combine all the assays
    sce <- sces[[1]]
    SummarizedExperiment::assays(sce) <- map(sces, function(sce) {
        SummarizedExperiment::assays(sce)[[1]]
    })

    sce
}

#' Converts a data frame into a single SCE
#'
#' @param i Suffix to be added to the column names, to make them unique
#' @param df The data frame to be converted
#' @param dir_prefix The path to the single cell experiment, minus the final segment
#' @param features The list of genes/rows of interest
#' @return A SingleCellExperiment object
#' @importFrom dplyr mutate
#' @importFrom HDF5Array loadHDF5SummarizedExperiment
#' @importFrom SummarizedExperiment colData<-
#' @importFrom tibble column_to_rownames
#' @importFrom utils head
#' @noRd
group_to_sce <- function(i, df, dir_prefix, features) {
    sce_path <- df$file_id_db |>
        head(1) |>
        file.path(
            dir_prefix,
            suffix = _
        )

    file.exists(sce_path) |>
        assert_that(
            msg = "Your cache does not contain a file you
                            attempted to query. Please provide the repository
                            parameter so that files can be synchronised from the
                            internet"
        )

    sce <- loadHDF5SummarizedExperiment(sce_path)
    # The cells we select here are those that are both available in the SCE
    # object, and requested for this particular file
    cells <- colnames(sce) |> intersect(df$.cell)
    # We need to make the cell names globally unique, which we can guarantee
    # by adding a suffix that is derived from file_id_db, which is the grouping
    # variable
    new_cellnames <- paste0(cells, "_", i)
    new_coldata <- df |>
        mutate(original_cell_id = .data$.cell, .cell = new_cellnames) |>
        column_to_rownames(".cell") |>
        as("DataFrame")

    features |>
        is.null() |>
        {
            `if`
        }(
            sce[, cells], {
                # Optionally subset the genes
                genes <- rownames(sce) |> intersect(features)
                sce[genes, cells]
        }) |>
        `colnames<-`(new_cellnames) |>
        `colData<-`(value = new_coldata)
}

#' Synchronises one or more remote assays with a local copy
#'
#' @param url A character vector of length one. The base HTTP URL from which to
#'   obtain the files.
#' @param cache_dir A character vector of length one. The local filepath to
#'   synchronise files to.
#' @param subdirs A character vector of subdirectories within the root URL to
#'   sync. These correspond to assays.
#' @param files A character vector containing one or more file_id_db entries
#' @returns A character vector consisting of file paths to all the newly
#'   downloaded files
#'
#' @return A character vector of files that have been downloaded
#' @importFrom purrr pmap_chr transpose
#' @importFrom httr modify_url GET write_disk stop_for_status parse_url
#' @importFrom dplyr tibble transmute filter full_join
#' @importFrom glue glue
#' @importFrom assertthat assert_that
#' @importFrom cli cli_alert_success cli_alert_info cli_abort
#' @noRd
#'
sync_assay_files <- function(
    url = parse_url(REMOTE_URL),
    cache_dir,
    subdirs,
    files
) {
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
            full_url = paste0(
                url$path,
                "/",
                .data$subdir,
                "/",
                .data$sample_id,
                "/",
                .data$filename
            ) |> map(~ modify_url(url, path = .)),

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
            # Don't bother downloading files that don't exist TODO: use some
            # kind of hashing to check if the remote file has changed, and
            # proceed with the download if it has. However this is low
            # importance as the repository is not likely to change often
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
#' @importFrom tools R_user_dir
#' @importFrom utils packageName
#' @noRd
#'
get_default_cache_dir <- function() {
    packageName() |>
        R_user_dir(
            "cache"
        ) |>
        normalizePath()
}

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

#' Downloads a parquet database of the Human Cell Atlas metadata to a local 
#' cache, and then opens it as a data frame. It can then be filtered and 
#' passed into [get_SingleCellExperiment()] 
#' to obtain a [`SingleCellExperiment`](SingleCellExperiment::SingleCellExperiment-class)
#'
#' @param remote_url Optional character vector of length 1. An HTTP URL pointing
#'   to the location of the parquet database.
#' @param cache_directory Optional character vector of length 1. A file path on
#'   your local system to a directory (not a file) that will be used to store
#'   metadata.parquet
#' @return A lazy data.frame subclass containing the metadata. You can interact
#'   with this object using most standard dplyr functions. For string matching,
#'   it is recommended that you use `stringr::str_like` to filter character
#'   columns, as `stringr::str_match` will not work.
#' @export
#' @examples
#' library(dplyr)
#' filtered_metadata <- get_metadata() |>
#'     filter(
#'         ethnicity == "African" &
#'             assay %LIKE% "%10x%" &
#'             tissue == "lung parenchyma" &
#'             cell_type %LIKE% "%CD4%"
#'     )
#'
#' @importFrom DBI dbConnect
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl
#' @importFrom httr progress
#' @importFrom cli cli_alert_info
#' @importFrom utils untar
get_metadata <- function(
    remote_url = "https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/metadata/metadata.0.2.2.parquet",
    cache_directory = get_default_cache_dir()
) {
    db_path <- file.path(cache_directory, "metadata.parquet")
    sync_remote_file(
        remote_url,
        db_path,
        progress(type = "down", con = stderr())
    )
    table <- duckdb() |>
        dbConnect(drv = _, read_only = TRUE) |>
        tbl(db_path)
}
