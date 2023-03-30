# Functions that relate to downloading count data into SingleCellExperiments

# We need to load utils now so it can be used at the top level
#' @include utils.R
# This is a hack to force Seurat packages to be loaded, and also to
# satisfy R CMD check. We don't need to attach them at all.
#' @importFrom Seurat as.SingleCellExperiment
NULL

# Maps user provided assay names to their corresponding paths in the repository
assay_map <- c(
    counts = "original",
    cpm = "cpm"
)

#' Base URL pointing to the count data at the current version
COUNTS_URL <- single_line_str(
    "https://object-store.rc.nectar.org.au/v1/
    AUTH_06d6e008e3e642da99d806ba3ea629c5/cellxgene-0.2.1-hdf5"
)
#' Current version of the counts. This will be incremented when a newer
#' version is released
COUNTS_VERSION <- "0.2.1"

#' Gets a SingleCellExperiment from curated metadata
#'
#' Given a data frame of Curated Atlas metadata obtained from [get_metadata()],
#' returns a [`SingleCellExperiment::SingleCellExperiment-class`] object
#' corresponding to the samples in that data frame
#'
#' @param data A data frame containing, at minimum, a `sample_` column, which
#'   corresponds to a single cell sample ID. This can be obtained from the
#'   [get_metadata()] function.
#' @param assays A character vector whose elements must be either "counts"
#'   and/or "cpm", representing the corresponding assay(s) you want to request.
#'   By default only the count assay is downloaded. If you are interested in
#'   comparing a limited amount of genes, the "cpm" assay is more appropriate.
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
#' @importFrom SingleCellExperiment SingleCellExperiment combineCols
#' @importFrom SummarizedExperiment colData assayNames<-
#' @importFrom httr parse_url
#' @importFrom assertthat assert_that has_name
#' @importFrom cli cli_alert_success cli_alert_info
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @importFrom S4Vectors DataFrame
#' @export
get_SingleCellExperiment <- function(
    data,
    assays = "counts",
    cache_directory = get_default_cache_dir(),
    repository = COUNTS_URL,
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
    has_name(raw_data, c("cell_", "file_id_db")) |> assert_that()

    versioned_cache_directory <- file.path(cache_directory, COUNTS_VERSION)
    versioned_cache_directory |> dir.create(
        showWarnings = FALSE,
        recursive = TRUE
    )

    subdirs <- assay_map[assays]

    # The repository is optional. If not provided we load only from the cache
    if (!is.null(repository)) {
        cli_alert_info("Synchronising files")
        parsed_repo <- parse_url(repository)
        parsed_repo$scheme |>
            `%in%`(c("http", "https")) |>
            assert_that()

        files_to_read <-
            raw_data |>
            pull(.data$file_id_db) |>
            unique() |>
            as.character() |>
            sync_assay_files(
                url = parsed_repo,
                cache_dir = versioned_cache_directory,
                files = _,
                subdirs = subdirs
            )
    }

    cli_alert_info("Reading files.")
    sces <- subdirs |>
        imap(function(current_subdir, current_assay) {
            # Build up an SCE for each assay
            dir_prefix <- file.path(
                versioned_cache_directory,
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
#' @param i Suffix to be added to the column names, to make them unique
#' @param df The data frame to be converted
#' @param dir_prefix The path to the single cell experiment, minus the final
#'   segment
#' @param features The list of genes/rows of interest
#' @return A SingleCellExperiment object
#' @importFrom dplyr mutate filter
#' @importFrom HDF5Array loadHDF5SummarizedExperiment
#' @importFrom SummarizedExperiment colData<-
#' @importFrom tibble column_to_rownames
#' @importFrom utils head
#' @importFrom cli cli_alert_warning cli_abort
#' @importFrom glue glue
#' @importFrom stringr str_replace_all
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
    cells <- colnames(sce) |> intersect(df$cell_)

    if (length(cells) < nrow(df)){
        str_replace_all(
            "Some cells were filtered out because of extremely low counts. The
            number of cells in the SingleCellExperiment will be less than the
            number of cells you have selected from the metadata."
        )
        df <- filter(df, .data$cell_ %in% cells)
    }
    else if (length(cells) > nrow(df)){
        cli_abort("This should never happen")
    }
    
    # Fix for https://github.com/tidyverse/dplyr/issues/6746
    force(i)
    
    new_coldata <- df |>
        # We need to make the cell names globally unique, which we can guarantee
        # by adding a suffix that is derived from file_id_db, which is the
        # grouping variable
        mutate(original_cell_id = .data$cell_, cell_ = glue("{cell_}_{i}")) |>
        column_to_rownames("cell_") |>
        as("DataFrame")
    
    `if`(
            is.null(features),
            sce[, new_coldata$original_cell_id],
            {
                # Optionally subset the genes
                genes <- rownames(sce) |> intersect(features)
                sce[genes, new_coldata$original_cell_id]
            }
    ) |>
        `colnames<-`(new_coldata$cell_) |>
        `colData<-`(value = new_coldata)
}

#' Synchronises one or more remote assays with a local copy
#' @param url A character vector of length one. The base HTTP URL from which to
#'   obtain the files.
#' @param cache_dir A character vector of length one. The local filepath to
#'   synchronise files to.
#' @param subdirs A character vector of subdirectories within the root URL to
#'   sync. These correspond to assays.
#' @param files A character vector containing one or more file_id_db entries
#' @returns A character vector consisting of file paths to all the newly
#'   downloaded files
#' @return A character vector of files that have been downloaded
#' @importFrom purrr pmap_chr map_chr
#' @importFrom httr modify_url
#' @importFrom dplyr transmute filter
#' @noRd
#'
sync_assay_files <- function(
    url = parse_url(COUNTS_URL),
    cache_dir,
    subdirs,
    files
) {
    # Find every combination of file name, sample id, and assay, since each
    # will be a separate file we need to download
    files <- expand.grid(
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
            ) |> map_chr(~ modify_url(url, path = .)),

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
        )
    
    report_file_sizes(files$full_url)
    
    pmap_chr(files, function(full_url, output_dir, output_file) {
        sync_remote_file(full_url, output_file)
        output_file
    }, .progress = list(name = "Downloading files"))
}
