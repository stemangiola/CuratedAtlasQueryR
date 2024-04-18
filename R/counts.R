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
#' @noRd
COUNTS_URL <- single_line_str(
    "https://object-store.rc.nectar.org.au/v1/
    AUTH_06d6e008e3e642da99d806ba3ea629c5/cellxgene-0.2.1-hdf5"
)
#' Current version of the counts. This will be incremented when a newer
#' version is released
#' @noRd
COUNTS_VERSION <- "0.2.1"

#' @inherit get_single_cell_experiment
#' @inheritDotParams get_single_cell_experiment
#' @importFrom cli cli_alert_warning
#' @export
get_SingleCellExperiment <- function(...){
    single_line_str("This function name is deprecated. 
                    Please use `get_single_cell_experiment()` instead") |>
        cli_alert_warning()
    
    get_single_cell_experiment(...)
}

#' Gets a SingleCellExperiment from curated metadata
#'
#' Given a data frame of Curated Atlas metadata obtained from [get_metadata()],
#' returns a [`SingleCellExperiment::SingleCellExperiment-class`] object
#' corresponding to the samples in that data frame
#' @inheritParams param_validation
#' @returns A SingleCellExperiment object, with one assay for each value in the
#'   assays argument
#' @examples
#' meta <- get_metadata() |> head(2)
#' sce <- get_single_cell_experiment(meta)
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
#' @importFrom S4Vectors DataFrame
#' @export
get_single_cell_experiment <- function(
    data,
    assays = "counts",
    cache_directory = get_default_cache_dir(),
    repository = COUNTS_URL,
    features = NULL
) {
  
  validated <- param_validation(data, assays, cache_directory, repository, features)
  
  # Extract variables from validation
  raw_data <- validated$raw_data
  versioned_cache_directory <- validated$versioned_cache_directory
  subdirs <- validated$subdirs
  
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
      
      sce_list <- raw_data |>
        dplyr::group_by(.data$file_id_db) |>
        # Load each file and attach metadata
        dplyr::summarise(sces = list(
          group_to_sce(
            dplyr::cur_group_id(),
            dplyr::cur_data_all(),
            dir_prefix,
            features
          )
        )) |>
        dplyr::pull(sces)
      
      # Check whether genes in a list of SCEs overlap, use gene intersection if overlap
      commonGenes <- sce_list |> check_gene_overlap()
      sce_list <- map(sce_list, function(sce) {
        sce[commonGenes,]
      }) |>
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
#' Gets a Summarized Experiment from curated metadata
#'
#' Given a data frame of Curated Atlas metadata obtained from [get_metadata()],
#' returns a [`SummarizedExperiment::SummarizedExperiment-class`] object
#' corresponding to the samples in that data frame
#' @inheritParams param_validation
#' @examples
#' \dontrun{
#' meta <- get_metadata() |> 
#'         filter(sample_ %in% c('068502277538ef5559154b543167fefa',
#'         '07605571f51519f71da03704f056fa43'))
#' sme <- get_pseudobulk(meta, repository=NULL,cache_directory = "~/projects/caq/pseudobulk/0.2.1")
#' }
#'
#' @importFrom dplyr pull filter as_tibble inner_join collect
#' @importFrom tibble column_to_rownames
#' @importFrom purrr reduce map map_int imap keep
#' @importFrom BiocGenerics cbind
#' @importFrom glue glue
#' @importFrom SummarizedExperiment colData assayNames<-
#' @importFrom httr parse_url
#' @importFrom assertthat assert_that has_name
#' @importFrom cli cli_alert_success cli_alert_info
#' @importFrom rlang .data
#' @importFrom S4Vectors DataFrame
#' @export
get_pseudobulk <- function(
    data,
    assays = "counts",
    cache_directory = get_default_cache_dir(),
    repository = NULL, # cloud container not ready yet
    features = NULL
) {
  
  validated <- param_validation(data, assays, cache_directory, repository, features)
  
  # Extract variables from validated parameters
  raw_data <- validated$raw_data
  versioned_cache_directory <- validated$versioned_cache_directory
  subdirs <- validated$subdirs
  
  # The repository is optional. If not provided we load only from the cache
  if (!is.null(repository)) {
    cli_alert_info("Synchronising files")
    parsed_repo <- parse_url(repository)
    parsed_repo$scheme |>
      `%in%`(c("http", "https")) |>
      assert_that()
    
    files_to_read <-
      raw_data |>
      pull(.data$sample_) |>
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
  smes <- subdirs |>
    imap(function(current_subdir, current_assay) {
      # Build up a Summarized Experiment for each assay
      dir_prefix <- file.path(
        versioned_cache_directory,
        current_subdir
      )
      
      sme_list <- raw_data |>
        dplyr::group_by(.data$sample_) |>
        # Load each file and attach metadata
        dplyr::summarise(smes = list(
          group_to_sme(
            dplyr::cur_group_id(),
            dplyr::cur_data_all(),
            dir_prefix,
            features
          )
        )) |>
        dplyr::pull(smes)
      
      # Check whether genes/rows in a list of SummarizedExperiment overlap, use gene intersection if does
      commonGenes <- sme_list |> check_gene_overlap()
      sme_list <- map(sme_list, function(sme) {
        sme[commonGenes,]
      }) |>
        do.call(cbind, args = _)
    })
  
  cli_alert_info("Compiling Summarized Experiment.")
  # Combine all the assays
  
  sme <- smes[[1]]
  SummarizedExperiment::assays(sme) <- map(smes, function(sme) {
    SummarizedExperiment::assays(sme)[[1]]
  })
  
  sme
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
            msg = "Your cache does not contain a file {sce_path} you
                            attempted to query. Please provide the repository
                            parameter so that files can be synchronised from the
                            internet" |> glue()
        )

    sce <- loadHDF5SummarizedExperiment(sce_path)
    # The cells we select here are those that are both available in the SCE
    # object, and requested for this particular file
    cells <- colnames(sce) |> intersect(df$cell_)

    if (length(cells) < nrow(df)){
        single_line_str(
            "The number of cells in the SingleCellExperiment will be less than the
            number of cells you have selected from the metadata.
            Are cell IDs duplicated? Or, do cell IDs correspond to the counts file?
            "
        ) |> cli_alert_warning()
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

#' Converts a data frame into a Summarized Experiment
#' @inheritParams group_to_sce
#' @return A Summarized Experiment object
#' @importFrom utils head
#' @noRd
group_to_sme <- function(i, df, dir_prefix, features) {
  sme_path <- df$sample_ |>
    head(1) |>
    paste0(".rds") |>
    file.path(
      dir_prefix,
      suffix=_
    )
  
  file.exists(sme_path) |>
    assert_that(
      msg = "Your cache does not contain a file {sme_path} you
                            attempted to query. Please provide the repository
                            parameter so that files can be synchronised from the
                            internet" |> glue()
    )
  
  sme <- readRDS(sme_path)
  force(i)
  
  sme
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

#' Checks whether genes in a list of SingleCellExperiment or SummarizedExperiment objects overlap
#' @param obj_list A list of SingleCellExperiment or SummarizedExperiment objects
#' @return A character vector of genes intersection across objects
#' @importFrom purrr map reduce
#' @importFrom cli cli_alert_warning
#' @noRd
check_gene_overlap <- function(obj_list) {
  gene_lists <- map(obj_list, rownames)
  common_genes <- reduce(gene_lists, intersect)
  if (any(lengths(gene_lists) != length(common_genes))) {
    single_line_str(
      "CuratedAtlasQuery says: Not all genes completely overlap across the provided objects.
      Counts are generated by genes intersection"
    ) |> cli_alert_warning()
  }
  
  common_genes
}

#' Validate parameters for Summarized Experiment analysis
#' @param data A data frame containing, at minimum, `cell_`, `file_id_db`, `sample_` column, which
#'   correspond to a single cell ID, file subdivision for internal use, and a single cell sample ID. 
#'   They can be obtained from the [get_metadata()] function.
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
#' @return A list of elements:
#' \itemize{
#'  \item{raw_data}{Data after being converted to an in-memory data frame }
#'  \item{versioned_cache_directory}{The path to the cache directory}
#'  \item{subdirs}{Vector of subdirectory names from the `assays` input} 
#' }
#' @importFrom dplyr collect
#' @importFrom assertthat assert_that has_name
#' @importFrom cli cli_alert_info
#' @importFrom rlang .data
#' @keywords internal
param_validation <- function(data,
                             assays,
                             cache_directory,
                             repository,
                             features
                             ) {
  # Parameter validation 
  assays %in% names(assay_map) |>
    all() |>
    assert_that(msg = 'assays must be a character vector containing "counts" and/or
            "cpm"')
  assert_that(
    !anyDuplicated(assays),
    inherits(cache_directory, "character"),
    is.null(repository) || is.character(repository),
    is.null(features) || is.character(features)
  )
  
  # Data parameter validation (last, because it's slower)
  ## Evaluate the promise now so that we get a sensible error message
  force(data)
  ## We have to convert to an in-memory table here, or some of the dplyr
  ## operations will fail when passed a database connection
  cli_alert_info("Realising metadata.")
  raw_data <- collect(data)
  assert_that(inherits(raw_data, "tbl"),
              has_name(raw_data, c("sample_", "cell_", "file_id_db")))
  
  versioned_cache_directory <- cache_directory
  versioned_cache_directory |> dir.create(showWarnings = FALSE,
                                          recursive = TRUE)
  
  subdirs <- assay_map[assays]
  
  list(raw_data = raw_data, versioned_cache_directory = versioned_cache_directory, subdirs = subdirs)
}

