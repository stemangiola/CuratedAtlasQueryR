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

#' Base URL pointing to the pseudobulk counts at the current version
#' @noRd
pseudobulk_url <- single_line_str(
  "https://object-store.rc.nectar.org.au/v1/
  AUTH_06d6e008e3e642da99d806ba3ea629c5/pseudobulk-0.1.0"
)

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
#' @inheritDotParams get_summarized_experiment
#' @examples
#' meta <- get_metadata() |> head(2)
#' sce <- get_single_cell_experiment(meta)
#' @export
get_single_cell_experiment <- function(...){
  get_summarized_experiment(..., type = "sce")
}

#' Gets a Pseudobulk from curated metadata
#' 
#' @inheritDotParams get_summarized_experiment
#' @examples
#' \dontrun{
#' meta <- get_metadata() |> filter(tissue_harmonised == "lung")
#' pseudobulk <- meta |> get_pseudobulk()
#' }
#' @export
get_pseudobulk <- function(...) {
  get_summarized_experiment(..., type = "pseudobulk")
}


#' Gets a Summarized Experiment from curated metadata
#' 
#' Given a data frame of Curated Atlas metadata obtained from [get_metadata()],
#' returns a [`SummarizedExperiment::SummarizedExperiment-class`] object
#' corresponding to the samples in that data frame
#' @inheritParams param_validation
#' @param type A character vector of length one. Either sce for Single Cell Experiment,
#' or pseudobulk.
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
get_summarized_experiment <- function(
    data,
    assays = "counts",
    cache_directory = get_default_cache_dir(),
    repository = COUNTS_URL,
    features = NULL,
    type = "sce"
) {
  repository <- if (type == "sce") COUNTS_URL else pseudobulk_url 
  validated <- param_validation(data, assays, cache_directory, repository, features)
  
  # Extract variables from validation
  raw_data <- validated$raw_data
  versioned_cache_directory <- validated$versioned_cache_directory
  subdirs <- validated$subdirs
  
  file_id_col <- if (type == "sce") "file_id_db" else "file_id"
  
  # The repository is optional. If not provided we load only from the cache
  if (!is.null(repository)) {
    cli_alert_info("Synchronising files")
    parsed_repo <- parse_url(repository)
    parsed_repo$scheme |>
      `%in%`(c("http", "https")) |>
      assert_that()
    
    files_to_read <-
      raw_data |>
      pull(file_id_col) |>
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
  experiments <- subdirs |>
    imap(function(current_subdir, current_assay) {
      # Build up an experiment for each assay
      dir_prefix <- file.path(
        versioned_cache_directory,
        current_subdir
      )
      
      experiment_list <- raw_data |>
        dplyr::group_by(.data[[file_id_col]]) |>
        dplyr::summarise(experiments = list(
          group_to_sme(
            dplyr::cur_group_id(),
            dplyr::cur_data_all(),
            dir_prefix,
            features,
            type
          )
        )) |>
        dplyr::pull(experiments)
      
      commonGenes <- experiment_list |> check_gene_overlap()
      experiment_list <- map(experiment_list, function(exp) {
        exp[commonGenes,]
      }) |>
        do.call(cbind, args = _)
    })
  
  cli_alert_info("Compiling Experiment.")
  # Combine all the assays
  experiment <- experiments[[1]]
  SummarizedExperiment::assays(experiment) <- map(experiments, function(exp) {
    SummarizedExperiment::assays(exp)[[1]]
  })
  
  experiment
}

#' Converts a data frame into a Summarized Experiment
#' @param i Suffix to be added to the column names, to make them unique
#' @param df The data frame to be converted
#' @param dir_prefix The path to the single cell experiment, minus the final
#'   segment
#' @param features The list of genes/rows of interest
#' @return A Summarized Experiment object
#' @importFrom dplyr mutate filter
#' @importFrom HDF5Array loadHDF5SummarizedExperiment
#' @importFrom SummarizedExperiment colData<-
#' @importFrom tibble column_to_rownames
#' @importFrom utils head
#' @importFrom cli cli_alert_warning cli_abort
#' @importFrom glue glue
#' @importFrom stringr str_replace_all
#' @noRd
group_to_sme <- function(i, df, dir_prefix, features, type = "sce") {
  # Set file name based on type
  file_id_col <- if (type == "sce") "file_id_db" else "file_id"
  experiment_path <- df[[file_id_col]] |>
    head(1) |>
    file.path(
      dir_prefix,
      suffix=_
    )
  
  # Check if file exists
  file.exists(experiment_path) |>
    assert_that(
      msg = "Your cache does not contain the file {experiment_path} you
                            attempted to query. Please provide the repository
                            parameter so that files can be synchronised from the
                            internet" |> glue()
    )
  
  # Load experiment
  experiment <- loadHDF5SummarizedExperiment(experiment_path)
  
  # Fix for https://github.com/tidyverse/dplyr/issues/6746
  force(i)
  
  if (type == "sce") {
    # Process specific to SCE
    cells <- colnames(experiment) |> intersect(df$cell_)
    
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
    
    new_coldata <- df |>
      mutate(original_cell_id = .data$cell_, cell_ = glue("{cell_}_{i}")) |>
      column_to_rownames("cell_") |>
      as("DataFrame")
    
    experiment <- `if`(
      is.null(features),
      experiment[, new_coldata$original_cell_id],
      {
        # Optionally subset the genes
        genes <- rownames(experiment) |> intersect(features)
        experiment[genes, new_coldata$original_cell_id]
      }
    ) |>
      `colnames<-`(new_coldata$cell_) |>
      `colData<-`(value = new_coldata)
  }
  else if (type == "pseudobulk") {
    # remove cell-level annotations
    cell_level_anno <- c("cell_", "cell_type", "confidence_class", "file_id_db",
                         "cell_annotation_blueprint_singler",
                         "cell_annotation_monaco_singler", 
                         "cell_annotation_azimuth_l2")
    
    new_coldata <- df |> select(-dplyr::all_of(cell_level_anno)) |> distinct() |>
      mutate(sample_identifier = glue("{sample_}___{cell_type_harmonised}"),
             original_sample_id = .data$sample_identifier) |>
      column_to_rownames("original_sample_id")

    experiment <- `if`(
      is.null(features),
      experiment[, new_coldata$sample_identifier] |>
        `colnames<-`(new_coldata$sample_identifier) |>
        `colData<-`(value = DataFrame(new_coldata))
    )
    experiment
  }
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

#' Checks whether genes in a list of SummarizedExperiment objects overlap
#' @param obj_list A list of SummarizedExperiment objects
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
#' @param data A data frame containing, at minimum, `cell_`, `file_id_db`, `file_id` column, which
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
              has_name(raw_data, c("file_id", "cell_", "file_id_db")))
  
  versioned_cache_directory <- cache_directory
  versioned_cache_directory |> dir.create(showWarnings = FALSE,
                                          recursive = TRUE)
  
  subdirs <- assay_map[assays]
  list(raw_data = raw_data, versioned_cache_directory = versioned_cache_directory, 
       subdirs = subdirs)
}

