#' Checks importing criteria for new atlas 
#' Generating counts per million from provided raw counts
#'
#' @param metadata_tbl A tibble
#' @param cache_dir Optional character vector of length 1. A file path on
#'   your local system to a directory (not a file) that will be used to store
#'   `metadata.parquet`
#' @param counts_path_df A data.frame RDS of original counts defined by users
#' @export
#' @return A metadata.parquet with version, and a directory stores counts per million in cache directory
#' @examples
#' example_metadata <- get_metadata() |> head(3) |> as_tibble()
#' import_metadata_counts(metadata_tbl = example_metadata,
#'                        cache_dir = get_default_cache_dir(),
#'                        counts_path_df = "/Users/shen.m/projects/caq/my_counts_df.rds")
#'
#' @importFrom assertthat assert_that
#' @importFrom checkmate check_tibble check_directory_exists check_set_equal check_true check_character check_subset check_file_exists
#' @importFrom dplyr tbl
#' @importFrom cli cli_alert_info
#' @importFrom glue glue
#' @importFrom arrow write_parquet
#' @importFrom fs dir_copy
#' @importFrom purrr walk2
import_metadata_counts <- function(metadata_tbl,
                                   cache_dir = CuratedAtlasQueryR:::get_default_cache_dir(),
                                   counts_path_df) {
  # Load counts_path dataframe from RDS
  counts_path <- readRDS(counts_path_df)
  
  # Convert to tibble if metadata_tbl is not a tibble 
  metadata_tbl <- metadata_tbl |> as_tibble()
  check_directory_exists(cache_dir)
  check_directory_exists(counts_path)
  
  
  
  # create original and cpm folders in cache directory if not exist (in order to append new counts to existing ones)
  if (!dir.exists(file.path(cache_dir, "original"))) {
    dir.create(cache_dir, "original", recursive = TRUE)
  }
  
  if (!dir.exists(file.path(cache_dir, "cpm"))) {
    dir.create(cache_dir, "cpm", recursive = TRUE)
  }
  
  # check count H5 directory name not included in the cache directory original
  all(!dir(file.path(counts_path)) %in% dir(file.path(cache_dir, 'original'))) |>
    check_true() |>
    assert_that(msg = 'Count H5 directory name should not duplicate with that in cache directory')
  
  # if checkpoints above pass, generate cpm
  cli_alert_info("Generating cpm from {.path {counts_path$file_path}}. ")
  
  counts_path <-
    counts_path |> mutate(
      original_path = file.path(original_dir, basename(file_id)),
      cpm_path = file.path(cache_dir, "cpm", basename(file_id))
    )
  
  purrr::walk2(counts_path$file_path, counts_path$cpm_path, ~{
    get_counts_per_million(
      input_file_rds = .x,
      output_file = .y
    )
    dir_copy(.x |> dirname(), file.path(cache_dir, "original"))
    data = readRDS(.x)
    saveHDF5SummarizedExperiment(data, counts_path$original_path, replace=TRUE)
  })
  
  cli_alert_info("cpm are generated in {.path {cache_dir}/cpm}. ")
  
  sce <- metadata_tbl |>
    get_single_cell_experiment(cache_directory = cache_dir)
  se <- metadata_tbl |> get_seurat(cache_directory = cache_dir)
  
  counts_data <- assays(sce)$counts |>
    as_tibble()
  assert_that(inherits(sce, "SingleCellExperiment"),
              msg = "SingleCellExperiment Object is not created from metadata.")
  assert_that(!'^ENS' %in% rownames(se[['originalexp']]),
              msg = "Gene names cannot contain Ensembl IDs.")
  assert_that(all(counts_data >= 0),
              msg = "Counts for SingleCellExperiment cannot be negative.")
  
  # check metadata sample file ID match the count file ID
  all(metadata_tbl |> pull(file_id_db) %in% dir(file.path(cache_dir, "original"))) |> 
    assert_that(msg = "The metadata sample file ID and the count file ID does not match")
  
  # check the number of sub directories in original match cpm 
  check_set_equal(length(list.dirs(file.path(cache_dir, "original"))), length(list.dirs(file.path(cache_dir, "cpm"))))
  
  # check the metadata contains cell_, file_id_db, sample_ with correct types
  check_true("cell_" %in% names(metadata_tbl))
  check_true("file_id_db" %in% names(metadata_tbl)) 
  check_true("sample_" %in% names(metadata_tbl))
  metadata_tbl |> select(cell_, file_id_db, sample_) |> sapply(class) |> check_character()
  
  # check cell_ values in metadata_tbl is unique
  (anyDuplicated(metadata_tbl$cell_) == 0 ) |> assert_that(msg = "cell_ in the metadata must be unique")
  
  # check cell_ values are not duplicated when join with parquet
  cells <- get_metadata() |> select(cell_) |> as_tibble()
  (!any(metadata_tbl$cell_ %in% cells$cell_)) |> assert_that(msg = "cell_ in the metadata should not duplicate with that exists in API")
  
  # check age_days is either -99 or greater than 365
  assert_that(any(metadata_tbl$age_days==-99 | metadata_tbl$age_days> 365),
              msg = "age_days should be either -99 for unknown or greater than 365")
  
  # check sex capitalisation then convert to lower case 
  metadata_tbl <- metadata_tbl |> mutate(sex = tolower(sex))
  metadata_tbl |> distinct(sex) |> pull() |> check_subset(c('female','male','unknown'))
  
  # convert metadata_tbl to parquet if above checkpoints pass
  write_parquet(metadata_tbl, file.path(cache_dir, glue("metadata.parquet")))
  
  
  
}



