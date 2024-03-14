#' Checks importing criteria for new atlas 
#' Generats counts per million from provided raw counts
#'
#' @param sce_rds A SingleCellExperiment object Rds
#' @param cache_dir Optional character vector of length 1. A file path on
#'   your local system to a directory (not a file) that will be used to store
#'   `metadata.parquet`
#' @export
#' @return A metadata.parquet, and a counts per million directory in provided cache directory
#' @examples
#' import_metadata_counts(sce_rds = "~/projects/caq/import_api_pipelines/12eb5fe25994253c1d320ca590a6e999/se.rds",
#'                        cache_directory = "~/projects/caq/cache_for_testing")
#'
#' @importFrom assertthat assert_that
#' @importFrom checkmate check_tibble check_directory_exists check_set_equal check_true check_character check_subset check_file_exists
#' @importFrom dplyr tbl
#' @importFrom cli cli_alert_info
#' @importFrom glue glue
#' @importFrom arrow write_parquet
#' @importFrom fs dir_copy
#' @importFrom purrr walk2
#' @importFrom openssl md5
import_metadata_counts <- function(
  sce_rds = data,  
  cache_dir = get_default_cache_dir()) {
  
  original_dir <- file.path(cache_dir, "original")
  
  # Load sce, identify metadata and counts matrix
  sce_obj <- readRDS(sce_rds)
  metadata_tbl <- metadata(sce_obj)$data
  counts_matrix <- sce@assays@data$X
  
  # Convert to tibble if metadata_tbl is not a tibble
  metadata_tbl <- metadata_tbl |> as_tibble()
  
  # create file_id_db from dataset_id
  metadata_tbl <- metadata_tbl |> mutate(file_id_db = dataset_id |> md5() |> as.character())
  
  # create original and cpm folders in cache directory if not exist (in order to append new counts to existing ones)
  if (!dir.exists(original_dir)) {
    dir.create(cache_dir, "original", recursive = TRUE)
  }
  
  if (!dir.exists(file.path(cache_dir, "cpm"))) {
    dir.create(cache_dir, "cpm", recursive = TRUE)
  }
  
  # existing metadata genes
  genes <- get_metadata() |> head(1) |> 
    get_single_cell_experiment(cache_directory = cache_dir) |> rownames()
  
  # check count H5 directory name not included in the cache directory original
  all(!metadata_tbl$file_id_db %in% dir(original_dir)) |>
    check_true() |>
    assert_that(msg = "Count H5 directory name should not duplicate with that in cache directory")

  counts_data <- sce_obj@assays@data$X
  assert_that(inherits(sce_obj, "SingleCellExperiment"),
              msg = "SingleCellExperiment Object is not created from metadata.")
  assert_that(!"^ENS" %in% genes,
              msg = "Gene names cannot contain Ensembl IDs.")
  assert_that(all(counts_data >= 0),
              msg = "Counts for SingleCellExperiment cannot be negative.")
  
  # check the metadata contains cell_, file_id_db, sample_ with correct types
  check_true("cell_" %in% names(metadata_tbl))
  check_true("file_id_db" %in% names(metadata_tbl)) 
  metadata_tbl |> select(cell_, file_id_db) |> sapply(class) |> check_character()
  
  # check cell_ values in metadata_tbl is unique
  (anyDuplicated(metadata_tbl$cell_) == 0 ) |> assert_that(msg = "cell_ in the metadata must be unique")
  
  # check cell_ values are not duplicated when join with parquet
  cells <- get_metadata() |> select(cell_) |> as_tibble()
  (!any(metadata_tbl$cell_ %in% cells$cell_)) |> assert_that(msg = "cell_ in the metadata should not duplicate with that exists in API")
  
  # check age_days is either -99 or greater than 365
  if (any(colnames(metadata_tbl) == "age_days")) {
    assert_that(all(metadata_tbl$age_days==-99 | metadata_tbl$age_days> 365),
                msg = "age_days should be either -99 for unknown or greater than 365")
  }
  
  # check sex capitalisation then convert to lower case 
  if (any(colnames(metadata_tbl) == "sex")) {
    metadata_tbl <- metadata_tbl |> mutate(sex = tolower(sex))
    metadata_tbl |> distinct(sex) |> pull() |> check_subset(c("female","male","unknown"))
  }
  
  counts_path <-
    metadata_tbl |> select(file_id_db) |> mutate(
      file_path = sce_rds,
      original_path = file.path(original_dir, basename(file_id_db)),
      cpm_path = file.path(cache_dir, "cpm", basename(file_id_db))
    )
  
  # if checkpoints above pass, generate cpm
  cli_alert_info("Generating cpm from {.path {metadata_tbl$file_id_db}}. ")
  
  # Generate cpm from counts
  get_counts_per_million(input_file_rds = sce_rds, output_file = counts_path$cpm_path)
  dir_copy(sce_rds |> dirname(), counts_path$original_path)
  rds_data = readRDS(sce_rds)
  # check whether new data has the same gene set as existing metadata
  (rownames(rds_data) |> length() == genes |> length()) |>
    assert_that(msg = "CuratedAtlasQuery reports:
                  The number of genes in the existing metadata count does not match the number of genes specified")
  saveHDF5SummarizedExperiment(rds_data, counts_path$original_path, replace=TRUE)
  
  cli_alert_info("cpm are generated in {.path {counts_path$cpm_path}}. ")
  
  # check metadata sample file ID match the count file ID in cache directory
  all(metadata_tbl |> pull(file_id_db) %in% dir(original_dir)) |> 
    assert_that(msg = "The metadata sample file ID and the count file ID does not match")
  
  # convert metadata_tbl to parquet if above checkpoints pass
  write_parquet(metadata_tbl, file.path(cache_dir, glue("metadata.parquet")))
}



