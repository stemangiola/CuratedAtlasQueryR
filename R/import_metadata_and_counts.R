

#' Checks importing criteria for new atlas 
#' Calculating counts per million from raw counts
#'
#' @param metadata_tbl A tibble
#' @param cache_dir Optional character vector of length 1. A file path on
#'   your local system to a directory (not a file) that will be used to store
#'   `metadata.parquet`
#' @param counts_path A character vector of original counts path
#' @param version A character vector of metadata version. Default is 0.2.3 
#' @export
#' @return A metadata.parquet with version, and a directory stores counts per million in cache directory
#' @examples
#' example_metadata <- get_metadata() |> head(3) |> as_tibble()
#' import_metadata_counts(metadata_tbl = example_metadata,
#'                        cache_dir = get_default_cache_dir(),
#'                        counts_path = "/var/folders/ls/99n281zx4bbd73kllmc1rc0h0005lj/T//RtmpUjWDqt/original",
#'                        counts_to_cpm_file = "/var/folders/ls/99n281zx4bbd73kllmc1rc0h0005lj/T//RtmpUjWDqt/counts_per_million.R"
#'
#' @importFrom assertthat assert_that
#' @importFrom checkmate check_tibble check_directory_exists check_set_equal check_true check_character check_subset check_file_exists
#' @importFrom dplyr tbl
#' @importFrom cli cli_alert_info
#' @importFrom glue glue
#' @importFrom arrow write_parquet
#' @importFrom fs dir_copy
import_metadata_counts <- function(metadata_tbl,
                                   cache_dir = CuratedAtlasQueryR:::get_default_cache_dir(),
                                   version = "0.2.3",
                                   counts_path) {
  # Convert to tibble if metadata_tbl is not a tibble 
  metadata_tbl <- metadata_tbl |> as_tibble()
  check_tibble(metadata_tbl)
  check_directory_exists(cache_dir)
  check_directory_exists(counts_path)
  dir_copy(counts_path, cache_dir)
  file_copy(counts_to_cpm_file, cache_dir)
  check_directory_exists(file.path(cache_dir, "original"))
  sce <- metadata_tbl |> 
    get_single_cell_experiment()
  se <- metadata_tbl |> 
    get_seurat()

  counts_data <- assays(sce)$counts |> 
    as_tibble()
  assert_that(inherits(sce, "SingleCellExperiment"),
              msg = "SingleCellExperiment Object is not created from metadata.")
  assert_that(!'^ENS' %in% rownames(se[['originalexp']]),
              msg = "Gene names cannot contain Ensembl IDs.")
  assert_that(all(counts_data >= 0), 
              msg = "Counts for SingleCellExperiment cannot be negative.")
  all(metadata_tbl |> pull(file_id_db) %in% dir(file.path(cache_dir, "original"))) |> 
    assert_that(msg = "The metadata sample file ID and the count file ID does not match")
  write_parquet(metadata_tbl, file.path(cache_dir, glue("metadata.{version}.parquet")))
  
  # check the number of sub directories in original match cpm 
  check_set_equal(length(list.dirs(file.path(cache_dir, "original"))), length(list.dirs(file.path(cache_dir, "cpm"))))
  
  # check the metadata contains cell_, file_id_db, sample_ with correct types
  check_true("cell_" %in% names(metadata_tbl))
  check_true("file_id_db" %in% names(metadata_tbl)) 
  check_true("sample_" %in% names(metadata_tbl))
  metadata_tbl |> select(cell_, file_id_db, sample_) |> sapply(class) |> check_character()
  
  # check age_days is either -99 or greater than 365
  assert_that(any(metadata_tbl$age_days==-99 | metadata_tbl$age_days> 365),
              msg = "age_days should be either -99 for unknown or greater than 365")
  
  # check sex capitalisation then convert to lower case 
  metadata_tbl <- metadata_tbl |> mutate(sex = tolower(sex))
  metadata_tbl |> distinct(sex) |> pull() |> check_subset(c('female','male','unknown'))
  
  # if checkpoints above pass, generate cpm
  cli_alert_info("Generating cpm from {.path {cache_dir}/original}. ")
  
  input_dir <- file.path(cache_dir, "original")
  output_dir <- file.path(cache_dir, "cpm")
  
  for (subdir in list.files(input_dir, full.names = TRUE)) {
    file_name <- basename(subdir)
    count_path <- glue("{input_dir}/{file_name}/")
    cpm_path <- glue("{output_dir}/{file_name}/")
    # Iterate get_counts_per_million function to generate cpm for each raw count
    get_counts_per_million(input_file = count_path,
                           output_file = cpm_path)
  }
  
  cli_alert_info("cpm are generated in {.path {cache_dir}/cpm}. ")
  
}



