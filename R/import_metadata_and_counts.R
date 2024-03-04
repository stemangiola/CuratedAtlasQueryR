# Maps user provided assay names to their corresponding paths in the repository
assay_map <- c(counts = "original",
               cpm = "cpm")

#' Checks metadata meets importing criteria. 
#' Returns 
#' Calculating counts per million from raw counts
#' Gets the Curated Atlas metadata as a data frame.
#'
#'
#' @param meta_input
#' @param cache_dir Optional character vector of length 1. A file path on
#'   your local system to a directory (not a file) that will be used to store
#'   `metadata.parquet`
#' @param meta_output
#' @return A directory stores counts per million
#' @export
#' @examples
#' example_metadata <- get_metadata() |> head(3) |> as_tibble()
#' check_new_metadata <- import_metadata_counts(meta_input = example_metadata,
#'                                              meta_output = "example_metadata",
#'                                              cache_dir = tempdir())
#'
#' @importFrom assertthat assertthat
#' @importFrom checkmate check_tibble check_directory_exists check_set_equal check_true check_character check_subset check_file_exists
#' @importFrom dplyr tbl
#' @importFrom cli cli_alert_info
#' @importFrom glue glue
#' @importFrom arrow read_parquet


import_metadata_counts <- function(meta_input,
                                   cache_dir = CuratedAtlasQueryR:::get_default_cache_dir(),
                                   meta_output,
                                   version = "0.2.3") {
  # A few checks from here
  counts <- assay_map["counts"]
  cpm <- assay_map['cpm']
  check_tibble(meta_input)
  check_directory_exists(cache_dir)
  check_directory_exists(file.path(cache_dir, counts))
  sce <- meta_input |> 
    get_single_cell_experiment()
  se <- meta_input |> 
    get_seurat()
  
  counts_data <- assays(sce)$counts |> 
    as_tibble()
  assert_that(inherits(sce, "SingleCellExperiment"),
              msg = "SingleCellExperiment Object is not created from metadata.")
  assert_that(!'^ENS' %in% rownames(se[['originalexp']]),
              msg = "Gene names cannot contain Ensembl IDs.")
  assert_that(all(counts_data >= 0), 
              msg = "Counts for SingleCellExperiment cannot be negative.")
  all(meta_input |> pull(file_id_db) %in% dir(file.path(cache_dir, counts))) |> 
    assert_that(msg = "The metadata sample file ID and the count file ID does not match")
  write_parquet(meta_input, file.path(cache_dir, glue("{meta_output}.{version}.parquet")))
  
  # check the number of sub directories in original match cpm 
  check_set_equal(length(list.dirs(file.path(cache_dir, counts))), length(list.dirs(file.path(cache_dir, cpm))))
  
  # check the metadata contains cell_, file_id_db, sample_ with correct types
  check_true("cell_" %in% names(meta_input))
  check_true("file_id_db" %in% names(meta_input)) 
  check_true("sample_" %in% names(meta_input))
  meta_input |> select(cell_, file_id_db, sample_) |> sapply(class) |> check_character()
  
  # check age_days is either -99 or greater than 365
  assert_that(any(meta_input$age_days==-99 | meta_input$age_days> 365),
              msg = "age_days should be either -99 for unknown or greater than 365")
  
  # check sex capitalisation then convert to lower case 
  meta_input <- meta_input |> mutate(sex = tolower(sex))
  meta_input |> distinct(sex) |> pull() |> check_subset(c('female','male','unknown'))
  
  # if checkpoints above pass, generate cpm
  cli_alert_info("Generating cpm from {.path {cache_dir}/{counts}}. ")
  check_file_exists("count_per_millions.R")
  
  input_dir <- file.path(cache_dir, counts)
  output_dir <- file.path(cache_dir, cpm)
  
  for (subdir in list.files(input_dir, full.names = TRUE)) {
    file_name <- basename(subdir)
    count_path <- glue("{input_dir}/{file_name}/")
    cpm_path <- glue("{output_dir}/{file_name}/")
    # Iterate count_per_millions.R to generate cpm for each raw count
    command <-
      glue("Rscript ./count_per_millions.R {count_path} {cpm_path}" )
    system(command)
  }
  
  cli_alert_info("cpm are generated in {.path {cache_dir}/{cpm}}. ")
  
}



