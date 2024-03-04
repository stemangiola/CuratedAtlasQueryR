library(CuratedAtlasQueryR)
library(dplyr)
library(SingleCellExperiment)
library(assertthat)
library(checkmate)
library(arrow)

# Maps user provided assay names to their corresponding paths in the repository
assay_map <- c(counts = "original",
               cpm = "cpm")

cols <- get_metadata(remote_url = get_database_url("fibrosis.0.2.3.parquet")) |>
  select(1:22) |> head(10) |> as_tibble()

col_map <- tibble(
  colnames = names(cols),
  class = sapply(cols, class)
)

import_metadata_counts <- function(meta_input,
                                   cache_dir = CuratedAtlasQueryR:::get_default_cache_dir(),
                                   meta_output,
                                   version = "0.2.3") {
  # A few checks from here
  counts <- assay_map["counts"]
  cpm <- assay_map['cpm']
  checkTibble(meta_input)
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
  arrow::write_parquet(meta_input, file.path(cache_dir, glue("{meta_output}.{version}.parquet")))
  
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



