library(CuratedAtlasQueryR)
library(dplyr)
library(SingleCellExperiment)
library(assertthat)
library(checkmate)
library(arrow)

# import_metadata_and_counts(metadata_tbl, list_of_count_matriced)
# 
# 
# will accept
# 
# - single_cell_expeirment object, with certain characteristics (we need to decide the checks to do)
# 
# - tibble of metadata, with certain characteristics (we need to decide the checks to do)
# 
# what will do

# - the metadata has the minimum set of columns, with the right types (we could use https://github.com/mllg/checkmate)
# 
# - check that the cache does not include multiple version of the same metadata_file



# Maps user provided assay names to their corresponding paths in the repository
assay_map <- c(counts = "original",
               cpm = "cpm")

import_metadata_counts <- function(meta_input,
                                   cache_dir = CuratedAtlasQueryR:::get_default_cache_dir(),
                                   # it will be solved once the pull request is submitted and merged
                                   meta_output,
                                   version = "0.2.3") {
  # A few checks from here
  counts <- assay_map["counts"]
  cpm <- assay_map['cpm']
  checkTibble(meta_input)
  check_directory_exists(cache_dir) # should we add the constraint, if exist, delete? if not, continue?
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
  
  # check the number of sub directories in original match cpm 
  check_set_equal(length(list.dirs(file.path(cache_dir, counts))), length(list.dirs(file.path(cache_dir, cpm))))
  
}



