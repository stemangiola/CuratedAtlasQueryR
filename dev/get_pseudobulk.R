library(dplyr)
library(cellxgenedp)
library(tidyverse)
#library(tidySingleCellExperiment)
library(stringr)
library(scMerge)
library(glue)
library(DelayedArray)
library(HDF5Array)
library(openssl)
library(stringr)
library(CuratedAtlasQueryR)
library(purrr)

library(dbplyr)
library(DBI)
library(duckdb)
library(tidySingleCellExperiment)
library(tidySummarizedExperiment)
library(magrittr)

metadata_DB = "/vast/projects/cellxgene_curated/metadata_annotated_0.2.3.parquet"
root_directory = "/vast/projects/cellxgene_curated"

metadata =
  duckdb() |>
  dbConnect(drv = _, read_only = TRUE) |>
  tbl(metadata_DB)
# 
# CREATE MAKEFILE
data_directory = glue("{root_directory}/0.2.1/original")
output_directory = glue("{root_directory}/pseudobulk_0.2.3")
output_directory |> dir.create( showWarnings = FALSE, recursive = TRUE)
tab = "\t"
metadata |>
  distinct(file_id_db) |>
  as_tibble() |>
  mutate(
    input_file_path = glue("{data_directory}/{file_id_db}/assays.h5") |> as.character(),
    output_file_path = glue("{output_directory}/{file_id_db}.rds") |> as.character()
  ) |>
	mutate(commands = pmap(
	  list(input_file_path, file_id_db, output_file_path),
		~
			c(
				glue("CATEGORY=pseudobulk\nMEMORY=30000\nCORES=1\nWALL_TIME=30000"),
				glue(
					"{..3}:{..1} {metadata_DB}\n{tab}Rscript get_pseudobulk.R {..1} {..2} {..3}"
				)
			)
	))  |>
	pull(commands) |>
	unlist() |>
	write_lines(glue("~/PostDoc/CuratedAtlasQueryR/dev/get_pseudobulk.makeflow"))


# Read arguments
args = commandArgs(trailingOnly = TRUE)
input_file = args[[1]]
file_id_db_ = args[[2]]
output_file = args[[3]]



  metadata  = 
  metadata |>
  filter(file_id_db==file_id_db_) |>
  filter(confidence_class < 4) 
  
  if(metadata |> as_tibble() |> nrow() |> equals(0)){
    NULL |>
    saveRDS(output_file)
  } else {
    metadata |>
    get_single_cell_experiment(
      cache_directory = "/vast/projects/cellxgene_curated"
    ) |>
    aggregate_cells(c(sample_, cell_type_harmonised)) |>
    saveRDS(output_file)
}

