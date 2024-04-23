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
  tbl(metadata_DB) |>

  filter(confidence_class < 4) |>

  # Attach lineage
  left_join(
    read_csv("~/PostDoc/immuneHealthyBodyMap/metadata_cell_type.csv") |>
      replace_na(list(lineage_1 = "other_non_immune")),
    copy=TRUE
  ) |>
  mutate(is_immune = if_else(lineage_1 == "immune", "immune", "other")) |>

  # Fix typo
  mutate(tissue_harmonised = tissue_harmonised |> str_replace("plcenta", "placenta"))


# # CREATE MAKEFILE
# data_directory = glue("{root_directory}/0.2.1/original")
# output_directory = glue("/home/users/allstaff/mangiola.s/PostDoc/immuneHealthyBodyMap/pseudobulk_0.2.3.4")
# script_directory = "/home/users/allstaff/mangiola.s/PostDoc/CuratedAtlasQueryR/dev"
# 
# output_directory |> dir.create( showWarnings = FALSE, recursive = TRUE)
# tab = "\t"
# metadata |>
#   distinct(file_id, tissue_harmonised, cell_type_harmonised, is_immune) |>
#   as_tibble() |>
#   mutate(
#     tissue_harmonised = tissue_harmonised |> str_replace_all(" ", "__"),
#     cell_type_harmonised = cell_type_harmonised |> str_replace_all(" ", "__")
#   ) |>
#   mutate(
#     file_name = glue("{tissue_harmonised}____{cell_type_harmonised}____{is_immune}____{file_id}.rds") |>
#       as.character()|>
#       str_replace_all("/", "__"),
#     output_file_path =
#       glue("{output_directory}/{file_name}") |>
#       as.character()
#   ) |>
# 	mutate(commands = pmap(
# 	  list(output_file_path, file_id, tissue_harmonised, cell_type_harmonised, is_immune),
# 		~
# 			c(
# 				glue("CATEGORY=pseudobulk\nMEMORY=30000\nCORES=1\nWALL_TIME=30000"),
# 				glue(
# 					"{..1}:{metadata_DB}\n{tab}Rscript {script_directory}/get_pseudobulk.R {..2} {..3} {..4} {..5} {..1}"
# 				)
# 			)
# 	))  |>
# 	pull(commands) |>
# 	unlist() |>
# 	write_lines(glue(glue("{output_directory}/get_pseudobulk.makeflow")))


# Read arguments
args = commandArgs(trailingOnly = TRUE)
file_id = args[[1]]
tissue_harmonised = args[[2]] |> str_replace_all("____", " ")
cell_type_harmonised = args[[3]] |> str_replace_all("____", " ")
is_immune = args[[4]] 
output_file = args[[5]] 

metadata  = 
  metadata |>
  filter(
    file_id==!!file_id & 
      tissue_harmonised == !!tissue_harmonised & 
      cell_type_harmonised == !!cell_type_harmonised &
      is_immune == !!is_immune
  ) |>
    get_single_cell_experiment(
      cache_directory = "/vast/projects/cellxgene_curated"
    ) |>
    aggregate_cells(c(sample_, cell_type_harmonised)) |>
    saveRDS(output_file)

