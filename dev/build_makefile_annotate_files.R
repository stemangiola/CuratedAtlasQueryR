
library(zellkonverter)
library(Seurat)
library(SingleCellExperiment) # load early to avoid masking dplyr::count()
library(tidySingleCellExperiment)
library(dplyr)
library(cellxgenedp)
library(tidyverse)
#library(tidySingleCellExperiment)
library(stringr)
library(scMerge)
library(glue)
library(DelayedArray)
library(HDF5Array)
library(tidyseurat)
library(celldex)
library(SingleR)
library(glmGamPoi)
source("utility.R")


# # CREATE MAKEFILE
tab = "\t"
root_directory = "/vast/scratch/users/mangiola.s/human_cell_atlas"
annotated_data_directory = glue("{root_directory}/annotated_data")
light_data_directory = glue("{root_directory}/splitted_light_data")
metadata = glue("{root_directory}/metadata.rds")
cell_type_df = "metadata_cell_type.csv"

light_file_paths = dir(light_data_directory, full.names = TRUE)
.sample = basename(light_file_paths) |> tools::file_path_sans_ext()
annotated_file_paths = glue("{annotated_data_directory}/{.sample}.rds")

c(
	glue("CATEGORY=light_data\nMEMORY=30024\nCORES=1\nWALL_TIME=10000"),
	glue("{annotated_file_paths}:{light_file_paths} {metadata} {cell_type_df}\n{tab}Rscript annotate_files.R {light_file_paths} {metadata} {cell_type_df} {annotated_file_paths}")
)  |>
	write_lines(glue("annotate_files.makeflow"))
