library(zellkonverter)
library(SingleCellExperiment) # load early to avoid masking dplyr::count()
library(dplyr)
library(cellxgenedp)
library(tidyverse)
#library(tidySingleCellExperiment)
library(stringr)
library(scMerge)
library(glue)
library(DelayedArray)
library(HDF5Array)
# 
# # CREATE MAKEFILE
# tab = "\t"
# root_directory = "/vast/scratch/users/mangiola.s/human_cell_atlas"
# metadata_directory = glue("{root_directory}/metadata")
# raw_data_directory = glue("{root_directory}/raw_data")
# input_files = dir(raw_data_directory, full.names=TRUE) |> basename()
# 
# c(
# 	glue("CATEGORY=metadata\nMEMORY=80024\nCORES=1\nWALL_TIME=10000"),
# 	glue("{metadata_directory}/{input_files}:{raw_data_directory}/{input_files}\n{tab}Rscript get_metadata.R {raw_data_directory}/{input_files} {metadata_directory}/{input_files}")
# )  |>
# 	write_lines(glue("get_metadata.makeflow"))


# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file = args[[1]]
output_file = args[[2]]

# Create directory
output_file |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

# Read metadata
col_data = 
	readH5AD(input_file,	reader = "R",	use_hdf5 = TRUE	) |> 
	colData() 

if( ncol(col_data) == 0 )
	col_data = 
	readH5AD(input_file,	use_hdf5 = TRUE	) |> 
	colData() 

col_data |> 
 	as_tibble() |>
 	mutate(file_id = basename(input_file) |> tools::file_path_sans_ext()) |> 
 	mutate(file_path = input_file) |> 
	mutate_if(is.character, as.factor) |> 
	saveRDS(output_file)
	
