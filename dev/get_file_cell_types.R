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

# # CREATE MAKEFILE
# tab = "\t"
# root_directory = "/vast/scratch/users/mangiola.s/human_cell_atlas"
# file_cell_types_directory = glue("{root_directory}/file_cell_types")
# raw_data_directory = glue("{root_directory}/raw_data")
# input_files_path = dir(raw_data_directory, full.names=TRUE)
# input_files = input_files_path |> basename()
# output_files = input_files |> str_replace("H5AD$", "rds")
# output_files_path = glue("{file_cell_types_directory}/{output_files}")
# 
# c(
# 	glue("CATEGORY=file_cell_types\nMEMORY=80024\nCORES=1\nWALL_TIME=10000"),
# 	glue("{output_files_path}:{input_files_path}\n{tab}Rscript get_file_cell_types.R {input_files_path} {output_files_path}")
# )  |>
# 	write_lines(glue("get_file_cell_types.makeflow"))



# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file = args[[1]]
output_file = args[[2]]

# Create directory
output_file |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

# Read file_cell_types
col_data = 
	readH5AD(input_file,	reader = "R",	use_hdf5 = TRUE	) |> 
	colData() 

if( ncol(col_data) == 0 )
	col_data = 
	readH5AD(input_file,	use_hdf5 = TRUE	) |> 
	colData() 

col_data |> 
 	as_tibble(rownames = ".cell") |>
 	mutate(file_id = basename(input_file) |> tools::file_path_sans_ext()) |> 
 	distinct(file_id, cell_type) |> 
	mutate_if(is.character, as.factor) |> 
	saveRDS(output_file)
	
