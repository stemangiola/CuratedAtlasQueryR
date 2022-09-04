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
# metadata_directory = glue("{root_directory}/metadata")
# raw_data_directory = glue("{root_directory}/raw_data")
# files_metadata = glue("{root_directory}/files_metadata.rds")
# 
# input_files_path = dir(raw_data_directory, full.names=TRUE)
# input_files = input_files_path |> basename()
# output_files = input_files |> str_replace("H5AD$", "rds")
# output_files_path = glue("{metadata_directory}/{output_files}")
# metadata_path = glue("{root_directory}/metadata.rds")
# c(
# 	glue("CATEGORY=get_metadata\nMEMORY=80024\nCORES=1\nWALL_TIME=10000"),
# 	glue("{output_files_path}:{input_files_path}\n{tab}Rscript get_metadata.R {input_files_path} {output_files_path}"),
# 	glue("CATEGORY=merge_metadata\nMEMORY=80024\nCORES=1\nWALL_TIME=10000"),
# 	glue("{metadata_path}:{paste(output_files_path, collapse = \" \")} {files_metadata}\n{tab}Rscript merge_metadata.R {paste(output_files_path, collapse = \" \")} {files_metadata} {metadata_path}")
# )  |>
# 	write_lines(glue("get_metadata.makeflow"))

source("utility.R")

root_directory = "/vast/scratch/users/mangiola.s/human_cell_atlas"
splitted_light_data_directory = glue("{root_directory}/splitted_light_data")

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

# Process
if( ncol(col_data) == 0 )
	col_data = 
	readH5AD(input_file,	use_hdf5 = TRUE	) |> 
	colData() 

col_data |> 
 	as_tibble(rownames = ".cell") |>
	
	# Link file IDs
	mutate(file_id = basename(input_file) |> tools::file_path_sans_ext()) |> 
	mutate(split_file_id = glue("{file_id}___{clean_cell_type(cell_type)}") |> as.character()) |> 
 	mutate(file_path = glue("{splitted_light_data_directory}/{split_file_id}.H5AD" |> as.character())) |> 
	
	# Make cell names unique
	mutate(.cell = glue("{.cell}_{file_id}")) |> 
	
	mutate_if(is.character, as.factor) |> 
	saveRDS(output_file)
	
