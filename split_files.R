library(zellkonverter)
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

source("utility.R")

# 
# 
# # CREATE MAKEFILE
# tab = "\t"
# root_directory = "/vast/scratch/users/mangiola.s/human_cell_atlas"
# raw_data_directory = glue("{root_directory}/raw_data")
# splitted_data_directory = glue("{root_directory}/splitted_data")
# 
# file_cell_types_directory = glue("{root_directory}/file_cell_types")
# input_files_path = dir(file_cell_types_directory, full.names = TRUE)
# file_cell_type_table =
# 	input_files_path |>
# 	map_dfr(readRDS) |>
# 	mutate(cell_type = cell_type |> 
# 				 	clean_cell_type()
# 				) |>
# 	mutate(output_file = glue("{file_id}___{cell_type}")) |>
# 	mutate(
# 		input_file_path = glue("{raw_data_directory}/{file_id}.H5AD") |> as.character(),
# 		output_file_path = glue("{splitted_data_directory}/{output_file}.H5AD" |> as.character())
# 	) |>
# 	with_groups(input_file_path, ~ .x |> summarise(output_files_path = paste(output_file_path, collapse=" ")))
# 
# c(
# 	glue("CATEGORY=split_data\nMEMORY=60024\nCORES=1\nWALL_TIME=10000"),
# 	glue("{file_cell_type_table$output_files_path}:{file_cell_type_table$input_file_path}\n{tab}Rscript split_files.R {file_cell_type_table$input_file_path} {file_cell_type_table$output_files_path}")
# )  |>
# 	write_lines(glue("split_files.makeflow"))



# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file = args[[1]]
output_files_path = args[2:length(args)]

output_files_df = 
	output_files_path |> 
	enframe(value = "output_file_path") |> 
	mutate(file_id = basename(output_file_path) |> tools::file_path_sans_ext()) |> 
	separate(file_id, c("file_id", "cell_type"), sep="___") |> 
	select(cell_type, output_file_path)
	
# Create directory
output_files_df$output_file_path[1] |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

# Read file_cell_types
data = readH5AD(input_file,	reader = "R",	use_hdf5 = TRUE	) 

if(ncol(colData(data)) ==0 )
	data = readH5AD(input_file,	use_hdf5 = TRUE	) 


data |> 
	
	# Mutate cell types prettify
	mutate(cell_type = cell_type |> 
				 	clean_cell_type()
			) |> 

	# add output files
	left_join(output_files_df, by ="cell_type") |> 

	nest(data = -output_file_path) |> 
	
	mutate(saved = map2(data, output_file_path, ~ {
		h5 = .x
		colData(h5) = h5 |> colData() |> droplevels()
		
		# Deal with writeH5AD error for 1-line files
		if(ncol(h5)==1){
			h5 = cbind(h5, h5)
			colnames(h5) = colnames(h5)[1] |> c(glue("{colnames(h5)[1]}_FAKE_CELL_DELETE"))
		}

		h5 |>	writeH5AD(.y)
	})) 
	
