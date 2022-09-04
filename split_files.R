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
# CREATE MAKEFILE
# tab = "\t"
# root_directory = "/vast/scratch/users/mangiola.s/human_cell_atlas"
# raw_data_directory = glue("{root_directory}/raw_data")
# splitted_data_directory = glue("{root_directory}/splitted_data")
# 
# 
# metadata_path = glue('{root_directory}/metadata.rds')
# 

file_cell_type_table |>
	# tibble(split_file_paths, splitted_file_paths) |>
	mutate(Mb = map_dbl(input_file_path, ~ (file.info(.x)$size /1e6) |> as.integer() )) |>
	mutate(memory = pmax(Mb * 10, 10000)) |>
	rowid_to_column() |>
	mutate(commands = pmap(list(output_files_path, input_file_path,  memory, rowid), ~
												 	c(
												 		glue("CATEGORY=split_data{..4}\nMEMORY={..3}\nCORES=1\nWALL_TIME=10000"),
												 		glue("{..1}:{..2}\n{tab}Rscript split_files.R {..2} {..1}")
												 	)
	))  |>
	pull(commands) |>
	unlist() |>
	write_lines(glue("split_files.makeflow"))

# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file = args[[1]]
output_files_path = args[2:length(args)]

# If I have one cell type I don't need to do anything
if(length(output_files_path) == 1){
	system(glue("cp {input_file} {output_files_path[[1]]}"))
} else{
	
	# Otherwise 
	output_files_df = 
		output_files_path |> 
		enframe(value = "output_file_path") |> 
		mutate(file_id = basename(output_file_path) |> tools::file_path_sans_ext()) |> 
		separate(file_id, c("file_id", "cell_type"), sep="___") |> 
		select(cell_type, output_file_path)
	
	# Create directory
	output_files_df$output_file_path[1] |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)
	
	# Read file_cell_types
	data = readH5AD(input_file,	use_hdf5 = TRUE	) 
	
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
}

	
