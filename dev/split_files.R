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
#
# metadata_path = glue('{root_directory}/metadata.rds')
# file_cell_types_directory = glue("{root_directory}/file_cell_types")
# input_files_path = dir(file_cell_types_directory, full.names = TRUE)
#
# metadata = readRDS(metadata_path)
#
# file_cell_type_table =
# 	metadata |>
# 	distinct(.sample, file_id) |>
# 	mutate(
# 		input_file_path = glue("{raw_data_directory}/{file_id}.H5AD") |> as.character(),
# 		output_file_path = glue("{splitted_data_directory}/{.sample}.H5AD" |> as.character())
# 	) |>
# 	with_groups(input_file_path, ~ .x |> summarise(output_files_path = paste(output_file_path, collapse=" ")))
#
# file_cell_type_table |>
# 	# tibble(split_file_paths, splitted_file_paths) |>
# 	mutate(Mb = map_dbl(input_file_path, ~ (file.info(.x)$size /1e6) |> as.integer() )) |>
# 	mutate(memory = pmax(Mb * 10, 10000)) |>
# 	rowid_to_column() |>
# 	mutate(commands = pmap(list(output_files_path, input_file_path,  memory, rowid), ~
# 												 	c(
# 												 		glue("CATEGORY=split_data{..4}\nMEMORY={max(40000, ..3)}\nCORES=1\nWALL_TIME=30000"),
# 												 		glue("{..1}:{..2} {metadata_path}\n{tab}Rscript split_files.R {..2} {metadata_path} {splitted_data_directory}")
# 												 	)
# 	))  |>
# 	pull(commands) |>
# 	unlist() |>
# 	write_lines(glue("split_files.makeflow"))

# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file = args[[1]]
metadata_file = args[[2]]
directory_out = args[[3]]

metadata = readRDS(metadata_file)
file_id = basename(input_file) |> tools::file_path_sans_ext() |> str_split("___") %>% .[[1]] %>% .[1]
#directory_out = output_files_path[[1]] |> dirname()

# # If I have one cell type I don't need to do anything
# if(length(output_files_path) == 1){
# 	system(glue("cp {input_file} {output_files_path[[1]]}"))
# } else{



# Create directory
directory_out |> dir.create( showWarnings = FALSE, recursive = TRUE)

data =	readH5AD(input_file,	use_hdf5 = TRUE	)

colnames(data) = glue("{colnames(data)}_{file_id}")

data |>
	left_join(
		select(metadata, .cell, .sample),
		by=".cell"
	) |>

	nest(data = -.sample) |>

	mutate(saved = map2(data, .sample, ~ {
		print(.y)

		h5 = .x
		colData(h5) = h5 |> colData() |> droplevels()

		# Deal with writeH5AD error for 1-line files
		if(ncol(h5)==1){
			h5 = cbind(h5, h5)
			colnames(h5) = colnames(h5)[1] |> c(glue("{colnames(h5)[1]}_FAKE_CELL_DELETE"))
		}

		h5 |>	writeH5AD(glue("{directory_out}/{.y}.H5AD"))
	})) |>
	select(.sample)
# }


