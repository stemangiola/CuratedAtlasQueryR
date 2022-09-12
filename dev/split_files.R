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
options(scipen = 999)
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
# #
# metadata = readRDS(metadata_path)
#
# metadata |>
# 	distinct(.sample, file_id) |>
# 	mutate(
# 		input_file_path = glue("{raw_data_directory}/{file_id}.H5AD") |> as.character(),
# 		output_file_path = glue("{splitted_data_directory}/{.sample}" |> as.character())
# 	) |>
#
# 	mutate(Mb = map_dbl(input_file_path, ~ (file.info(.x)$size / 1e6) |> as.integer())) |>
# 	mutate(memory = pmax(Mb * 2, 30000)) |>
# 	mutate(
# 		memory = case_when(
# 			input_file_path %in% c(
# 				"/vast/scratch/users/mangiola.s/human_cell_atlas/raw_data/b50b15f1-bf19-4775-ab89-02512ec941a6.H5AD",
# 				"/vast/scratch/users/mangiola.s/human_cell_atlas/raw_data/56e0359f-ee8d-4ba5-a51d-159a183643e5.H5AD",
# 				"/vast/scratch/users/mangiola.s/human_cell_atlas/raw_data/51f114ae-232a-4550-a910-934e175db814.H5AD",
# 				"/vast/scratch/users/mangiola.s/human_cell_atlas/raw_data/21ca95b3-776b-4fa2-9956-09a07c0e5224.H5AD"
# 			) ~ 160000,
# 			input_file_path ==	"/vast/scratch/users/mangiola.s/human_cell_atlas/raw_data/08247324-7ab7-45c5-8bd6-6c22676761ed.H5AD" ~ 200000,
#
# 			TRUE ~ memory
# 		)
# 	) |>
# 	rowid_to_column() |>
# 	mutate(commands = pmap(
# 		list(output_file_path, input_file_path,  memory, rowid),
# 		~
# 			c(
# 				glue("CATEGORY=split_data{..4}\nMEMORY={..3}\nCORES=1\nWALL_TIME=30000"),
# 				glue(
# 					"{..1}:{..2} {metadata_path}\n{tab}Rscript split_files.R {..2} {metadata_path} {..1}"
# 				)
# 			)
# 	))  |>
# 	pull(commands) |>
# 	unlist() |>
# 	write_lines(glue("split_files.makeflow"))

# Read arguments
args = commandArgs(trailingOnly = TRUE)
input_file = args[[1]]
metadata_file = args[[2]]
output_file = args[[3]]

file_id = basename(input_file) |> tools::file_path_sans_ext() |> str_split("___") %>% .[[1]] %>% .[1]
.sample = basename(output_file) |> tools::file_path_sans_ext()

metadata =
	readRDS(glue("{dirname(input_file)}/../metadata/{file_id}.rds")) |>
	select(.cell, .sample)


# Create directory
output_file |> dir.create(showWarnings = FALSE, recursive = TRUE)

data = readH5AD(input_file,	use_hdf5 = TRUE)

colnames(data) = metadata |> pull(.cell)

h5 =
	data |>
	left_join(metadata) |>
	filter(.sample == !!.sample)

# if (ncol(h5) == 1) {
# 	h5 = cbind(h5, h5)
# 	colnames(h5) = colnames(h5)[1] |> c(glue("{colnames(h5)[1]}_FAKE_CELL_DELETE"))
# }

h5 |>	saveHDF5SummarizedExperiment(output_file, replace=TRUE, verbose=TRUE)
