library(zellkonverter)
library(SingleCellExperiment) # load early to avoid masking dplyr::count()
library(tidySingleCellExperiment)
library(dplyr)
# library(cellxgenedp)
library(tidyverse)
#library(tidySingleCellExperiment)
library(stringr)
# library(scMerge)
library(glue)
library(DelayedArray)
library(HDF5Array)
library(stringr)
library(HCAquery)
#
# #
# # source("utility.R")
# options(scipen = 999)
# #
#
# # CREATE MAKEFILE
# tab = "\t"
# root_directory = "/vast/projects/RCP/human_cell_atlas"
# raw_data_directory = glue("{root_directory}/raw_data")
# splitted_data_directory = glue("{root_directory}/splitted_data_0.2")
#
# metadata_directory = glue('{root_directory}/metadata_0.2')
# metadata_path = glue('{root_directory}/metadata_0.2.rds')
# file_cell_types_directory = glue("{root_directory}/file_cell_types")
# input_files_path = dir(file_cell_types_directory, full.names = TRUE)
# #
#
# readRDS(metadata_path) |>
# 	distinct(.sample, file_id) |>
# 	mutate(
# 		input_file_path = glue("{raw_data_directory}/{file_id}.H5AD") |> as.character(),
# 		output_file_path = glue("{splitted_data_directory}/{.sample}" |> as.character())
# 	) |>
#
# 	mutate(Mb = map_dbl(input_file_path, ~
# 												( (file.info(glue("{.x}"))$size /1e6) |> as.integer() )
# 	)) |>
# 	mutate(memory = pmax(Mb * 8, 10000)) |>
#
# 	mutate(memory = case_when(
# 		input_file_path %in% c(
# 			"/vast/projects/RCP/human_cell_atlas//raw_data/51f114ae-232a-4550-a910-934e175db814.H5AD"
# 		)~ 160000,
# 		TRUE ~ memory
# 	)) |>
#
# 	mutate(
# 		memory = case_when(
# 			input_file_path %in% c(
# 		"/vast/projects/RCP/human_cell_atlas/raw_data/b50b15f1-bf19-4775-ab89-02512ec941a6.H5AD",
# 				"/vast/projects/RCP/human_cell_atlas/raw_data/56e0359f-ee8d-4ba5-a51d-159a183643e5.H5AD",
# 				"/vast/projects/RCP/human_cell_atlas/raw_data/51f114ae-232a-4550-a910-934e175db814.H5AD",
# 				"/vast/projects/RCP/human_cell_atlas/raw_data/21ca95b3-776b-4fa2-9956-09a07c0e5224.H5AD",
# 		"/vast/projects/RCP/human_cell_atlas/raw_data/327927c7-c365-423c-9ebc-07acb09a0c1a.H5AD",
# 		"/vast/projects/RCP/human_cell_atlas/raw_data/5500774a-6ebe-4ddf-adce-90302b7cd007.H5AD",
# 		"/vast/projects/RCP/human_cell_atlas/raw_data/56ed3d2a-a8cf-4c60-a184-f4e3e4af5176.H5AD"
#
# 			) ~ 160000,
# 			input_file_path ==	"/vast/projects/RCP/human_cell_atlas/raw_data/08247324-7ab7-45c5-8bd6-6c22676761ed.H5AD" ~ 160000,
#
# 			TRUE ~ memory
# 		)
# 	) |>
# 	mutate(metadata_file_id = glue("{metadata_directory}/{file_id}.rds")) |>
# 	rowid_to_column() |>
# 	mutate(commands = pmap(
# 		list(output_file_path, input_file_path,  memory, rowid, metadata_file_id),
# 		~
# 			c(
# 				glue("CATEGORY=split_data{..4}\nMEMORY={..3}\nCORES=1\nWALL_TIME=30000"),
# 				glue(
# 					"{..1}:{..2} {..5}\n{tab}Rscript split_files.R {..2} {..5} {..1}"
# 				)
# 			)
# 	))  |>
# 	pull(commands) |>
# 	unlist() |>
# 	write_lines(glue("~/PostDoc/HCAquery/dev/split_files.makeflow"))

# Read arguments
args = commandArgs(trailingOnly = TRUE)
input_file = args[[1]]
metadata_file = args[[2]]
output_file = args[[3]]

file_id = basename(input_file) |> tools::file_path_sans_ext() |> str_split("___") %>% .[[1]] %>% .[1]
.sample = basename(output_file) |> tools::file_path_sans_ext()


# Create directory
output_file |> dir.create(showWarnings = FALSE, recursive = TRUE)

data = readH5AD(input_file,	use_hdf5 = TRUE, reader = "R")
data_for_rows_and_columns = readH5AD(input_file,	use_hdf5 = TRUE, skip_assays = TRUE)
colnames(data) = data_for_rows_and_columns |> colnames()
rownames(data) = data_for_rows_and_columns|> rownames()
rm(data_for_rows_and_columns)
gc()

# Check row names and colnames
if(colnames(data) |> is.null() | rownames(data) |> is.null()) stop()
if(!rownames(data) |> sample(100, replace = TRUE) |> str_detect("^ENSG") |> all()) stop()

h5 =
	data |>
	inner_join(

		# Metadata
		readRDS(metadata_file) |>
			filter(.sample == !!.sample) |>
			select(.cell) |>
			as_tibble()
	)

# if (ncol(h5) == 1) {
# 	h5 = cbind(h5, h5)
# 	colnames(h5) = colnames(h5)[1] |> c(glue("{colnames(h5)[1]}_FAKE_CELL_DELETE"))
# }

h5 |>	saveHDF5SummarizedExperiment(output_file, replace=TRUE, verbose=TRUE, )
