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
library(openssl)

source("utility.R")
options(scipen = 999)
#

# # CREATE MAKEFILE
# tab = "\t"
# root_directory = "/vast/scratch/users/mangiola.s/human_cell_atlas"
# raw_data_directory = glue("{root_directory}/raw_data")
# splitted_DB2_data_directory = glue("{root_directory}/splitted_DB2_data")
#
#
# metadata_path = glue('{root_directory}/metadata.rds')
# file_cell_types_directory = glue("{root_directory}/file_cell_types")
# input_files_path = dir(file_cell_types_directory, full.names = TRUE)
# #
# ## metadata = readRDS(metadata_path)
#
#
# # CREATE MAKEFILE
# tab = "\t"
# root_directory = "/vast/scratch/users/mangiola.s/human_cell_atlas"
# splitted_light_data_directory = glue("{root_directory}/splitted_light_data")
# DB_data_directory = glue("{root_directory}/splitted_DB_data")
# gene_names = glue("{root_directory}/gene_names.rds")
# files_metadata = glue("{root_directory}/files_metadata.rds")
# metadata_path = glue('{root_directory}/metadata.rds')
#
#
# get_metadata() |>
# 	distinct(.sample, file_id, cell_type) |>
# 	as_tibble() |>
# 	unite("file_id_db", c(.sample, cell_type), remove = FALSE) |>
# 	mutate(file_id_db = file_id_db |> md5() |> as.character()) |>
#
# 	unite("file_id_db2", c(file_id, cell_type), remove = FALSE) |>
# 	mutate(file_id_db2 = file_id_db2 |> md5() |> as.character()) |>
#
# 	mutate(input_file_path = glue("{DB_data_directory}/{file_id_db}") |> as.character()) |>
#
# 	## mutate(Mb = map_dbl(input_file_path, ~ (file.info(glue("{.x}/se.rds"))$size /1e6) |> as.integer() )) |>
# 	mutate(Mb = 1) |>
# 	group_by(file_id_db2) |>
#
# 	summarise(
# 		input_file_path = paste(input_file_path, collapse=" "),
# 		memory = pmax(Mb * 10) |> max(10000)
# 	) |>
#
# 	mutate(		output_file_path = glue("{splitted_DB2_data_directory}/{file_id_db2}" |> as.character())) |>
#
# 	rowid_to_column() |>
# 		mutate(commands = pmap(
# 			list(output_file_path, input_file_path,  memory, rowid),
# 			~
# 				c(
# 					glue("CATEGORY=split_data{..4}\nMEMORY={..3}\nCORES=1\nWALL_TIME=30000"),
# 					glue(
# 						"{..1}:{..2} {metadata_path}\n{tab}Rscript DB2_files.R {..2} {metadata_path} {..1}"
# 					)
# 				)
# 		))  |>
# 		pull(commands) |>
# 	unlist() |>
# 	write_lines(glue("dev/DB2_files.makeflow"))



# Read arguments
args = commandArgs(trailingOnly = TRUE)
input_file = args[1:(length(args)-2)] [[1]]
metadata_file = args[[length(args)-1]]
output_file = args[[length(args)]]

output_file |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

data =
	input_file |>
	map(loadHDF5SummarizedExperiment) %>%
	do.call(cbind, .)

# Select just the X assay
sce = SingleCellExperiment(list(X = data@assays@data$X))
rownames(sce) = rownames(data)
colnames(sce) = colnames(data)

rm(data)
gc()

sce |>	saveHDF5SummarizedExperiment(output_file, replace=TRUE)

