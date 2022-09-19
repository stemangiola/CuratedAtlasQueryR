
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

metadata_df = readRDS(metadata)

metadata_df |>
	distinct(.sample, file_id) |>
	filter(.sample %in% samples_to_use) |>
	mutate(
		input_file_path = glue("{light_data_directory}/{.sample}") |> as.character(),
		output_file_path = glue("{annotated_data_directory}/{.sample}.rds" |> as.character())
	) |>

	mutate(Mb = map_dbl(input_file_path, ~ (file.info(glue("{.x}/se.rds"))$size /1e6) |> as.integer() )) |>
	mutate(memory = pmax(Mb * 20 + 20000, 40000)) |>

	# mutate(memory = case_when(
	# 	output_file_path %in% c(
	# 		"/vast/scratch/users/mangiola.s/human_cell_atlas/annotated_data/305a8bd9b8e529a967feb5f73cc8c4df.rds" ,
	# 		"/vast/scratch/users/mangiola.s/human_cell_atlas/annotated_data/087c2093be040a404c9685af1ecb3c65.rds",
	# 		"/vast/scratch/users/mangiola.s/human_cell_atlas/annotated_data/86a6d20305d912e98318ad4d1d5d1814.rds",
	# 		"/vast/scratch/users/mangiola.s/human_cell_atlas/annotated_data/829b99a569ec9ebb5fdd1b0b29208aaf.rds",
	# 		"/vast/scratch/users/mangiola.s/human_cell_atlas/annotated_data/176f8892a21bec1bd7bdbc4181af75ed.rds",
	# 		"/vast/scratch/users/mangiola.s/human_cell_atlas/annotated_data/23c822334c194bceb576a9ccb1db5929.rds",
	# 		"/vast/scratch/users/mangiola.s/human_cell_atlas/annotated_data/3c20ba18525fb5e0b41cb8ea189b5d33.rds",
	# 		"/vast/scratch/users/mangiola.s/human_cell_atlas/annotated_data/522dde7ab389d65b265d4cd598576f31.rds",
	# 		"/vast/scratch/users/mangiola.s/human_cell_atlas/annotated_data/2cf3bb4ffbb2024a9ca04baec073ae14.rds",
	# 		"/vast/scratch/users/mangiola.s/human_cell_atlas/annotated_data/c8ff7c63b3152a25c338cc279b31ab07.rds",
	# 		"/vast/scratch/users/mangiola.s/human_cell_atlas/annotated_data/5be263dbc1384b3cec21c5d3c580f838.rds",
	# 		"/vast/scratch/users/mangiola.s/human_cell_atlas/annotated_data/024d53b702b1846a476cabe5d691f992.rds",
	# 		"/vast/scratch/users/mangiola.s/human_cell_atlas/annotated_data/9da244f06591fa49e5649c65ed3b0934.rds",
	# 		)~ 160000,
	# 	TRUE ~ memory
	# )) |>

	rowid_to_column() |>
	mutate(commands = pmap(list(output_file_path, input_file_path,  memory, rowid, file_id), ~
												 	c(
												 		glue("CATEGORY=light_data{..4}\nMEMORY={..3}\nCORES=1\nWALL_TIME=30000"),
												 		glue("{..1}:{..2} {metadata} {cell_type_df}\n{tab}Rscript annotate_files.R {..2} {metadata} {cell_type_df} {..1}")
												 	)
	))  |>
	pull(commands) |>
	unlist() |>
	write_lines(glue("annotate_files.makeflow"))


# c(
# 	glue("CATEGORY=light_data\nMEMORY=30024\nCORES=1\nWALL_TIME=10000"),
# 	glue("{annotated_file_paths}:{light_file_paths} {metadata} {cell_type_df}\n{tab}Rscript annotate_files.R {light_file_paths} {metadata} {cell_type_df} {annotated_file_paths}")
# )  |>
# 	write_lines(glue("annotate_files.makeflow"))
