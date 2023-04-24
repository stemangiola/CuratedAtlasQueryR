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
# library(CuratedAtlasQueryR)
library(openssl)


# CREATE MAKEFILE
tab = "\t"
root_directory = "/vast/projects/cellxgene_curated"
splitted_light_data_directory = "/vast/projects/RCP/human_cell_atlas/splitted_light_data" #glue("{root_directory}/splitted_light_data")
DB_data_directory = glue("{root_directory}/splitted_DB_data")
gene_names = glue("{root_directory}/gene_names.rds")
files_metadata = glue("{root_directory}/files_metadata.rds")
metadata_path = glue('{root_directory}/metadata.rds')


get_metadata() |>
	distinct(.sample, cell_type) |>
	as_tibble() |>
	unite("file_id_db", c(.sample, cell_type), remove = FALSE) |>
	mutate(file_id_db = file_id_db |> md5()) |>

	mutate(
		input_file_path = glue("{splitted_light_data_directory}/{.sample}") |> as.character(),
		output_file_path = glue("{DB_data_directory}/{file_id_db}" |> as.character())
	) |>
	nest(data = -input_file_path) |>

	mutate(Mb = map_dbl(input_file_path, ~ (file.info(glue("{.x}/se.rds"))$size /1e6) |> as.integer() )) |>
	mutate(memory = pmax(Mb * 10, 10000)) |>
	unnest(data) |>
	rowid_to_column() |>
	mutate(commands = pmap(list(output_file_path, input_file_path,  memory, rowid, file_id_db), ~
												 	c(
												 		glue("CATEGORY=DB_data{..4}\nMEMORY={..3}\nCORES=1\nWALL_TIME=10000"),
												 		glue("{..1}:{..2} {gene_names} {files_metadata}\n{tab}Rscript DB_files.R {..2} {gene_names} {files_metadata} {..5} {..1}")
												 	)
												 	))  |>
	pull(commands) |>
	unlist() |>
	write_lines(glue("dev/DB_files.makeflow"))



# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file = args[[1]]
gene_names = args[[2]]
files_metadata = args[[3]]
file_id_db = args[[4]]
output_file = args[[5]]

# Create directory
output_file |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)
sample = basename(input_file) |> tools::file_path_sans_ext()

data =
	loadHDF5SummarizedExperiment(input_file	) |>
	inner_join(
		get_metadata() |>
			filter(.sample == !!sample) |>
			select(.cell, .sample, cell_type) |>
			as_tibble() |>
			unite("file_id_db", c(.sample, cell_type), remove = FALSE) |>
			mutate(file_id_db = file_id_db |> md5()) |>
			filter(file_id_db == !!file_id_db) |>
			select(.cell) |>
			distinct()
	)

# Select just the X assay
sce = SingleCellExperiment(list(X = data@assays@data$X))
rownames(sce) = rownames(data)
colnames(sce) = colnames(data)

sce |>	saveHDF5SummarizedExperiment(output_file, replace=TRUE)

