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




# # CREATE MAKEFILE
# tab = "\t"
# root_directory = "/vast/scratch/users/mangiola.s/human_cell_atlas"
# splitted_data_directory = glue("{root_directory}/splitted_data")
# light_data_directory = glue("{root_directory}/splitted_light_data")
# gene_names = glue("{root_directory}/gene_names.rds")
# files_metadata = glue("{root_directory}/files_metadata.rds")
# metadata_path = glue('{root_directory}/metadata.rds')
#
# metadata = readRDS(metadata_path)
#
#
# metadata |>
# 	distinct(.sample, file_id) |>
# 	mutate(
# 		input_file_path = glue("{splitted_data_directory}/{.sample}") |> as.character(),
# 		output_file_path = glue("{light_data_directory}/{.sample}" |> as.character())
# 	) |>
#
# 	mutate(Mb = map_dbl(input_file_path, ~ (file.info(glue("{.x}/se.rds"))$size /1e6) |> as.integer() )) |>
# 	mutate(memory = pmax(Mb * 10, 10000)) |>
# 	rowid_to_column() |>
# 	mutate(commands = pmap(list(output_file_path, input_file_path,  memory, rowid, file_id), ~
# 												 	c(
# 												 		glue("CATEGORY=light_data{..4}\nMEMORY={..3}\nCORES=1\nWALL_TIME=10000"),
# 												 		glue("{..1}:{..2} {gene_names} {files_metadata}\n{tab}Rscript light_files.R {..2} {gene_names} {files_metadata} {..5} {..1}")
# 												 	)
# 												 	))  |>
# 	pull(commands) |>
# 	unlist() |>
# 	write_lines(glue("light_files.makeflow"))





# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file = args[[1]]
gene_names = args[[2]]
files_metadata = args[[3]]
file_id = args[[4]]
output_file = args[[5]]


# Create directory
output_file |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

# Read file_cell_types
# data = readH5AD(input_file, reader = "R",	use_hdf5 = TRUE	)
# # Read file_cell_types
# if(ncol(colData(data)) >0 & "Cell" %in% colnames(colData(data))){
#
#
# 	colnames(data) = colData(data)$Cell
# 	data@colData = data@colData[,!colnames(data@colData) %in% "Cell"]
#
# } else
data = loadHDF5SummarizedExperiment(input_file	)

rownames(data) = rowData(data)$feature_name
rowData(data) = NULL
colData(data) = NULL
reducedDims(data) = NULL
data@int_colData$colPairs = data@int_colData$colPairs[,0]


# Complete gene set
missing_genes = readRDS(gene_names)  |> setdiff(rownames(data))

missing_matrix =
	HDF5RealizationSink(c(length(missing_genes),ncol(data)), as.sparse = TRUE) |>
	as("DelayedArray")

rownames(missing_matrix) = missing_genes
colnames(missing_matrix) = colnames(data)

missing_sce = SingleCellExperiment(list(X=missing_matrix),  colData=colData(data))
missing_sce@int_colData = data@int_colData


# Select just the X assay
data@assays@data = data@assays@data |> as.list() %>% .[1] |> SimpleList()


# Make cell name unique
data = data |> rbind(missing_sce	)

transformation =
	readRDS(files_metadata) |>
	filter(file_id == !!file_id) |>
	mutate(transformation = case_when(
		x_normalization |> str_detect("log|ln|Log") ~ "log",
		x_normalization |> str_detect("square root") ~ "square_root",
		TRUE ~ "none"
	)) |>
	pull(transformation) |>
	unique()

data@assays@data$X =
	data@assays@data$X |>
	when(
		transformation=="log" ~ {
			.x = expm1(.)
			#type(.x) = "integer"
			.x
		},
		transformation=="square_root" ~ {
			.x = (.)^2
			#type(.x) = "integer"
			.x
		},
		~ {
			.x = (.)
			#type(.x) = "integer"
			.x
		}
	)


data |>	saveHDF5SummarizedExperiment(output_file, replace=TRUE)

