library(zellkonverter)
library(SingleCellExperiment)
library(tidyverse)
library(purrr)
library(glue)

# Read arguments
args = commandArgs(trailingOnly=TRUE)
root_directory = args[[1]]
raw_data_directory = glue("{root_directory}/raw_data")

# Read gene names
gene_names = 
	dir(raw_data_directory, full.names = TRUE) |> 
	map(	~ rowData(readH5AD(.x, reader = "R", use_hdf5 = TRUE))$feature_name
	)  |> 
	unlist() |> 
	unique() |> 
	saveRDS(glue("{root_directory}/gene_names.rds"))