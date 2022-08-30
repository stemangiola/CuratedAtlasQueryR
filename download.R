suppressPackageStartupMessages({
	library(zellkonverter)
	library(SingleCellExperiment) # load early to avoid masking dplyr::count()
	library(dplyr)
	library(cellxgenedp)
})

library(tidyverse)
#library(tidySingleCellExperiment)
library(stringr)
library(scMerge)
library(glue)
library(DelayedArray)
library(HDF5Array)

db <- db()

# Arguments
args = commandArgs(trailingOnly=TRUE)
root_directory = args[[1]]

files_metadata = 
	datasets(db) |>
	left_join(
		files(db) |> filter(filetype=="H5AD"),
		by = "dataset_id"
	) 

files_metadata |> saveRDS(glue("{root_directory}/files_metadata.rds"))

files_metadata |>

	# Get organism list and filter human
	mutate(organism_name = map_chr(organism, ~ .x |> map(~.x$label) |> paste(collapse=", ") )) |>
	filter(organism_name |> str_detect("Homo sapiens")) |>

	# Download
	files_download(dry.run = FALSE, cache_path = "{root_directory}/raw_data/") |>

	# Save file list
	saveRDS("{root_directory}/file_location.rds")

