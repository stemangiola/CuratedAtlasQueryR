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
library(openssl)

# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file_paths = args[1:(length(args)-2)]
files_metadata = args[[length(args)-1]]
output_file = args[[length(args)]]

# Read metadata
# metadata =
# 	input_file_paths|>
# 	map(readRDS)

# Get which column don't have too many NAs
common_colnames =
	input_file_paths |>
	imap_dfr(
		~ .x %>%
			{print(.y); (.)} |>
			readRDS() |>
			colnames() |>
			as_tibble()
	) |>
	dplyr::count(value) |>
	mutate(n_datasets = length(input_file_paths)) |>
	filter(n > (n_datasets / 2)) |>
	pull(value)

print(common_colnames)

# Get all metadata
input_file_paths  |>

	# Select core columns
	imap(~ .x %>%
			 	{print(.y); (.)} |>
			 	readRDS() |>
			 	select(.cell, .sample, .sample_name, one_of(common_colnames)) |>
			 	mutate_if(is.factor, as.character)
					) |>
	bind_rows() |>

		unite("file_id_db", c(.sample, cell_type), remove = FALSE) |>
		mutate(file_id_db = file_id_db |> md5()) |>

	# Curate tissue
	left_join(
		read_csv("dev/tissue_label_curated.csv"),
		by="tissue"
	) |>

	#mutate_if(is.character, as.factor) |>
	# Add files metadata
	left_join(readRDS(files_metadata) |> select_if(function(x) !is.list(x)), by="file_id") |>

	saveRDS(output_file)



