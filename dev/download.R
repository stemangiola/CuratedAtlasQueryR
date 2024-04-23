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
root_directory = "/vast/projects/cellxgene_curated/raw_data_Apr_2024"

files_metadata = 
	datasets(db) |>
	left_join(
		files(db) |> filter(filetype=="H5AD"),
		by = "dataset_id"
	) 

files_metadata |> saveRDS(glue("{root_directory}/files_metadata.rds"))


test = 
  files_metadata |> 
  slice(1:50) |> 
  nest(data = c(dataset_id, dataset_version_id, filetype, url)) |> 
  mutate(has_donor_id = map_lgl(
    data,
    ~ {
      h5_path = .x |> files_download(dry.run = FALSE)
      has_donor_id = 
        h5_path |> 
        readH5AD(use_hdf5 = TRUE	) |> 
        colData() |> 
        as_tibble() |> 
        select(any_of("donor_id")) |> 
        ncol() >
        0
      file.remove(h5_path)
      has_donor_id
    }
  )) |> 
  unnest(data) |> 
  select(dataset_version_id, has_donor_id)

files_metadata |>

	# Get organism list and filter human
	mutate(organism_name = map_chr(organism, ~ .x |> map(~.x$label) |> paste(collapse=", ") )) |>
	filter(organism_name |> str_detect("Homo sapiens")) |>

	# Download
	files_download(dry.run = FALSE, cache_path = "{root_directory}/raw_data/") |>

	# Save file list
	saveRDS("{root_directory}/file_location.rds")

