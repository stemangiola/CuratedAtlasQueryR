library(zellkonverter)
library(SingleCellExperiment)
library(tidyverse)
library(purrr)
library(glue)
library(HCAquery)

# Read arguments
args = commandArgs(trailingOnly=TRUE)
root_directory = "/vast/projects/RCP/human_cell_atlas" # args[[1]]
metadata_sql = glue("{root_directory}/metadata_annotated.sqlite")
raw_data_directory = glue("{root_directory}/splitted_data_0.2")



samples =
	# get_metadata(metadata_sql) |>
	readRDS("/vast/projects/RCP/human_cell_atlas/metadata_annotated.rds") |>
	distinct(file_id, .sample) |>
	group_by(file_id) |>
	slice(1) |>
	pull(.sample)

# Read gene names
dir(raw_data_directory, full.names = TRUE) |>
	str_subset(samples  |> str_c(collapse = "|")) |>
	imap(	~ {
		print(.y)

		AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
													keys = rownames(loadHDF5SummarizedExperiment(.x)),
													keytype = "ENSEMBL",
													column = "SYMBOL",
													multiVals = "first"
		)
	}
	)  |>
	unlist() |>
	unique() |>
	saveRDS(glue("{root_directory}/gene_names.rds"))
