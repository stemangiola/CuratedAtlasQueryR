library(zellkonverter)
library(SingleCellExperiment)
library(tidyverse)
library(purrr)
library(glue)
library(CuratedAtlasQueryR)
library(HDF5Array)

library(dbplyr)
library(DBI)
library(duckdb)


# Read arguments
args = commandArgs(trailingOnly=TRUE)
metadata_DB = "/vast/projects/cellxgene_curated/metadata_annotated_0.2.3.parquet"
root_directory = "/vast/projects/cellxgene_curated" # args[[1]]
raw_data_directory = glue("{root_directory}/splitted_data_0.2")

samples =
  duckdb() |>
  dbConnect(drv = _, read_only = TRUE) |>
  tbl(metadata_DB) |>
	distinct(file_id, sample_) |>
  as_tibble() |> 
	group_by(file_id) |>
	slice(1) |>
	pull(sample_)

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
