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

source("utility.R")



# # CREATE MAKEFILE
# tab = "\t"
# root_directory = "/vast/scratch/users/mangiola.s/human_cell_atlas"
# annotated_data_directory = glue("{root_directory}/annotated_data")
# light_data_directory = glue("{root_directory}/splitted_light_data")
# metadata = glue("{root_directory}/metadata.rds")
# cell_type_df = "metadata_cell_type.csv"
#
# light_file_paths = dir(light_data_directory, full.names = TRUE)
# file_id = basename(light_file_paths) |> tools::file_path_sans_ext()
# annotated_file_paths = glue("{annotated_data_directory}/{file_id}.rds")
#
# c(
# 	glue("CATEGORY=light_data\nMEMORY=20024\nCORES=1\nWALL_TIME=10000"),
# 	glue("{annotated_file_paths}:{light_file_paths} {metadata} {cell_type_df}\n{tab}Rscript annotate_files.R {light_file_paths} {metadata} {cell_type_df} {annotated_file_paths}")
# )  |>
# 	write_lines(glue("annotate_files.makeflow"))



# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file = args[[1]]
metadata = args[[2]]
cell_type_df = args[[3]]
output_file = args[[4]]

file_id = basename(input_file) |> tools::file_path_sans_ext() |> str_split("___") %>% .[[1]] %>% .[1]

# Create directory
output_file |>  dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

# Read file_cell_types
data =
	readH5AD(input_file,	use_hdf5 = TRUE	)  |>

	# add lineage 1
	left_join(readRDS(metadata) |> distinct(.cell, .sample, cell_type)) |>
	left_join(read_csv(cell_type_df)) |>
	filter(lineage_1 == "immune")

if(ncol(data) <= 5){
		tibble(.cell = character(), lineage_2 = character(), lineage_3 = character()) |>

		# Save
		saveRDS(output_file)
} else {

	# CHANGE REFERENCE
	# reference_azimuth=
	# 	readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/reference_azimuth.rds") |>
	# 	RunUMAP(dims = 1:30, spread = 0.5,min.dist  = 0.01, n.neighbors = 10L, return.model=TRUE, umap.method = 'uwot')

	reference_azimuth = readRDS("reference_azimuth_NEW_UWOT.rds")

	data@assays@data$X = data@assays@data$X |>  as("dgCMatrix")

	# Normalise
	data =
		data |>

		# Convert
		as.Seurat(counts = "X",  data = NULL) |>
		nest(data = -.sample) |>
		filter(map_int(data, ~ ncol(.x))>5) |>
		mutate(data = map(
		data,
		~ {
browser()
			x = .x

			VariableFeatures(x) = reference_azimuth |> VariableFeatures()

			x	|>

				# Normalise RNA - not informed by smartly selected variable genes
				SCTransform(assay="originalexp") |>
				RunPCA(approx=FALSE) |>
				RunUMAP(dims = 1:30, spread = 0.5,min.dist  = 0.01, n.neighbors = 10L, return.model=TRUE, umap.method = 'uwot')
			#
			anchors <- FindTransferAnchors(
				reference = reference_azimuth,
				query = x,
				normalization.method = "SCT",
				reference.reduction = "pca",
				dims = 1:30
			)

			MapQuery(
				anchorset = anchors,
				query = x,
				reference = reference_azimuth ,
				refdata = list(
					celltype.l1 = "celltype.l1",
					celltype.l2 = "celltype.l2"
				),
				reference.reduction = "pca",
				reduction.model = "umap"
			)
		}
	)) |>
		unnest(data) |>

		# Just select essential information
		as_tibble() |>
		select(.cell, lineage_2 = predicted.celltype.l1, lineage_3 = predicted.celltype.l2, contains("refUMAP")) |>

		# Save
		saveRDS(output_file)
}

