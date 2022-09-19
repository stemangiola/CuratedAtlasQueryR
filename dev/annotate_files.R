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




# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file = args[[1]]
metadata = args[[2]]
cell_type_df = args[[3]]
output_file = args[[4]]

# Create directory
output_file |>  dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

# Read file_cell_types
data =
	loadHDF5SummarizedExperiment(input_file	)  |>

	# add lineage 1
	left_join(readRDS(metadata) |> distinct(.cell, .sample, cell_type)) |>
	left_join(read_csv(cell_type_df)) |>
	filter(lineage_1 == "immune")

if(ncol(data) <= 30){
	tibble(.cell = character()) |>

		# Save
		saveRDS(output_file)
} else {

	# CHANGE REFERENCE
	# reference_azimuth=
	# 	readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/reference_azimuth.rds") |>
	# 	RunUMAP(dims = 1:30, spread = 0.5,min.dist  = 0.01, n.neighbors = 10L, return.model=TRUE, umap.method = 'uwot')

	reference_azimuth = readRDS("reference_azimuth_NEW_UWOT.rds")

	# Selectonly interesting genes
	data = data[rownames(data) %in% VariableFeatures(reference_azimuth),]

	# Convert to Seurat matrix
	data@assays@data$X = data@assays@data$X |>  as("dgCMatrix")
	data =
		data |>

		# Convert
		as.Seurat(counts = "X",  data = NULL)

	# If I have negative values
	if((data@assays$originalexp@counts |> as.matrix() |> min()) < 0)
		data@assays$originalexp@counts[data@assays$originalexp@counts<0] <- 0

	VariableFeatures(data) = reference_azimuth |> VariableFeatures()

	# data = data	|>
	#
	# 	# Normalise RNA - not informed by smartly selected variable genes
	# 	# _Originally posted by @ChristophH in https://github.com/satijalab/seurat/issues/3618#issuecomment-719492054_
	# 	SCTransform(assay="originalexp", method = 'glmGamPoi') |>
	# 	RunPCA(approx=FALSE) |>
	# 	RunUMAP(dims = 1:30, spread = 0.5, min.dist  = 0.01, n.neighbors = 10L, return.model=TRUE, umap.method = 'uwot')
	#
	anchors <- FindTransferAnchors(
		reference = reference_azimuth,
		query = data,
		normalization.method = "SCT",
		reference.reduction = "pca",
		dims = 1:30
	)

	data =
		tryCatch(
		expr = {
			MapQuery(
				anchorset = anchors,
				query = data,
				reference = reference_azimuth ,
				refdata = list(
					celltype.l1 = "celltype.l1",
					celltype.l2 = "celltype.l2"
				),
				reference.reduction = "pca",
				reduction.model = "umap"
			)
		},
		error = function(e){
			print(e)
			data
		}
	)


	blueprint <- BlueprintEncodeData()

	annotation <-
		data |>
		as.SingleCellExperiment() |>
		SingleR(ref = blueprint, assay.type.test=1,
						labels = blueprint$label.fine)

	data |>
		left_join(
			annotation  |>
				as_tibble(rownames=".cell") |>
				select(.cell, blueprint_singler = first.labels)
		) |>

	# Just select essential information
	as_tibble() |>
		select(.cell, one_of("predicted.celltype.l1", "predicted.celltype.l2"), blueprint_singler, contains("refUMAP")) |>

		# Save
		saveRDS(output_file)
}

