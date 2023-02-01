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
library(HCAquery)
library(BiocParallel)


# Read arguments
args = commandArgs(trailingOnly=TRUE)
sample_files = args[1:(length(args)-6)]
file_id_column = args[[length(args)-5]]
light_data_directory = args[[length(args)-4]]
file_for_annotation_workflow = args[[length(args)-3]]
cell_type_df = args[[length(args)-2]]
output_dir = args[[length(args)-1]]
memory_Mb = args[[length(args)]] |> as.integer()

library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize = (memory_Mb - 10000)  * 1024^2)

print(file_id_column )
			print(light_data_directory )
						print(file_for_annotation_workflow )
									print(cell_type_df )
												print(output_dir)
# Create directory
output_dir |> dir.create( showWarnings = FALSE, recursive = TRUE)

reference_azimuth = readRDS("reference_azimuth_NEW_UWOT.rds")

metadata_few_columns = readRDS(file_for_annotation_workflow)

data =
	sample_files |>
	enframe(value = "input_file") |>
	mutate(file_id = file_id_column) |>
	mutate(sce = map2(
		input_file, file_id,
		~ {
			.x |>
				loadHDF5SummarizedExperiment() |>
			#sce[rownames(sce) %in% VariableFeatures(reference_azimuth),] |>

				# Filter Immune
				left_join(metadata_few_columns |> filter(file_id == .y) |> select(-file_id)) |>
				left_join(read_csv(cell_type_df)) |>
				filter(lineage_1 == "immune")
	})) |>

	unnest_single_cell_experiment(sce)

	# pull(sce) %>%
	# do.call(cbind, .)

if(ncol(data) <= 30){

	data |>
		distinct(file_id, .sample) |>
		mutate(
			saved = map(.sample, ~ 	tibble(.cell = character()) |> saveRDS(glue("{output_dir}/{.x}.rds")))
			)

} else {

	# CHANGE REFERENCE
	# reference_azimuth=
	# 	readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/reference_azimuth.rds") |>
	# 	RunUMAP(dims = 1:30, spread = 0.5,min.dist  = 0.01, n.neighbors = 10L, return.model=TRUE, umap.method = 'uwot')


	data_seurat = data

	# Selectonly interesting genes
	data_seurat = data_seurat[rownames(data_seurat) %in% VariableFeatures(reference_azimuth),]

	# Convert to Seurat matrix
	data_seurat@assays@data$X = data_seurat@assays@data$X |>  as("dgCMatrix")
	data_seurat =
		data_seurat |>

		# Convert
		as.Seurat(counts = "X",  data = NULL)

	# If I have negative values
	if((data_seurat@assays$originalexp@counts |> as.matrix() |> min()) < 0)
		data_seurat@assays$originalexp@counts[data_seurat@assays$originalexp@counts<0] <- 0

	VariableFeatures(data_seurat) = reference_azimuth |> VariableFeatures()

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
		query = data_seurat,
		normalization.method = "SCT",
		reference.reduction = "pca",
		dims = 1:30
	)

	data_seurat =
		tryCatch(
		expr = {
			MapQuery(
				anchorset = anchors,
				query = data_seurat,
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
			data_seurat
		}
	) |>
		as_tibble()

	gc()

	blueprint <- BlueprintEncodeData()

	library(scuttle)

	annotation_blueprint <-
		data |>
		logNormCounts(assay.type = "X") |>
		SingleR(ref = blueprint, assay.type.test=1,
						labels = blueprint$label.fine,
						BPPARAM=MulticoreParam(4)
					)

	rm(blueprint)
	gc()

	MonacoImmuneData = MonacoImmuneData()

	annotation_monaco <-
		data |>
		logNormCounts(assay.type = "X") |>
		SingleR(ref = MonacoImmuneData, assay.type.test=1,
						labels = MonacoImmuneData$label.fine,
						BPPARAM=MulticoreParam(4)
						)

	rm(data)
	gc()

	data_seurat |>
		left_join(
			annotation_blueprint  |>
				as_tibble(rownames=".cell") |>
				select(.cell, blueprint_singler = labels)
		) |>
		left_join(
			annotation_monaco  |>
				as_tibble(rownames=".cell") |>
				select(.cell, monaco_singler = labels)
		) |>

		# Just select essential information
		as_tibble() |>
		select(.cell, .sample, one_of("predicted.celltype.l1", "predicted.celltype.l2"), blueprint_singler, monaco_singler, contains("refUMAP")) |>

		# Save
		nest(data = -.sample) |>
		mutate(saved = map2(
			data, .sample,
			~ .x |> saveRDS(glue("{output_dir}/{.y}.rds"))
		))


}

