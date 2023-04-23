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
library(CuratedAtlasQueryR)
library(BiocParallel)
library(scuttle)

# Read arguments
args = commandArgs(trailingOnly=TRUE)
sample_files = args[1:(length(args)-6)]
file_id_column = args[[length(args)-5]]
light_data_directory = args[[length(args)-4]]
file_for_annotation_workflow = args[[length(args)-3]]
cell_type_df = args[[length(args)-2]]
output_dir = args[[length(args)-1]]
memory_Mb = args[[length(args)]] |> as.integer()

# library(future)
# plan("multicore", workers = 4)
# options(future.globals.maxSize = (memory_Mb - 10000)  * 1024^2)

# print(file_id_column )
# 			print(light_data_directory )
# 						print(file_for_annotation_workflow )
# 									print(cell_type_df )
# 												print(output_dir)

# Create directory
output_dir |> dir.create( showWarnings = FALSE, recursive = TRUE)

metadata_few_columns = readRDS(file_for_annotation_workflow)


# SINGLER
blueprint <- BlueprintEncodeData()
MonacoImmuneData = MonacoImmuneData()

print("Start SingleR")

data_singler =
	sample_files |>
	enframe(value = "input_file") |>
	mutate(file_id = file_id_column) |>
	mutate(singler = map2(
		input_file, file_id,
		~ {
print(.x)
			sce =
				.x |>
				loadHDF5SummarizedExperiment() |>
				#sce[rownames(sce) %in% VariableFeatures(reference_azimuth),] |>

				# Filter Immune
				left_join(metadata_few_columns |> filter(file_id == .y) |> select(-file_id)) |>
				left_join(read_csv(cell_type_df)) |>
				filter(lineage_1 == "immune") |>

				logNormCounts(assay.type = "X")

			# If no cell kill
			if(ncol(sce)==0)
				return(tibble(.cell = character()))

			# If I have negative values
			if((	sce@assays@data$X[,1:min(100, ncol(sce))] |>  min()) < 0)
				sce@assays@data$X[	sce@assays@data$X<0] <- 0

			if(ncol(sce)==1){
				sce = cbind(sce, sce)
				colnames(sce)[2]= "dummy___"
			}

			left_join(
				sce |>
					SingleR(
						ref = blueprint,
						assay.type.test= 1,
						labels = blueprint$label.fine
					)  |>
					as_tibble(rownames=".cell") |>
					nest(blueprint_scores = starts_with("score")) |>
					select(-one_of("delta.next"),- pruned.labels) |>
					rename( blueprint_singler = labels),

				sce |>
					SingleR(
						ref = MonacoImmuneData,
						assay.type.test= 1,
						labels = MonacoImmuneData$label.fine
					)  |>
					as_tibble(rownames=".cell") |>

					nest(monaco_scores = starts_with("score")) |>
					select(-delta.next,- pruned.labels) |>
					rename( monaco_singler = labels)
			) |>
				filter(.cell!="dummy___")

		})) |>

	unnest(singler) |>
	select(-name)

# If not immune cells
if(nrow(data_singler) == 0){

	sample_files |>
		enframe(value = "file") |>
		mutate(.sample = file |>  basename() |> tools::file_path_sans_ext()) |>
		mutate(
			saved = map(.sample, ~ 	tibble(.cell = character()) |> saveRDS(glue("{output_dir}/{.x}.rds")))
		)

} else if(nrow(data_singler) <= 30){

# If too little immune cells


	data_singler |>
		mutate(.sample = input_file |>  basename() |> tools::file_path_sans_ext()) |>
		nest(data = -c(.sample, file_id)) |>
		mutate(
			saved = map(.sample, ~ 	tibble(.cell = character()) |> saveRDS(glue("{output_dir}/{.x}.rds")))
		)


} else{

print("Start Seurat")

reference_azimuth = readRDS("reference_azimuth_NEW_UWOT.rds")


data_seurat =
	sample_files |>

	enframe(value = "input_file") |>
	mutate(file_id = file_id_column) |>
	mutate(seu = map2(
		input_file, file_id,
		~ {
			sce =
				.x |>
				loadHDF5SummarizedExperiment() |>
				#sce[rownames(sce) %in% VariableFeatures(reference_azimuth),] |>

				# Filter Immune
				left_join(metadata_few_columns |> filter(file_id == .y) |> select(-file_id)) |>
				left_join(read_csv(cell_type_df)) |>
				filter(lineage_1 == "immune") %>%

				.[rownames(.) %in% VariableFeatures(reference_azimuth),]

			sce@assays@data$X = sce@assays@data$X |>  as("dgCMatrix")

			sce |>
				as.Seurat(counts = "X",  data = NULL)
		})) |>

	unnest_seurat(seu)


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
		as_tibble() |>
		select(.cell, .sample, one_of("predicted.celltype.l1", "predicted.celltype.l2"), contains("refUMAP"))


	data_seurat |>

		left_join(	data_singler	, by = join_by(.cell)	) |>

		# Just select essential information

		# Save
		nest(data = -c(.sample, file_id)) |>
		mutate(saved = map2(
			data, .sample,
			~ .x |> saveRDS(glue("{output_dir}/{.y}.rds"))
		))


}

