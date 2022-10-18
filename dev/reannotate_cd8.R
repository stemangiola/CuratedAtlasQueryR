library(scuttle)
library(celldex)
library(SingleR)
library(BiocParallel)
library(HCAquery)
library(SingleCellExperiment)
library(tidySingleCellExperiment)
library(digest)
library(openssl)
library(stringr)
library(HDF5Array)
monocyte_reference =celldex::MonacoImmuneData()


readRDS("/vast/projects/RCP/human_cell_atlas/metadata.rds") |>
	filter(.cell %in% (
		readRDS("dev/curated_annotation.rds") |>
			filter(cell_type_harmonised=="cd8 t") |>
			pull(.cell)
	)) |>
	unite("file_id_db", c(file_id, cell_type), remove = FALSE) |>
	mutate(file_id_db = file_id_db |> md5() |> as.character()) |>

	get_SingleCellExperiment("/vast/projects/RCP/human_cell_atlas/splitted_DB2_data") |>

	# Annotate
	logNormCounts(assay.type = "X") |>

	SingleR(
		ref = monocyte_reference,
		assay.type.test=1,
		labels = monocyte_reference$label.fine
	) |>
	saveRDS("reannotate_cd8.rds")
