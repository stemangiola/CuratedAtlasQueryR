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

clean_cell_types = function(.x){
	.x |>
		str_remove_all("\\+") |>
		str_remove_all("cells") |>
		str_remove_all("cell") |>
		str_remove_all("blast") |>
		str_remove_all("-") |>
		str_trim()
}

annotation_harmonised =
	dir("/vast/scratch/users/mangiola.s/human_cell_atlas/annotated_data/", full.names = TRUE) |>
	map_dfr(readRDS) |>

	# Format
	mutate(across(c(predicted.celltype.l1, predicted.celltype.l2, blueprint_singler ),	tolower	)) |>
	mutate(across(c(predicted.celltype.l1, predicted.celltype.l2, blueprint_singler ),	clean_cell_types	))



annotation_harmonised =
	annotation_harmonised |>
	mutate(strong_evidence = predicted.celltype.l2 == blueprint_singler) |>
	mutate(strong_evidence = case_when(
		strong_evidence == TRUE ~ strong_evidence,
		predicted.celltype.l2 == "cd14 mono" & blueprint_singler == "monocytes" ~ TRUE,
		predicted.celltype.l2 == "b naive" & blueprint_singler == "naive b" ~ TRUE,
		predicted.celltype.l2 == "cd16 mono" & blueprint_singler == "monocytes" ~ TRUE,
		predicted.celltype.l2 == "plasma" & blueprint_singler == "plasma" ~ TRUE,
		predicted.celltype.l2 == "b memory" & blueprint_singler == "class-switched memory b" ~ TRUE,
		predicted.celltype.l2 == "nk_cd56bright" & blueprint_singler == "nk" ~ TRUE,
		predicted.celltype.l2 == "treg" & blueprint_singler == "tregs" ~ TRUE,
		predicted.celltype.l2 == "b memory" & blueprint_singler == "memory b" ~ TRUE,
		predicted.celltype.l2 == "nk proliferating" & blueprint_singler == "nk" ~ TRUE,
		predicted.celltype.l2 == "cd8 tem" & blueprint_singler == "cd8 t" ~ TRUE,
		predicted.celltype.l2 == "cdc2" & blueprint_singler == "dc" ~ TRUE,
		predicted.celltype.l2 == "cd4 naive" & blueprint_singler == "cd4 t" ~ TRUE,
		predicted.celltype.l2 == "cd4 tcm" & blueprint_singler == "cd4 t" ~ TRUE,
		predicted.celltype.l2 == "cd8 tcm" & blueprint_singler == "cd8 t" ~ TRUE,
		predicted.celltype.l2 == "cd4 tem" & blueprint_singler == "cd4 t" ~ TRUE,
		predicted.celltype.l2 == "sdc" & blueprint_singler == "dc" ~ TRUE,
		predicted.celltype.l2 == "mait" & blueprint_singler == "cd8 t" ~ TRUE,
		TRUE ~ FALSE
	))

# annotation_harmonised |>
# 	filter(!strong_evidence) |>
# 	filter(!predicted.celltype.l2 %in% c("eryth", "doublet" )) |>
# 	count(predicted.celltype.l2, blueprint_singler, strong_evidence) |>
# 	arrange(desc(n)) |>
# 	slice(-c(1:200)) |>
# 	print(n=99)

job::job({
	annotation_harmonised |>  saveRDS("/vast/scratch/users/mangiola.s/human_cell_atlas/annotation_harmonised.rds")

})

library(HCAquery)
metadata_df = get_metadata()

# Integrate with metadata

annotation =
	metadata_df |>
	select(.cell, cell_type, file_id) |>
	as_tibble() |>
	left_join(read_csv("dev/metadata_cell_type.csv"),  by = "cell_type") |>
	left_join(annotation_harmonised, by = ".cell") |>

	# Clen cell types
	mutate(cell_type = cell_type |> clean_cell_types())

annotation |>
	filter(lineage_1=="immune") |>
	count(cell_type, predicted.celltype.l2, blueprint_singler, strong_evidence) |>
	arrange(!strong_evidence, desc(n)) |>
	write_csv("dev/annotation_confirm.csv")
