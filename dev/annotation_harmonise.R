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
library(stringr)
library(purrr)

# source("utility.R")

clean_cell_types = function(.x){
	.x |>
		str_remove_all("\\+") |>
		str_remove_all("cells") |>
		str_remove_all("cell") |>
		str_remove_all("blast") |>
		str_remove_all("-") |>
		str_trim()
}

# annotation_harmonised =
# 	dir("/vast/scratch/users/mangiola.s/human_cell_atlas/annotated_data/", full.names = TRUE) |>
# 	map_dfr(readRDS) |>
#
# 	# Format
# 	mutate(across(c(predicted.celltype.l1, predicted.celltype.l2, blueprint_singler ),	tolower	)) |>
# 	mutate(across(c(predicted.celltype.l1, predicted.celltype.l2, blueprint_singler ),	clean_cell_types	))
#
#
#
# annotation_harmonised =
# 	annotation_harmonised |>
# 	mutate(strong_evidence = predicted.celltype.l2 == blueprint_singler) |>
# 	mutate(strong_evidence = case_when(
# 		strong_evidence == TRUE ~ strong_evidence,
# 		predicted.celltype.l2 == "cd14 mono" & blueprint_singler == "monocytes" ~ TRUE,
# 		predicted.celltype.l2 == "b naive" & blueprint_singler == "naive b" ~ TRUE,
# 		predicted.celltype.l2 == "cd16 mono" & blueprint_singler == "monocytes" ~ TRUE,
# 		predicted.celltype.l2 == "plasma" & blueprint_singler == "plasma" ~ TRUE,
# 		predicted.celltype.l2 == "b memory" & blueprint_singler == "class-switched memory b" ~ TRUE,
# 		predicted.celltype.l2 == "nk_cd56bright" & blueprint_singler == "nk" ~ TRUE,
# 		predicted.celltype.l2 == "treg" & blueprint_singler == "tregs" ~ TRUE,
# 		predicted.celltype.l2 == "b memory" & blueprint_singler == "memory b" ~ TRUE,
# 		predicted.celltype.l2 == "nk proliferating" & blueprint_singler == "nk" ~ TRUE,
# 		predicted.celltype.l2 == "cd8 tem" & blueprint_singler == "cd8 t" ~ TRUE,
# 		predicted.celltype.l2 == "cdc2" & blueprint_singler == "dc" ~ TRUE,
# 		predicted.celltype.l2 == "cd4 naive" & blueprint_singler == "cd4 t" ~ TRUE,
# 		predicted.celltype.l2 == "cd4 tcm" & blueprint_singler == "cd4 t" ~ TRUE,
# 		predicted.celltype.l2 == "cd8 tcm" & blueprint_singler == "cd8 t" ~ TRUE,
# 		predicted.celltype.l2 == "cd4 tem" & blueprint_singler == "cd4 t" ~ TRUE,
# 		predicted.celltype.l2 == "sdc" & blueprint_singler == "dc" ~ TRUE,
# 		predicted.celltype.l2 == "mait" & blueprint_singler == "cd8 t" ~ TRUE,
# 		TRUE ~ FALSE
# 	))

# annotation_harmonised |>
# 	filter(!strong_evidence) |>
# 	filter(!predicted.celltype.l2 %in% c("eryth", "doublet" )) |>
# 	count(predicted.celltype.l2, blueprint_singler, strong_evidence) |>
# 	arrange(desc(n)) |>
# 	slice(-c(1:200)) |>
# 	print(n=99)

# job::job({
# 	annotation_harmonised |>  saveRDS("/vast/scratch/users/mangiola.s/human_cell_atlas/annotation_harmonised.rds")
#
# })

annotation_harmonised = readRDS("/vast/scratch/users/mangiola.s/human_cell_atlas/annotation_harmonised.rds")

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

# annotation |>
# 	filter(lineage_1=="immune") |>
# 	count(cell_type, predicted.celltype.l2, blueprint_singler, strong_evidence) |>
# 	arrange(!strong_evidence, desc(n)) |>
# 	write_csv("dev/annotation_confirm.csv")


annotation_crated_confirmed =
	read_csv("dev/annotation_confirm_manually_curated.csv") |>
	filter(!is.na(azhimut_confirmed) | !is.na(blueprint_confirmed)) |>
	filter(azhimut_confirmed + blueprint_confirmed > 0) |>

	# Format
	mutate(cell_type_harmonised = case_when(
		azhimut_confirmed ~ predicted.celltype.l2,
		blueprint_confirmed ~ blueprint_singler
	)) |>

	mutate(confidence_class = 1)



# To avoid immune cell annotation if very contrasting evidence
blueprint_definitely_non_immune = c(   "astrocytes" , "chondrocytes"  , "endothelial"  ,  "epithelial" ,  "fibros"  ,  "keratinocytes" ,    "melanocytes"  , "mesangial"  ,  "mv endothelial",   "myocytes" ,  "neurons"  ,  "pericytes" ,  "preadipocytes" , "skeletal muscle"  ,  "smooth muscle"      )

annotation_crated_UNconfirmed =

	# Read
	read_csv("dev/annotation_confirm_manually_curated.csv") |>
	filter(is.na(azhimut_confirmed) | (azhimut_confirmed + blueprint_confirmed) == 0) |>

	# Annotate
	mutate(cell_type_original = cell_type) |>
	## mutate(cell_type = cell_type |> clean_cell_types()) |>
	mutate(cell_type = cell_type |> tolower()) |>
	mutate(cell_type = cell_type |> str_remove_all(",")) |>
	mutate(cell_type = cell_type |> str_remove("alphabeta")) |>
	mutate(cell_type = cell_type |> str_remove_all("positive")) |>
	mutate(cell_type = cell_type |> str_replace("cd4  t", "cd4")) |>
	mutate(cell_type = cell_type |> str_replace("regulatory t", "treg")) |>
	mutate(cell_type = cell_type |> str_remove("thymusderived")) |>
	mutate(cell_type = cell_type |> str_remove("human")) |>
	mutate(cell_type = cell_type |> str_remove("igg ")) |>
	mutate(cell_type = cell_type |> str_remove("igm ")) |>
	mutate(cell_type = cell_type |> str_remove("iga ")) |>
	mutate(cell_type = cell_type |> str_remove("group [0-9]")) |>
	mutate(cell_type = cell_type |> str_remove("common")) |>
	mutate(cell_type = cell_type |> str_remove("cd45ro")) |>
	mutate(cell_type = cell_type |> str_remove("type i")) |>
	mutate(cell_type = cell_type |> str_remove("germinal center")) |>
	mutate(cell_type = cell_type |> str_remove("iggnegative")) |>
	mutate(cell_type = cell_type |> str_remove("terminally differentiated")) |>



	mutate(cell_type = if_else(cell_type |> str_detect("macrophage"), "macrophage", cell_type) ) |>
	mutate(cell_type = if_else(cell_type |> str_detect(" treg"), "treg", cell_type) ) |>
	mutate(cell_type = if_else(cell_type |> str_detect(" dendritic"), "dendritic", cell_type) ) |>
	mutate(cell_type = if_else(cell_type |> str_detect(" thelper"), "thelper", cell_type) ) |>
	mutate(cell_type = if_else(cell_type |> str_detect("thelper "), "thelper", cell_type) ) |>
	mutate(cell_type = if_else(cell_type |> str_detect("gammadelta"), "tgd", cell_type) ) |>
	mutate(cell_type = if_else(cell_type |> str_detect("natural killer"), "nk", cell_type) ) |>


	mutate(cell_type = cell_type |> str_replace_all("  ", " ")) |>


	mutate(cell_type = cell_type |> str_replace("myeloid leukocyte", "myeloid")) |>
	mutate(cell_type = cell_type |> str_replace("effector memory", "tem")) |>
	mutate(cell_type = cell_type |> str_replace("effector", "tem")) |>
	mutate(cell_type = cell_type |> str_replace_all("cd8 t", "cd8")) |>
	mutate(cell_type = cell_type |> str_replace("central memory", "tcm")) |>
	mutate(cell_type = cell_type |> str_replace("gammadelta t", "gdt")) |>
	mutate(cell_type = cell_type |> str_replace("nonclassical monocyte", "cd16 monocyte")) |>
	mutate(cell_type = cell_type |> str_replace("classical monocyte", "cd14 monocyte")) |>
	mutate(cell_type = cell_type |> str_replace("follicular b", "b")) |>
	mutate(cell_type = cell_type |> str_replace("unswitched memory", "memory")) |>

	mutate(cell_type = cell_type |> str_trim()) |>

	mutate(cell_type_harmonised = "") |>

	# Classify strong evidence
	mutate(blueprint_confirmed = if_else(cell_type |> str_detect("cd8 cytokine secreting tem t") & blueprint_singler == "nk", T, blueprint_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type |> str_detect("cd8 cytotoxic t") & blueprint_singler == "nk",  T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type |> str_detect("cd8alphaalpha intraepithelial t") & predicted.celltype.l2 == "cd8 tem" & blueprint_singler == "cd8 tem", T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type |> str_detect("mature t") & strong_evidence & predicted.celltype.l2  |> str_detect("tem|tcm"), T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type |> str_detect("myeloid") & strong_evidence & predicted.celltype.l2  == "cd16 mono", T, azhimut_confirmed) ) |>

	# Classify weak evidence
	mutate(azhimut_confirmed = if_else(cell_type %in% c("b", "B") & predicted.celltype.l2   == "b memory" & blueprint_singler == "classswitched memory b", T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type %in% c("b", "B") & predicted.celltype.l2   %in% c("b memory", "b intermediate", "b naive", "plasma") & !blueprint_singler %in% c("classswitched memory b", "memory b", "naive b"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type %in% c("b", "B") & !predicted.celltype.l2   %in% c("b memory", "b intermediate", "b naive") & blueprint_singler %in% c("classswitched memory b", "memory b", "naive b", "plasma"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type == "activated cd4" & predicted.celltype.l2  %in% c("cd4 tcm", "cd4 tem", "tregs"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type == "activated cd4" & blueprint_singler  %in% c("cd4 tcm", "cd4 tem", "tregs"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type == "activated cd8" & predicted.celltype.l2  %in% c("cd8 tcm", "cd8 tem"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type == "activated cd8" & blueprint_singler  %in% c("cd8 tcm", "cd8 tem"), T, blueprint_confirmed) ) |>

	mutate(azhimut_confirmed = if_else(cell_type == "cd14 cd16 monocyte" & predicted.celltype.l2  %in% c("cd14 mono", "cd16 mono"), T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type == "cd14 cd16negative classical monocyte" & predicted.celltype.l2  %in% c("cd14 mono"), T, azhimut_confirmed) ) |>
	mutate(cell_type_harmonised = if_else(cell_type == "cd14 cd16negative classical monocyte" & blueprint_singler  %in% c("monocytes"), "cd14 mono", cell_type_harmonised) ) |>
	mutate(azhimut_confirmed = if_else(cell_type == "cd14 monocyte" & predicted.celltype.l2  %in% c("cd14 mono"), T, azhimut_confirmed) ) |>
	mutate(cell_type_harmonised = if_else(cell_type == "cd14 monocyte" & blueprint_singler  %in% c("monocytes"), "cd14 mono", cell_type_harmonised) ) |>
	mutate(azhimut_confirmed = if_else(cell_type == "cd14low cd16 monocyte" & predicted.celltype.l2  %in% c("cd16 mono"), T, azhimut_confirmed) ) |>
	mutate(cell_type_harmonised = if_else(cell_type == "cd14low cd16 monocyte" & blueprint_singler  %in% c("monocytes"), "cd16 mono", cell_type_harmonised) ) |>
	mutate(azhimut_confirmed = if_else(cell_type == "cd16 monocyte" & predicted.celltype.l2  %in% c("cd16 mono"), T, azhimut_confirmed) ) |>
	mutate(cell_type_harmonised = if_else(cell_type == "cd16 monocyte" & blueprint_singler  %in% c("monocytes"), "cd16 mono", cell_type_harmonised) ) |>


	mutate(blueprint_confirmed = if_else(cell_type == "monocyte" & blueprint_singler  |> str_detect("monocyte|macrophage") & !predicted.celltype.l2 |> str_detect(" mono"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type == "monocyte" & predicted.celltype.l2 |> str_detect(" mono"), T, azhimut_confirmed) ) |>


	mutate(azhimut_confirmed = if_else(cell_type == "cd4" & predicted.celltype.l2 |> str_detect("cd4|treg") & !blueprint_singler  |> str_detect("cd4"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type == "cd4" & !predicted.celltype.l2 |> str_detect("cd4") & blueprint_singler  |> str_detect("cd4|treg"), T, blueprint_confirmed) ) |>

	mutate(azhimut_confirmed = if_else(cell_type == "cd8" & predicted.celltype.l2 |> str_detect("cd8") & !blueprint_singler  |> str_detect("cd8"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type == "cd8" & !predicted.celltype.l2 |> str_detect("cd8") & blueprint_singler  |> str_detect("cd8"), T, blueprint_confirmed) ) |>


	mutate(azhimut_confirmed = if_else(cell_type == "memory t" & predicted.celltype.l2 |> str_detect("tem|tcm") & !blueprint_singler  |> str_detect("tem|tcm"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type == "memory t" & !predicted.celltype.l2 |> str_detect("tem|tcm") & blueprint_singler  |> str_detect("tem|tcm"), T, blueprint_confirmed) ) |>


	mutate(azhimut_confirmed = if_else(cell_type == "cd8alphaalpha intraepithelial t" & predicted.celltype.l2 |> str_detect("cd8") & !blueprint_singler  |> str_detect("cd8"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type == "cd8alphaalpha intraepithelial t" & !predicted.celltype.l2 |> str_detect("cd8") & blueprint_singler  |> str_detect("cd8"), T, blueprint_confirmed) ) |>

	mutate(azhimut_confirmed = if_else(cell_type == "cd8hymocyte" & predicted.celltype.l2 |> str_detect("cd8") & !blueprint_singler  |> str_detect("cd8"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type == "cd8hymocyte" & !predicted.celltype.l2 |> str_detect("cd8") & blueprint_singler  |> str_detect("cd8"), T, blueprint_confirmed) ) |>

	mutate(azhimut_confirmed = if_else(cell_type  |> str_detect("memory b") & predicted.celltype.l2 =="b memory", T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type  |> str_detect("memory b") & blueprint_singler |> str_detect("memory b"), T, blueprint_confirmed) ) |>

	mutate(azhimut_confirmed = if_else(cell_type == "immature b" & predicted.celltype.l2 =="b naive", T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type == "immature b" & blueprint_singler |> str_detect("naive b"), T, blueprint_confirmed) ) |>

	mutate(azhimut_confirmed = if_else(cell_type == "mature b" & predicted.celltype.l2 %in% c("b memory", "b intermediate"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type == "mature b" & blueprint_singler |> str_detect("memory b"), T, blueprint_confirmed) ) |>

	mutate(azhimut_confirmed = if_else(cell_type == "naive b" & predicted.celltype.l2 %in% c("b naive"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type == "naive b" & blueprint_singler |> str_detect("naive b"), T, blueprint_confirmed) ) |>

	mutate(azhimut_confirmed = if_else(cell_type == "transitional stage b" & predicted.celltype.l2 %in% c("b intermediate"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type == "transitional stage b" & blueprint_singler |> str_detect("naive b") & !predicted.celltype.l2 %in% c("b intermediate"), T, blueprint_confirmed) ) |>

	mutate(azhimut_confirmed = if_else(cell_type == "memory b" & predicted.celltype.l2 %in% c("b intermediate"), T, azhimut_confirmed) ) |>


	mutate(azhimut_confirmed = if_else(cell_type %in% c( "precursor b", "prob") & predicted.celltype.l2 %in% c("b naive") & !blueprint_singler %in% c("clp","hcs", "mpp"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type %in% c( "precursor b", "prob") & blueprint_singler |> str_detect("naive b") & predicted.celltype.l2 %in% c("hspc"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type %in% c( "precursor b", "prob") & predicted.celltype.l2 %in% c("hspc"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type %in% c( "precursor b", "prob") & blueprint_singler %in% c("clp","hcs", "mpp"), T, blueprint_confirmed) ) |>

	mutate(azhimut_confirmed = if_else(cell_type %in% c( "plasma") & predicted.celltype.l2 %in% c("b intermediate", "b memory") & blueprint_singler %in% c("classswitched memory b","memory b"), T, azhimut_confirmed) ) |>


	mutate(azhimut_confirmed = case_when(
		cell_type %in% c("cd4 cytotoxic t", "cd4 helper t") & predicted.celltype.l2 == "cd4 ctl" & blueprint_singler != "cd4 tcm" ~ T,
		cell_type %in% c("cd4 cytotoxic t", "cd4 helper t") & predicted.celltype.l2 == "cd4 tem" & blueprint_singler != "cd4 tcm" ~ T,
		TRUE ~ azhimut_confirmed
	) ) |>
	mutate(blueprint_confirmed = case_when(
		cell_type %in% c("cd4 cytotoxic t", "cd4 helper t") & blueprint_singler == "cd4 tem" & predicted.celltype.l2 != "cd4 tcm" ~ T,
		cell_type %in% c("cd4 cytotoxic t", "cd4 helper t") & blueprint_singler == "cd4 t" & predicted.celltype.l2 != "cd4 tcm" ~ T,
		TRUE ~ blueprint_confirmed
	) ) |>

	mutate(azhimut_confirmed = if_else(cell_type == "cd4hymocyte" & predicted.celltype.l2 |> str_detect("cd4|treg") & !blueprint_singler  |> str_detect("cd4"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type == "cd4hymocyte" & !predicted.celltype.l2 |> str_detect("cd4") & blueprint_singler  |> str_detect("cd4|treg"), T, blueprint_confirmed) ) |>

	mutate(azhimut_confirmed = case_when(
		cell_type %in% c("cd8 memory t") & predicted.celltype.l2 == "cd8 tem" & blueprint_singler != "cd8 tcm" ~ T,
		cell_type %in% c("cd8 memory t") & predicted.celltype.l2 == "cd8 tcm" & blueprint_singler != "cd8 tem" ~ T,
		TRUE ~ azhimut_confirmed
	) ) |>
	mutate(blueprint_confirmed = case_when(
		cell_type %in% c("cd8 memory t") & predicted.celltype.l2 != "cd8 tem" & blueprint_singler == "cd8 tcm" ~ T,
		cell_type %in% c("cd8 memory t") & predicted.celltype.l2 != "cd8 tcm" & blueprint_singler == "cd8 tem" ~ T,
		TRUE ~ blueprint_confirmed
	) ) |>

	mutate(azhimut_confirmed = case_when(
		cell_type %in% c("cd4 memory t") & predicted.celltype.l2 == "cd4 tem" & blueprint_singler != "cd8 tcm" ~ T,
		cell_type %in% c("cd4 memory t") & predicted.celltype.l2 == "cd4 tcm" & blueprint_singler != "cd8 tem" ~ T,
		TRUE ~ azhimut_confirmed
	) ) |>
	mutate(blueprint_confirmed = case_when(
		cell_type %in% c("cd4 memory t") & predicted.celltype.l2 != "cd4 tem" & blueprint_singler == "cd4 tcm" ~ T,
		cell_type %in% c("cd4 memory t") & predicted.celltype.l2 != "cd4 tcm" & blueprint_singler == "cd4 tem" ~ T,
		TRUE ~ blueprint_confirmed
	) ) |>

	mutate(azhimut_confirmed = if_else(cell_type %in% c( "t") & blueprint_singler =="cd8 t" & predicted.celltype.l2 |> str_detect("cd8"), T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type %in% c( "t") & blueprint_singler =="cd4 t" & predicted.celltype.l2 |> str_detect("cd4|treg"), T, azhimut_confirmed) ) |>

	mutate(blueprint_confirmed = if_else(cell_type %in% c( "treg") & blueprint_singler %in% c("tregs"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type %in% c( "treg") & predicted.celltype.l2 == "treg", T, azhimut_confirmed) ) |>


	mutate(blueprint_confirmed = if_else(cell_type %in% c( "tcm cd4") & blueprint_singler %in% c("cd4 tcm"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type %in% c( "tcm cd4") & predicted.celltype.l2 == "cd4 tcm", T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type %in% c( "tcm cd8") & blueprint_singler %in% c("cd8 tcm"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type %in% c( "tcm cd8") & predicted.celltype.l2 == "cd8 tcm", T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type %in% c( "tem cd4") & blueprint_singler %in% c("cd4 tem"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type %in% c( "tem cd4") & predicted.celltype.l2 == "cd4 tem", T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type %in% c( "tem cd8") & blueprint_singler %in% c("cd8 tem"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type %in% c( "tem cd8") & predicted.celltype.l2 == "cd8 tem", T, azhimut_confirmed) ) |>

	mutate(azhimut_confirmed = if_else(cell_type %in% c( "tgd") & predicted.celltype.l2 == "gdt", T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type %in% c( "activated cd4") & predicted.celltype.l2 == "cd4 proliferating", T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type %in% c( "activated cd8") & predicted.celltype.l2 == "cd8 proliferating", T, azhimut_confirmed) ) |>



	mutate(azhimut_confirmed = if_else(cell_type %in% c("naive cd4", "naive t") & predicted.celltype.l2 %in% c("cd4 naive"), T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type %in% c("naive cd8", "naive t") & predicted.celltype.l2 %in% c("cd8 naive"), T, azhimut_confirmed) ) |>

	mutate(azhimut_confirmed = if_else(cell_type %in% c( "prot") & predicted.celltype.l2 %in% c("cd4 naive") & !blueprint_singler |> str_detect("clp|hcs|mpp|cd8"), T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type %in% c( "prot") & predicted.celltype.l2 %in% c("cd8 naive") & !blueprint_singler |> str_detect("clp|hcs|mpp|cd4"), T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type %in% c( "prot") & predicted.celltype.l2 %in% c("hspc"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type %in% c( "prot") & blueprint_singler %in% c("clp","hcs", "mpp"), T, blueprint_confirmed) ) |>

	mutate(azhimut_confirmed = if_else(cell_type == "dendritic" & predicted.celltype.l2 %in% c("asdc", "cdc2", "cdc1", "pdc"), T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type == "double negative t regulatory" & predicted.celltype.l2 == "dnt", T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type %in% c( "early t lineage precursor", "immature innate lymphoid") & blueprint_singler %in% c("clp","hcs", "mpp"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type %in% c( "early t lineage precursor", "immature innate lymphoid") & predicted.celltype.l2 == "hspc" & blueprint_singler != "clp", T, azhimut_confirmed) ) |>

	mutate(blueprint_confirmed = if_else(cell_type %in% c("ilc1", "ilc2", "innate lymphoid") & blueprint_singler == "nk", T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type %in% c("ilc1", "ilc2", "innate lymphoid") & predicted.celltype.l2 %in% c( "nk", "ilc", "nk proliferating"), T, azhimut_confirmed) ) |>

	mutate(blueprint_confirmed = if_else(cell_type %in% c( "immature t") & blueprint_singler %in% c("naive t"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type %in% c( "immature t") & predicted.celltype.l2 == "t naive", T, azhimut_confirmed) ) |>

	mutate(cell_type_harmonised = if_else(cell_type == "fraction a prepro b", "naive b", cell_type_harmonised))  |>
	mutate(blueprint_confirmed = if_else(cell_type == "granulocyte" & blueprint_singler %in% c("eosinophils", "neutrophils"), T, blueprint_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type %in% c("immature neutrophil", "neutrophil") & blueprint_singler %in% c( "neutrophils"), T, blueprint_confirmed) ) |>

	mutate(blueprint_confirmed = if_else(cell_type |> str_detect("megakaryocyte") & blueprint_singler |> str_detect("megakaryocyte"), T, blueprint_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type |> str_detect("macrophage") & blueprint_singler |> str_detect("macrophage"), T, blueprint_confirmed) ) |>

	mutate(blueprint_confirmed = if_else(cell_type %in% c( "nk") & blueprint_singler %in% c("nk"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type %in% c( "nk") & predicted.celltype.l2 %in% c("nk", "nk proliferating", "nk_cd56bright", "ilc"), T, azhimut_confirmed) ) |>


	# If identical force
	mutate(azhimut_confirmed = if_else(cell_type == predicted.celltype.l2 , T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type == blueprint_singler , T, blueprint_confirmed) ) |>

	# Perogenitor
	mutate(azhimut_confirmed = if_else(cell_type  |> str_detect("progenitor|hematopoietic|precursor") & predicted.celltype.l2  == "hspc", T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type  |> str_detect("progenitor|hematopoietic|precursor") & blueprint_singler %in% c("clp","hcs", "mpp"), T, blueprint_confirmed) ) |>

	# Mast cells
	mutate(cell_type_harmonised = if_else(cell_type == "mast", "mast", cell_type_harmonised))  |>


	# Visualise
	#distinct(cell_type, predicted.celltype.l2, blueprint_singler, strong_evidence, azhimut_confirmed, blueprint_confirmed) |>
	arrange(!strong_evidence, cell_type) |>

	# set cell names
	mutate(cell_type_harmonised = case_when(
		cell_type_harmonised == "" & azhimut_confirmed ~ predicted.celltype.l2,
		cell_type_harmonised == "" & blueprint_confirmed ~ blueprint_singler,
		TRUE ~ cell_type_harmonised
	)) |>

	# Add NA
	mutate(cell_type_harmonised = case_when(cell_type_harmonised != "" ~ cell_type_harmonised)) |>

	# Add unannotated cells because datasets were too small
	mutate(cell_type_harmonised = case_when(
		is.na(cell_type_harmonised) & cell_type  |> str_detect("progenitor|hematopoietic|stem|precursor") ~ "stem",

		is.na(cell_type_harmonised) & cell_type == "cd14 monocyte" ~ "cd14 mono",
		is.na(cell_type_harmonised) & cell_type == "cd16 monocyte" ~ "cd16 mono",
		is.na(cell_type_harmonised) & cell_type %in% c("cd4 cytotoxic t", "tem cd4") ~ "cd4 tem",
		is.na(cell_type_harmonised) & cell_type %in% c("cd8 cytotoxic t", "tem cd8") ~ "cd8 tem",
		is.na(cell_type_harmonised) & cell_type == "macrophage" ~ "macrophage",
		is.na(cell_type_harmonised) & cell_type %in% c("mature b", "memory b", "plasma", "transitional stage b") ~ "b memory",
		is.na(cell_type_harmonised) & cell_type == "mucosal invariant t" ~ "mait",
		is.na(cell_type_harmonised) & cell_type == "naive b" ~ "b naive",
		is.na(cell_type_harmonised) & cell_type == "nk" ~ "nk",
		is.na(cell_type_harmonised) & cell_type == "naive cd4" ~"cd4 naive",
		is.na(cell_type_harmonised) & cell_type == "naive cd8" ~"cd8 naive",
		is.na(cell_type_harmonised) & cell_type == "treg" ~ "treg",
		is.na(cell_type_harmonised) & cell_type == "tgd" ~ "tgd",
		TRUE ~ cell_type_harmonised
	)) |>

	mutate(confidence_class = case_when(
		!is.na(cell_type_harmonised) & strong_evidence ~ 2,
		!is.na(cell_type_harmonised) & !strong_evidence ~ 3
	)) |>



	# Lowest grade annotation UNreliable
	mutate(cell_type_harmonised = case_when(

		# Get origincal annotation
		is.na(cell_type_harmonised) & cell_type %in% c("neutrophil", "granulocyte") ~ cell_type,
		is.na(cell_type_harmonised) & cell_type %in% c("conventional dendritic", "dendritic") ~ "cdc",
		is.na(cell_type_harmonised) & cell_type %in% c("classical monocyte") ~ "cd14 mono",

		# Get Seurat annotation
		is.na(cell_type_harmonised) & predicted.celltype.l2 != "eryth" & !is.na(predicted.celltype.l2) ~ predicted.celltype.l2,
		is.na(cell_type_harmonised) & !blueprint_singler %in% c(
			"astrocytes", "smooth muscle", "preadipocytes", "mesangial", "myocytes",
			"doublet", "melanocytes", "chondrocytes", "mv endothelial", "fibros",
			"neurons", "keratinocytes", "endothelial", "epithelial", "skeletal muscle", "pericytes", "erythrocytes", "adipocytes"
			) & !is.na(blueprint_singler) ~ blueprint_singler,
		TRUE ~ cell_type_harmonised

	)) |>

	# Lowest grade annotation UNreliable
	mutate(cell_type_harmonised = case_when(

		# Get origincal annotation
		!cell_type_harmonised %in% c("doublet", "platelet") ~ cell_type_harmonised

	)) |>

	mutate(confidence_class = case_when(
		is.na(confidence_class) & !is.na(cell_type_harmonised) ~ 4,
		TRUE ~ confidence_class
	)) |>


	# Rename
	rename(cell_type_temporary_shortened = cell_type) |>
	rename(cell_type = cell_type_original)

# Another passage

# annotated_samples = annotation_crated_UNconfirmed |> filter(!is.na(cell_type_harmonised)) |>  distinct( cell_type, .sample, file_id)
#
# annotation_crated_UNconfirmed |>
# 	filter(is.na(cell_type_harmonised))  |>
# 	count(cell_type ,    cell_type_harmonised ,predicted.celltype.l2 ,blueprint_singler) |>
# 	arrange(desc(n)) |>
# 	print(n=99)


annotation_all =
	annotation_crated_confirmed |>
	bind_rows(
		annotation_crated_UNconfirmed
	) |>

	# I have multiple confidence_class per combination of labels
	distinct() |>
	with_groups(c(cell_type, predicted.celltype.l2, blueprint_singler), ~ .x |> arrange(confidence_class) |> slice(1)) |>

	# Simplify after harmonisation
	mutate(cell_type_harmonised =	case_when(
		cell_type_harmonised %in% c("b memory", "b intermediate", "plasma", "classswitched memory b", "memory b" ) ~ "b memory",
		cell_type_harmonised %in% c("b naive", "naive b") ~ "b naive",
		cell_type_harmonised %in% c("nk_cd56bright", "nk", "nk proliferating", "ilc") ~ "nk",
		cell_type_harmonised %in% c("mpp", "clp", "hspc", "mep", "cmp", "hsc", "gmp") ~ "stem",
		cell_type_harmonised %in% c("macrophages",  "macrophages m1", "macrophages m2") ~ "macrophage",
		cell_type_harmonised %in% c("treg",  "tregs") ~ "treg",
		cell_type_harmonised %in% c("cd8 proliferating",  "cd8 tem") ~ "cd8 tem",
		cell_type_harmonised %in% c("cd4 proliferating",  "cd4 tem") ~ "cd4 tem",
		cell_type_harmonised %in% c("eosinophils",  "neutrophils", "granulocyte") ~ "granulocyte",
		cell_type_harmonised %in% c("cdc",  "cdc1", "cdc2", "dc") ~ "cdc",

		TRUE ~ cell_type_harmonised
	)) |>
	select(cell_type, cell_type_harmonised, predicted.celltype.l2, blueprint_singler, confidence_class) |>
	distinct()


curated_annotation =
	annotation |>
	filter(lineage_1=="immune") |>
	select(.cell, cell_type, predicted.celltype.l2, blueprint_singler) |>
	left_join(
		annotation_all ,
		by = c("cell_type", "predicted.celltype.l2", "blueprint_singler")
	) |>
	select(.cell, cell_type, cell_type_harmonised, confidence_class, cell_annotation_azimuth_l2 = predicted.celltype.l2, cell_annotation_blueprint_singler = blueprint_singler)



# Reannotation of generic cell types
reannotate_cd4 <-
	readRDS("dev/reannotate_cd4.rds")$scores |>
	as_tibble(rownames = ".cell") |>
	select("Central memory CD8 T cells", "Effector memory CD8 T cells"  , "Naive CD8 T cells" ) |>
	mutate(.cell = rownames(readRDS("dev/reannotate_cd4.rds"))) |>
	pivot_longer(
		c(`Central memory CD8 T cells`, `Effector memory CD8 T cells`  , `Naive CD8 T cells` ),
		names_to = "cell_type_Monaco", values_to = "score"
		) |>
	mutate(cell_type_Monaco = cell_type_Monaco |> str_replace_all("CD8", "CD4")) |>
	with_groups(.cell, ~ .x |> arrange(desc(score)) |> slice(1)) |>
	mutate(cell_type_Monaco = case_when(
		cell_type_Monaco == "Effector memory CD4 T cells" ~ "cd4 tem",
		cell_type_Monaco == "Central memory CD4 T cells" ~ "cd4 tcm",
		cell_type_Monaco == "Naive CD4 T cells" ~ "cd4 naive"
	)) |>
	mutate(cell_type_harmonised = "cd4 t")

reannotate_cd8 <-
	readRDS("dev/reannotate_cd8.rds")$scores |>
	as_tibble(rownames = ".cell") |>
	select("Central memory CD8 T cells", "Effector memory CD8 T cells"  , "Naive CD8 T cells" ) |>
	mutate(.cell = rownames(readRDS("dev/reannotate_cd8.rds"))) |>
	pivot_longer(
		c(`Central memory CD8 T cells`, `Effector memory CD8 T cells`  , `Naive CD8 T cells` ),
		names_to = "cell_type_Monaco", values_to = "score"
	) |>
	with_groups(.cell, ~ .x |> arrange(desc(score)) |> slice(1)) |>
	mutate(cell_type_Monaco = case_when(
		cell_type_Monaco == "Effector memory CD8 T cells" ~ "cd8 tem",
		cell_type_Monaco == "Central memory CD8 T cells" ~ "cd8 tcm",
		cell_type_Monaco == "Naive CD8 T cells" ~ "cd8 naive"
	)) |>
	mutate(cell_type_harmonised = "cd8 t")

reannotate_monocytes <-
	readRDS("dev/reannotate_monocytes.rds")$scores |>
	as_tibble(rownames = ".cell") |>
	select("Non classical monocytes", "Classical monocytes"   ) |>
	mutate(.cell = rownames(readRDS("dev/reannotate_monocytes.rds"))) |>
	pivot_longer(
		c(`Non classical monocytes`, `Classical monocytes`   ),
		names_to = "cell_type_Monaco", values_to = "score"
	) |>
	with_groups(.cell, ~ .x |> arrange(desc(score)) |> slice(1)) |>
	mutate(cell_type_Monaco = case_when(
		cell_type_Monaco == "Non classical monocytes" ~ "cd16 mono",
		cell_type_Monaco == "Classical monocytes" ~ "cd14 mono"
	)) |>
	mutate(cell_type_harmonised = "monocytes")

library(glue)


curated_annotation =

	# Fix cell ID
	get_metadata() |>
	select(.cell, .sample) |>
	as_tibble() |>
	mutate(.cell_combined = glue("{.cell}_{.sample}")) |>

	# Add cell type
	inner_join(
		curated_annotation,
		by=c(".cell_combined" = ".cell")
	) |>
	left_join(
		reannotate_cd4 |>
			bind_rows(reannotate_cd8) |>
			bind_rows(reannotate_monocytes)
	) |>
	mutate(cell_type_harmonised = if_else(
		!is.na(cell_type_Monaco),
		cell_type_Monaco,
		cell_type_harmonised
	)) |>
	select(-.cell_combined) |>
	select(-cell_type_Monaco, -score) |> 
  
  # Replace NA
  mutate(cell_type_harmonised = if_else(is.na(cell_type_harmonised), "immune_unclassified", cell_type_harmonised)) |>
  
  # Add non immune
  select(-cell_type) |> 
  full_join(
    get_metadata("dev/metadata.SQLite") |> 
      select(.cell, .sample) |> 
      as_tibble()
  ) |> 
  mutate(cell_type_harmonised = if_else(is.na(cell_type_harmonised), "non_immune", cell_type_harmonised))

# Save
job::job({
	curated_annotation |>
    mutate(across(contains("cell_"), as.factor)) |> 
	saveRDS("dev/curated_annotation.rds")
})


cell_metadata_with_harmonised_annotation =
	curated_annotation |>
  mutate(.cell = .cell |> str_remove(.sample) |> str_remove("_$")) |> 
	left_join(
		get_metadata("dev/metadata.SQLite") |>
			select(.cell, .sample, file_id, file_id_db, tissue, disease, is_primary_data.x, is_primary_data.y, name) |>
		  left_join(read_csv("dev/tissue_label_curated.csv"), copy=TRUE) |> 
		  as_tibble(),
		by=c(".cell", ".sample")
	) |> 
  distinct()  |> 
  
  # Drop secondary data often cell type subsets
  filter(is_primary_data.x==TRUE)

# Filter samples that do not have immune
cell_metadata_with_harmonised_annotation = 
  cell_metadata_with_harmonised_annotation |> 
  nest(data = -c(.sample, tissue_harmonised)) |> 
  filter(map_int(data, ~ .x |> filter(cell_type_harmonised != "non_immune") |> nrow()) > 0) |> 
  unnest(data)

# Tissue with no immune
cell_metadata_with_harmonised_annotation = 
  cell_metadata_with_harmonised_annotation |> 
  filter(disease == "normal") |> 

  # Filter tissues
  nest(data = -c(cell_type_harmonised, tissue_harmonised)) |> 
  add_count(tissue_harmonised, name = "n_cell_type_in_tissue") |> 
  filter(n_cell_type_in_tissue>=18) |> 
  add_count(cell_type_harmonised, name = "n_tissue_in_cell_type") |> 
  filter(n_tissue_in_cell_type>=26) |> 
  unnest(data)


cell_metadata_with_harmonised_annotation |> 
  saveRDS("dev/cell_metadata_with_harmonised_annotation.rds")



annotated_samples =
	cell_metadata_with_harmonised_annotation |> filter(!is.na(cell_type_harmonised)) |>  distinct( cell_type, .sample, file_id)

# Cell types that most need attention
cell_metadata_with_harmonised_annotation |> anti_join(annotated_samples) |> select(contains("cell_")) |> count(cell_type ,    cell_type_harmonised ,cell_annotation_azimuth_l2 ,cell_annotation_blueprint_singler) |> arrange(desc(n)) |> print(n=99)

# How many samples miss annotation
cell_metadata_with_harmonised_annotation |> anti_join(annotated_samples)  |>  distinct(cell_type,  cell_type_harmonised, .sample) |>  distinct( cell_type, .sample) |> count(cell_type)  |> arrange(desc(n))

# Histo of annotation
cell_metadata_with_harmonised_annotation |> filter(!is.na(cell_type_harmonised)) |>  distinct( cell_type_harmonised, .sample) |> count(.sample) |> pull(n) |> hist(breaks=30)

# NEEDED CHANGES - ANIMAL CELLS AND HEMATOPOIETIC SHOULD BE EVALUATED AS IMMUNE FOR REANNNOTATION?

# Tissue with no immune

temp_tissue = 
  cell_metadata_with_harmonised_annotation |> 
  filter(cell_type_harmonised=="pdc") |> 
  pull(tissue_harmonised)  
  
cell_metadata_with_harmonised_annotation |> 
  filter(!tissue_harmonised %in% temp_tissue) |> 
  distinct(tissue_harmonised) |> 
  

cell_metadata_with_harmonised_annotation |> 
  filter(disease == "normal") |> 
  filter(!is.na(cell_type_harmonised)) |>  
  distinct( cell_type_harmonised, tissue_harmonised) |>
  count(tissue_harmonised) |> arrange(n) |> print(n=99)


