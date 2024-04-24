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

metadata_file = "/vast/projects/cellxgene_curated//metadata_0.2.rds"
file_curated_annotation_merged = "~/PostDoc/CuratedAtlasQueryR/dev/cell_type_curated_annotation_0.2.3.rds"
file_metadata_annotated = "/vast/projects/cellxgene_curated/metadata_annotated_0.2.3.rds"
annotation_directory = "/vast/projects/cellxgene_curated//annotated_data_0.2/"

# metadata_file = "/vast/projects/cellxgene_curated//metadata.rds"
# file_curated_annotation_merged = "~/PostDoc/CuratedAtlasQueryR/dev/cell_type_curated_annotation.rds"
# file_metadata_annotated = "/vast/projects/cellxgene_curated//metadata_annotated.rds"
# annotation_directory = "/vast/projects/cellxgene_curated//annotated_data_0.1/"


annotation_harmonised =
	dir(annotation_directory, full.names = TRUE) |>
	enframe(value="file") |>
	tidyr::extract(  file,".sample", "/([a-z0-9]+)\\.rds", remove = F) |>
	mutate(data = map(file, ~ .x |> readRDS() |> select(-contains("score")) )) |>
	unnest(data) |>

	# Format
	mutate(across(c(predicted.celltype.l1, predicted.celltype.l2, blueprint_singler, monaco_singler, ),	tolower	)) |>
	mutate(across(c(predicted.celltype.l1, predicted.celltype.l2, blueprint_singler, monaco_singler, ),	clean_cell_types	)) |>

	# Format
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



job::job({
	annotation_harmonised |>  saveRDS("~/PostDoc/CuratedAtlasQueryR/dev/annotated_data_0.2_temp_table.rds")
})

clean_cell_types_deeper = function(x){
  x |> 
    # Annotate
    mutate(cell_type_clean = cell_type_clean |> tolower()) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove_all(",")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("alphabeta")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove_all("positive")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("cd4  t", "cd4")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("regulatory t", "treg")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("thymusderived")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("human")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("igg ")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("igm ")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("iga ")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("group [0-9]")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("common")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("cd45ro")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("type i")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("germinal center")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("iggnegative")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("terminally differentiated")) |>
    
    mutate(cell_type_clean = if_else(cell_type_clean |> str_detect("macrophage"), "macrophage", cell_type_clean) ) |>
    mutate(cell_type_clean = if_else(cell_type_clean == "mononuclear phagocyte", "macrophage", cell_type_clean) ) |>
    
    mutate(cell_type_clean = if_else(cell_type_clean |> str_detect(" treg"), "treg", cell_type_clean) ) |>
    mutate(cell_type_clean = if_else(cell_type_clean |> str_detect(" dendritic"), "dendritic", cell_type_clean) ) |>
    mutate(cell_type_clean = if_else(cell_type_clean |> str_detect(" thelper"), "thelper", cell_type_clean) ) |>
    mutate(cell_type_clean = if_else(cell_type_clean |> str_detect("thelper "), "thelper", cell_type_clean) ) |>
    mutate(cell_type_clean = if_else(cell_type_clean |> str_detect("gammadelta"), "tgd", cell_type_clean) ) |>
    mutate(cell_type_clean = if_else(cell_type_clean |> str_detect("natural killer"), "nk", cell_type_clean) ) |>
    
    
    mutate(cell_type_clean = cell_type_clean |> str_replace_all("  ", " ")) |>
    
    
    mutate(cell_type_clean = cell_type_clean |> str_replace("myeloid leukocyte", "myeloid")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("effector memory", "tem")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("effector", "tem")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace_all("cd8 t", "cd8")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("central memory", "tcm")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("gammadelta t", "gdt")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("nonclassical monocyte", "cd16 monocyte")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("classical monocyte", "cd14 monocyte")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("follicular b", "b")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("unswitched memory", "memory")) |>
    
    mutate(cell_type_clean = cell_type_clean |> str_trim()) 
}


# annotation_harmonised = readRDS("~/PostDoc/CuratedAtlasQueryR/dev/annotated_data_0.2_temp_table.rds")

# library(CuratedAtlasQueryR)
metadata_df = readRDS(metadata_file)

# Integrate with metadata

annotation =
	metadata_df |>
	select(.cell, cell_type, file_id, .sample) |>
	as_tibble() |>
	left_join(read_csv("~/PostDoc/CuratedAtlasQueryR/dev/metadata_cell_type.csv"),  by = "cell_type") |>
	left_join(annotation_harmonised, by = c(".cell", ".sample")) |>

	# Clen cell types
	mutate(cell_type_clean = cell_type |> clean_cell_types())

# annotation |>
# 	filter(lineage_1=="immune") |>
# 	count(cell_type, predicted.celltype.l2, blueprint_singler, strong_evidence) |>
# 	arrange(!strong_evidence, desc(n)) |>
# 	write_csv("~/PostDoc/CuratedAtlasQueryR/dev/annotation_confirm.csv")


annotation_crated_confirmed =
	read_csv("~/PostDoc/CuratedAtlasQueryR/dev/annotation_confirm_manually_curated.csv") |>

	# TEMPORARY
	rename(cell_type_clean = cell_type) |>

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
	read_csv("~/PostDoc/CuratedAtlasQueryR/dev/annotation_confirm_manually_curated.csv") |>

	# TEMPORARY
	rename(cell_type_clean = cell_type) |>

	filter(is.na(azhimut_confirmed) | (azhimut_confirmed + blueprint_confirmed) == 0) |>

  clean_cell_types_deeper() |> 

	mutate(cell_type_harmonised = "") |>

	# Classify strong evidence
	mutate(blueprint_confirmed = if_else(cell_type_clean |> str_detect("cd8 cytokine secreting tem t") & blueprint_singler == "nk", T, blueprint_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean |> str_detect("cd8 cytotoxic t") & blueprint_singler == "nk",  T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean |> str_detect("cd8alphaalpha intraepithelial t") & predicted.celltype.l2 == "cd8 tem" & blueprint_singler == "cd8 tem", T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean |> str_detect("mature t") & strong_evidence & predicted.celltype.l2  |> str_detect("tem|tcm"), T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean |> str_detect("myeloid") & strong_evidence & predicted.celltype.l2  == "cd16 mono", T, azhimut_confirmed) ) |>

	# Classify weak evidence
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c("b", "B") & predicted.celltype.l2   == "b memory" & blueprint_singler == "classswitched memory b", T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c("b", "B") & predicted.celltype.l2   %in% c("b memory", "b intermediate", "b naive", "plasma") & !blueprint_singler %in% c("classswitched memory b", "memory b", "naive b"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean %in% c("b", "B") & !predicted.celltype.l2   %in% c("b memory", "b intermediate", "b naive") & blueprint_singler %in% c("classswitched memory b", "memory b", "naive b", "plasma"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean == "activated cd4" & predicted.celltype.l2  %in% c("cd4 tcm", "cd4 tem", "tregs"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean == "activated cd4" & blueprint_singler  %in% c("cd4 tcm", "cd4 tem", "tregs"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean == "activated cd8" & predicted.celltype.l2  %in% c("cd8 tcm", "cd8 tem"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean == "activated cd8" & blueprint_singler  %in% c("cd8 tcm", "cd8 tem"), T, blueprint_confirmed) ) |>

	# Monocyte macrophage
	mutate(azhimut_confirmed = if_else(cell_type_clean == "cd14 cd16 monocyte" & predicted.celltype.l2  %in% c("cd14 mono", "cd16 mono"), T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean == "cd14 cd16negative classical monocyte" & predicted.celltype.l2  %in% c("cd14 mono"), T, azhimut_confirmed) ) |>
	mutate(cell_type_harmonised = if_else(cell_type_clean == "cd14 cd16negative classical monocyte" & blueprint_singler  %in% c("monocytes"), "cd14 mono", cell_type_harmonised) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean == "cd14 monocyte" & predicted.celltype.l2  %in% c("cd14 mono"), T, azhimut_confirmed) ) |>
	mutate(cell_type_harmonised = if_else(cell_type_clean == "cd14 monocyte" & blueprint_singler  %in% c("monocytes"), "cd14 mono", cell_type_harmonised) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean == "cd14low cd16 monocyte" & predicted.celltype.l2  %in% c("cd16 mono"), T, azhimut_confirmed) ) |>
	mutate(cell_type_harmonised = if_else(cell_type_clean == "cd14low cd16 monocyte" & blueprint_singler  %in% c("monocytes"), "cd16 mono", cell_type_harmonised) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean == "cd16 monocyte" & predicted.celltype.l2  %in% c("cd16 mono"), T, azhimut_confirmed) ) |>
	mutate(cell_type_harmonised = if_else(cell_type_clean == "cd16 monocyte" & blueprint_singler  %in% c("monocytes"), "cd16 mono", cell_type_harmonised) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean == "monocyte" & blueprint_singler  |> str_detect("monocyte|macrophage") & !predicted.celltype.l2 |> str_detect(" mono"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean == "monocyte" & predicted.celltype.l2 |> str_detect(" mono"), T, azhimut_confirmed) ) |>


	mutate(azhimut_confirmed = if_else(cell_type_clean == "cd4" & predicted.celltype.l2 |> str_detect("cd4|treg") & !blueprint_singler  |> str_detect("cd4"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean == "cd4" & !predicted.celltype.l2 |> str_detect("cd4") & blueprint_singler  |> str_detect("cd4|treg"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean == "cd8" & predicted.celltype.l2 |> str_detect("cd8") & !blueprint_singler  |> str_detect("cd8"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean == "cd8" & !predicted.celltype.l2 |> str_detect("cd8") & blueprint_singler  |> str_detect("cd8"), T, blueprint_confirmed) ) |>


	mutate(azhimut_confirmed = if_else(cell_type_clean == "memory t" & predicted.celltype.l2 |> str_detect("tem|tcm") & !blueprint_singler  |> str_detect("tem|tcm"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean == "memory t" & !predicted.celltype.l2 |> str_detect("tem|tcm") & blueprint_singler  |> str_detect("tem|tcm"), T, blueprint_confirmed) ) |>


	mutate(azhimut_confirmed = if_else(cell_type_clean == "cd8alphaalpha intraepithelial t" & predicted.celltype.l2 |> str_detect("cd8") & !blueprint_singler  |> str_detect("cd8"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean == "cd8alphaalpha intraepithelial t" & !predicted.celltype.l2 |> str_detect("cd8") & blueprint_singler  |> str_detect("cd8"), T, blueprint_confirmed) ) |>

	mutate(azhimut_confirmed = if_else(cell_type_clean == "cd8hymocyte" & predicted.celltype.l2 |> str_detect("cd8") & !blueprint_singler  |> str_detect("cd8"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean == "cd8hymocyte" & !predicted.celltype.l2 |> str_detect("cd8") & blueprint_singler  |> str_detect("cd8"), T, blueprint_confirmed) ) |>

	# B cells
	mutate(azhimut_confirmed = if_else(cell_type_clean  |> str_detect("memory b") & predicted.celltype.l2 =="b memory", T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean  |> str_detect("memory b") & blueprint_singler |> str_detect("memory b"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean == "immature b" & predicted.celltype.l2 =="b naive", T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean == "immature b" & blueprint_singler |> str_detect("naive b"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean == "mature b" & predicted.celltype.l2 %in% c("b memory", "b intermediate"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean == "mature b" & blueprint_singler |> str_detect("memory b"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean == "naive b" & predicted.celltype.l2 %in% c("b naive"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean == "naive b" & blueprint_singler |> str_detect("naive b"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean == "transitional stage b" & predicted.celltype.l2 %in% c("b intermediate"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean == "transitional stage b" & blueprint_singler |> str_detect("naive b") & !predicted.celltype.l2 %in% c("b intermediate"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean == "memory b" & predicted.celltype.l2 %in% c("b intermediate"), T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "precursor b", "prob") & predicted.celltype.l2 %in% c("b naive") & !blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "precursor b", "prob") & blueprint_singler |> str_detect("naive b") & predicted.celltype.l2 %in% c("hspc"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "precursor b", "prob") & predicted.celltype.l2 %in% c("hspc"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "precursor b", "prob") & blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), T, blueprint_confirmed) ) |>

	# Plasma cells
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "plasma") & predicted.celltype.l2 == "plasma" , T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "plasma") & predicted.celltype.l2 == "plasma" , T, blueprint_confirmed) ) |>

	mutate(azhimut_confirmed = case_when(
		cell_type_clean %in% c("cd4 cytotoxic t", "cd4 helper t") & predicted.celltype.l2 == "cd4 ctl" & blueprint_singler != "cd4 tcm" ~ T,
		cell_type_clean %in% c("cd4 cytotoxic t", "cd4 helper t") & predicted.celltype.l2 == "cd4 tem" & blueprint_singler != "cd4 tcm" ~ T,
		TRUE ~ azhimut_confirmed
	) ) |>
	mutate(blueprint_confirmed = case_when(
		cell_type_clean %in% c("cd4 cytotoxic t", "cd4 helper t") & blueprint_singler == "cd4 tem" & predicted.celltype.l2 != "cd4 tcm" ~ T,
		cell_type_clean %in% c("cd4 cytotoxic t", "cd4 helper t") & blueprint_singler == "cd4 t" & predicted.celltype.l2 != "cd4 tcm" ~ T,
		TRUE ~ blueprint_confirmed
	) ) |>

	mutate(azhimut_confirmed = if_else(cell_type_clean == "cd4hymocyte" & predicted.celltype.l2 |> str_detect("cd4|treg") & !blueprint_singler  |> str_detect("cd4"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean == "cd4hymocyte" & !predicted.celltype.l2 |> str_detect("cd4") & blueprint_singler  |> str_detect("cd4|treg"), T, blueprint_confirmed) ) |>

	mutate(azhimut_confirmed = case_when(
		cell_type_clean %in% c("cd8 memory t") & predicted.celltype.l2 == "cd8 tem" & blueprint_singler != "cd8 tcm" ~ T,
		cell_type_clean %in% c("cd8 memory t") & predicted.celltype.l2 == "cd8 tcm" & blueprint_singler != "cd8 tem" ~ T,
		TRUE ~ azhimut_confirmed
	) ) |>
	mutate(blueprint_confirmed = case_when(
		cell_type_clean %in% c("cd8 memory t") & predicted.celltype.l2 != "cd8 tem" & blueprint_singler == "cd8 tcm" ~ T,
		cell_type_clean %in% c("cd8 memory t") & predicted.celltype.l2 != "cd8 tcm" & blueprint_singler == "cd8 tem" ~ T,
		TRUE ~ blueprint_confirmed
	) ) |>

	mutate(azhimut_confirmed = case_when(
		cell_type_clean %in% c("cd4 memory t") & predicted.celltype.l2 == "cd4 tem" & blueprint_singler != "cd8 tcm" ~ T,
		cell_type_clean %in% c("cd4 memory t") & predicted.celltype.l2 == "cd4 tcm" & blueprint_singler != "cd8 tem" ~ T,
		TRUE ~ azhimut_confirmed
	) ) |>
	mutate(blueprint_confirmed = case_when(
		cell_type_clean %in% c("cd4 memory t") & predicted.celltype.l2 != "cd4 tem" & blueprint_singler == "cd4 tcm" ~ T,
		cell_type_clean %in% c("cd4 memory t") & predicted.celltype.l2 != "cd4 tcm" & blueprint_singler == "cd4 tem" ~ T,
		TRUE ~ blueprint_confirmed
	) ) |>

	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "t") & blueprint_singler =="cd8 t" & predicted.celltype.l2 |> str_detect("cd8"), T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "t") & blueprint_singler =="cd4 t" & predicted.celltype.l2 |> str_detect("cd4|treg"), T, azhimut_confirmed) ) |>

	mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "treg") & blueprint_singler %in% c("tregs"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "treg") & predicted.celltype.l2 == "treg", T, azhimut_confirmed) ) |>


	mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "tcm cd4") & blueprint_singler %in% c("cd4 tcm"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "tcm cd4") & predicted.celltype.l2 == "cd4 tcm", T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "tcm cd8") & blueprint_singler %in% c("cd8 tcm"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "tcm cd8") & predicted.celltype.l2 == "cd8 tcm", T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "tem cd4") & blueprint_singler %in% c("cd4 tem"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "tem cd4") & predicted.celltype.l2 == "cd4 tem", T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "tem cd8") & blueprint_singler %in% c("cd8 tem"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "tem cd8") & predicted.celltype.l2 == "cd8 tem", T, azhimut_confirmed) ) |>

	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "tgd") & predicted.celltype.l2 == "gdt", T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "activated cd4") & predicted.celltype.l2 == "cd4 proliferating", T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "activated cd8") & predicted.celltype.l2 == "cd8 proliferating", T, azhimut_confirmed) ) |>



	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c("naive cd4", "naive t") & predicted.celltype.l2 %in% c("cd4 naive"), T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c("naive cd8", "naive t") & predicted.celltype.l2 %in% c("cd8 naive"), T, azhimut_confirmed) ) |>

	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "prot") & predicted.celltype.l2 %in% c("cd4 naive") & !blueprint_singler |> str_detect("clp|hcs|mpp|cd8"), T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "prot") & predicted.celltype.l2 %in% c("cd8 naive") & !blueprint_singler |> str_detect("clp|hcs|mpp|cd4"), T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "prot") & predicted.celltype.l2 %in% c("hspc"), T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "prot") & blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), T, blueprint_confirmed) ) |>

	mutate(azhimut_confirmed = if_else(cell_type_clean == "dendritic" & predicted.celltype.l2 %in% c("asdc", "cdc2", "cdc1", "pdc"), T, azhimut_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean == "double negative t regulatory" & predicted.celltype.l2 == "dnt", T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "early t lineage precursor", "immature innate lymphoid") & blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "early t lineage precursor", "immature innate lymphoid") & predicted.celltype.l2 == "hspc" & blueprint_singler != "clp", T, azhimut_confirmed) ) |>

	mutate(blueprint_confirmed = if_else(cell_type_clean %in% c("ilc1", "ilc2", "innate lymphoid") & blueprint_singler == "nk", T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c("ilc1", "ilc2", "innate lymphoid") & predicted.celltype.l2 %in% c( "nk", "ilc", "nk proliferating"), T, azhimut_confirmed) ) |>

	mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "immature t") & blueprint_singler %in% c("naive t"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "immature t") & predicted.celltype.l2 == "t naive", T, azhimut_confirmed) ) |>

	mutate(cell_type_harmonised = if_else(cell_type_clean == "fraction a prepro b", "naive b", cell_type_harmonised))  |>
	mutate(blueprint_confirmed = if_else(cell_type_clean == "granulocyte" & blueprint_singler %in% c("eosinophils", "neutrophils"), T, blueprint_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean %in% c("immature neutrophil", "neutrophil") & blueprint_singler %in% c( "neutrophils"), T, blueprint_confirmed) ) |>

	mutate(blueprint_confirmed = if_else(cell_type_clean |> str_detect("megakaryocyte") & blueprint_singler |> str_detect("megakaryocyte"), T, blueprint_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean |> str_detect("macrophage") & blueprint_singler |> str_detect("macrophage"), T, blueprint_confirmed) ) |>

	mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "nk") & blueprint_singler %in% c("nk"), T, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "nk") & predicted.celltype.l2 %in% c("nk", "nk proliferating", "nk_cd56bright", "ilc"), T, azhimut_confirmed) ) |>


	# If identical force
	mutate(azhimut_confirmed = if_else(cell_type_clean == predicted.celltype.l2 , T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean == blueprint_singler , T, blueprint_confirmed) ) |>

	# Perogenitor
	mutate(azhimut_confirmed = if_else(cell_type_clean  |> str_detect("progenitor|hematopoietic|precursor") & predicted.celltype.l2  == "hspc", T, azhimut_confirmed) ) |>
	mutate(blueprint_confirmed = if_else(cell_type_clean  |> str_detect("progenitor|hematopoietic|precursor") & blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), T, blueprint_confirmed) ) |>

	# Generic original annotation and stem for new annotations
		mutate(azhimut_confirmed = if_else(
			cell_type_clean  %in% c("T cell", "myeloid cell", "leukocyte", "myeloid leukocyte", "B cell") &
				predicted.celltype.l2  == "hspc" &
				blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), T, azhimut_confirmed) ) |>

	# Omit mature for stem
	mutate(blueprint_confirmed = if_else(cell_type_clean  |> str_detect("mature") & blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), F, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean  |> str_detect("mature") & predicted.celltype.l2  == "hspc", F, azhimut_confirmed) ) |>

	# Omit megacariocyte for stem
	mutate(blueprint_confirmed = if_else(cell_type_clean  == "megakaryocyte" & blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), F, blueprint_confirmed) ) |>
	mutate(azhimut_confirmed = if_else(cell_type_clean  == "megakaryocyte" & predicted.celltype.l2  == "hspc", F, azhimut_confirmed) ) |>

	# Mast cells
	mutate(cell_type_harmonised = if_else(cell_type_clean == "mast", "mast", cell_type_harmonised))  |>


	# Visualise
	#distinct(cell_type_clean, predicted.celltype.l2, blueprint_singler, strong_evidence, azhimut_confirmed, blueprint_confirmed) |>
	arrange(!strong_evidence, cell_type_clean) |>

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
		is.na(cell_type_harmonised) & cell_type_clean  |> str_detect("progenitor|hematopoietic|stem|precursor") ~ "stem",

		is.na(cell_type_harmonised) & cell_type_clean == "cd14 monocyte" ~ "cd14 mono",
		is.na(cell_type_harmonised) & cell_type_clean == "cd16 monocyte" ~ "cd16 mono",
		is.na(cell_type_harmonised) & cell_type_clean %in% c("cd4 cytotoxic t", "tem cd4") ~ "cd4 tem",
		is.na(cell_type_harmonised) & cell_type_clean %in% c("cd8 cytotoxic t", "tem cd8") ~ "cd8 tem",
		is.na(cell_type_harmonised) & cell_type_clean |> str_detect("macrophage") ~ "macrophage",
		is.na(cell_type_harmonised) & cell_type_clean %in% c("mature b", "memory b", "transitional stage b") ~ "b memory",
		is.na(cell_type_harmonised) & cell_type_clean == "mucosal invariant t" ~ "mait",
		is.na(cell_type_harmonised) & cell_type_clean == "naive b" ~ "b naive",
		is.na(cell_type_harmonised) & cell_type_clean == "nk" ~ "nk",
		is.na(cell_type_harmonised) & cell_type_clean == "naive cd4" ~"cd4 naive",
		is.na(cell_type_harmonised) & cell_type_clean == "naive cd8" ~"cd8 naive",
		is.na(cell_type_harmonised) & cell_type_clean == "treg" ~ "treg",
		is.na(cell_type_harmonised) & cell_type_clean == "tgd" ~ "tgd",
		TRUE ~ cell_type_harmonised
	)) |>

	mutate(confidence_class = case_when(
		!is.na(cell_type_harmonised) & strong_evidence ~ 2,
		!is.na(cell_type_harmonised) & !strong_evidence ~ 3
	)) |>

	# Lowest grade annotation UNreliable
	mutate(cell_type_harmonised = case_when(

		# Get origincal annotation
		is.na(cell_type_harmonised) & cell_type_clean %in% c("neutrophil", "granulocyte") ~ cell_type_clean,
		is.na(cell_type_harmonised) & cell_type_clean %in% c("conventional dendritic", "dendritic") ~ "cdc",
		is.na(cell_type_harmonised) & cell_type_clean %in% c("classical monocyte") ~ "cd14 mono",

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
	))

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
  clean_cell_types_deeper() |> 
	bind_rows(
		annotation_crated_UNconfirmed
	) |>

	# I have multiple confidence_class per combination of labels
	distinct() |>
	with_groups(c(cell_type_clean, predicted.celltype.l2, blueprint_singler), ~ .x |> arrange(confidence_class) |> slice(1)) |>

	# Simplify after harmonisation
	mutate(cell_type_harmonised =	case_when(
		cell_type_harmonised %in% c("b memory", "b intermediate", "classswitched memory b", "memory b" ) ~ "b memory",
		cell_type_harmonised %in% c("b naive", "naive b") ~ "b naive",
		cell_type_harmonised %in% c("nk_cd56bright", "nk", "nk proliferating", "ilc") ~ "ilc",
		cell_type_harmonised %in% c("mpp", "clp", "hspc", "mep", "cmp", "hsc", "gmp") ~ "stem",
		cell_type_harmonised %in% c("macrophages",  "macrophages m1", "macrophages m2") ~ "macrophage",
		cell_type_harmonised %in% c("treg",  "tregs") ~ "treg",
		cell_type_harmonised %in% c("gdt",  "tgd") ~ "tgd",
		cell_type_harmonised %in% c("cd8 proliferating",  "cd8 tem") ~ "cd8 tem",
		cell_type_harmonised %in% c("cd4 proliferating",  "cd4 tem") ~ "cd4 tem",
		cell_type_harmonised %in% c("eosinophils",  "neutrophils", "granulocyte", "neutrophil") ~ "granulocyte",
		cell_type_harmonised %in% c("cdc",  "cdc1", "cdc2", "dc") ~ "cdc",

		TRUE ~ cell_type_harmonised
	)) |>
	dplyr::select(cell_type_clean, cell_type_harmonised, predicted.celltype.l2, blueprint_singler, confidence_class) |>
	distinct()



curated_annotation =
	annotation |>
  clean_cell_types_deeper() |> 
	filter(lineage_1=="immune") |>
	dplyr::select(
		.cell, .sample, cell_type, cell_type_clean, predicted.celltype.l2, blueprint_singler, monaco_singler) |>
	left_join(
		annotation_all ,
		by = c("cell_type_clean", "predicted.celltype.l2", "blueprint_singler")
	) |>
	dplyr::select(
		.cell, .sample, cell_type, cell_type_harmonised, confidence_class,
		cell_annotation_azimuth_l2 = predicted.celltype.l2, cell_annotation_blueprint_singler = blueprint_singler,
		cell_annotation_monaco_singler = monaco_singler
	) |>

# Reannotation of generic cell types
	mutate(cell_type_harmonised = case_when(
		cell_type_harmonised=="cd4 t" & cell_annotation_monaco_singler |> str_detect("effector memory") ~ "cd4 tem",
		cell_type_harmonised=="cd4 t" & cell_annotation_monaco_singler |> str_detect("mait") ~ "mait",
		cell_type_harmonised=="cd4 t" & cell_annotation_monaco_singler |> str_detect("central memory") ~ "cd4 tcm",
		cell_type_harmonised=="cd4 t" & cell_annotation_monaco_singler |> str_detect("naive") ~ "cd4 naive",
		cell_type_harmonised=="cd8 t" & cell_annotation_monaco_singler |> str_detect("effector memory") ~ "cd8 tem",
		cell_type_harmonised=="cd8 t" & cell_annotation_monaco_singler |> str_detect("central memory") ~ "cd8 tcm",
		cell_type_harmonised=="cd8 t" & cell_annotation_monaco_singler |> str_detect("naive") ~ "cd8 naive",
		cell_type_harmonised=="monocytes" & cell_annotation_monaco_singler |> str_detect("non classical") ~ "cd16 mono",
		cell_type == "nonclassical monocyte" & cell_type_harmonised=="monocytes" & cell_annotation_monaco_singler =="intermediate monocytes"   ~ "cd16 mono",
		cell_type_harmonised=="monocytes" & cell_annotation_monaco_singler |> str_detect("^classical") ~ "cd14 mono",
		cell_type == "classical monocyte" & cell_type_harmonised=="monocytes" & cell_annotation_monaco_singler =="intermediate monocytes"   ~ "cd14 mono",
		cell_type_harmonised=="monocytes" & cell_annotation_monaco_singler =="myeloid dendritic" & str_detect(cell_annotation_azimuth_l2, "cdc")   ~ "cdc",


		TRUE ~ cell_type_harmonised
	)) |>

	# Change CD4 classification for version 0.2.1
	mutate(confidence_class = if_else(
		cell_type_harmonised |> str_detect("cd4|mait|treg|tgd") & cell_annotation_monaco_singler %in% c("terminal effector cd4 t", "naive cd4 t", "th2", "th17", "t regulatory", "follicular helper t", "th1/th17", "th1", "nonvd2 gd t", "vd2 gd t"),
		3,
		confidence_class
	)) |>

	# Change CD4 classification for version 0.2.1
	mutate(cell_type_harmonised = if_else(
		cell_type_harmonised |> str_detect("cd4|mait|treg|tgd") & cell_annotation_monaco_singler %in% c("terminal effector cd4 t", "naive cd4 t", "th2", "th17", "t regulatory", "follicular helper t", "th1/th17", "th1", "nonvd2 gd t", "vd2 gd t"),
		cell_annotation_monaco_singler,
		cell_type_harmonised
	)) |>


	mutate(cell_type_harmonised = cell_type_harmonised |>
				 	str_replace("naive cd4 t", "cd4 naive") |>
				 	str_replace("th2", "cd4 th2") |>
				 	str_replace("^th17$", "cd4 th17") |>
				 	str_replace("t regulatory", "treg") |>
				 	str_replace("follicular helper t", "cd4 fh") |>
				 	str_replace("th1/th17", "cd4 th1/th17") |>
				 	str_replace("^th1$", "cd4 th1") |>
				 	str_replace("nonvd2 gd t", "tgd") |>
				 str_replace("vd2 gd t", "tgd")
	) |>

	# add immune_unclassified
	mutate(cell_type_harmonised = if_else(cell_type_harmonised == "monocytes", "immune_unclassified", cell_type_harmonised)) |>
	mutate(cell_type_harmonised = if_else(is.na(cell_type_harmonised), "immune_unclassified", cell_type_harmonised)) |>
	mutate(confidence_class = if_else(is.na(confidence_class), 5, confidence_class)) |>

	# drop uncommon cells
	mutate(cell_type_harmonised = if_else(cell_type_harmonised %in% c("cd4 t", "cd8 t", "asdc", "cd4 ctl"), "immune_unclassified", cell_type_harmonised))


# Further rescue of unannotated cells, manually

# curated_annotation |>
# 	filter(cell_type_harmonised == "immune_unclassified") |>
# 	count(cell_type   ,       cell_type_harmonised ,confidence_class ,cell_annotation_azimuth_l2 ,cell_annotation_blueprint_singler ,cell_annotation_monaco_singler) |>
# 	arrange(desc(n)) |>
# 	write_csv("curated_annotation_still_unannotated_0.2.csv")


curated_annotation =
	curated_annotation |>
	left_join(
		read_csv("~/PostDoc/CuratedAtlasQueryR/dev/curated_annotation_still_unannotated_0.2_manually_labelled.csv") |>
			select(cell_type, cell_type_harmonised_manually_curated = cell_type_harmonised, confidence_class_manually_curated = confidence_class, everything()),
		by = join_by(cell_type, cell_annotation_azimuth_l2, cell_annotation_blueprint_singler, cell_annotation_monaco_singler)
	) |>
	mutate(
		confidence_class = if_else(cell_type_harmonised == "immune_unclassified", confidence_class_manually_curated, confidence_class),
		cell_type_harmonised = if_else(cell_type_harmonised == "immune_unclassified", cell_type_harmonised_manually_curated, cell_type_harmonised),
	) |>
	select(-contains("manually_curated"), -n) |>

	# drop uncommon cells
	mutate(cell_type_harmonised = if_else(cell_type_harmonised %in% c("cd4 tcm", "cd4 tem"), "immune_unclassified", cell_type_harmonised))



	# # Recover confidence class == 4

	# curated_annotation |>
	# 	filter(confidence_class==4) |>
	# 	count(cell_type   ,       cell_type_harmonised ,confidence_class ,cell_annotation_azimuth_l2 ,cell_annotation_blueprint_singler ,cell_annotation_monaco_singler) |>
	# 	arrange(desc(n)) |>
	# 	write_csv("curated_annotation_still_unannotated_0.2_confidence_class_4.csv")

curated_annotation =
	curated_annotation |>
	left_join(
		read_csv("~/PostDoc/CuratedAtlasQueryR/dev/curated_annotation_still_unannotated_0.2_confidence_class_4_manually_labelled.csv") |>
			select(confidence_class_manually_curated = confidence_class, everything()),
		by = join_by(cell_type, cell_type_harmonised, cell_annotation_azimuth_l2, cell_annotation_blueprint_singler, cell_annotation_monaco_singler)
	) |>
	mutate(
		confidence_class = if_else(confidence_class == 4 & !is.na(confidence_class_manually_curated), confidence_class_manually_curated, confidence_class)
	) |>
	select(-contains("manually_curated"), -n)

# Correct fishy stem cell labelling
# If stem for the study's annotation and blueprint is non-immune it is probably wrong, 
# even because the heart has too many progenitor/stem
curated_annotation =
  curated_annotation |>
  mutate(confidence_class = case_when(
    cell_type_harmonised == "stem" & cell_annotation_blueprint_singler %in% c(
      "skeletal muscle", "adipocytes", "epithelial", "smooth muscle", "chondrocytes", "endothelial"
    ) ~ 5,
    TRUE ~ confidence_class
  ))


curated_annotation_merged =

	# Fix cell ID
	metadata_df |>
	dplyr::select(.cell, .sample, cell_type) |>
	as_tibble() |>

	# Add cell type
	left_join(curated_annotation |> dplyr::select(-cell_type), by = c(".cell", ".sample")) |>

  # Add non immune
	mutate(cell_type_harmonised = if_else(is.na(cell_type_harmonised), "non_immune", cell_type_harmonised)) |>
	mutate(confidence_class = if_else(is.na(confidence_class) & cell_type_harmonised == "non_immune", 1, confidence_class)) |>

	# For some unknown reason
	distinct()


curated_annotation_merged |>

	# Save
	saveRDS(file_curated_annotation_merged)

metadata_annotated =
	curated_annotation_merged |>

	# merge with the rest of metadata
	left_join(
		metadata_df |>
			as_tibble(),
		by=c(".cell", ".sample", "cell_type")
	)

# Replace `.` with `_` for all column names as it can create difficoulties for MySQL and Python
colnames(metadata_annotated) = colnames(metadata_annotated) |> str_replace_all("\\.", "_")
metadata_annotated = metadata_annotated |> rename(cell_ = `_cell`, sample_ = `_sample`)

harmonise_names_non_immune = function(metadata){
  
  # grep("\\bfibroblast\\b", c("myofibroblast", "fibroblast", "other fibroblast"))
  
  metadata$cell_type_harmonised =  metadata$cell_type
  
  table(metadata$cell_type_harmonised[grepl("\\bepithelial\\b", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("\\bepithelial\\b", metadata$cell_type_harmonised),
                                          "epithelial_cell",  ## does not include myoepithelial.
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("\\bfibroblast\\b", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("\\bfibroblast\\b", metadata$cell_type_harmonised),
                                          "fibroblast",  ## does not include myofibroblast.
                                          metadata$cell_type_harmonised)
  
  
  table(metadata$cell_type_harmonised[grepl("\\bendothelial\\b", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("\\bendothelial\\b", metadata$cell_type_harmonised),
                                          "endothelial_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("Mueller", metadata$cell_type_harmonised)])
  table(metadata$cell_type_harmonised[grepl("Muller", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(metadata$cell_type_harmonised=="Mueller cell" | metadata$cell_type_harmonised=="Muller cell",
                                          "Muller_cell",  ## Merge Mueller cells and Muller cells.
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("\\bneuron\\b", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("\\bneuron\\b", metadata$cell_type_harmonised),
                                          "neuron",  ## does not include interneuron and neuronal receptor cell.
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("amplifying cell", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("amplifying cell", metadata$cell_type_harmonised),
                                          "amplifying_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("stem cell", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("stem cell", metadata$cell_type_harmonised),
                                          "stem_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("progenitor cell", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("progenitor cell", metadata$cell_type_harmonised),
                                          "progenitor_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("acinar cell", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("acinar cell", metadata$cell_type_harmonised),
                                          "acinar_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("goblet cell", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("goblet cell", metadata$cell_type_harmonised),
                                          "goblet_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("thymocyte", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("thymocyte", metadata$cell_type_harmonised),
                                          "thymocyte",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("urothelial", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("urothelial", metadata$cell_type_harmonised),
                                          "urothelial_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("fat", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("fat", metadata$cell_type_harmonised),
                                          "fat_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("pneumocyte", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("pneumocyte", metadata$cell_type_harmonised),
                                          "pneumocyte",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("pneumocyte", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("pneumocyte", metadata$cell_type_harmonised),
                                          "pneumocyte",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("mesothelial", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("mesothelial", metadata$cell_type_harmonised),
                                          "mesothelial_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("enteroendocrine", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("enteroendocrine", metadata$cell_type_harmonised),
                                          "enteroendocrine_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("enterocyte", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("enterocyte", metadata$cell_type_harmonised),
                                          "enterocyte",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("basal", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("basal", metadata$cell_type_harmonised),
                                          "basal_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("stromal", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("stromal", metadata$cell_type_harmonised),
                                          "stromal_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("retina", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("retina", metadata$cell_type_harmonised),
                                          "retinal_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("ciliated", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("ciliated", metadata$cell_type_harmonised),
                                          "ciliated_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("pericyte", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("pericyte", metadata$cell_type_harmonised),
                                          "pericyte_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("trophoblast", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("trophoblast", metadata$cell_type_harmonised),
                                          "trophoblast",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("brush", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("brush", metadata$cell_type_harmonised),
                                          "brush_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("serous", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("serous", metadata$cell_type_harmonised),
                                          "serous cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("hepatocyte", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("hepatocyte", metadata$cell_type_harmonised),
                                          "hepatocyte",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("melanocyte", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("melanocyte", metadata$cell_type_harmonised),
                                          "melanocyte",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("myocyte", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("myocyte", metadata$cell_type_harmonised),
                                          "myocyte",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("promyelocyte", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("promyelocyte", metadata$cell_type_harmonised),
                                          "promyelocyte",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("cholangiocyte", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("cholangiocyte", metadata$cell_type_harmonised),
                                          "cholangiocyte",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("melanocyte", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("melanocyte", metadata$cell_type_harmonised),
                                          "melanocyte",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("muscle", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("myoblast", metadata$cell_type_harmonised),
                                          "myoblast", ## Discussed with Stefano on Teams on 16/12/2022.
                                          metadata$cell_type_harmonised)
  metadata$cell_type_harmonised <- ifelse(grepl("satellite", metadata$cell_type_harmonised),
                                          "satellite_cell", ## Discussed with Stefano on Teams on 16/12/2022.
                                          metadata$cell_type_harmonised)
  metadata$cell_type_harmonised <- ifelse(grepl("muscle", metadata$cell_type_harmonised),
                                          "muscle_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("progenitor", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("progenitor", metadata$cell_type_harmonised),
                                          "progenitor_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("progenitor", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("progenitor", metadata$cell_type_harmonised),
                                          "progenitor_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("erythrocyte", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("erythrocyte", metadata$cell_type_harmonised),
                                          "erythrocyte",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("myoepithelial", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("myoepithelial", metadata$cell_type_harmonised),
                                          "myoepithelial_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("myofibroblast", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("myofibroblast", metadata$cell_type_harmonised),
                                          "myofibroblast_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("pancreatic", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("pancreatic", metadata$cell_type_harmonised),
                                          "pancreatic_cell",  ## Discussed with Stefano on Teams on 19/12/2022.
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("renal", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("renal", metadata$cell_type_harmonised),
                                          "renal_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("epidermal", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("epidermal", metadata$cell_type_harmonised),  ## Includes epidermal cell and epidermal Langerhans cell.
                                          "epidermal_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("cortical", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("cortical", metadata$cell_type_harmonised),
                                          "cortical_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("interstitial", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("interstitial", metadata$cell_type_harmonised),
                                          "interstitial_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("neuroendocrine", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("neuroendocrine", metadata$cell_type_harmonised),
                                          "neuroendocrine_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("granular", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("granular", metadata$cell_type_harmonised),
                                          "granular_cell", ## Discussed with Stefano on Teams on 19/12/2022.
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("kidney", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("kidney", metadata$cell_type_harmonised),
                                          "kidney_cell",  ## Discussed with Stefano on Teams on 19/12/2022.
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("paneth", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("paneth", metadata$cell_type_harmonised),
                                          "paneth_cell",
                                          metadata$cell_type_harmonised)
  
  table(metadata$cell_type_harmonised[grepl("bipolar", metadata$cell_type_harmonised)])
  metadata$cell_type_harmonised <- ifelse(grepl("bipolar", metadata$cell_type_harmonised),
                                          "bipolar_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- gsub(" " , "_", metadata$cell_type_harmonised)
  
  
  
  
  table(metadata$cell_type_harmonised[grepl("glial", metadata$cell_type_harmonised)])
  ##  glial cell, microglial cell, radial glial cell
  ## https://www.simplypsychology.org/glial-cells.html#:~:text=Glial%20cells%20are%20a%20general,that%20keep%20the%20brain%20functioning.
  
  table(metadata$cell_type_harmonised[grepl("hematopoietic", metadata$cell_type_harmonised)])
  ## hematopoietic cell, hematopoietic precursor cell
  ## precursor cell = stem_cell?
  
  table(metadata$cell_type_harmonised[grepl("papillary", metadata$cell_type_harmonised)])
  ## papillary tips cell = renal_cell ?
  
  
  metadata
}

dictionary_connie_non_immune = 
  metadata_annotated |> 
  filter(cell_type_harmonised == "non_immune") |> 
  distinct(cell_type) |> 
  harmonise_names_non_immune() |> 
  rename(cell_type_harmonised_non_immune = cell_type_harmonised )

metadata_annotated = 
  metadata_annotated |> 
  left_join(dictionary_connie_non_immune) |> 
  mutate(cell_type_harmonised = if_else(cell_type_harmonised=="non_immune", cell_type_harmonised_non_immune, cell_type_harmonised)) |> 
  select(-cell_type_harmonised_non_immune)

# Save
metadata_annotated |>
	saveRDS(file_metadata_annotated)

# # Build SQLite
# library(RSQLite)
# library(DBI)
# library(dplyr)
# 
# sqllite_output_path = tools::file_path_sans_ext(file_metadata_annotated)
# con <- dbConnect(SQLite(), dbname=glue("{sqllite_output_path}.sqlite"))
# dbWriteTable(con, "metadata", metadata_annotated)
# dbDisconnect(con)

# Build Parquet
arrow::write_parquet(metadata_annotated, glue("{tools::file_path_sans_ext(file_metadata_annotated)}.parquet"))

# # Filter samples that do not have immune
# cell_metadata_with_harmonised_annotation =
#   cell_metadata_with_harmonised_annotation |>
#   nest(data = -c(.sample, tissue_harmonised)) |>
#   filter(map_int(data, ~ .x |> filter(cell_type_harmonised != "non_immune") |> nrow()) > 0) |>
#   unnest(data)
#
# # Tissue with no immune
# cell_metadata_with_harmonised_annotation =
#   cell_metadata_with_harmonised_annotation |>
#   filter(disease == "normal") |>
#
#   # Filter tissues
#   nest(data = -c(cell_type_harmonised, tissue_harmonised)) |>
#   add_count(tissue_harmonised, name = "n_cell_type_in_tissue") |>
#   filter(n_cell_type_in_tissue>=18) |>
#   add_count(cell_type_harmonised, name = "n_tissue_in_cell_type") |>
#   filter(n_tissue_in_cell_type>=26) |>
#   unnest(data)



#
# annotated_samples =
# 	cell_metadata_with_harmonised_annotation |> filter(!is.na(cell_type_harmonised)) |>  distinct( cell_type, .sample, file_id)
#
# # Cell types that most need attention
# cell_metadata_with_harmonised_annotation |> anti_join(annotated_samples) |> select(contains("cell_")) |> count(cell_type ,    cell_type_harmonised ,cell_annotation_azimuth_l2 ,cell_annotation_blueprint_singler) |> arrange(desc(n)) |> print(n=99)
#
# # How many samples miss annotation
# cell_metadata_with_harmonised_annotation |> anti_join(annotated_samples)  |>  distinct(cell_type,  cell_type_harmonised, .sample) |>  distinct( cell_type, .sample) |> count(cell_type)  |> arrange(desc(n))
#
# # Histo of annotation
# cell_metadata_with_harmonised_annotation |> filter(!is.na(cell_type_harmonised)) |>  distinct( cell_type_harmonised, .sample) |> count(.sample) |> pull(n) |> hist(breaks=30)

# NEEDED CHANGES - ANIMAL CELLS AND HEMATOPOIETIC SHOULD BE EVALUATED AS IMMUNE FOR REANNNOTATION?

# Tissue with no immune

# temp_tissue =
#   cell_metadata_with_harmonised_annotation |>
#   filter(cell_type_harmonised=="pdc") |>
#   pull(tissue_harmonised)
#
# cell_metadata_with_harmonised_annotation |>
#   filter(!tissue_harmonised %in% temp_tissue) |>
#   distinct(tissue_harmonised) |>
#
#
# cell_metadata_with_harmonised_annotation |>
#   filter(disease == "normal") |>
#   filter(!is.na(cell_type_harmonised)) |>
#   distinct( cell_type_harmonised, tissue_harmonised) |>
#   count(tissue_harmonised) |> arrange(n) |> print(n=99)


