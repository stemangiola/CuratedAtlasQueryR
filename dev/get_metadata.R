library(zellkonverter)
library(SingleCellExperiment) # load early to avoid masking dplyr::count()
library(dplyr)
library(cellxgenedp)
library(tidyverse)
#library(tidySingleCellExperiment)
library(stringr)
library(scMerge)
library(glue)
library(DelayedArray)
library(HDF5Array)
library(openssl)

# # CREATE MAKEFILE
# tab = "\t"
# root_directory = "/vast/projects/cellxgene_curated"
# metadata_directory = glue("{root_directory}/metadata_0.2")
# raw_data_directory = glue("{root_directory}/raw_data")
# files_metadata = glue("{root_directory}/files_metadata.rds")
#
# input_files_path = dir(raw_data_directory, full.names=TRUE)
# input_files = input_files_path |> basename()
#
# files_thta_require_more_memory = c(
# 	"5c64f247-5b7c-4842-b290-65c722a65952",
# 	"56e0359f-ee8d-4ba5-a51d-159a183643e5",
# 	"51f114ae-232a-4550-a910-934e175db814",
# 	"21ca95b3-776b-4fa2-9956-09a07c0e5224"
# )
#
#
# output_files = input_files |> str_replace("H5AD$", "rds")
#
#
# output_files_path = glue("{metadata_directory}/{output_files}")
# metadata_path = glue("{root_directory}/metadata_0.2.rds")
#
# # output_files_path = output_files_path |> str_replace(root_directory, my_root_directory)
#
# heavy_file_pattern = "5c64f247-5b7c-4842-b290-65c722a65952|56e0359f-ee8d-4ba5-a51d-159a183643e5|51f114ae-232a-4550-a910-934e175db814|21ca95b3-776b-4fa2-9956-09a07c0e5224"
# input_files_path_heavy = input_files_path |> str_subset(heavy_file_pattern)
# input_files_path_light = input_files_path |> str_subset(heavy_file_pattern, negate = TRUE)
# output_files_path_heavy = output_files_path |> str_subset(heavy_file_pattern)
# output_files_path_light = output_files_path |> str_subset(heavy_file_pattern, negate = TRUE)
#
# c(
# 	glue("CATEGORY=get_metadata_light\nMEMORY=80024\nCORES=1\nWALL_TIME=10000"),
# 	glue("{output_files_path_light}:{input_files_path_light}\n{tab}Rscript get_metadata.R {input_files_path_light} {output_files_path_light}"),
#
# 	glue("CATEGORY=get_metadata_heavy\nMEMORY=200024\nCORES=1\nWALL_TIME=10000"),
# 	glue("{output_files_path_heavy}:{input_files_path_heavy}\n{tab}Rscript get_metadata.R {input_files_path_heavy} {output_files_path_heavy}"),
#
# 	glue("CATEGORY=merge_metadata\nMEMORY=80024\nCORES=1\nWALL_TIME=10000"),
# 	glue("{metadata_path}:{paste(output_files_path, collapse = \" \")} {files_metadata}\n{tab}Rscript merge_metadata.R {paste(output_files_path, collapse = \" \")} {files_metadata} {metadata_path}")
# )  |>
# 	write_lines(glue("~/PostDoc/CuratedAtlasQueryR/dev/get_metadata.makeflow"))

source("utility.R")

root_directory = "/vast/projects/RCP/human_cell_atlas"
splitted_light_data_directory = glue("{root_directory}/splitted_light_data")

# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file = args[[1]]
input_file = "/home/users/allstaff/mangiola.s/.cache/R/cellxgenedp/curation/v1/c480e527-5725-4699-bd8a-e09535b23ba8.h5ad"
output_file = args[[2]]

# Create directory
output_file |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

# Read metadata
data = readH5AD(input_file,	use_hdf5 = TRUE	)

col_data = data |> colData()

rm(data)
gc()

col_data |>
 	as_tibble(rownames = ".cell") |>

	# Link file IDs
	mutate(file_id = basename(input_file) |> tools::file_path_sans_ext()) %>%

	# Clean


 	# Sort sample ID
 	# Fix some sample id missing
 	when(unique(.$file_id) %in% c(
 		"11b86bc3-6d4d-4e28-903a-0361ea8f6bdf",
 		"492b0613-ff5b-4fca-a585-503fc4102e4f",
 		"11b86bc3-6d4d-4e28-903a-0361ea8f6bdf",
 		"0e8f9ce4-46e5-434e-9ca0-e769d1dd27ea"
 	) ~ mutate(., PatientID = glue("{sample} {replicate} {time_point} {target}") |> as.character()) , ~ (.)) %>%
 	when(unique(.$file_id) %in% c(
 		"0273924c-0387-4f44-98c5-2292dbaab11e",
 		"a16bec18-5c9f-40ad-8169-12c5199c7506",
 		"556bb449-bbef-43d3-9487-87031fc0decb"
 	) ~ mutate(., PatientID = glue("{Collection.ID} {Genotype} {Location}")|> as.character()) , ~ (.)) %>%
 	when(unique(.$file_id) %in% c(
 		"b83afdc1-baa1-42c0-bd5b-cb607084757d"
 	) ~ mutate(., PatientID = glue("{sex} {development_stage} {disease}")|> as.character()) , ~ (.)) %>%
 	when(unique(.$file_id) %in% c(
 		"3fe53a40-38ff-4f25-b33b-e4d60f2289ef",
 		"5c1cc788-2645-45fb-b1d9-2f43d368bba8"
 	) ~ mutate(., PatientID = glue("{Batch} {Fetus_id} {Development_day} {sex} {tissue} {disease}")|> as.character()) , ~ (.)) |>



 	# select(
 	# 	.cell,
 	# 	# Select most common columns across datasets
 	# 	one_of(common_colnames),
 	#
 	# 	# select columns for sample ID
 	# 	contains("sample", ignore.case = T),
 	# 	contains("donor", ignore.case = T),
 # 	contains("uuid", ignore.case = T),
 # 	contains("patient", ignore.case = T),
 # 	contains("specimen", ignore.case = T),
 # 	contains("library", ignore.case = T),
 # 	one_of("index"),
 # 	one_of("barcode"),
 # 	one_of("catalog_number"),
 # 	one_of("dataset_id"),
 # 	one_of("batch") ,
 # 	one_of("development_stage"),
 # 	one_of("Fetus_id"),
 # 	one_of("individual"),
 # 	one_of("orig.ident"),
 # 	one_of("Name"),
 # 	one_of("tissue"),
 # 	file_id
 # ) |>
 mutate_if(is.factor, as.character) |>
 	type_convert(guess_integer = TRUE) |>
 	mutate_if(is.integer, as.character) |>

 	# Convert types
 	when("donor_id" %in% colnames(.) ~ mutate(., donor_id = donor_id |> as.character() ), ~(.)) |>
 	when("Cluster" %in% colnames(.) ~ mutate(., Cluster = Cluster |> as.character() ), ~(.)) |>
 	when("cluster_id" %in% colnames(.) ~ mutate(., cluster_id = cluster_id |> as.character() ), ~(.)) |>
 	when("Batch" %in% colnames(.) ~ mutate(., Batch = Batch |> as.character() ), ~(.)) |>
 	when("batch" %in% colnames(.) ~ mutate(., batch = batch |> as.character() ), ~(.)) |>
 	when("age" %in% colnames(.) ~ mutate(., age = age |> as.character() ), ~(.)) |>
 	when("BMI" %in% colnames(.) ~ mutate(., BMI = BMI |> as.character() ), ~(.)) |>
 	when("donor_BMI" %in% colnames(.) ~ mutate(., donor_BMI = donor_BMI |> as.character() ), ~(.)) |>
 	when("author_cell_type" %in% colnames(.) ~ mutate(., author_cell_type = author_cell_type |> as.character() ), ~(.)) |>
 	when("time_point" %in% colnames(.) ~ mutate(., time_point = time_point |> as.character() ), ~(.)) |>
 	when("cluster" %in% colnames(.) ~ mutate(., cluster = cluster |> as.character() ), ~(.)) |>
 	when("ClusterID" %in% colnames(.) ~ mutate(., ClusterID = ClusterID |> as.character() ), ~(.)) |>
 	when("Stage" %in% colnames(.) ~ mutate(., Stage = Stage |> as.character() ), ~(.)) |>
 	when("individual" %in% colnames(.) ~ mutate(., individual = individual |> as.character() ), ~(.)) |>
 	when("recurrent_cluster" %in% colnames(.) ~ mutate(., recurrent_cluster = recurrent_cluster |> as.character() ), ~(.)) |>
 	when("PatientID" %in% colnames(.) ~ mutate(., PatientID = PatientID |> as.character() ), ~(.)) |>
 	when("PMI" %in% colnames(.) ~ mutate(., PMI = PMI |> as.character() ), ~(.)) |>
 	when("n_genes" %in% colnames(.) ~ mutate(., n_genes = n_genes |> as.numeric() ), ~(.)) |>
 	when("n_counts" %in% colnames(.) ~ mutate(., n_counts = n_counts |> as.numeric() ), ~(.)) |>
 	when("n_genes_by_counts" %in% colnames(.) ~ mutate(., n_genes_by_counts = n_genes_by_counts |> as.numeric() ), ~(.)) |>
 	when("nUMI" %in% colnames(.) ~ mutate(., nUMI = nUMI |> as.numeric() ), ~(.)) |>
 	when("percent.cortex" %in% colnames(.) ~ mutate(., percent.cortex = percent.cortex |> as.character() ), ~(.)) |>
 	when("percent.medulla" %in% colnames(.) ~ mutate(., percent.medulla = percent.medulla |> as.character() ), ~(.)) |>
 	when("Age" %in% colnames(.) ~ mutate(., Age = Age |> as.numeric() ), ~(.)) |>
 	when("nCount_RNA" %in% colnames(.) ~ mutate(., nCount_RNA = nCount_RNA |> as.numeric() ), ~(.)) |>
 	when("is_primary_data" %in% colnames(.) ~ mutate(., is_primary_data = is_primary_data |> as.character() ), ~(.)) |>

 	#mutate(across(contains("cluster", ignore.case = TRUE), ~ as.character)) |>
 	select(-one_of('PCW')) %>%

 	# Sort sample ID. It works but not elegant.
 	# Based on observation of strangely behaving datasets, where sample ID is not clear
 	when("sampleID" %in% colnames(.) & !"PatientID" %in% colnames(.) ~
 			 	mutate(., PatientID = as.character(sampleID )) |>  select(-sampleID), ~(.)) %>%
 	when("Patient" %in% colnames(.) ~ mutate(., Sample = NA |> as.character()), ~(.)) |>
 	mutate(sample_placeholder = NA |> as.character()) %>%
 	when(unique(.$file_id)=="e40591e7-0e5a-4bef-9b60-7015abe5b17f" ~ mutate(., sample_placeholder = glue("{batch} {development_stage}") |> as.character()), ~ (.)) %>%
 	when(unique(.$file_id)=="39b6cc45-8c5c-4f7b-944c-58f66da5efb1" ~ mutate(., sample_placeholder =sample_id), ~ (.))  %>%
 	when(unique(.$file_id)=="443d6a0e-dbcb-4002-8af0-628e7d4a18fa" ~ mutate(., sample_placeholder =sample_id), ~ (.))  %>%
 	when(unique(.$file_id)=="a91f075b-52d5-4aa3-8ecc-86c4763a49b3" ~ mutate(., sample_placeholder =sample), ~ (.))  %>%
 	when(unique(.$file_id)=="0af763e1-0e2f-4de6-9563-5abb0ad2b01e" ~ mutate(., sample_placeholder ="only_one_culture"), ~ (.))  %>%
 	when(unique(.$file_id)=="5c64f247-5b7c-4842-b290-65c722a65952" ~ mutate(., sample_placeholder ="only_one_culture"), ~ (.))  %>%
 	when(unique(.$file_id)=="d6f92754-e178-4202-b86f-0f430e965d72" ~ mutate(., sample_placeholder =orig.ident), ~ (.))  %>%
 	when(unique(.$file_id)=="c790ef7a-1523-4627-8603-d6a02f8f4877" ~ mutate(., sample_placeholder =orig.ident), ~ (.))  %>%
 	when(unique(.$file_id)=="1e81a742-e457-4fc6-9c39-c55189ec9dc2" ~ mutate(., sample_placeholder =orig.ident), ~ (.))  %>%
 	when(unique(.$file_id)=="351ef284-b59e-43a5-83ba-0eb907dc282c" ~ mutate(., sample_placeholder =orig.ident), ~ (.))  %>%
	when(unique(.$file_id)=="f498030e-246c-4376-87e3-90b28c7efb00" ~ mutate(., sample_placeholder =Name), ~ (.))  %>%

	# These are the datasets with too few cells per inferred samples, therefore simplifying
	when(unique(.$file_id)=="e3a56e00-8417-4d82-9d35-3fab3aac12f2" ~ mutate(., SpecimenID =NA), ~ (.))  %>%
	when(unique(.$file_id)=="17b34e42-bbd2-494b-bf32-b9229344a3f6" ~ mutate(., Sample =NA), ~ (.))  %>%

	# Fix huge samples for plate experiments
	extract(.cell, "experiment___", "(^expr?[0-9]+)", remove = F) |>
	mutate(experiment___ = if_else(file_id=="3fe53a40-38ff-4f25-b33b-e4d60f2289ef", experiment___, "")) |>

 	# Empirically infer samples from many characteristics
 	unite(".sample_name", one_of(
 		"sample_placeholder",
 		"Sample",
 		"SampleID",
 		"sample_uuid",
 		"Sample_ID",
 		"scRNASeq_sample_ID",
 		"Sample_Tag",
 		"Sample.ID",
 		"sample_names",
 		"Short_Sample",
 		"Sample.ID.short",
 		"Sample.name",
 		"patient",
 		"Donor.ID",
 		"donor_id",
 		"donor",
 		"PatientID",
 		"donor_uuid",
 		"library_uuid",
 		"suspension_uuid",
 		"Patient",
 		"tissue_section_uuid",
 		"DonorID",
 		"specimen",
 		"SpecimenID",
 		"Fetus_id",
 		"individual",
 		"tissue",
 		"development_stage",
 		"assay",
 		"experiment___",
 		"disease"
 	), na.rm = TRUE, sep = "___", remove = F) |>


 	# Add sample hash
 	mutate(.sample = md5(glue("{.sample_name}{file_id}"))) |>

 	# # Make cell unique
 	# mutate(.cell = glue("{.cell}_{.sample}")) |>

 	# make lighter
 	mutate_if(is.character, as.factor) |>

	saveRDS(output_file)




# # TROUBLESHOOT FOR DECIDING SAMPLE ID - DO NOT DELETE
# x=metadata |>
# 	map_dfr(
# 		~ .x |>
#
# 			# Fix some sample id missing
# 			when(unique(.$file_id) %in% c(
# 				"11b86bc3-6d4d-4e28-903a-0361ea8f6bdf",
# 				"492b0613-ff5b-4fca-a585-503fc4102e4f",
# 				"11b86bc3-6d4d-4e28-903a-0361ea8f6bdf",
# 				"0e8f9ce4-46e5-434e-9ca0-e769d1dd27ea"
# 			) ~ .x |> mutate(PatientID = glue("{sample} {replicate} {time_point} {target}") |> as.character()) , ~ (.)) |>
# 			when(unique(.$file_id) %in% c(
# 				"0273924c-0387-4f44-98c5-2292dbaab11e",
# 				"a16bec18-5c9f-40ad-8169-12c5199c7506",
# 				"556bb449-bbef-43d3-9487-87031fc0decb"
# 			) ~ .x |> mutate(PatientID = glue("{Collection.ID} {Genotype} {Location}")|> as.character()) , ~ (.)) |>
# 			when(unique(.$file_id) %in% c(
# 				"b83afdc1-baa1-42c0-bd5b-cb607084757d"
# 			) ~ .x |> mutate(PatientID = glue("{sex} {development_stage} {disease}")|> as.character()) , ~ (.)) |>
# 			when(unique(.$file_id) %in% c(
# 				"3fe53a40-38ff-4f25-b33b-e4d60f2289ef",
# 				"5c1cc788-2645-45fb-b1d9-2f43d368bba8"
# 			) ~ .x |> mutate(PatientID = glue("{Batch} {Fetus_id} {Development_day} {sex} {tissue} {disease}")|> as.character()) , ~ (.)) |>
#
# 			# Subset
# 			select(
# 				contains("sample", ignore.case = T),
# 				contains("donor", ignore.case = T),
# 				contains("uuid", ignore.case = T),
# 				contains("patient", ignore.case = T),
# 				contains("specimen", ignore.case = T),
# 				contains("library", ignore.case = T),
# 				one_of("index"),
# 				one_of("barcode"),
# 				one_of("catalog_number"),
# 				one_of("dataset_id"),
# 				one_of("batch") ,
# 				one_of("development_stage"),
# 				one_of("Fetus_id"),
# 				one_of("individual"),
# 				one_of("orig.ident"),
# 				one_of("Name"),
# 				one_of("tissue"),
# 				file_id
#
# 			) |>
# 			when("batch" %in% colnames(.) ~ mutate(., batch = batch |> as.character() ), ~(.))
# 	) |>
#
# 	mutate_if(is.factor, as.character)
#
# y=x |>
# 	mutate(PatientID = if_else(is.na(PatientID), sampleID, PatientID)) |>
# 	select(-sampleID) |>
# 	mutate(Sample = case_when(is.na(Patient) ~ Sample)) |>
# 	replace_na(list(PatientID = "")) |>
# 	mutate(PatientID = case_when(
# 		file_id=="e40591e7-0e5a-4bef-9b60-7015abe5b17f" ~ glue("{batch} {development_stage}") |> as.character(),
# 		file_id=="39b6cc45-8c5c-4f7b-944c-58f66da5efb1" ~ sample_id,
# 		file_id=="443d6a0e-dbcb-4002-8af0-628e7d4a18fa" ~ sample_id,
# 		file_id=="a91f075b-52d5-4aa3-8ecc-86c4763a49b3" ~ sample,
# 		file_id=="0af763e1-0e2f-4de6-9563-5abb0ad2b01e" ~ "only_one_culture",
# 		file_id=="5c64f247-5b7c-4842-b290-65c722a65952" ~ "only_one_culture",
# 		file_id=="d6f92754-e178-4202-b86f-0f430e965d72" ~ orig.ident,
# 		file_id=="c790ef7a-1523-4627-8603-d6a02f8f4877" ~ orig.ident,
# 		file_id=="1e81a742-e457-4fc6-9c39-c55189ec9dc2" ~ orig.ident,
# 		file_id=="351ef284-b59e-43a5-83ba-0eb907dc282c" ~ orig.ident,
# 		file_id=="f498030e-246c-4376-87e3-90b28c7efb00" ~ Name,
# 		T ~ PatientID
# 	)) |>
#
# 	unite(".sample_name", c(
# 		Sample, SampleID, sample_uuid, Sample_ID,
# 		scRNASeq_sample_ID, Sample_Tag, Sample.ID, sample_names,
# 		Short_Sample, Sample.ID.short, Sample.name,
# 		patient, Donor.ID, donor_id, donor, PatientID, donor_uuid,
# 		library_uuid, suspension_uuid, Patient, tissue_section_uuid, DonorID, specimen, SpecimenID, Fetus_id, individual, tissue
# 	), na.rm = TRUE, sep = "___", remove = F)
#
# y |>
# 	dplyr::count(file_id, .sample_name ) |>
# 	arrange(desc(n)) |>
# 	View()

