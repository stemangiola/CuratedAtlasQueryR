suppressPackageStartupMessages({
	library(zellkonverter)
	library(SingleCellExperiment) # load early to avoid masking dplyr::count()
	library(dplyr)
	library(cellxgenedp)
})

library(tidyverse)
#library(tidySingleCellExperiment)
library(stringr)
library(scMerge)
library(glue)
library(DelayedArray)
library(HDF5Array)

# # CREATE MAKEFILE
# metadata_directory = glue("{root_directory}/metadata")
# raw_data_directory = glue("{root_directory}/raw_data")
# input_files = dir(raw_data_directory, full.names=TRUE) |> basename()
# 
# c(
# 	glue("CATEGORY=metadata\nMEMORY=30024\nCORES=1\nWALL_TIME=10000"),
# 	glue("{metadata_directory}/{input_files}:{raw_data_directory}/{input_files}\n{tab}Rscript metadata.R {raw_data_directory}/{input_files} {metadata_directory}/{input_files}")
# )  |> 
# 	write_lines(glue("metadata.makeflow"))


# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file = args[[1]]
output_file = args[[2]]

# local_file = readRDS("/vast/scratch/users/mangiola.s/human_cell_atlas/file_location.rds")
local_file = dir(glue("{root_directory}/raw_data"), full.names = TRUE)

# Read metadata
metadata = 
	local_file %>% 
	map(	~ readH5AD(
		.x, 
		#reader = "R", 
		use_hdf5 = TRUE
	) |> 
		colData() |> 
			 	as_tibble() |>
			 	mutate(file_id = basename(.x) |> tools::file_path_sans_ext()) |> 
			 	mutate(file_path = .x)
	)

# Get which column don't have too many NAs
common_colnames = 
	metadata |> 
	map_dfr(
		~ colnames(.x) |> as_tibble() |> mutate(dataset_id = unique(.x$dataset_id))
	) |>  
	count(value) |> 
	mutate(n_datasets = length(local_file)) |> 
	filter(n > (n_datasets / 2)) |>
	pull(value)

# Get all metadata
metadata  |> 
	
	# Clean
	map( ~ 
			 	.x |> 
			 	select(one_of(common_colnames)) |> 
			 	mutate_if(is.factor, as.character) |> 
			 	type_convert(guess_integer = TRUE) |> 
			 	mutate_if(is.integer, as.character) |> 
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
			 	#mutate(across(contains("cluster", ignore.case = TRUE), ~ as.character)) |> 
			 	select(-one_of('PCW'))
			 
			 
			) |> 
	
	# Integrate
	bind_rows()  |> 
	mutate_if(is.character, as.factor) |> 
	
	# Add files metadata
	left_join(readRDS(glue("{root_directory}/files_metadata.rds"))) |> 
	
	saveRDS(output_path_metadata)
	

# 
# metadata = 
# 	metadata |> 
# 	left_join(read_csv("metadata_cell_type.csv"))






x=metadata_ |> 
	map_dfr(
		~ .x |>
			select(
				contains("sample", ignore.case = T), 
				contains("donor", ignore.case = T), 
				contains("uuid", ignore.case = T), 
				contains("patient", ignore.case = T),
				contains("specimen", ignore.case = T), 
				contains("library", ignore.case = T), 
				one_of("index"),
				one_of("barcode"), 
				one_of("catalog_number"),
				one_of("dataset_id"),
				one_of("batch") ,
				one_of("development_stage"),
				one_of("Fetus_id"),
				one_of("individual"),
				one_of("orig.ident")
			) |>
			when("batch" %in% colnames(.) ~ mutate(., batch = batch |> as.character() ), ~(.)) 
	) |> 

	mutate_if(is.factor, as.character)

# x |> 
# 	unite("sample___", c(
# 		Sample, sampleID, sample_uuid, SampleID, Sample_ID, 
# 		scRNASeq_sample_ID, Sample_Tag, Sample.ID, sample_names, 
# 		Short_Sample, Sample.ID.short, Sample.name,
# 		Specimen.ID, patient, Donor.ID, donor_id, donor, PatientID, donor_uuid, 
# 		library_uuid, suspension_uuid, Patient, tissue_section_uuid, DonorID, specimen, SpecimenID, Fetus_id, individual
# 	), na.rm = TRUE, sep = "___") |> 
# 	mutate_if(is.factor, as.character) 
	


y=x |> 
	mutate(PatientID = if_else(is.na(PatientID), sampleID, PatientID)) |> 
	select(-sampleID) |> 
	mutate(Sample = case_when(is.na(Patient) ~ Sample)) |> 
	replace_na(list(PatientID = "")) |> 
	mutate(PatientID = case_when(
		dataset_id=="e40591e7-0e5a-4bef-9b60-7015abe5b17f" ~ glue("{batch} {development_stage}") |> as.character(),
		dataset_id=="39b6cc45-8c5c-4f7b-944c-58f66da5efb1" ~ sample_id,
		dataset_id=="443d6a0e-dbcb-4002-8af0-628e7d4a18fa" ~ sample_id,
		dataset_id=="a91f075b-52d5-4aa3-8ecc-86c4763a49b3" ~ sample,
		dataset_id=="0af763e1-0e2f-4de6-9563-5abb0ad2b01e" ~ "only_one_culture",
		dataset_id=="d6f92754-e178-4202-b86f-0f430e965d72" ~ orig.ident,
		T ~ PatientID
	)) |> 
	
	unite("sample___", c(
		Sample, SampleID, sample_uuid, Sample_ID, 
		scRNASeq_sample_ID, Sample_Tag, Sample.ID, sample_names, 
		Short_Sample, Sample.ID.short, Sample.name,
		patient, Donor.ID, donor_id, donor, PatientID, donor_uuid, 
		library_uuid, suspension_uuid, Patient, tissue_section_uuid, DonorID, specimen, SpecimenID, Fetus_id, individual
	), na.rm = TRUE, sep = "___", remove = F) 

y |> 
	count(dataset_id, sample___ ) |> 
	arrange(desc(n)) |> 
	View()

y |> filter(sample___ == "M1TX_210218_161_A01___H21.33.001")