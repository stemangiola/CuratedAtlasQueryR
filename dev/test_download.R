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

# db <- db()

# datasets(db) |> 
# 	distinct(organism) |> 
# 	unnest_wider(organism) |> 
# 	unnest_wider(organism)


# datasets(db) |> 
# 	left_join(
# 		files(db) |> filter(filetype=="H5AD"), 
# 		by = "dataset_id"
# 	) |> 
# 	
# 	# Get organism list and filter human
# 	mutate(organism_name = map_chr(organism, ~ .x |> map(~.x$label) |> paste(collapse=", ") )) |>  
# 	filter(organism_name |> str_detect("Homo sapiens")) |> 
# 	
# 	# Download
# 	files_download(dry.run = FALSE, cache_path = "/vast/scratch/users/mangiola.s/human_cell_atlas/raw_data/") |> 
# 	
# 	# Save file list
# 	saveRDS("/vast/scratch/users/mangiola.s/human_cell_atlas/file_location.rds")

local_file = readRDS("/vast/scratch/users/mangiola.s/human_cell_atlas/file_location.rds")

# Get all genes
gene_names = 
	local_file |> 
	map(	~ rowData(readH5AD(.x, reader = "R", use_hdf5 = TRUE))$feature_name
					)  |> 
	unlist() |> 
	unique() 

matedata = readRDS("metadata.rds") 


all_data = 
	local_file |> 
	#head(2) |> 
	map(
		~ {
			cat(".")
			h5 = readH5AD(.x, reader = "R", use_hdf5 = TRUE) 
			rownames(h5) = rowData(h5)$feature_name
			rowData(h5) = NULL
			colData(h5) = NULL
				# colData(h5) |> 
				# as_tibble()  |> 
				# bind_rows(metadata |> slice(0)) |> 
				# select(colnames(metadata)) |> 
				# type_convert(guess_integer = TRUE) |> 
				# DataFrame()
			
			reducedDims(h5) = NULL
			
			# Complete gene set
			missing_genes = gene_names  |> setdiff(rownames(h5))

			missing_matrix = HDF5RealizationSink(c(length(missing_genes),ncol(h5)), as.sparse = TRUE) |> 
			as("DelayedArray")
			
			rownames(missing_matrix) = missing_genes
			colnames(missing_matrix) = colnames(h5)
			
			missing_sce = SingleCellExperiment(list(X=missing_matrix),  colData=colData(h5))
			missing_sce@int_colData = h5@int_colData
			
			
			# Select just the X assay
			h5@assays@data = h5@assays@data |> as.list() %>% .[1] |> SimpleList()
			
			# Make cell name unique
			h5 = h5 |> rbind(missing_sce	)
			h5 = h5 |> tidySingleCellExperiment::mutate(dataset_id = basename(.x) |> tools::file_path_sans_ext())
			
			if(is.null(colnames(h5))) colnames(h5) = seq_len(ncol(h5)) |> as.character()
			colnames(h5) = glue("{colnames(h5)}_{basename(.x) |> tools::file_path_sans_ext()}")
			h5
	
		} 
	)  %>%
	do.call(cbind, .)

# Save on disk
all_data |> 
	mutate(save_by = cell_type |> str_replace_all(" ", "_") |> str_replace_all(",", "_")) |> 
	nest(data = -save_by) |> 
	mutate(saved = map2(
		data, save_by
		~ {
			print(.y)
			writeH5AD(.x, glue("{.y}.hdf5"))
		}
	))





x = readH5AD("adipocyte_of_epicardial_fat_of_left_ventricle.hdf5", use_hdf5 = TRUE)



xdelay |> 
saveHDF5SummarizedExperiment("delay")
delay = loadHDF5SummarizedExperiment("delay")
