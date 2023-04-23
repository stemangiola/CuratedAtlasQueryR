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


# # # CREATE MAKEFILE
# tab = "\t"
# root_directory = "/vast/projects/cellxgene_curated"
# split_data_directory = glue("{root_directory}/splitted_DB2_data_0.2")
# scaled_data_directory = glue("{root_directory}/splitted_DB2_data_scaled_0.2.1")
# 
# dir(split_data_directory) |>
# 	map( ~ glue("{scaled_data_directory}/{.x}:{split_data_directory}/{.x}\n{tab}Rscript scale_files.R {split_data_directory}/{.x} {scaled_data_directory}/{.x}")
# ) |>
# 	prepend(glue("CATEGORY=scale_data\nMEMORY=10000\nCORES=2\nWALL_TIME=30000")) |>
# 	unlist()  |>
# 	write_lines(glue("~/PostDoc/CuratedAtlasQueryR/dev/scale_files.makeflow"))





# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file = args[[1]]
output_file = args[[2]]

# Create directory
output_file |>  dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

# Read file_cell_types
data = loadHDF5SummarizedExperiment(input_file	)

# Avoid completely empty cells
col_sums = colSums(data@assays@data$X)
which_to_select = which(col_sums >0 & col_sums < Inf)

sce = SingleCellExperiment(list(counts_per_million = scuttle::calculateCPM(data[,which_to_select ,drop=FALSE ], assay.type = "X")))
rownames(sce) = rownames(data[,which_to_select  ])
colnames(sce) = colnames(data[,which_to_select  ])

# Avoid scaling zeros
which_to_select_not = setdiff(1:ncol(data), which_to_select)
sce_zero = SingleCellExperiment(list(counts_per_million = data@assays@data$X[,which_to_select_not ,drop=FALSE ]))
rownames(sce_zero) = rownames(data[,which_to_select_not  ])
colnames(sce_zero) = colnames(data[,which_to_select_not ])

sce = sce |> cbind(sce_zero)

sce = sce[,colnames(data)]



rm(data)
gc()

sce |>	saveHDF5SummarizedExperiment(output_file, replace=TRUE)
