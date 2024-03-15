#' Generating counts per million from SingleCellExperiment object
#'
#' @param input_sce_obj A SingleCellExperiment object read from RDS
#' @param output_dir A character vector of counts per million path
#' @param hd5_file_dir A character vector of HDF5 file path after converting from RDS
#' @export
#' @return A directory stores counts per million 
#' @examples
#' input_sce_obj <- sample_sce_obj
#' output_dir <- "~/projects/caq/cache_for_testing/cpm/8df700ab083ab215e618fe732e2f3dfe"
#' hd5_file_dir <- "~/projects/caq/cache_for_testing/original/8df700ab083ab215e618fe732e2f3dfe"
#' get_counts_per_million(input_sce_obj, output_dir, hd5_file_dir)
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom dplyr tbl
#' @importFrom glue glue
#' @importFrom HDF5Array loadHDF5SummarizedExperiment saveHDF5SummarizedExperiment
#' @importFrom scuttle calculateCPM
#' @importFrom DelayedArray realize
#' @importFrom purrr map

get_counts_per_million <- function(input_sce_obj, output_dir, hd5_file_dir) {
  
  # Create directories
  output_dir |>  dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)
  output_dir |> dir.create( showWarnings = FALSE, recursive = TRUE)
  hd5_file_dir |> dir.create(showWarnings = FALSE, recursive = TRUE)
  
  saveRDS(input_sce_obj, file.path(hd5_file_dir, "se.rds"))
  
  data = input_sce_obj
  
  # Assign HDF5 matrix for counts data
  hdf5_file <- file.path(hd5_file_dir, "assays.h5")
  hdf5_counts <- writeHDF5Array(assays(data)$X, filepath = hdf5_file, name = "assays.h5")
  assays(data, withDimnames = FALSE)[["X"]] <- hdf5_counts
  
  # Avoid completely empty cells
  col_sums = colSums(as.matrix(data@assays@data$X))
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
  
  # Check if there is a memory issue 
  assays(sce) <- assays(sce) |> purrr::map(realize)
  
  sce |>	saveRDS(file.path(output_dir,"se.rds"))
  sce |> saveHDF5SummarizedExperiment(output_dir, replace = TRUE)
} 

