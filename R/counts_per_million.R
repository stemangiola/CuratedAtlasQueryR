#' Generating counts per million from SingleCellExperiment object
#'
#' @param input_sce_obj A SingleCellExperiment object read from RDS
#' @param output_dir A character vector of counts per million path
#' @param hd5_file_dir A character vector of HDF5 file path after converting from RDS
#' @return A directory stores counts per million 
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom dplyr tbl
#' @importFrom glue glue
#' @importFrom HDF5Array saveHDF5SummarizedExperiment
#' @importFrom purrr map
#' @noRd
get_counts_per_million <- function(input_sce_obj, output_dir, hd5_file_dir) {
  
  # Create directories
  output_dir |>  dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

  # Save SCE to the cache directory counts folder
  input_sce_obj |> saveHDF5SummarizedExperiment(hd5_file_dir)
  
  data <- input_sce_obj
  
  # Avoid completely empty cells
  col_sums <- colSums(as.matrix(data@assays@data$X))
  which_to_select <- which(col_sums >0 & col_sums < Inf)
  
  sce <- SingleCellExperiment(list(counts_per_million = scuttle::calculateCPM(data[,which_to_select ,drop=FALSE ], assay.type = "X")))
  rownames(sce) <- rownames(data[,which_to_select  ])
  colnames(sce) <- colnames(data[,which_to_select  ])
  
  # Avoid scaling zeros
  which_to_select_not <- setdiff(1:ncol(data), which_to_select)
  sce_zero <- SingleCellExperiment(list(counts_per_million = data@assays@data$X[,which_to_select_not ,drop=FALSE ]))
  rownames(sce_zero) <- rownames(data[,which_to_select_not  ])
  colnames(sce_zero) <- colnames(data[,which_to_select_not ])
  
  sce <- sce |> cbind(sce_zero)
  
  sce <- sce[,colnames(data)]

  rm(data)
  gc()
  
  # Check if there is a memory issue 
  SummarizedExperiment::assays(sce) <- SummarizedExperiment::assays(sce) |> map(DelayedArray::realize)
  
  sce |> saveHDF5SummarizedExperiment(output_dir)
} 

