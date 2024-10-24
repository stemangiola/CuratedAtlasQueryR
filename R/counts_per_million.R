#' Generating counts per million from a SingleCellExperiment object
#'
#' @param input_sce_obj A SingleCellExperiment object read from RDS
#' @param output_dir A character vector of counts per million path
#' @param hd5_file_dir A character vector of HDF5 file path after converting from RDS
#' @return A directory stores counts per million 
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assay assays assays<-
#' @importFrom HDF5Array saveHDF5SummarizedExperiment
#' @importFrom purrr map
#' @noRd
get_counts_per_million <- function(input_sce_obj, output_dir, hd5_file_dir) {

  # Save SCE to the cache directory counts folder
  input_sce_obj |> saveHDF5SummarizedExperiment(hd5_file_dir, replace=TRUE)
  
  data <- input_sce_obj
  
  # Avoid completely empty cells
  col_sums <- colSums(as.matrix(assay(data)))
  selected_cols <- which(col_sums >0 & col_sums < Inf)
  
  assay_name <- assays(data) |> names()
  sce <- SingleCellExperiment(list(counts_per_million = scuttle::calculateCPM(data[,selected_cols ,drop=FALSE ], assay.type = assay_name)))
  rownames(sce) <- rownames(data[,selected_cols  ])
  colnames(sce) <- colnames(data[,selected_cols  ])
  
  # Avoid scaling zeros
  sce_zero <- SingleCellExperiment(list(counts_per_million = assay(data)[, !selected_cols ,drop=FALSE ]))
  rownames(sce_zero) <- rownames(data[, !selected_cols  ])
  colnames(sce_zero) <- colnames(data[, !selected_cols ])
  
  sce <- sce |> cbind(sce_zero)
  
  sce <- sce[,colnames(data)]

  rm(data)
  
  # Check if there is a memory issue 
  assays(sce) <- assays(sce) |> map(DelayedArray::realize)
  
  sce |> saveHDF5SummarizedExperiment(output_dir, replace = TRUE)
} 

