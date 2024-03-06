#' Generating counts per million
#'
#' @param input_file A character vector of raw counts path
#' @param output_file A character vector of counts per million path
#' @export
#' @return A directory stores counts per million 
#' @examples
#' input_file <- "/var/folders/ls/99n281zx4bbd73kllmc1rc0h0005lj/T/RtmpUjWDqt/original/12eb5fe25994253c1d320ca590a6e681"
#' output_file <- "/var/folders/ls/99n281zx4bbd73kllmc1rc0h0005lj/T/RtmpUjWDqt/cpm/12eb5fe25994253c1d320ca590a6e681"
#' get_counts_per_million(input_file, output_file)
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom dplyr tbl
#' @importFrom glue glue
#' @importFrom HDF5Array loadHDF5SummarizedExperiment
#' @importFrom scuttle calculateCPM

get_counts_per_million <- function(input_file, output_file) {
  
  # Create directory
  output_file |>  dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)
  
  # Read file_cell_types
  data = loadHDF5SummarizedExperiment(input_file)
  
  # Avoid completely empty cells
  col_sums = colSums(data@assays@data$X)
  which_to_select = which(col_sums >0 & col_sums < Inf)
  
  sce = SingleCellExperiment(list(counts_per_million = calculateCPM(data[,which_to_select ,drop=FALSE ], assay.type = "X")))
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
} 

