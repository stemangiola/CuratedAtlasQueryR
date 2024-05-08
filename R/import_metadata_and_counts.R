#' Import and process metadata and counts for a SingleCellExperiment object
#'
#' @param sce_obj A SingleCellExperiment object from RDS
#' @param cache_dir Optional character vector of length 1. A file path on
#'   your local system to a directory (not a file) that will be used to store
#'   `metadata.parquet`
#' @export
#' @return A metadata.parquet strip from the SingleCellExperiment object. 
#' Directories store counts and counts per million in the provided cache directory.
#' @importFrom checkmate check_true check_character check_subset assert
#' @importFrom dplyr select distinct pull
#' @importFrom cli cli_alert_info
#' @importFrom rlang .data
#' @importFrom SingleCellExperiment reducedDims rowData reducedDims<-
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment assay
#' @importFrom stringr str_detect
#' @examples
#' data(sample_sce_obj)
#' import_metadata_counts(sample_sce_obj,
#'                        cache_dir = get_default_cache_dir())
import_metadata_counts <- function(
    sce_obj,  
    cache_dir = get_default_cache_dir()
  ) {
  original_dir <- file.path(cache_dir, "original")
  
  # Identify metadata and counts matrix
  metadata_tbl <- metadata(sce_obj)$data
  counts_matrix <- assay(sce_obj)
  
  # Identify whether genes in SingleCellxExperiment object are in ensembl nomenclature
  genes <- rowData(sce_obj) |> rownames()
  assert(sce_obj |> inherits( "SingleCellExperiment"),
              "sce_obj is not identified as SingleCellExperiment object.")
  assert(!str_detect(genes, "^ENSG%") |> all(), 
              "Gene names in SingleCellExperiment object cannot contain Ensembl IDs.")
  assert(all(counts_matrix >= 0),
              "Counts for SingleCellExperiment cannot be negative.")
  
  # Convert to tibble if not provided
  metadata_tbl <- metadata_tbl |> as_tibble()
  
  # Create file_id_db from dataset_id
  metadata_tbl <-
    metadata_tbl |> mutate(file_id_db = .data$dataset_id |> openssl::md5() |> as.character())
  metadata(sce_obj)$data <- metadata_tbl
  
  # Remove existing reducedDim slot to enable get_SCE API functionality 
  if (length(names(reducedDims(sce_obj))) >0 ) {
    reducedDims(sce_obj) <- NULL
  }
  
  # Create original and cpm folders in the cache directory if not exist
  if (!dir.exists(original_dir)) {
    cache_dir |> file.path("original") |> dir.create(recursive = TRUE)
  }
  
  if (!dir.exists(file.path(cache_dir, "cpm"))) {
    cache_dir |> file.path("cpm") |> dir.create(recursive = TRUE)
  }
  
  # Check whether count H5 directory has been generated
  all(!metadata_tbl$file_id_db %in% dir(original_dir)) |>
    check_true() |>
    assert("The filename for count assay (file_id_db) already exists in the cache directory.")
  
  # Check the metadata contains cell_, file_id_db, sample_ with correct types
  check_true("cell_" %in% colnames(metadata_tbl))
  check_true("file_id_db" %in% names(metadata_tbl)) 
  pull(metadata_tbl, .data$cell_) |> class() |> check_character()
  select(metadata_tbl, .data$file_id_db) |> class() |> check_character()
  
  # Check cell_ values in metadata_tbl is unique
  (anyDuplicated(metadata_tbl$cell_) == 0 ) |> assert("Cell names (cell_) in the metadata must be unique.")
  
  # Check cell_ values are not duplicated when join with parquet
  cells <- select(get_metadata(cache_directory = cache_dir), .data$cell_) |> as_tibble()
  (!any(metadata_tbl$cell_ %in% cells$cell_)) |> 
    assert("Cell names (cell_) should not clash with cells that already exist in the atlas.")
  
  # Check age_days is either -99 or greater than 365
  if (any(colnames(metadata_tbl) == "age_days")) {
    assert(all(metadata_tbl$age_days==-99 | metadata_tbl$age_days> 365),
                "age_days should be either -99 for unknown or greater than 365.")
  }
  
  # Check sex capitalisation then convert to lower case 
  if (any(colnames(metadata_tbl) == "sex")) {
    metadata_tbl <- metadata_tbl |> mutate(sex = tolower(.data$sex))
    distinct(metadata_tbl, .data$sex) |> pull(.data$sex) |> check_subset(c("female","male","unknown"))
  }
  counts_path <-
    select(metadata_tbl, .data$file_id_db) |> mutate(
      original_path = file.path(original_dir, basename(.data$file_id_db)),
      cpm_path = file.path(cache_dir, "cpm", basename(.data$file_id_db))
    ) |>
    distinct()
  
  # Generate cpm from counts
  cli_alert_info("Generating cpm from {.path {counts_path$file_id_db}}. ")
  get_counts_per_million(input_sce_obj = sce_obj, output_dir = counts_path$cpm_path, hd5_file_dir = counts_path$original_path)
  saveHDF5SummarizedExperiment(sce_obj, counts_path$original_path, replace=TRUE)
  cli_alert_info("cpm are generated in {.path {counts_path$cpm_path}}. ")
  
  # check metadata sample file ID match the count file ID in cache directory
  all(metadata_tbl |> pull(.data$file_id_db) %in% dir(original_dir)) |> 
    assert("The filename for count assay, which matches the file_id_db column in the metadata, already exists in the cache directory.")
  
  # convert metadata_tbl to parquet if above checkpoints pass
  arrow::write_parquet(metadata_tbl, file.path(cache_dir, "metadata.parquet"))
}

