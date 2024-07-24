#' Import and process metadata and counts for a SingleCellExperiment object
#'
#' @param sce_obj A SingleCellExperiment object from RDS, the metadata slot of which
#' must contain `cell_` and `dataset_id`
#' @param cache_dir Optional character vector of length 1. A file path on
#'   your local system to a directory (not a file) that will be used to store
#'   `metadata.parquet`
#' @param pseudobulk Optional character. Set to TRUE for generating and importing pseudobulk,
#' the metadata slot of which must contain `file_id`, `cell_type_harmonised` and `sample_`
#' 
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
#' @importFrom tidySingleCellExperiment aggregate_cells
#' @importFrom tidybulk quantile_normalise_abundance
#' @examples
#' data(sample_sce_obj)
#' import_one_sce(sample_sce_obj,
#'                cache_dir = get_default_cache_dir())
#' @references Mangiola, S., M. Milton, N. Ranathunga, C. S. N. Li-Wai-Suen, 
#'   A. Odainic, E. Yang, W. Hutchison et al. "A multi-organ map of the human 
#'   immune system across age, sex and ethnicity." bioRxiv (2023): 2023-06.
#'   doi:10.1101/2023.06.08.542671.
#' @source [Mangiola et al.,2023](https://www.biorxiv.org/content/10.1101/2023.06.08.542671v3)
import_one_sce <- function(
    sce_obj,  
    cache_dir = get_default_cache_dir(),
    pseudobulk = FALSE
  ) {
  original_dir <- file.path(cache_dir, "original")
  
  # Identify metadata and counts matrix
  metadata_tbl <- metadata(sce_obj)$data
  counts_matrix <- assay(sce_obj)
  
  # Check if the input contains only one dataset_id
  (metadata_tbl |> distinct(dataset_id) |> nrow() == 1) |>
    check_true() |> 
    assert("sce_obj metadata should contain one dataset_id at each time.")
  
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
  file_id_db <- metadata_tbl$dataset_id |> unique() |> openssl::md5() |> as.character()
  metadata_tbl <-
    metadata_tbl |> mutate(file_id_db = file_id_db)
  metadata(sce_obj)$data <- metadata_tbl
  
  # Remove existing reducedDim slot to enable get_SCE API functionality 
  if (length(names(reducedDims(sce_obj))) >0 ) {
    reducedDims(sce_obj) <- NULL
  }
  
    # Pseudobulk checkpoint 
  pseudobulk_sample <- c("sample_", "cell_type_harmonised")
  if (isTRUE(pseudobulk)) {
    assert(
      all(pseudobulk_sample %in% (colData(sce_obj) |> colnames()) ),
      "Character in pseudobulk_sample must contain sample_ and cell_type_harmonised columns 
      in the SingleCellExperiment column metadata")
    
    assert(c(pseudobulk_sample, "file_id") %in% (names(metadata_tbl)) |> all() ,
           "SingleCellExperiment metadata must at least contain sample_, cell_type_harmonised,
           file_id for pseudobulk generation"
    ) }
  
  # Create original and cpm folders in the cache directory if not exist
  if (!dir.exists(original_dir)) {
    cache_dir |> file.path("original") |> dir.create(recursive = TRUE)
  }
  
  if (!dir.exists(file.path(cache_dir, "cpm"))) {
    cache_dir |> file.path("cpm") |> dir.create(recursive = TRUE)
  }
  
  # Check whether count H5 directory has been generated
  all(!file_id_db %in% dir(original_dir)) |>
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
  
  original_path <- file.path(original_dir, basename(file_id_db))
  cpm_path <- file.path(cache_dir, "cpm", basename(file_id_db))
  
  # Generate cpm from counts
  cli_alert_info("Generating cpm from {file_id_db}. ")
  get_counts_per_million(input_sce_obj = sce_obj, output_dir = cpm_path, hd5_file_dir = original_path)
  saveHDF5SummarizedExperiment(sce_obj, original_path, replace=TRUE)
  cli_alert_info("cpm are generated in {.path {cpm_path}}. ")
  
  # check metadata sample file ID match the count file ID in cache directory
  all(metadata_tbl |> pull(.data$file_id_db) %in% dir(original_dir)) |> 
    assert("The filename for count assay, which matches the file_id_db column in 
           the metadata, already exists in the cache directory.")
  
  # convert metadata_tbl to parquet if above checkpoints pass
  arrow::write_parquet(metadata_tbl, file.path(cache_dir, "metadata.parquet"))
  
  # generate pseudobulk counts and quantile_normalised counts
  if (isTRUE(pseudobulk)) {
    file_id <- metadata_tbl$file_id |> unique() |> as.character()
      
    cli_alert_info("Generating pseudobulk counts from {file_id}. ")
    pseudobulk_counts <- sce_obj |> aggregate_cells(c(sample_, cell_type_harmonised)) 
    
    normalised_counts_best_distribution <- assay(pseudobulk_counts, "counts") |> as.matrix() |>
      preprocessCore::normalize.quantiles.determine.target()
    
    normalised_counts <- pseudobulk_counts |> quantile_normalise_abundance(
      method="preprocesscore_normalize_quantiles_use_target",
      target_distribution = normalised_counts_best_distribution
    )
    
    assay(normalised_counts, "counts") <- NULL
    names(assays(normalised_counts)) <- "quantile_normalised"
    
    if (!dir.exists(file.path(cache_dir, "pseudobulk/original"))) {
      cache_dir |> file.path("pseudobulk/original") |> dir.create(recursive = TRUE)
    }
    
    if (!dir.exists(file.path(cache_dir, "pseudobulk/quantile_normalised"))) {
      cache_dir |> file.path("pseudobulk/quantile_normalised") |> dir.create(recursive = TRUE)
    }
    
    path <- file.path(cache_dir, "pseudobulk")
    pseudobulk_counts_path <- file.path(path, "original", basename(file_id))
    pseudobulk_qnorm_path <- file.path(path, "quantile_normalised", basename(file_id))
    
    saveHDF5SummarizedExperiment(pseudobulk_counts, pseudobulk_counts_path )
    saveHDF5SummarizedExperiment(normalised_counts, pseudobulk_qnorm_path )
    
    cli_alert_info("pseudobulk are generated in {.path {path}}. ")
  }
  
}


