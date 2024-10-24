library(glue)
library(dplyr)
library(arrow)
library(SingleCellExperiment)
library(tibble)
library(SummarizedExperiment)
library(glue)
library(purrr)
library(zellkonverter)
library(tidyr)
library(ggplot2)
library(plotly)
library(targets)
library(stringr)
library(CuratedAtlasQueryR)
library(fs)
library(HPCell)
library(crew.cluster)
directory = "/home/users/allstaff/shen.m/scratch/Census_rerun/split_h5ad_based_on_sample_id/"
sample_anndata <- dir(glue("{directory}"), full.names = T)
downloaded_samples_tbl <- read_parquet("/home/users/allstaff/shen.m/scratch/Census_rerun/census_samples_to_download_groups.parquet")
downloaded_samples_tbl <- downloaded_samples_tbl |>
  rename(cell_number = list_length) |> 
  mutate(cell_number = cell_number |> as.integer(),
         file_name = glue("{directory}{sample_2}.h5ad") |> as.character(), 
         tier = case_when(
           cell_number < 500 ~ "tier_1", cell_number >= 500 &
             cell_number < 1000 ~ "tier_2", cell_number >= 1000 &
             cell_number < 10000 ~ "tier_3", cell_number >= 10000 ~ "tier_4"
         ))

result_directory = "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024"
sample_meta <- tar_read(metadata_dataset_id_common_sample_columns, store = glue("{result_directory}/_targets"))
sample_tbl = downloaded_samples_tbl |> left_join(get_metadata() |> select(dataset_id, contains("norm")) |>
                                 distinct() |> filter(!is.na(x_normalization)) |>
                                 as_tibble(), by = "dataset_id")


sample_tbl <- sample_tbl |> left_join(sample_meta, by = "dataset_id") |> distinct(file_name, tier, cell_number, dataset_id, sample_2,
                                                                                  x_normalization, x_approximate_distribution) |>
  mutate(transform_method = case_when(str_like(x_normalization, "C%") ~ "log",
                                      x_normalization == "none" ~ "log",
                                      x_normalization == "normalized" ~ "log",
                                      is.na(x_normalization) & is.na(x_approximate_distribution) ~ "log",
                                      is.na(x_normalization) & x_approximate_distribution == "NORMAL" ~ "NORMAL",
                                      is.na(x_normalization) & x_approximate_distribution == "COUNT" ~ "COUNT",
                                      str_like(x_normalization, "%canpy%") ~ "log1p",
                                      TRUE ~ x_normalization)) |>
  
  mutate(method_to_apply =  case_when(transform_method %in% c("log","LogNormalization","LogNormalize","log-normalization") ~ "exp",
                                      is.na(x_normalization) & is.na(x_approximate_distribution) ~ "exp",
                                      str_like(transform_method, "Counts%") ~ "exp",
                                      str_like(transform_method, "%log2%") ~ "exp",
                                      transform_method %in% c("log1p", "log1p, base e", "Scanpy",
                                                              "scanpy.api.pp.normalize_per_cell method, scaling factor 10000") ~ "expm1",
                                      transform_method == "log1p, base 2" ~ "expm1",
                                      transform_method == "NORMAL" ~ "exp",
                                      transform_method == "COUNT" ~ "identity",
                                      is.na(transform_method) ~ "identity"
  ) ) |>
  mutate(comment = case_when(str_like(x_normalization, "Counts%")  ~ "a checkpoint for max value of Assay must <= 50",
                             is.na(x_normalization) & is.na(x_approximate_distribution) ~ "round negative value to 0",
                             x_normalization == "normalized" ~ "round negative value to 0"
  ))

sample_tbl <- sample_tbl |> mutate(transformation_function = map(
    method_to_apply,
    ~ ( function(data) {
      assay_name <- data@assays |> names() |> magrittr::extract2(1)
      counts <- assay(data, assay_name)
      density_est <- counts |> as.matrix() |> density()
      mode_value <- density_est$x[which.max(density_est$y)]
      if (mode_value < 0 )  counts <- counts + abs(mode_value)
      
      # Scale max counts to 20 to avoid any downstream failure
      if (.x != "identity" && (max(counts) > 20)) {
        scale_factor = 20 / max(counts)
        counts <- counts * scale_factor}
      
      counts <- transform_method(counts)
      # round counts to avoid potential substraction error due to different digits print out
      counts <- counts |> round(5)
      majority_gene_counts = names(which.max(table(as.vector(counts)))) |> as.numeric()
      if (majority_gene_counts != 0) {
        counts <- counts - majority_gene_counts
      }
      
      # Avoid downstream failures negative counts
      if((counts[,1:min(10000, ncol(counts))] |> min()) < 0)
        counts[counts < 0] <- 0
      
      # Assign counts back to data
      assay(data, assay_name) <- counts
      
      col_sums <- colSums(counts)
      # Drop all zero cells
      data <- data[, col_sums > 0]
      
      # Avoid downstream binding error
      rowData(data) = NULL

      data
      
    }) |>
      # Meta programming, replacing the transformation programmatically
      substitute( env = list(transform_method = as.name(.x))) |>
      # Evaluate back to a working function
      eval()
  ))

#sample_tbl |> saveRDS("~/scratch/Census_rerun/sample_tbl_input_for_hpcell.rds")
sample_tbl <- readRDS("~/scratch/Census_rerun/sample_tbl_input_for_hpcell.rds")

# Set the parent directory where the subdirectories will be created
# parent_dir <- "~/scratch/Census_rerun/"
# 
# # Directory names to create
# dir_names <- paste0("run", 1:25)
# 
# # Full paths of the directories
# full_dir_paths <- file.path(parent_dir, dir_names)
# 
# # Create each directory if it does not exist
# for (dir_path in full_dir_paths) {
#   if (!dir_exists(dir_path)) {
#     dir_create(dir_path)
#   }
# }

# Run 1000 samples per run. Save log and result in the corresponding store
store = "~/scratch/Census_rerun/run3/"
setwd(glue("{store}"))
sliced_sample_tbl = sample_tbl |> slice(2001:3000) |> select(file_name, tier, cell_number, dataset_id,
                                                          sample_2, transformation_function)

# Enable sample_names.rds to store sample names for the input
sample_names <- sliced_sample_tbl |> pull(file_name)  |> set_names(sliced_sample_tbl |> pull(sample_2))

sample_names |>
  initialise_hpc(
    gene_nomenclature = "ensembl",
    data_container_type = "anndata",
    store = store,
    tier = sliced_sample_tbl |> pull(tier),
    computing_resources = list(
      crew_controller_slurm(
        name = "tier_1",
        script_lines = "#SBATCH --mem 35G",
        slurm_cpus_per_task = 1,
        workers = 200,
        tasks_max = 1,
        verbose = T
      ),
      
      crew_controller_slurm(
        name = "tier_2",
        script_lines = "#SBATCH --mem 60G",
        slurm_cpus_per_task = 1,
        workers = 50,
        tasks_max = 1,
        verbose = T
      ),
      crew_controller_slurm(
        name = "tier_3",
        script_lines = "#SBATCH --mem 90G",
        slurm_cpus_per_task = 1,
        workers = 25,
        tasks_max = 1,
        verbose = T
      ),
      crew_controller_slurm(
        name = "tier_4",
        script_lines = "#SBATCH --mem 100G",
        slurm_cpus_per_task = 1,
        workers = 14,
        tasks_max = 1,
        verbose = T
      )
    )
    
  ) |>
  tranform_assay(fx = sliced_sample_tbl |>
                   pull(transformation_function),
                 target_output = "sce_transformed") |>
  
  # Remove empty outliers based on RNA count threshold per cell
  remove_empty_threshold(target_input = "sce_transformed", RNA_feature_threshold = 200) |>
  
  # Remove dead cells
  remove_dead_scuttle(target_input = "sce_transformed") |>
  
  # Score cell cycle
  score_cell_cycle_seurat(target_input = "sce_transformed") |>
  
  # Remove doublets
  remove_doublets_scDblFinder(target_input = "sce_transformed") |>
  
  # Annotation
  annotate_cell_type(target_input = "sce_transformed", azimuth_reference = "pbmcref") |>
  
  normalise_abundance_seurat_SCT(
    factors_to_regress = c("subsets_Mito_percent", "subsets_Ribo_percent", "G2M.Score"),
    target_input = "sce_transformed"
  )


