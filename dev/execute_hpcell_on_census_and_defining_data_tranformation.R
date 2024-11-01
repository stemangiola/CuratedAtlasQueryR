library(glue)
library(dplyr)
library(tibble)
library(glue)
library(purrr)
library(targets)
library(stringr)
library(HPCell)
library(crew.cluster)
library(arrow)
library(CuratedAtlasQueryR)
directory = "/vast/scratch/users/shen.m/Census_rerun/split_h5ad_based_on_sample_id/"
sample_anndata <- dir(glue("{directory}"), full.names = T)
downloaded_samples_tbl <- read_parquet("/vast/scratch/users/shen.m/Census_rerun/census_samples_to_download_groups.parquet")
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

# setwd(glue("{store}"))
sliced_sample_tbl = 
  sample_tbl |> 
  select(file_name, tier, cell_number, dataset_id, sample_2, method_to_apply)

# sliced_sample_tbl =
#   sliced_sample_tbl |>
#   slice(1:20)

# Enable sample_names.rds to store sample names for the input
sample_names <-
  sliced_sample_tbl |> 
  pull(file_name) |> 
  str_replace("/home/users/allstaff/shen.m/scratch", "/vast/scratch/users/shen.m") |> 
  set_names(sliced_sample_tbl |> pull(sample_2))
tiers = sliced_sample_tbl |> pull(tier)
functions = sliced_sample_tbl |>  pull(method_to_apply)
# rm(sliced_sample_tbl)
# gc()

job::job({
  
  library(HPCell)
  
  sample_names |>
    initialise_hpc(
      gene_nomenclature = "ensembl",
      data_container_type = "anndata",
      store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/census_hpcell_oct_2024/target_store",
      # store = "/vast/projects/cellxgene_curated/census_hpcell_oct_2024/target_store_1_20",
      tier = tiers,
      computing_resources = list(
        crew_controller_slurm(
          name = "tier_1", 
          script_lines = "#SBATCH --mem 8G",
          slurm_cpus_per_task = 1, 
          workers = 300, 
          tasks_max = 10,
          verbose = T,
          launch_max = 5
        ),
        
        crew_controller_slurm(
          name = "tier_2",
          script_lines = "#SBATCH --mem 10G",
          slurm_cpus_per_task = 1,
          workers = 200,
          tasks_max = 10,
          verbose = T,
          launch_max = 5
        ),
        crew_controller_slurm(
          name = "tier_3",
          script_lines = "#SBATCH --mem 15G",
          slurm_cpus_per_task = 1,
          workers = 100,
          tasks_max = 10,
          verbose = T,
          launch_max = 5
        ),
        crew_controller_slurm(
          name = "tier_4",
          script_lines = "#SBATCH --mem 70G",
          slurm_cpus_per_task = 1,
          workers = 100,
          tasks_max = 10,
          verbose = T,
          launch_max = 5
        )
      ),
      verbosity = "summary",
      # debug_step = "annotation_tbl_tier_4",
      update = "never", 
      error = "continue",
      garbage_collection = 100, 
      workspace_on_error = TRUE
      
    ) |> 
    tranform_assay(fx = functions, target_output = "sce_transformed") |>
    
    # # Remove empty outliers based on RNA count threshold per cell
    remove_empty_threshold(target_input = "sce_transformed", RNA_feature_threshold = 200) |>
    
    # Annotation
    annotate_cell_type(target_input = "sce_transformed", azimuth_reference = "pbmcref") |> 
    
    print()
  
  
})



tar_meta(starts_with("annotation_tbl"), store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/census_hpcell_oct_2024/target_store") |> 
  filter(!error |> is.na()) |> arrange(desc(time)) |> select(error, name)

# I have to check the input of this NULL target 
# annotation_tbl_tier_1_ecac957542df0c20

tar_workspace(annotation_tbl_tier_4_8fcdb6edc543d3ea, store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/census_hpcell_oct_2024/target_store")
annotation_label_transfer(sce_transformed_tier_4, empty_droplets_tbl = empty_tbl_tier_4, reference_azimuth = "pbmcref", feature_nomenclature = gene_nomenclature)

tar_read_raw("annotation_tbl_tier_1_ecac957542df0c20", store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/census_hpcell_oct_2024/target_store")


# tar_invalidate(starts_with("annotation_tbl_tier_"), store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/census_hpcell_oct_2024/target_store")
