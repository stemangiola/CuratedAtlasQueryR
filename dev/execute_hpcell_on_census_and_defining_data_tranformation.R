library(dplyr)
library(tibble)
library(glue)
library(purrr)
library(stringr)
library(HPCell)
library(arrow)
library(CuratedAtlasQueryR)
directory = "/vast/scratch/users/shen.m/Census_rerun/split_h5ad_based_on_sample_id/"
sample_anndata <- dir(glue("{directory}"), full.names = T)
downloaded_samples_tbl <- read_parquet("/vast/scratch/users/shen.m/Census_rerun/census_samples_to_download_groups.parquet")
downloaded_samples_tbl <- downloaded_samples_tbl |>
  dplyr::rename(cell_number = list_length) |>
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
      store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/census_hpcell_oct_2024/target_store/",
      gene_nomenclature = "ensembl",
      data_container_type = "anndata",
      tier = tiers,
      computing_resources = list(
        crew_controller_slurm(
          name = "tier_1", 
          script_lines = "#SBATCH --mem 8G",
          slurm_cpus_per_task = 1, 
          workers = 300, 
          tasks_max = 10,
          verbose = T,
          launch_max = 10, 
          seconds_idle = 30
        ),
        
        crew_controller_slurm(
          name = "tier_2",
          script_lines = "#SBATCH --mem 10G",
          slurm_cpus_per_task = 1,
          workers = 200,
          tasks_max = 10,
          verbose = T,
          launch_max = 10, 
          seconds_idle = 30
        ),
        crew_controller_slurm(
          name = "tier_3",
          script_lines = "#SBATCH --mem 15G",
          slurm_cpus_per_task = 1,
          workers = 100,
          tasks_max = 10,
          verbose = T,
          launch_max = 5, 
          seconds_idle = 30
        ),
        crew_controller_slurm(
          name = "tier_4",
          script_lines = "#SBATCH --mem 200G", # There should be a tier 5 with 12Gb for sample of 40K cells
          slurm_cpus_per_task = 1,
          workers = 100,
          tasks_max = 10,
          verbose = T,
          launch_max = 5, 
          seconds_idle = 30
        )
      ),
      verbosity = "summary",
      # debug_step = "annotation_tbl_tier_4",
      update = "never", 
      error = "continue",
      garbage_collection = 100, 
      workspace_on_error = TRUE
      
    ) |> 
    transform_assay(fx = functions, target_output = "sce_transformed") |>
    
    # # Remove empty outliers based on RNA count threshold per cell
    remove_empty_threshold(target_input = "sce_transformed", RNA_feature_threshold = 200) |>
    
    # Annotation
    annotate_cell_type(target_input = "sce_transformed", azimuth_reference = "pbmcref") |> 
    
    print()
  
  
})



tar_meta(starts_with("annotation_tbl_"), store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/census_hpcell_oct_2024/target_store") |> 
  filter(!data |> is.na()) |> arrange(desc(time)) |> select(error, name)

# I have to check the input of this NULL target 
# annotation_tbl_tier_1_ecac957542df0c20

tar_workspace(annotation_tbl_tier_4_a95d2334c8388111, store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/census_hpcell_oct_2024/target_store")
annotation_label_transfer(sce_transformed_tier_4, empty_droplets_tbl = empty_tbl_tier_4, reference_azimuth = "pbmcref", feature_nomenclature = gene_nomenclature)


#' Pipeline for Lightening Annotations in High-Performance Computing Environment
#' 
#' This pipeline is designed to read, process, and "lighten" large annotation tables in an HPC environment.
#' It uses the `targets` package for reproducibility and `crew` for efficient job scheduling on a Slurm cluster.
#' The `lighten_annotation` function selects and processes specific columns from large tables to reduce memory usage.
#' 
#' The pipeline consists of:
#' - **Crew Controllers**: Four tiers of Slurm controllers with varying memory allocations to optimize resource usage.
#' - **Targets**:
#'   - `my_store`: Defines the path to the target storage directory, ensuring all targets use the correct storage location.
#'   - `target_name`: Retrieves metadata to identify branch targets for annotation.
#'   - `annotation_tbl_light`: Applies `lighten_annotation` to process each target name, optimally running with `tier_1` resources.
#' 
#' @libraries:
#'   - `dplyr`, `magrittr`, `tibble`, `targets`, `tarchetypes` for data manipulation and pipeline structure.
#'   - `crew`, `crew.cluster` for parallel computation and cluster scheduling in an HPC environment.
#' 
#' @options:
#'   - Memory settings, garbage collection frequency, and error handling are set to handle large data efficiently.
#'   - The `cue` option is set to `never` for forced target updates if needed.
#'   - `controller` is a group of Slurm controllers to manage computation across memory tiers.
#' 
#' @function `lighten_annotation`: Processes each annotation table target, unnesting and selecting specific columns to reduce data size.
#'
#' @example Usage:
#'   The pipeline script is saved as `/vast/scratch/users/mangiola.s/lighten_annotation_tbl_target.R` and can be run using `tar_make()`.
tar_script({
  library(dplyr)
  library(magrittr)
  library(tibble)
  library(targets)
  library(tarchetypes)
  library(crew)
  library(crew.cluster)
  tar_option_set(
    memory = "transient", 
    garbage_collection = 100, 
    storage = "worker", 
    retrieval = "worker", 
    error = "continue", 
    #debug = "annotation_tbl_light", 
    cue = tar_cue(mode = "never"), 
    controller = crew_controller_group(
      list(
        crew_controller_slurm(
          name = "tier_1", 
          script_lines = "#SBATCH --mem 8G",
          slurm_cpus_per_task = 1, 
          workers = 200, 
          tasks_max = 10,
          verbose = T,
          launch_max = 5, 
          seconds_idle = 30
        ),
        
        crew_controller_slurm(
          name = "tier_2",
          script_lines = "#SBATCH --mem 10G",
          slurm_cpus_per_task = 1,
          workers = 200,
          tasks_max = 10,
          verbose = T,
          launch_max = 5, 
          seconds_idle = 30
        ),
        crew_controller_slurm(
          name = "tier_3",
          script_lines = "#SBATCH --mem 15G",
          slurm_cpus_per_task = 1,
          workers = 200,
          tasks_max = 10,
          verbose = T,
          launch_max = 5, 
          seconds_idle = 30
        ),
        crew_controller_slurm(
          name = "tier_4",
          script_lines = "#SBATCH --mem 50G",
          slurm_cpus_per_task = 1,
          workers = 30,
          tasks_max = 10,
          verbose = T,
          launch_max = 5, 
          seconds_idle = 30
        )
      )
    ), 
    trust_object_timestamps = TRUE, 
    workspaces = "annotation_tbl_light_ffcd3d5a64bedf1f"
  )
  
  lighten_annotation = function(target_name, my_store ){
    annotation_tbl = tar_read_raw( target_name,  store = my_store )
    if(annotation_tbl |> is.null()) { 
      warning("this annotation is null -> ", target_name)
      return(NULL) 
      }
    
    annotation_tbl |> 
      unnest(blueprint_scores_fine) |> 
      select(.cell, blueprint_first.labels.fine, monaco_first.labels.fine, any_of("azimuth_predicted.celltype.l2"), monaco_scores_fine, contains("macro"), contains("CD4") ) |> 
      unnest(monaco_scores_fine) |> 
      select(.cell, blueprint_first.labels.fine, monaco_first.labels.fine, any_of("azimuth_predicted.celltype.l2"), contains("macro") , contains("CD4"), contains("helper"), contains("Th")) |> 
      rename(cell_ = .cell)
  }
  
  list(
    
    # The input DO NOT DELETE
    tar_target(my_store, "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/census_hpcell_oct_2024/target_store", deployment = "main"),
    
    tar_target(
      target_name,
      tar_meta(
        starts_with("annotation_tbl_tier_"), 
        store = my_store) |> 
        filter(type=="branch") |> 
        pull(name),
      deployment = "main"
    )    ,
    
    tar_target(
      annotation_tbl_light,
      lighten_annotation(target_name, my_store),
      packages = c("dplyr", "tidyr"),
      pattern = map(target_name),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "tier_1")
      )
    )
  )
  
  
}, script = "/vast/scratch/users/mangiola.s/lighten_annotation_tbl_target.R", ask = FALSE)

job::job({
  
  tar_make(
    script = "/vast/scratch/users/mangiola.s/lighten_annotation_tbl_target.R", 
    store = "/vast/scratch/users/mangiola.s/lighten_annotation_tbl_target", 
    reporter = "summary"
  )
  
})

# Sample metadata
library(arrow)
library(dplyr)
library(duckdb)

write_parquet_to_parquet = function(data_tbl, output_parquet, compression = "gzip") {
  
  # Establish connection to DuckDB in-memory database
  con_write <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  # Register `data_tbl` within the DuckDB connection (this doesn't load it into memory)
  duckdb::duckdb_register(con_write, "data_tbl_view", data_tbl)
  
  # Use DuckDB's COPY command to write `data_tbl` directly to Parquet with compression
  copy_query <- paste0("
  COPY data_tbl_view TO '", output_parquet, "' (FORMAT PARQUET, COMPRESSION '", compression, "');
  ")
  
  # Execute the COPY command
  dbExecute(con_write, copy_query)
  
  # Unregister the temporary view
  duckdb::duckdb_unregister(con_write, "data_tbl_view")
  
  # Disconnect from the database
  dbDisconnect(con_write, shutdown = TRUE)
}


# Write annotation light
cell_metadata <- 
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata.parquet')")
  ) |>
  mutate(cell_ = paste0(cell_, "___", dataset_id)) |> 
  select(cell_, observation_joinid, contains("cell_type"), dataset_id,  self_reported_ethnicity, tissue, donor_id,  sample_id, is_primary_data, assay)


cell_annotation = 
  tar_read(annotation_tbl_light, store = "/vast/scratch/users/mangiola.s/lighten_annotation_tbl_target") |> 
  rename(
    blueprint_first_labels_fine = blueprint_first.labels.fine, 
    monaco_first_labels_fine = monaco_first.labels.fine, 
    azimuth_predicted_celltype_l2 = azimuth_predicted.celltype.l2
  ) 

empty_droplet = 
  tar_read(empty_tbl_tier_1, store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/census_hpcell_oct_2024/target_store") |> 
  c(tar_read(empty_tbl_tier_2, store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/census_hpcell_oct_2024/target_store")) |> 
  c(tar_read(empty_tbl_tier_3, store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/census_hpcell_oct_2024/target_store")) |> 
  c(tar_read(empty_tbl_tier_4, store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/census_hpcell_oct_2024/target_store")) |> 
  bind_rows() |> 
  rename(cell_ = .cell)
  

cell_metadata |> 
  left_join(empty_droplet, copy=TRUE) |>  
  left_join(cell_annotation, copy=TRUE) |> 
  write_parquet_to_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_annotation.parquet")

system("~/bin/rclone copy /vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_annotation.parquet box_adelaide:/Mangiola_ImmuneAtlas/reannotation_consensus/")

tar_workspace(
  "annotation_tbl_light_ffcd3d5a64bedf1f", 
  script = "/vast/scratch/users/mangiola.s/lighten_annotation_tbl_target.R", 
  store = "/vast/scratch/users/mangiola.s/lighten_annotation_tbl_target", 
)

# tar_workspaces(annotation_tbl_light_c8078b8175604dd3,  store = "/vast/scratch/users/mangiola.s/lighten_annotation_tbl_target")

# tar_invalidate(starts_with("annotation_tbl_tier_"), store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/census_hpcell_oct_2024/target_store")
