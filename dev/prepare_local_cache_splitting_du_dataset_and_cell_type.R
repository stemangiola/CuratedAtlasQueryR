
# Group samples by dataset_id, cell_type

# This script sets up a robust and scalable data processing pipeline for single-cell RNA sequencing (scRNA-seq) datasets using the targets package in R, which facilitates reproducible and efficient workflows. Specifically, the code orchestrates the ingestion and preprocessing of multiple SingleCellExperiment objects corresponding to different datasets (dataset_id) and targets (target_name). It leverages high-performance computing resources through the crew package, configuring multiple SLURM-based controllers (tier_1 to tier_4) to handle varying computational loads efficiently.
# 
# The pipeline performs several key steps:
#   
#   1.	Data Retrieval: It reads raw SingleCellExperiment objects for each target, ensuring that only successfully loaded data proceeds further.
# 2.	Normalization: Calculates Counts Per Million (CPM) for each cell to normalize gene expression levels across cells and samples.
# 3.	Data Aggregation: Groups the data by dataset_id and tar_group, then combines the SingleCellExperiment objects within each group into a single object, effectively consolidating the data for each dataset.
# 4.	Metadata Integration: Joins additional metadata, such as cell types, by connecting to a DuckDB database and fetching relevant information from a Parquet file. This enriches the single-cell data with essential annotations.
# 5.	Cell Type Segmentation: Splits the combined SingleCellExperiment objects into separate objects based on cell_type, facilitating downstream analyses that are specific to each cell type.
# 6.	Data Saving with Error Handling: Generates unique identifiers for each cell type within a dataset and saves both the raw counts and CPM-normalized data to specified directories. It includes special handling for cases where a cell type has only one cell, duplicating the data to prevent errors during the saving process.
# 
# By integrating targets, crew, and various data manipulation packages (dplyr, tidyverse, SingleCellExperiment), this script ensures that large-scale scRNA-seq data processing is efficient, reproducible, and capable of leveraging parallel computing resources. It is designed to handle edge cases gracefully and provides a clear framework for preprocessing scRNA-seq data, which is essential for subsequent analyses such as clustering, differential expression, and cell type identification.

library(targets)
library(tidyverse)
store_file_cellNexus = "/vast/scratch/users/mangiola.s/targets_prepare_database_split_datasets_chunked"




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
                 storage = "worker", retrieval = "worker", error = "continue", 
                 #  debug = "dataset_id_sce_ce393fc1e85f2cbc", 
                 cue = tar_cue(mode = "never"), 
                  workspace_on_error = TRUE,
                 controller = crew_controller_group(
                   list(
                     crew_controller_slurm(
                       name = "tier_1", 
                       script_lines = "#SBATCH --mem 8G",
                       slurm_cpus_per_task = 1, 
                       workers = 300, 
                       tasks_max = 10,
                       verbose = T,
                       launch_max = 5, 
                       seconds_idle = 30
                     ),
                     
                     crew_controller_slurm(
                       name = "tier_2",
                       script_lines = "#SBATCH --mem 10G",
                       slurm_cpus_per_task = 1,
                       workers = 300,
                       tasks_max = 10,
                       verbose = T,
                       launch_max = 5, 
                       seconds_idle = 30
                     ),
                     crew_controller_slurm(
                       name = "tier_3",
                       script_lines = "#SBATCH --mem 30G",
                       slurm_cpus_per_task = 1,
                       workers = 300,
                       tasks_max = 10,
                       verbose = T,
                       launch_max = 5, 
                       seconds_idle = 30
                     ),
                     crew_controller_slurm(
                       name = "tier_4",
                       script_lines = "#SBATCH --mem 80G",
                       slurm_cpus_per_task = 1,
                       workers = 40,
                       tasks_max = 10,
                       verbose = T,
                       launch_max = 5, 
                       seconds_idle = 30
                     )
                   )
                 ), 
                 trust_object_timestamps = TRUE
  )
  

  save_anndata = function(target_name_grouped_by_dataset_id, my_store, cache_directory, use_future = FALSE){
    
    if(use_future) my_map = future_map2
    else my_map = map2
    
    my_dataset_id = unique(target_name_grouped_by_dataset_id$dataset_id) 
    
    my_consensus = 
      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_annotation_new.parquet')")
      ) |> 
      filter(dataset_id == my_dataset_id) |> 
      # mutate(cell_ = paste(cell_, dataset_id, sep="___")) |> # ALREAD DONE
      
      # Lighten it up 
      select(cell_, dataset_id, cell_type_unified_ensemble) 
    
    
    # rename(
    #   cell_type_unified_ensemble = cell_type_unified_ensemble, 
    #   cell_type_immune_confidence = confidence_score,
    #   cell_type_immune_independent_annotation_consensus = data_driven_consensus
    # ) 
    
    
    #---------------
    ### WHICH RELEASE
    #---------------
  
    cache_directory_counts = glue("{cache_directory}/original")
    cache_directory_cpm = glue("{cache_directory}/cpm")
    
    dir.create(cache_directory_counts, recursive = TRUE, showWarnings = FALSE)
    dir.create(cache_directory_cpm, recursive = TRUE, showWarnings = FALSE)
    
    # Begin processing the data pipeline with the initial dataset 'target_name_grouped_by_dataset_id'
   sce_df = 
     target_name_grouped_by_dataset_id |> 
      
      # Step 1: Read raw data for each 'target_name' and store it in a new column 'sce'
      mutate(
        sce = map(
          target_name,
          ~ tar_read_raw(.x, store = my_store)  # Read the raw SingleCellExperiment object
        )
      ) |>
      
      # This should not be needed, but there are some data sets with zero cells 
      filter(!map_lgl(sce, is.null))
    
   if(nrow(sce_df) == 0) return(NULL)
   
   # plan(multisession, workers = 20)
   
   sce_df = 
     sce_df |> 
      
      # Step 2: Remove any entries where 'sce' is NULL (i.e., failed to read data)
      filter(!map_lgl(sce, is.null)  ) |>
      
      # Step 3: Calculate CPM (Counts Per Million) for each 'sce' object
      mutate(
        sce = map(
          sce,
          ~ { 
            # Calculate CPM using 'calculateCPM' on the 'X' assay and store it in a new assay 'cpm'
            assay(.x, "cpm") = calculateCPM(.x, assay.type = "X") 
            .x  # Return the modified 'sce' object
          }
        )
      ) |> 
      
      # Step 4: Group the data by 'dataset_id' and 'tar_group' for further summarization
      group_by(dataset_id, tar_group, chunk) |>
      
      # FORCEFULLY drop all but counts and metadata 
      # int_colData(.x) = DataFrame(row.names = colnames(.x))
      # Creates error
      # THIS SHOULD HAVE BEEN DONE IN THE TRANFORM HPCell
      mutate(sce = map(sce, ~  SingleCellExperiment(assay = assays(.x), colData = colData(.x)) )) |> 
      
      # Step 5: Combine all 'sce' objects within each group into a single 'sce' object
      summarise(
        sce = list(   do.call(cbind, args = sce) )
      ) |>
      
      # Step 6: Join additional metadata to each 'sce' object
      mutate(
        sce = map(
          sce,
          ~ .x |> 
            left_join(my_consensus, by =c(".cell" = "cell_", "dataset_id"), copy = TRUE)  # Merge consensus into 'sce'
        )
      ) |>
      
      # Step 7: Split each 'sce' object based on 'cell_type' and create a nested tibble
      mutate(
        sce = map(
          sce,
          ~ .x |> 
            HPCell:::splitColData(colData(.x)$cell_type_unified_ensemble) |>  # Split 'sce' by 'cell_type'
            enframe(name = "cell_type_unified_ensemble", value = "sce")  # Convert to tibble with 'cell_type' and 'sce' columns
        )
      ) |>
      
      # Step 8: Unnest the list of 'sce' objects to have one row per 'cell_type'
      unnest(sce) 
   
   
   
   sce_df = 
     sce_df |>
      
      # Step 9: Generate a unique file identifier for each 'cell_type' within a dataset
      unite("file_id_cellNexus", c(dataset_id, cell_type_unified_ensemble), sep = "___", remove=FALSE) |> 
      mutate(file_id_cellNexus = paste0(file_id_cellNexus, ".h5ad")) |> 
      mutate(file_id_cellNexus = map_chr(file_id_cellNexus, digest) ) |> 
      unite("file_id_cellNexus", c(file_id_cellNexus, chunk), sep = "___chunk") 
   
   
   # process_sce <- function(x, y, cache_directory_counts, my_assay) {
   #   tryCatch({
   #     # Check if the 'sce' has only one cell (column)
   #     if (ncol(assay(x, my_assay)) == 1) {
   #       # Duplicate the assay to prevent saving errors due to single-column matrices
   #       my_assay <- cbind(assay(x, my_assay), assay(x, my_assay))
   #       
   #       # Rename the second column to distinguish it
   #       colnames(my_assay)[2] <- paste0("DUMMY", "___", colnames(my_assay)[2])
   #     } else {
   #       # Use the assay as is if more than one cell
   #       my_assay <- assay(x, my_assay)
   #     }
   #     
   #     # Create a new SingleCellExperiment object with the adjusted counts assay
   #     new_sce <- SingleCellExperiment(assays = list(counts = my_assay))
   #     
   #     # Save the experiment data to the specified counts cache directory
   #     save_experiment_data(new_sce, glue("{cache_directory_counts}/{y}"))
   #     
   #     return(TRUE)  # Indicate successful saving
   #     
   #   }, error = function(e) {
   #     message("Error in process_sce: ", e$message)
   #     return(NA)  # Return NA to indicate an error occurred
   #   })
   # }
   # 
   # # Run `mclapply` in parallel across `y` and `z`
   # results <- mclapply(
   #   seq_along(y),
   #   function(i) process_sce(y[[i]], z[[i]], cache_directory_counts, "X"),
   #   mc.cores = 20
   # )
   # 
   # results <- mclapply(
   #   seq_along(y),
   #   function(i) process_sce(y[[i]], z[[i]], cache_directory_cpm, "cpm"),
   #   mc.cores = 20
   # )
   # 
   # sce_df |> select(-sce)
   
   sce_df|>

      # Step 10: Save the raw counts data for each 'sce' object
      mutate(
        saved_raw_counts = map2(
          sce, file_id_cellNexus,
          ~ {

            # Check if the 'sce' has only one cell (column)
            if(ncol(assay(.x)) == 1) {

              # Duplicate the assay to prevent saving errors due to single-column matrices
              my_assay = cbind(assay(.x), assay(.x))

              # Rename the second column to distinguish it
              colnames(my_assay)[2] = paste0("DUMMY", "___", colnames(my_assay)[2])
            } else {
              # Use the assay as is if more than one cell
              my_assay = assay(.x)
            }

            # Create a new SingleCellExperiment object with the adjusted counts assay
            SingleCellExperiment(assays = list(counts = my_assay)) |>

              # Save the experiment data to the specified counts cache directory
              save_experiment_data(glue("{cache_directory_counts}/{.y}"))

            return(TRUE)  # Indicate successful saving
          }
        )
      ) |>

      # Step 11: Save the CPM data for each 'sce' object
      mutate(
        saved_cpm = map2(
          sce, file_id_cellNexus,
          ~ {

            # Check if the 'cpm' assay has only one cell
            if(ncol(assay(.x, "cpm")) == 1) {

              # Duplicate the 'cpm' assay to avoid single-column issues
              my_assay = cbind(assay(.x, "cpm"), assay(.x, "cpm"))

              # Rename the second column
              colnames(my_assay)[2] = paste0("DUMMY", "___", colnames(my_assay)[2])
            } else {
              # Use the 'cpm' assay as is
              my_assay = assay(.x, "cpm")
            }

            # Create a new SingleCellExperiment object with the adjusted CPM assay
            SingleCellExperiment(assays = list(cpm = my_assay)) |>

              # Save the experiment data to the specified CPM cache directory
              save_experiment_data(glue("{cache_directory_cpm}/{.y}"))

            return(TRUE)  # Indicate successful saving
          }
        )
      ) |>

      # Step 12: Remove the 'sce' column from the final output as it's no longer needed
      select(-sce)
  }
  
  cbind_sce_by_dataset_id = function(target_name_grouped_by_dataset_id, chunk_tbl){

    
    my_dataset_id = unique(target_name_grouped_by_dataset_id$dataset_id) 
    


    # Begin processing the data pipeline with the initial dataset 'target_name_grouped_by_dataset_id'
    sce_df = 
      target_name_grouped_by_dataset_id |> 
      
      # Step 1: Read raw data for each 'target_name' and store it in a new column 'sce'
      mutate(
        sce = map(
          target_name,
          ~ tar_read_raw(.x, store = my_store)  # Read the raw SingleCellExperiment object
        )
      ) |>
      
      # This should not be needed, but there are some data sets with zero cells 
      filter(!map_lgl(sce, is.null)) 
      
    
    if(nrow(sce_df) == 0) return(NULL)
    
    # plan(multisession, workers = 20)
    
    sce_df = 
      sce_df |> 
      
      # # Step 4: Group the data by 'dataset_id' and 'tar_group' for further summarization
      # group_by(dataset_id, tar_group, chunk) |>
      # 
      
      # FORCEFULLY drop all but counts and metadata 
      # int_colData(.x) = DataFrame(row.names = colnames(.x))
      # Creates error
      # THIS SHOULD HAVE BEEN DONE IN THE TRANFORM HPCell
      mutate(sce = map(sce, ~  SingleCellExperiment(assay = assays(.x), colData = colData(.x)) )) |> 
      
      # Step 5: Combine all 'sce' objects within each group into a single 'sce' object
      summarise(
        sce = 
          do.call(cbind, args = sce) |> 
          left_join(chunk_tbl) |> 
          
      ) |>
      
      
      # Step 7: Split each 'sce' object based on 'cell_type' and create a nested tibble
      mutate(
        sce = map(
          sce,
          ~ .x |> 
            HPCell:::splitColData(colData(.x)$cell_type_unified_ensemble) |>  # Split 'sce' by 'cell_type'
            enframe(name = "cell_type_unified_ensemble", value = "sce")  # Convert to tibble with 'cell_type' and 'sce' columns
        )
      ) |>
      
      # Step 8: Unnest the list of 'sce' objects to have one row per 'cell_type'
      unnest(sce) 
    

    
    sce_df|>
      
      # Step 10: Save the raw counts data for each 'sce' object
      mutate(
        saved_raw_counts = map2(
          sce, file_id_cellNexus,
          ~ {
            
            # Check if the 'sce' has only one cell (column)
            if(ncol(assay(.x)) == 1) {
              
              # Duplicate the assay to prevent saving errors due to single-column matrices
              my_assay = cbind(assay(.x), assay(.x))
              
              # Rename the second column to distinguish it
              colnames(my_assay)[2] = paste0("DUMMY", "___", colnames(my_assay)[2])
            } else {
              # Use the assay as is if more than one cell
              my_assay = assay(.x)
            }
            
            # Create a new SingleCellExperiment object with the adjusted counts assay
            sce = SingleCellExperiment(assays = list(counts = my_assay)) 
            
            sce <<<>>> convert assay to HDF5 sparse integers
              
              # Save the experiment data to the specified counts cache directory
              save_experiment_data(glue("{cache_directory_counts}/{.y}"))
            
            return(TRUE)  # Indicate successful saving
          }
        )
      ) |>
      
      # Step 11: Save the CPM data for each 'sce' object
      mutate(
        saved_cpm = map2(
          sce, file_id_cellNexus,
          ~ {
            
            # Check if the 'cpm' assay has only one cell
            if(ncol(assay(.x, "cpm")) == 1) {
              
              # Duplicate the 'cpm' assay to avoid single-column issues
              my_assay = cbind(assay(.x, "cpm"), assay(.x, "cpm"))
              
              # Rename the second column
              colnames(my_assay)[2] = paste0("DUMMY", "___", colnames(my_assay)[2])
            } else {
              # Use the 'cpm' assay as is
              my_assay = assay(.x, "cpm")
            }
            
            # Create a new SingleCellExperiment object with the adjusted CPM assay
            SingleCellExperiment(assays = list(cpm = my_assay)) |>
              
              # Save the experiment data to the specified CPM cache directory
              save_experiment_data(glue("{cache_directory_cpm}/{.y}"))
            
            return(TRUE)  # Indicate successful saving
          }
        )
      ) |>
      
      # Step 12: Remove the 'sce' column from the final output as it's no longer needed
      select(-sce)
  }
  
  
  get_dataset_id = function(target_name, my_store){
    sce = tar_read_raw(target_name, store = my_store)
    if(sce |> is.null()) return("sce_0_cells")
    sce |> pull(dataset_id) |> unique()
  }
  
  list(
    
    # The input DO NOT DELETE
    tar_target(my_store, "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/census_hpcell_oct_2024/target_store", deployment = "main"),
    #tar_target(cache_directory, "/vast/projects/cellxgene_curated/cellNexus/cellxgene/29_10_2024", deployment = "main"),
    tar_target(cache_directory, "/vast/scratch/users/mangiola.s/cellNexus/cellxgene/11_11_2024", deployment = "main"),
    
    tar_target(
      target_name,
      tar_meta(
        starts_with("sce_transformed_tier_"), 
        store = my_store) |> 
        filter(type=="branch") |> 
        pull(name),
      deployment = "main"
    ),
    tar_target(
      dataset_id,
      get_dataset_id(target_name, my_store),
      packages = "tidySingleCellExperiment",
      pattern = map(target_name),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "tier_1")
      )
    ),
    tar_target(
      chunk_tbl,
      {
        cell_annotation = 
          tbl(
            dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
            sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_annotation_new.parquet')")
          ) |>
          select(cell_, dataset_id, cell_type_unified_ensemble) 
        
        
          tbl(
            dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
            sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_annotation.parquet')")
          ) |> 
          left_join(cell_annotation, copy=TRUE) |>
          
          # Define chunks
          group_by(dataset_id, cell_type_unified_ensemble, sample_id) |>
          summarise(cell_count = n(), .groups = "drop") |>
          group_by(dataset_id, cell_type_unified_ensemble) |>
          dbplyr::window_order(desc(cell_count)) |>
          mutate(chunk = cumsum(cell_count) %/% 20000 + 1) |> # max 20K cells per sample
          ungroup() |> 
            as_tibble()
      }, deployment = "main", packages = c("tidyverse", "dbplyr", "arrow", "dplyr", "duckdb")
    ),
    
    tar_target(
      sce_dataset_id,
      chunk_tbl
    )
    

    tar_target(
      target_name_grouped_by_dataset_id,
      tibble(target_name = target_name, dataset_id = dataset_id) |>
        group_by(dataset_id) |>
        tar_group(),
      iteration = "group",
      resources = tar_resources(
        crew = tar_resources_crew(controller = "tier_3")
      )
    ),
    tar_target(
      dataset_id_sce,
      save_anndata(target_name_grouped_by_dataset_id, my_store, cache_directory),
      pattern = map(target_name_grouped_by_dataset_id),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "HPCell", "digest", "scater", "arrow", "dplyr", "duckdb", "furrr", "future", "parallel"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "tier_4")
      )
    )
  )
  
  
}, script = paste0(store_file_cellNexus, "_target_script.R"), ask = FALSE)

job::job({
  
  store_file_cellNexus = "/vast/scratch/users/mangiola.s/targets_prepare_database_split_datasets_chunked"
  tar_make(script = paste0(store_file_cellNexus, "_target_script.R"), store = store_file_cellNexus, 
           reporter = "summary")
  
})

tar_make(script = paste0(store_file_cellNexus, "_target_script.R"), store = store_file_cellNexus, callr_function = NULL)

tar_workspaces(store = store_file_cellNexus)

tar_read_raw(name = "dataset_id_sce_0b71a2a4781c9418", store = store_file_cellNexus)

tar_meta(store = store_file_cellNexus) |> 
  arrange(desc(time)) |>
  filter(!error |> is.na()) |> 
  select(name, error)

tar_workspace("dataset_id_sce_b0a3540d0d2d257d", store = store_file_cellNexus, script = paste0(store_file_cellNexus, "_target_script.R"))

# tar_invalidate(dataset_id_sce, store = store_file_cellNexus)

# # PREPARE METADATA FINAL
# file_id_tbl = 
#   tbl(
#   dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
#   sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata_cell_type_consensus_v1_0_0.parquet')")
# ) |> 
#   distinct(dataset_id, cell_type_unified_ensemble) |> 
#   as_tibble() |> 
#   unite("file_id_cellNexus", c(dataset_id, cell_type_unified_ensemble), sep = "___", remove=FALSE) |> 
#   mutate(file_id_cellNexus = map_chr(file_id_cellNexus, digest::digest) )  
#   
  

library(arrow)
library(dplyr)
library(duckdb)


# Get Dharmesh metadata consensus
system("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/reannotation_consensus/cell_annotation_new.parquet /vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/")




job::job({
  
  get_file_ids = function(cell_annotation, cell_type_consensus_parquet){
    
    cell_consensus = 
      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue::glue("SELECT * FROM read_parquet('{cell_type_consensus_parquet}')"))
      ) |>
      select(cell_, dataset_id, cell_type_unified_ensemble, cell_type_unified) 
    
    
    tbl(
      dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
      sql(glue::glue("SELECT * FROM read_parquet('{cell_annotation}')"))
    ) |> 
      left_join(cell_consensus, copy=TRUE) |>
      
      # # Make sure I cover cell type even if consensus of harmonisation is not present (it should be the vast minority)
      # mutate(temp_cell_type_label_for_file_id = if_else(cell_type_unified_ensemble |> is.na(), cell_type, cell_type_unified)) |> 
      # mutate(temp_cell_type_label_for_file_id = if_else(temp_cell_type_label_for_file_id |> is.na(), cell_type, temp_cell_type_label_for_file_id)) |> 
      # 
      # Define chunks
      group_by(dataset_id, cell_type_unified, sample_id) |>
      summarise(cell_count = n(), .groups = "drop") |>
      group_by(dataset_id, cell_type_unified) |>
      dbplyr::window_order(desc(cell_count)) |>
      mutate(chunk = cumsum(cell_count) %/% 20000 + 1) |> # max 20K cells per sample
      ungroup() |> 
      as_tibble() |> 
      
      mutate(file_id_cellNexus_single_cell = 
               glue::glue("{dataset_id}___{cell_type_unified}") |> 
               sapply(digest::digest) |> 
               paste0("___", chunk)
      )
  }
  
  
  get_file_ids(
    "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_annotation.parquet",
    "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_annotation_new.parquet"
  )  |> 
    write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/file_id_cellNexus_single_cell.parquet")
  
  gc()
   
  con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  # Create a view for cell_annotation in DuckDB
  dbExecute(con, "
  CREATE VIEW cell_metadata AS
  SELECT 
    CONCAT(cell_, '___', dataset_id) AS cell_,
    dataset_id,
    *
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata.parquet')
")
  
  # Create views for other tables
  dbExecute(con, "
  CREATE VIEW cell_consensus AS
  SELECT *
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_annotation_new.parquet')
")
  
  dbExecute(con, "
  CREATE VIEW cell_annotation AS
  SELECT cell_, blueprint_first_labels_fine, monaco_first_labels_fine, azimuth_predicted_celltype_l2
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/annotation_tbl_light.parquet')
")
  
  dbExecute(con, "
  CREATE VIEW empty_droplet_df AS
  SELECT cell_, dataset_id, empty_droplet
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_annotation.parquet')
")
  
  dbExecute(con, "
  CREATE VIEW file_id_cellNexus_single_cell AS
  SELECT dataset_id, cell_type_unified, sample_id, file_id_cellNexus_single_cell
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/file_id_cellNexus_single_cell.parquet')
")
  
  # Perform the left join and save to Parquet
  copy_query <- "
  COPY (
    SELECT *
    FROM cell_metadata

  LEFT JOIN cell_annotation
    ON cell_annotation.cell_ = cell_metadata.cell_

  LEFT JOIN cell_consensus
    ON cell_consensus.cell_ = cell_metadata.cell_
    AND cell_consensus.dataset_id = cell_metadata.dataset_id
    
  LEFT JOIN empty_droplet_df
    ON empty_droplet_df.cell_ = cell_metadata.cell_
    AND empty_droplet_df.dataset_id = cell_metadata.dataset_id

  LEFT JOIN file_id_cellNexus_single_cell
    ON file_id_cellNexus_single_cell.sample_id = cell_consensus.sample_id
    AND file_id_cellNexus_single_cell.dataset_id = cell_consensus.dataset_id
    AND file_id_cellNexus_single_cell.cell_type_unified = cell_consensus.cell_type_unified
      
  ) TO '/vast/projects/cellxgene_curated/cellNexus/cell_metadata_cell_type_consensus_v1_0_1.parquet'
  (FORMAT PARQUET, COMPRESSION 'gzip');
"
  
  # Execute the final query to write the result to a Parquet file
  dbExecute(con, copy_query)
  
  # Disconnect from the database
  dbDisconnect(con, shutdown = TRUE)
  
  system("~/bin/rclone copy /vast/projects/cellxgene_curated/cellNexus/cell_metadata_cell_type_consensus_v1_0_1.parquet box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")
  
  
#   #----------------
#   # Add File ID
#   #----------------
#   
#   file_id_tbl = 
#     tbl(
#       dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
#       sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata_cell_type_consensus_TEMP.parquet')")
#     ) |> 
#     distinct(dataset_id, cell_type_unified_ensemble) |> 
#     as_tibble() |> 
#     unite("file_id_cellNexus", c(dataset_id, cell_type_unified_ensemble), sep = "___", remove=FALSE) |> 
#     mutate(file_id_cellNexus = map_chr(file_id_cellNexus, digest::digest) )  |> 
#     write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata_file_id_TEMP.parquet")
#   
#   
#   con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
#   
#   dbExecute(con, "
#   CREATE VIEW cell_metadata_cell_type_consensus AS
#   SELECT *
#   FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata_cell_type_consensus_TEMP.parquet')
# ")
#   
#   dbExecute(con, "
#   CREATE VIEW cell_metadata_file_id AS
#   SELECT *
#   FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata_file_id_TEMP.parquet')
# ")
#   
#   # Perform the join and save the result to Parquet
#   join_query <- "
#   COPY (
#     SELECT *
#     FROM cell_metadata_cell_type_consensus
# 
#     LEFT JOIN cell_metadata_file_id
#       ON cell_metadata_cell_type_consensus.dataset_id = cell_metadata_file_id.dataset_id
#       AND cell_metadata_cell_type_consensus.cell_type_unified_ensemble = cell_metadata_file_id.cell_type_unified_ensemble
# 
#   ) TO '/vast/projects/cellxgene_curated/cellNexus/cell_metadata_cell_type_consensus_v1_0_0.parquet'
#   (FORMAT PARQUET, COMPRESSION 'gzip');
# "
#   
#   # Execute the join query
#   dbExecute(con, join_query)
#   
#   # Disconnect from the database
#   dbDisconnect(con, shutdown = TRUE)
#   
#   file.remove("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata_file_id_TEMP.parquet")
#   file.remove("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata_cell_type_consensus_TEMP.parquet")

  })



# process_h5ad_files("/vast/projects/cellxgene_curated/cellxgene/", "/vast/projects/cellxgene_curated/cellNexus/cellxgene")

##############
#  PLOTS     #
##############


cell_metadata = 
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/cellNexus/cell_metadata_cell_type_consensus_v1_0_1.parquet')")
  ) 
  
cell_metadata_v1_0_0 |> 
  filter(age_days > 365) |> 
  distinct(donor_id, tissue_groups) |> 
  ggplot(aes(fct_infreq(tissue_groups))) +
  geom_bar() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

cell_metadata_v1_0_0 |> 
  filter(age_days > 365) |> 
  distinct(sample_id, donor_id, tissue_groups) |> 
  as_tibble() |> 
  nest(data = -tissue_groups) |> 
  mutate(n_sample = map_int(data, ~ .x |> distinct(sample_id) |> nrow())) |> 
  mutate(n_donor = map_int(data, ~ .x |> distinct(donor_id) |> nrow())) |> 
  ggplot(aes(n_sample, n_donor)) +
  geom_point() +
  geom_text(aes(label = tissue_groups)) +
  scale_y_log10() +
  scale_x_log10() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
