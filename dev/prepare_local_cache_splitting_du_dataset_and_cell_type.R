
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

store_file_cellNexus = "/vast/scratch/users/mangiola.s/targets_prepare_database_split_datasets"

# Get Dharmesh metadata consensus
system("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/reannotation_consensus/data_driven_consensus_new.parquet /vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/ ")

# 

tar_script({
  library(dplyr)
  library(magrittr)
  library(tibble)
  library(targets)
  library(tarchetypes)
  library(crew)
  library(crew.cluster)
  tar_option_set(memory = "transient", garbage_collection = 100, 
                 storage = "worker", retrieval = "worker", error = "continue", 
                 #  debug = "dataset_id_sce_b5312463451d7ee3", 
                 cue = tar_cue(mode = "never"), 
                 controller = crew_controller_group(
                   list(
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
                       workers = 300,
                       tasks_max = 10,
                       verbose = T,
                       launch_max = 5
                     ),
                     crew_controller_slurm(
                       name = "tier_3",
                       script_lines = "#SBATCH --mem 15G",
                       slurm_cpus_per_task = 1,
                       workers = 300,
                       tasks_max = 10,
                       verbose = T,
                       launch_max = 5
                     ),
                     crew_controller_slurm(
                       name = "tier_4",
                       script_lines = "#SBATCH --mem 50G",
                       slurm_cpus_per_task = 1,
                       workers = 30,
                       tasks_max = 10,
                       verbose = T,
                       launch_max = 5
                     )
                   )
                 ), 
                 trust_object_timestamps = TRUE
  )
  
  
  save_anndata = function(target_name_grouped_by_dataset_id, my_store){
    
    my_dataset_id = unique(target_name_grouped_by_dataset_id$dataset_id) 
    
    my_consensus = 
      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/consensus_output_new_plus_non_immune_harmonisation.parquet')")
      ) |> 
      filter(dataset_id == my_dataset_id) |> 
      mutate(cell_ = paste(cell_, dataset_id, sep="___")) |> 
      
      # Lighten it up 
      select(cell_, dataset_id, cell_type_immune_consensus = reannotation_consensus) 
    
    
    # rename(
    #   cell_type_immune_consensus = reannotation_consensus, 
    #   cell_type_immune_confidence = confidence_score,
    #   cell_type_immune_independent_annotation_consensus = data_driven_consensus
    # ) 
    
    
    #################
    ### WHICH RELEASE
    #################
    
    release_date = "29_10_2024"
    
    cache_directory = glue("/vast/projects/cellxgene_curated/cellxgene/{release_date}")
    cache_directory_counts = glue("{cache_directory}/original")
    cache_directory_cpm = glue("{cache_directory}/cpm")
    
    dir.create(cache_directory_counts, recursive = TRUE, showWarnings = FALSE)
    dir.create(cache_directory_cpm, recursive = TRUE, showWarnings = FALSE)
    
    # Begin processing the data pipeline with the initial dataset 'target_name_grouped_by_dataset_id'
    target_name_grouped_by_dataset_id |> 
      
      # Step 1: Read raw data for each 'target_name' and store it in a new column 'sce'
      mutate(
        sce = map(
          target_name,
          ~ tar_read_raw(.x, store = my_store)  # Read the raw SingleCellExperiment object
        )
      ) |>
      
      # Step 2: Remove any entries where 'sce' is NULL (i.e., failed to read data)
      filter(
        !map_lgl(sce, is.null)  # Keep only the rows where 'sce' is not NULL
      ) |>
      
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
      group_by(dataset_id, tar_group) |>
      
      # Step 5: Combine all 'sce' objects within each group into a single 'sce' object
      summarise(
        sce = list(
          do.call(cbind, args = sce)  # Column-bind the list of 'sce' objects
        )
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
            HPCell:::splitColData(colData(.x)$cell_type_immune_consensus) |>  # Split 'sce' by 'cell_type'
            enframe(name = "cell_type_immune_consensus", value = "sce")  # Convert to tibble with 'cell_type' and 'sce' columns
        )
      ) |>
      
      # Step 8: Unnest the list of 'sce' objects to have one row per 'cell_type'
      unnest(sce) |>
      
      # Step 9: Generate a unique file identifier for each 'cell_type' within a dataset
      unite("file_id_cellNexus", c(dataset_id, cell_type_immune_consensus), sep = "___", remove=FALSE) |> 
      mutate(file_id_cellNexus = map_chr(file_id_cellNexus, digest) )  |>
      
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
  
  get_dataset_id = function(target_name, my_store){
    sce = tar_read_raw(target_name, store = my_store)
    if(sce |> is.null()) return("sce_0_cells")
    sce |> pull(dataset_id) |> unique()
  }
  
  list(
    
    # The input DO NOT DELETE
    tar_target(my_store, "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/census_hpcell_oct_2024/target_store", deployment = "main"),
    
    tar_target(
      target_name,
      tar_meta(
        starts_with("sce_transformed_tier_"), 
        store = my_store) |> 
        filter(type=="branch") |> 
        pull(name),
      deployment = "main"
    )    ,
    
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
      target_name_grouped_by_dataset_id,
      tibble(target_name = target_name, dataset_id = dataset_id) |>
        group_by(dataset_id) |>
        tar_group(),
      iteration = "group",
      resources = tar_resources(
        crew = tar_resources_crew(controller = "tier_4")
      )
    ),
    tar_target(
      dataset_id_sce,
      save_anndata(target_name_grouped_by_dataset_id, my_store),
      pattern = map(target_name_grouped_by_dataset_id),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "HPCell", "digest", "scater", "arrow", "dplyr", "duckdb"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "tier_3")
      )
    )
  )
  
  
}, script = paste0(store_file_cellNexus, "_target_script.R"), ask = FALSE)

job::job({
  store_file_cellNexus = "/vast/scratch/users/mangiola.s/targets_prepare_database_split_datasets"
  tar_make(script = paste0(store_file_cellNexus, "_target_script.R"), store = store_file_cellNexus, reporter = "verbose_positives")
  
})

tar_make(script = paste0(store_file_cellNexus, "_target_script.R"), store = store_file_cellNexus, callr_function = NULL)

tar_read(dataset_id_sce, store = store_file_cellNexus)

tar_read_raw(name = "dataset_id_2c94b61d138014f6", store = store_file_cellNexus)

tar_meta(starts_with("dataset_id_sce_"), store = store_file_cellNexus)

tar_invalidate(dataset_id, store = store_file_cellNexus)
