
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
    
    # This because f7c1c579-2dc0-47e2-ba19-8165c5a0e353 includes 13K samples
    # It affects only very few datasets
    sample_chunk_df = 
      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue::glue("SELECT * FROM read_parquet('{cell_annotation}')"))
      ) |> 
      # Define chunks
      count(dataset_id, sample_id, name = "cell_count") |>  # Ensure unique dataset_id and sample_id combinations
      distinct(dataset_id, sample_id, cell_count) |>  # Ensure unique dataset_id and sample_id combinations
      group_by(dataset_id) |> 
      dbplyr::window_order(cell_count) |>  # Ensure dataset_id order
      mutate(sample_index = row_number()) |>  # Create sequential index within each dataset
      mutate(sample_chunk = (sample_index - 1) %/% 100 + 1) |>  # Assign chunks (up to 1000 samples per chunk)
      mutate(cell_chunk = cumsum(cell_count) %/% 20000 + 1) |> # max 20K cells per sample
      ungroup() 
    
    tbl(
      dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
      sql(glue::glue("SELECT * FROM read_parquet('{cell_annotation}')"))
    ) |> 
      left_join(cell_consensus, copy=TRUE) |>
      left_join(sample_chunk_df |> select(dataset_id, sample_chunk, cell_chunk, sample_id), copy=TRUE) |> 
      # # Make sure I cover cell type even if consensus of harmonisation is not present (it should be the vast minority)
      # mutate(temp_cell_type_label_for_file_id = if_else(cell_type_unified_ensemble |> is.na(), cell_type, cell_type_unified)) |> 
      # mutate(temp_cell_type_label_for_file_id = if_else(temp_cell_type_label_for_file_id |> is.na(), cell_type, temp_cell_type_label_for_file_id)) |> 
      # 
      # Define chunks
      group_by(dataset_id, sample_chunk, cell_chunk, cell_type_unified, sample_id) |>
      summarise(cell_count = n(), .groups = "drop") |>
      group_by(dataset_id, sample_chunk, cell_chunk, cell_type_unified) |>
      dbplyr::window_order(desc(cell_count)) |>
      mutate(chunk = cumsum(cell_count) %/% 20000 + 1) |> # max 20K cells per sample
      ungroup() |> 
      as_tibble() |> 
      
      mutate(file_id_cellNexus_single_cell = 
               glue::glue("{dataset_id}___{sample_chunk}___{cell_chunk}___{cell_type_unified}") |> 
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
  SELECT dataset_id, sample_chunk, cell_chunk, cell_type_unified, sample_id, file_id_cellNexus_single_cell
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/file_id_cellNexus_single_cell.parquet')
")
  
  # Perform the left join and save to Parquet
  copy_query <- "
  COPY (
     SELECT 
        cell_metadata.cell_ AS cell_id, -- Rename cell_ to cell_id
        cell_metadata.*,              -- Include all other columns from cell_metadata
        cell_annotation.*,            -- Include all columns from cell_annotation
        cell_consensus.*,             -- Include all columns from cell_consensus
        empty_droplet_df.*,           -- Include all columns from empty_droplet_df
        file_id_cellNexus_single_cell.* -- Include all columns from file_id_cellNexus_single_cell
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
      
  ) TO '/vast/projects/cellxgene_curated/cellNexus/cell_metadata_cell_type_consensus_v1_0_3.parquet'
  (FORMAT PARQUET, COMPRESSION 'gzip');
"
  
  # Execute the final query to write the result to a Parquet file
  dbExecute(con, copy_query)
  
  # Disconnect from the database
  dbDisconnect(con, shutdown = TRUE)
  
  system("~/bin/rclone copy /vast/projects/cellxgene_curated/cellNexus/cell_metadata_cell_type_consensus_v1_0_3.parquet box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")
  
})


cell_metadata = 
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/cellNexus/cell_metadata_cell_type_consensus_v1_0_3.parquet')")
  ) 


library(targets)
library(tidyverse)
store_file_cellNexus = "/vast/scratch/users/mangiola.s/targets_prepare_database_split_datasets_chunked_3"

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
    # debug = "target_name_grouped_by_dataset_id", 
    cue = tar_cue(mode = "never"), 
    workspace_on_error = TRUE,
    controller = crew_controller_group(
      list(
        crew_controller_slurm(
          name = "tier_1", 
          script_lines = "#SBATCH --mem 8G",
          slurm_cpus_per_task = 1, 
          workers = 300, 
          tasks_max = 50,
          verbose = T,
          crashes_error = 5, 
          seconds_idle = 30
        ),
        
        crew_controller_slurm(
          name = "tier_2",
          script_lines = "#SBATCH --mem 10G",
          slurm_cpus_per_task = 1,
          workers = 300,
          tasks_max = 10,
          verbose = T,
          crashes_error = 5, 
          seconds_idle = 30
        ),
        crew_controller_slurm(
          name = "tier_3",
          script_lines = "#SBATCH --mem 30G",
          slurm_cpus_per_task = 1,
          workers = 200,
          tasks_max = 10,
          verbose = T,
          crashes_error = 5, 
          seconds_idle = 30
        ),
        crew_controller_slurm(
          name = "tier_4",
          script_lines = "#SBATCH --mem 150G",
          slurm_cpus_per_task = 20,
          workers = 50,
          tasks_max = 10,
          verbose = T,
          crashes_error = 5, 
          seconds_idle = 30
        ),
        crew_controller_slurm(
          name = "tier_5",
          script_lines = "#SBATCH --mem 400G",
          slurm_cpus_per_task = 1,
          workers = 2,
          tasks_max = 10,
          verbose = T,
          crashes_error = 5, 
          seconds_idle = 30
        )
      )
    ), 
    trust_object_timestamps = TRUE, 
    workspaces = "dataset_id_sce_52dbec3c15f98d66"
  )
  
  
  save_anndata = function(dataset_id_sce, cache_directory, use_future = FALSE){
    
    dir.create(cache_directory, showWarnings = FALSE, recursive = TRUE)
    
    # Parallelise
    cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
    bp <- MulticoreParam(workers = cores , progressbar = TRUE)  # Adjust the number of workers as needed
    
    dataset_id_sce |> 
      purrr::transpose() |> 
      bplapply(
        FUN = function(x) {
          
          .x = x[[2]]
          .y = x[[1]]
          
          # Check if the 'sce' has only one cell (column)
          if(ncol(assay(.x)) == 1) {
            
            # Duplicate the assay to prevent saving errors due to single-column matrices
            my_assay = cbind(assay(.x), assay(.x))
            # Rename the second column to distinguish it
            colnames(my_assay)[2] = paste0("DUMMY", "___", colnames(my_assay)[2])
            
            cd = colData(.x)
            cd = cd |> rbind(cd)
            rownames(cd)[2] = paste0("DUMMY", "___", rownames(cd)[2])
            
           
            
           .x =  SingleCellExperiment(assay = list( my_assay ) |> set_names(names(assays(.x))[1]), colData = cd) 
          } 
          
          
          # TEMPORARY FOR SOME REASON THE MIN COUNTS IS NOT 0 FOR SOME SAMPLES
          .x = HPCell:::check_if_assay_minimum_count_is_zero_and_correct_TEMPORARY(.x, assays(.x) |> names() |> _[1])
          
          
          # My attempt to save a integer, sparse, delayed matrix (with zellkonverter it is not possible to save integers)
          # .x |> assay() |> type() <- "integer"
          # .x |> saveHDF5SummarizedExperiment("~/temp", as.sparse = T, replace = T)
          
          # Save the experiment data to the specified counts cache directory
          .x |> save_experiment_data(glue("{cache_directory}/{.y}"))
          
          return(TRUE)  # Indicate successful saving
        },
        BPPARAM = bp  # Use the defined parallel backend
      )
    
    return("saved")
    
    # dataset_id_sce|>
    #   
    #   # Step 10: Save the raw counts data for each 'sce' object
    #   mutate(
    #     saved_raw_counts = map2(
    #       sce, cell_type_unified_ensemble,
    #       ~ {
    #         
    #         # Check if the 'sce' has only one cell (column)
    #         if(ncol(assay(.x)) == 1) {
    #           
    #           # Duplicate the assay to prevent saving errors due to single-column matrices
    #           my_assay = cbind(assay(.x), assay(.x))
    #           
    #           # Rename the second column to distinguish it
    #           colnames(my_assay)[2] = paste0("DUMMY", "___", colnames(my_assay)[2])
    #         } else {
    #           # Use the assay as is if more than one cell
    #           my_assay = assay(.x)
    #         }
    #         
    #         # TEMPORARY FOR SOME REASON THE MIN COUNTS IS NOT 0 FOR SOME SAMPLES
    #         .x = HPCell:::check_if_assay_minimum_count_is_zero_and_correct_TEMPORARY(.x, assays(.x) |> names() |> _[1])
    #         
    #         
    #         # My attempt to save a integer, sparse, delayed matrix (with zellkonverter it is not possible to save integers)
    #         # .x |> assay() |> type() <- "integer"
    #         # .x |> saveHDF5SummarizedExperiment("~/temp", as.sparse = T, replace = T)
    #         
    #         # Save the experiment data to the specified counts cache directory
    #         .x |> save_experiment_data(glue("{cache_directory_counts}/{.y}"))
    #         
    #         return(TRUE)  # Indicate successful saving
    #       }
    #     )
    #   ) |>
    #   
    #   # Step 11: Save the CPM data for each 'sce' object
    #   mutate(
    #     saved_cpm = map2(
    #       sce, file_id_cellNexus,
    #       ~ {
    #         
    #         # Check if the 'cpm' assay has only one cell
    #         if(ncol(assay(.x, "cpm")) == 1) {
    #           
    #           # Duplicate the 'cpm' assay to avoid single-column issues
    #           my_assay = cbind(assay(.x, "cpm"), assay(.x, "cpm"))
    #           
    #           # Rename the second column
    #           colnames(my_assay)[2] = paste0("DUMMY", "___", colnames(my_assay)[2])
    #         } else {
    #           # Use the 'cpm' assay as is
    #           my_assay = assay(.x, "cpm")
    #         }
    #         
    #         # Create a new SingleCellExperiment object with the adjusted CPM assay
    #         SingleCellExperiment(assays = list(cpm = my_assay)) |>
    #           
    #           # Save the experiment data to the specified CPM cache directory
    #           save_experiment_data(glue("{cache_directory_cpm}/{.y}"))
    #         
    #         return(TRUE)  # Indicate successful saving
    #       }
    #     )
    #   ) |>
    #   
    #   # Step 12: Remove the 'sce' column from the final output as it's no longer needed
    #   select(-sce)
  }
  
  cbind_sce_by_dataset_id = function(target_name_grouped_by_dataset_id, file_id_db_file, my_store){
    
    my_dataset_id = unique(target_name_grouped_by_dataset_id$dataset_id) 
    
    file_id_db = 
      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue("SELECT * FROM read_parquet('{file_id_db_file}')"))
      ) |> 
      filter(dataset_id == my_dataset_id) |>
      select(cell_id, sample_id, dataset_id, file_id_cellNexus_single_cell) |> 
      as_tibble()
    
    # Parallelise
    cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
    bp <- MulticoreParam(workers = cores , progressbar = TRUE)  # Adjust the number of workers as needed
    
    # Begin processing the data pipeline with the initial dataset 'target_name_grouped_by_dataset_id'
    sce_df = 
      target_name_grouped_by_dataset_id |> 
      
      # Step 1: Read raw data for each 'target_name' and store it in a new column 'sce'
      mutate(
        sce = bplapply(
          target_name,
          FUN = function(x) tar_read_raw(x, store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/census_hpcell_oct_2024/target_store"),  # Read the raw SingleCellExperiment object
          BPPARAM = bp  # Use the defined parallel backend
        )
      ) |>
      
      # This should not be needed, but there are some data sets with zero cells 
      filter(!map_lgl(sce, is.null)) 
    
  
    
    if(nrow(sce_df) == 0) {
      warning("this chunk has no rows for somereason.")
      return(NULL)
    }
    
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
      summarise( sce =  list(do.call(cbind, args = sce) ) ) |>
      
      mutate(sce = map(sce,
                       ~ { .x = 
                         .x  |> 
                         left_join(file_id_db, by = join_by(.cell==cell_id, dataset_id==dataset_id, sample_id==sample_id)) 
                       .x |> 
                         HPCell:::splitColData(colData(.x)$file_id_cellNexus_single_cell) |>  # Split 'sce' by 'cell_type'
                         enframe(name = "file_id_cellNexus_single_cell", value = "sce")  # Convert to tibble with 'cell_type' and 'sce' columns
                       })) |> 
      # Step 8: Unnest the list of 'sce' objects to have one row per 'cell_type'
      unnest(sce) 
    
    
  }
  
  get_dataset_id = function(target_name, my_store){
    sce = tar_read_raw(target_name, store = my_store)
    
    if(sce |> is.null()) return(tibble(sample_id = character(), dataset_id= character(), target_name= character()))
    
    sce |> distinct(sample_id, dataset_id) |> mutate(target_name = !!target_name)
  }
  
  create_chunks_for_reading_and_saving = function(dataset_id_sample_id, cell_metadata){
    dataset_id_sample_id |> 
      left_join(
        tbl(
          dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
          sql(glue("SELECT * FROM read_parquet('{cell_metadata}')"))
        )   |> 
          distinct(dataset_id, sample_id, sample_chunk, cell_chunk) |> 
          as_tibble(), 
        copy=T
      )
  }
  
  list(
    
    # The input DO NOT DELETE
    tar_target(my_store, "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/census_hpcell_oct_2024/target_store", deployment = "main"),
    tar_target(cache_directory, "/vast/scratch/users/mangiola.s/cellNexus/cellxgene/19_11_2024", deployment = "main"),
    tar_target(
      cell_metadata,
      "/vast/projects/cellxgene_curated/cellNexus/cell_metadata_cell_type_consensus_v1_0_3.parquet", 
      packages = c( "arrow","dplyr","duckdb")
      
    ),
    
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
      dataset_id_sample_id,
      get_dataset_id(target_name, my_store),
      packages = "tidySingleCellExperiment",
      pattern = map(target_name),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "tier_1")
      )
    ),

    
    tar_target(
      target_name_grouped_by_dataset_id,
       create_chunks_for_reading_and_saving(dataset_id_sample_id, cell_metadata) |> 
        group_by(dataset_id, sample_chunk, cell_chunk) |>
        tar_group(),
      iteration = "group",
      resources = tar_resources(
        crew = tar_resources_crew(controller = "tier_3")
      ), 
      packages = c("arrow", "duckdb", "dplyr", "glue")
      
    ),
    
    tar_target(
      dataset_id_sce,
      cbind_sce_by_dataset_id(target_name_grouped_by_dataset_id, cell_metadata, my_store = my_store),
      pattern = map(target_name_grouped_by_dataset_id),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "HPCell", "digest", "scater", "arrow", "dplyr", "duckdb",  "BiocParallel", "parallelly"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "tier_3")
      )
    ),
  

    tar_target(
      saved_dataset,
      save_anndata(dataset_id_sce, paste0(cache_directory, "single_cell/counts")),
      pattern = map(dataset_id_sce),
      packages = c("tidySingleCellExperiment", "SingleCellExperiment", "tidyverse", "glue", "digest", "HPCell", "digest", "scater", "arrow", "dplyr", "duckdb", "BiocParallel", "parallelly"),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "tier_4")
      )
    )
  )
  
  
  
}, script = paste0(store_file_cellNexus, "_target_script.R"), ask = FALSE)

job::job({
  
  tar_make(
    script = paste0(store_file_cellNexus, "_target_script.R"), 
    store = store_file_cellNexus, 
    reporter = "summary" #, callr_function = NULL
  )
  
})

tar_make(script = paste0(store_file_cellNexus, "_target_script.R"), store = store_file_cellNexus, callr_function = NULL)

tar_workspace(target_name_grouped_by_dataset_id, store = store_file_cellNexus, script = paste0(store_file_cellNexus, "_target_script.R"))

job::job({ dataset_id_sce = tar_read_raw("dataset_id_sce_ed698af04de55aef", store = store_file_cellNexus) })

tar_meta(store = store_file_cellNexus) |> 
  arrange(desc(time)) |>
  filter(!error |> is.na()) |> 
  select(name, error)

tar_workspace("saved_dataset_ed3c60d5343228ee", store = store_file_cellNexus, script = paste0(store_file_cellNexus, "_target_script.R"))

# tar_invalidate(target_name_grouped_by_dataset_id, store = store_file_cellNexus)
# tar_delete(dataset_id_sce, , store = store_file_cellNexus)

tar_read(dataset_id_sce, store = store_file_cellNexus)

