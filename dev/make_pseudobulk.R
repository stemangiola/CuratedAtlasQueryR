library(tidyverse)
library(targets)
library(glue)

library(arrow)
library(dplyr)
library(duckdb)
# Get input


result_directory = "/vast/scratch/users/mangiola.s/pseudobulk_cellNexus_1_0_0"

tar_script({
    
    #-----------------------#
    # Input
    #-----------------------#
    library(tidyverse)
    library(targets)
    library(tarchetypes)
    library(glue)
    library(qs)
    library(crew)
    library(crew.cluster)
    
    #-----------------------#
    # Packages
    #-----------------------#
    tar_option_set(
      packages = c(
        "CuratedAtlasQueryR", "stringr", "tibble", "tidySingleCellExperiment", "dplyr", "Matrix",
        "Seurat", "glue", "qs",  "purrr", "tidybulk", "tidySummarizedExperiment", "edgeR", "crew", "magrittr", "digest", "glmmSeq", "readr", "forcats"
      ),
      
      memory = "transient", 
      garbage_collection = 100, 
      storage = "worker", 
      retrieval = "worker", 
      error = "continue", 
      #  debug = "dataset_id_sce_ce393fc1e85f2cbc", 
      cue = tar_cue(mode = "never"), 
      workspace_on_error = TRUE,
      format = "qs",

      #-----------------------#
      # SLURM
      #-----------------------#
      controller = crew_controller_group(
        
        
        crew_controller_slurm(
          name = "slurm_1_5",
          slurm_memory_gigabytes_per_cpu = 5,
          slurm_cpus_per_task = 1,
          workers = 400,
          tasks_max = 20,
          verbose = T
        ),
        crew_controller_slurm(
          name = "slurm_1_10",
          slurm_memory_gigabytes_per_cpu = 10,
          slurm_cpus_per_task = 1,
          workers = 300,
          tasks_max = 20,
          verbose = T
        ),
        crew_controller_slurm(
          name = "slurm_1_20",
          slurm_memory_gigabytes_per_cpu = 20,
          slurm_cpus_per_task = 1,
          workers = 200,
          tasks_max = 20,
          verbose = T
        ),
        crew_controller_slurm(
          name = "slurm_1_30",
          slurm_memory_gigabytes_per_cpu = 30,
          slurm_cpus_per_task = 1,
          workers = 200,
          tasks_max = 10,
          verbose = T
        ),
        crew_controller_slurm(
          name = "slurm_1_40",
          slurm_memory_gigabytes_per_cpu = 40,
          slurm_cpus_per_task = 1,
          workers = 100,
          tasks_max = 10,
          verbose = T
        ),
        crew_controller_slurm(
          name = "slurm_1_80",
          slurm_memory_gigabytes_per_cpu = 80,
          slurm_cpus_per_task = 1,
          workers = 30,
          tasks_max = 1,
          verbose = T
        ),
        crew_controller_slurm(
          name = "slurm_1_150",
          slurm_memory_gigabytes_per_cpu = 150,
          slurm_cpus_per_task = 1,
          workers = 2,
          tasks_max = 1,
          verbose = T
        ),
        crew_controller_slurm(
          name = "slurm_2_400",
          slurm_memory_gigabytes_per_cpu = 400,
          slurm_cpus_per_task = 2,
          workers = 1,
          tasks_max = 1,
          verbose = T
        )
      ),
      # debug = "pseudobulk_file_id_quantile_normalised_2_7240767f602a5810",
      
      resources = tar_resources(crew = tar_resources_crew("slurm_2_20")) 
      #, # Set the target you want to debug.
      #
    )
    
    split_metadata = function(){
      
      # CAQ_directory = "/vast/projects/cellxgene_curated/cellNexus"
      
      my_metadata =
        tbl(
          dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
          sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/cellNexus/cell_metadata_cell_type_consensus_v1_0_0.parquet')")
        ) |> 
        # get_metadata(cache_directory = CAQ_directory) |>
        
        # # DROP FETAL SCI-SEQ FOR THE MOMENT
        # filter(collection_id != "c114c20f-1ef4-49a5-9c2e-d965787fb90c") |> 
        # 
      
        # this, because for curated Atlas, query, I don't need of the other meta data, 
        # I can add it later from the full meta data 
        select(sample_id, cell_, file_id_cellNexus, cell_type_consensus_harmonised) |>
        filter(cell_type_consensus_harmonised=="other", cell_type_consensus_harmonised |> is.na()) |> 
        as_tibble() |> 
        mutate(chunk = sample_id) |> 
        nest(data = -c(chunk, file_id_cellNexus)) |> 
        
        # THIS SHOULD BE REMOVED IN THE FUTURE
        filter(!file_id_cellNexus |> is.na()) |> 
        mutate(file_id_cellNexus = paste0(file_id_cellNexus, ".h5ad")) |> 
        
        mutate(number_of_cells = map_int(data, nrow)) 
    
    }
    
    get_sce = 	function(tissue_cell_type_metadata, cache.path) {
      
      dir.create(cache.path, recursive = TRUE, showWarnings = FALSE)
      
      tissue_cell_type_metadata |>
        mutate(data = map2(
          data, file_id_cellNexus,
          ~ {
            
            # THIS IS TEMPORARY BECAUSE I HAVE BRAIN DATASETS WITH 1M CELL THAT ARE HARD TO SAVE IN THE DB
            if(!file.exists(paste0(cache.path, "/original/", .y))) {
              warning(paste0(cache.path, "/original/", .y, " does not exist in the cache") )
              return(NULL)
            } 
            
            .x |>
            mutate(file_id_cellNexus = .y) |> 
            CuratedAtlasQueryR:::get_data_container(
              repository = NULL, 
              cache_directory = cache.path, 
              grouping_column = "file_id_cellNexus"
            )
          }
          
            # get_single_cell_experiment(cache_directory = cache.path) 
          # |> 
          #   
          #   # !!!!
          #   # this is needed, because, for some reason, the local data repository, still has ENSEMBL genes, 
          #   # while the cloud repository does not
          #   _[sce_rownames,] 
          
          # |>
          # 	mutate(sample_se =
          # 				 	
          # 				 	# I need to fix Curated CellAtlas with disease sample, duplication for 
          # 				 	# file_id_cellNexus=="cc3ff54f-7587-49ea-b197-1515b6d98c4c", cell_type_consensus_harmonised=="stromal_cell"
          # 				 	# for lung
          # 				 	glue("{sample_id}___{disease}") |>
          # 				 	str_replace_all(" ", "_") |>
          # 				 	str_replace_all("/", "__")
          # 	) 
        )) |> 
        
        # THIS IS TEMPORARY BECAUSE I HAVE BRAIN DATASETS WITH 1M CELL THAT ARE HARD TO SAVE IN THE DB
        filter(!map_lgl(data, is.null) ) |> 
        filter(map_int(data, ncol)>0)
      
    }
    
    get_pseudobulk = 	function(sce_df) {
      
      sce_df |>
      
        
        mutate(data = map(
          data, 
          ~ {
            .x = tidySingleCellExperiment::aggregate_cells(.x, .sample = c(sample_id, cell_type_consensus_harmonised)	)
            
            # Decrease size
            # We will reattach rowames later
            my_assay = assay(.x, "counts") |>  as("sparseMatrix")
            # rownames(my_assay) = NULL
            SummarizedExperiment(assays = list(counts = my_assay), colData = colData(.x))
            
          }
        ))
      
    }
    
    aggregate = function(se_df, se_rownames){
      
      print("Start aggregate")
      gc()
      
 
        se_df |>
        
        # Add columns and filter
        mutate(data = map2(
          data, file_id_cellNexus,
          ~ {
            
            # ## TEMPORARY FIX BECAUSE I FORGOT TO ADD ASSAY NAME
            # assays(.x) = .x |> assays() |> as.list() |>  set_names("counts")
            # #####
            
            # Add columns
            se = 
              .x |>
              mutate(file_id_cellNexus = .y) |>
              select(-any_of(c("file_id_cellNexus", ".cell", "original_cell_id")))
            
            
            # # Identify samples with many genes
            # sample_with_many_genes =
            #   se |>
            #   assay("counts") |>
            #   apply(2, function(x) (x>0) |> which() |> length()) |>
            #   enframe() |>
            #   mutate(name = as.character(name)) |>
            #   filter(value > 5000) |>
            #   pull(name)
            # se = se[,sample_with_many_genes, drop=FALSE]
            # 
            # # Filter samples with too few cells
            # se = se |> filter(.aggregated_cells > 10)
            
            se
            
          }, 
          .progress=TRUE
        )) |>
        
        nest(data = -file_id_cellNexus) |> 
        mutate(data = map(data,
                          ~ {
                            se = .x |> 
                              pull(data) |> 
                              do.call(cbind, args= _)
                            
                            # # Filter very rare gene-transcripts
                            # all_zeros = assay(se, "counts") |> rowSums() |> equals(0) 
                            # se = se[!all_zeros,]
                            # lower_than_total_counts = assay(se, "counts") |> rowSums() < 15 
                            # se = se[!lower_than_total_counts,]
                            
                            
                            # Make it sparse
                            se@assays@data$counts = se@assays@data$counts |> as("sparseMatrix")
                            
                            # Attach rownames back
                            #rownames(se) = se_rownames
                            
                            se
                          })
              )
      
    }
    
    
    #-----------------------#
    # Pipeline
    #-----------------------#
    list(
      tar_target(cache.path, "/vast/projects/cellxgene_curated/cellNexus/cellxgene/29_10_2024/", deployment = "main"),
      
      # Get rownames
      tar_target(
        sce_rownames,
        tbl(
          dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
          sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/cellNexus/cell_metadata_cell_type_consensus_v1_0_0.parquet')")
        ) |>
          head(1) |> 
          mutate(file_id_cellNexus = paste0(file_id_cellNexus, ".h5ad")) |> 
          CuratedAtlasQueryR:::get_data_container(
            repository = NULL, 
            cache_directory = cache.path, 
            grouping_column = "file_id_cellNexus"
          ) |> 
          rownames(),
        resources = tar_resources(crew = tar_resources_crew("slurm_1_20")), 
        packages = c("arrow", "dplyr",  "duckdb")
      ),
      
      # Do metadata
      tar_target(
        metadata_split,
        split_metadata(),
        resources = tar_resources(crew = tar_resources_crew("slurm_1_40")), 
        packages = c("arrow", "dplyr",  "duckdb", "tidyr", "dplyr", "purrr")
      ),

      # Get SCE SMALL
      tarchetypes::tar_group_by(
        metadata_split_SMALL,
        metadata_split[metadata_split$number_of_cells<2500,],
        chunk,
        resources = tar_resources(crew = tar_resources_crew("slurm_1_20"))
      ),

      tarchetypes::tar_group_by(
        metadata_split_BIG,
        metadata_split[metadata_split$number_of_cells>=2500,],
        chunk,
        resources = tar_resources(crew = tar_resources_crew("slurm_1_80"))
      ),

      # Get SCE SMALL
      tar_target(
        metadata_split_pseudobulk_SMALL,
        metadata_split_SMALL |>  get_sce(cache.path) |> get_pseudobulk(),
        pattern = map(metadata_split_SMALL),
        resources = tar_resources(crew = tar_resources_crew("slurm_1_30"))
      ) ,

      # Get SCE BIG
      tar_target(
        metadata_split_pseudobulk_BIG,
        metadata_split_BIG |>  get_sce(cache.path) |> get_pseudobulk(),
        pattern = map(metadata_split_BIG),
        resources = tar_resources(crew = tar_resources_crew("slurm_1_80"))
      ),

      # Group samples
      tarchetypes::tar_group_by(
        metadata_grouped_pseudobulk,
        metadata_split_pseudobulk_SMALL |> bind_rows(metadata_split_pseudobulk_BIG),
        file_id_cellNexus,
        resources = tar_resources(crew = tar_resources_crew("slurm_2_400"))
      )  ,

      # Aggregate
      tar_target(
        pseudobulk_grouped_by_file_id,
        metadata_grouped_pseudobulk |> aggregate(sce_rownames) ,
        pattern = map(metadata_grouped_pseudobulk),
        resources = tar_resources(crew = tar_resources_crew("slurm_1_20"))
      ) ,

      # Quantile normalisation taking a "good looking distribution"
      tar_target(
        the_best_target_distribution,
        pseudobulk_grouped_by_file_id |>
          filter(file_id_cellNexus== "829f7db55e2ec1c73c75d83f2b950c73.h5ad") |>
          pull(data) |>
          _[[1]] |>
          assay( "counts") |>
         as.matrix() |>
          preprocessCore::normalize.quantiles.determine.target() ,
        resources = tar_resources(crew = tar_resources_crew("slurm_1_150"))
      ),

      tar_target(
        pseudobulk_file_id_quantile_normalised,
        pseudobulk_grouped_by_file_id |> mutate(data = map(
          data,
          
          # I DONT KNOW WHY SOME ROWNAMES ARE NA. IT COULD DEPEND ON CELLXGENE?
          ~ .x[!rownames(.x) |> is.na(),,drop=FALSE] |>
            quantile_normalise_abundance(
              method="preprocesscore_normalize_quantiles_use_target",
              target_distribution = the_best_target_distribution
            ) |>
            select(-counts)
        )) ,
        pattern = map(pseudobulk_grouped_by_file_id),
        resources = tar_resources(crew = tar_resources_crew("slurm_1_20"))
      )
      
    )
    
    
  }, 
  ask = FALSE, 
  script = glue("{result_directory}/_targets.R")
)

job::job({
  
  tar_make(
    # callr_function = NULL,
    reporter = "summary",
    script = glue("{result_directory}/_targets.R"),
    store = glue("{result_directory}/_targets")
  )
  
})


tar_meta(store = glue("{result_directory}/_targets")) |> 
  filter(!error |> is.na()) |> 
  arrange(desc(time)) |> 
  select(name, error) |>  
  View()

library(SummarizedExperiment)

x = tar_read(metadata_split_pseudobulk_SMALL, store = glue("{result_directory}/_targets"))


x = tar_read_raw("metadata_split_pseudobulk_SMALL", store = glue("{result_directory}/_targets"), branches = 1)


y = tar_read_raw("metadata_grouped_pseudobulk", store = glue("{result_directory}/_targets"), branches = 1:500)
  
  x |> filter(chunk == "03d9eb203bc24f826a30525d2e4d155a") |> pull(data) |> _[[10]] |> rownames() |> head()

x = x |> mutate(
  data = map(data, ~ .x[!rownames(.x) |> is.na(),,drop=FALSE] )
) |> 
  mutate(rownames = map(data, rownames))

tar_read_raw("metadata_split_pseudobulk_SMALL_74d9c5660535fcb7", store = glue("{result_directory}/_targets"), branches = 1)

# tar_invalidate(metadata_grouped_pseudobulk, store = glue("{result_directory}/_targets"))

tar_workspace(
  pseudobulk_file_id_quantile_normalised_5c8a7878e1e03d6e, 
  script = glue("{result_directory}/_targets.R"),
  store = glue("{result_directory}/_targets")
)

tar_meta(metadata_grouped_pseudobulk, store = glue("{result_directory}/_targets") )
