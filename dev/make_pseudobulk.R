library(tidyverse)
library(targets)
library(glue)
# Get input


result_directory = "/vast/projects/cellxgene_curated/pseudobulk_0.2.3.4"

tar_script(
  {
    
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
      
      garbage_collection = TRUE,
      #trust_object_timestamps = TRUE,
      memory = "transient",
      storage = "worker",
      retrieval = "worker",
      #error = "continue",
      format = "qs",
      
      #-----------------------#
      # SLURM
      #-----------------------#
      controller = crew_controller_group(
        
        
        crew_controller_slurm(
          name = "slurm_1_5",
          slurm_memory_gigabytes_per_cpu = 5,
          slurm_cpus_per_task = 1,
          workers = 200,
          verbose = T
        ),
        crew_controller_slurm(
          name = "slurm_1_10",
          slurm_memory_gigabytes_per_cpu = 10,
          slurm_cpus_per_task = 1,
          workers = 100,
          verbose = T
        ),
        crew_controller_slurm(
          name = "slurm_1_20",
          slurm_memory_gigabytes_per_cpu = 20,
          slurm_cpus_per_task = 1,
          workers = 50,
          verbose = T
        ),
        crew_controller_slurm(
          name = "slurm_1_40",
          slurm_memory_gigabytes_per_cpu = 40,
          slurm_cpus_per_task = 1,
          workers = 30,
          verbose = T
        ),
        crew_controller_slurm(
          name = "slurm_1_80",
          slurm_memory_gigabytes_per_cpu = 80,
          slurm_cpus_per_task = 1,
          workers = 20,
          verbose = T
        ),
        crew_controller_slurm(
          name = "slurm_2_400",
          slurm_memory_gigabytes_per_cpu = 400,
          slurm_cpus_per_task = 2,
          workers = 1,
          verbose = T
        )
      ),
      resources = tar_resources(crew = tar_resources_crew("slurm_2_20")) ,
      debug = "pseudobulk_file_id_3e12a178"
      #, # Set the target you want to debug.
      #cue = tar_cue(mode = "never")		
    )
    
    split_metadata = function(){
      
      CAQ_directory = "/vast/projects/cellxgene_curated"
      
      my_metadata =
        get_metadata(cache_directory = CAQ_directory) |>
        distinct(file_id, sample_, cell_, file_id_db, cell_type_harmonised) |>
        as_tibble() |> 
        mutate(chunk = sample_) |> 
        nest(data = -c(chunk, file_id)) |> 
        mutate(number_of_cells = map_int(data, nrow))
      
    }
    
    get_sce = 	function(tissue_cell_type_metadata) {
      
      tissue_cell_type_metadata |>
        mutate(data = pmap(
          list(data),
          ~ ..1 |>
            get_single_cell_experiment(cache_directory = "/vast/projects/cellxgene_curated") 
          # |>
          # 	mutate(sample_se =
          # 				 	
          # 				 	# I need to fix Curated CellAtlas with disease sample, duplication for 
          # 				 	# file_id=="cc3ff54f-7587-49ea-b197-1515b6d98c4c", cell_type_harmonised=="stromal_cell"
          # 				 	# for lung
          # 				 	glue("{sample_}___{disease}") |>
          # 				 	str_replace_all(" ", "_") |>
          # 				 	str_replace_all("/", "__")
          # 	) 
        ))
      
    }
    
    get_pseudobulk = 	function(sce_df) {
      
      sce_df |>
        mutate(data = map(
          data, 
          ~ {
            .x = tidySingleCellExperiment::aggregate_cells(.x, .sample = c(sample_, cell_type_harmonised)	)
            
            # Decrease size
            # We will reattach rowames later
            my_assay = assay(.x, "counts") |>  as("sparseMatrix")
            rownames(my_assay) = NULL
            SummarizedExperiment(assays = my_assay, colData = colData(.x))
            
          }
        ))
      
    }
    
    aggregate = function(se_df, se_rownames){
      
      print("Start aggregate")
      gc()
      
      se = 
        se_df |>
        
        # Add columns and filter
        mutate(data = map2(
          data, file_id,
          ~ {
            
            ## TEMPORARY FIX BECAUSE I FORGOT TO ADD ASSAY NAME
            assays(.x) = .x |> assays() |> as.list() |>  set_names("counts")
            #####
            
            # Add columns
            se = 
              .x |>
              mutate(file_id = .y) |>
              select(-any_of(c("file_id_db", ".cell", "original_cell_id")))
            
            
            # Identify samples with many genes
            sample_with_many_genes =
              se |>
              assay("counts") |>
              apply(2, function(x) (x>0) |> which() |> length()) |>
              enframe() |>
              mutate(name = as.character(name)) |>
              filter(value > 5000) |>
              pull(name)
            se = se[,sample_with_many_genes, drop=FALSE]
            
            # Filter samples with too few cells
            se = se |> filter(.aggregated_cells > 10)
            
            se
            
          }, 
          .progress=TRUE
        )) |>
        
        nest(data = -file_id) |> 
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
                            rownames(se) = se_rownames
                            
                            se
                          }))
      
    }
    
    
    #-----------------------#
    # Pipeline
    #-----------------------#
    list(
      
      # # Do metadata
      # tarchetypes::tar_group_by(
      #   chunks_df,
      #   define_chunks(),
      #   file_id, sample_chunk, 
      #   resources = tar_resources(crew = tar_resources_crew("slurm_1_40"))
      # ),
      
      # Get rownames
      tar_target(
        sce_rownames,
        get_metadata() |> head(1) |> get_single_cell_experiment() |> rownames(),
        resources = tar_resources(crew = tar_resources_crew("slurm_1_20"))
      ),
      
      # Do metadata
      tar_target(
        metadata_split,
        split_metadata(),
        resources = tar_resources(crew = tar_resources_crew("slurm_1_40"))
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
        resources = tar_resources(crew = tar_resources_crew("slurm_1_20"))
      ),
      
      # Get SCE SMALL
      tar_target(
        metadata_split_pseudobulk_SMALL,
        metadata_split_SMALL |>  get_sce() |> get_pseudobulk(),
        pattern = map(metadata_split_SMALL),
        iteration = "group",
        resources = tar_resources(crew = tar_resources_crew("slurm_1_20"))
      ),
      
      # Get SCE SMALL
      tar_target(
        metadata_split_pseudobulk_BIG,
        metadata_split_BIG |>  get_sce() |> get_pseudobulk(),
        pattern = map(metadata_split_BIG),
        iteration = "group",
        resources = tar_resources(crew = tar_resources_crew("slurm_1_80"))
      ),
      
      # Group samples
      tarchetypes::tar_group_by(
        metadata_grouped_pseudobulk,
        metadata_split_pseudobulk_SMALL |> bind_rows(metadata_split_pseudobulk_BIG),
        file_id,
        resources = tar_resources(crew = tar_resources_crew("slurm_2_400")), 
        deployment = "main"
      ),
      
      # Aggregate
      tar_target(
        pseudobulk_file_id,
        metadata_grouped_pseudobulk |> aggregate(sce_rownames) ,
        pattern = map(metadata_grouped_pseudobulk),
        iteration = "group",
        resources = tar_resources(crew = tar_resources_crew("slurm_1_20"))
      )
      
    )
    
    
  }, 
  ask = FALSE, 
  script = glue("{result_directory}/_targets.R")
)


tar_make(
  #callr_function = NULL,
  # reporter = "summary",
  script = glue("{result_directory}/_targets.R"),
  store = result_directory
)



# tar_read(pseudobulk_file_id, store = "/vast/projects/cellxgene_curated/pseudobulk_0.2.3.4")

