library(tidyverse)
library(targets)
library(glue)

library(arrow)
library(dplyr)
library(duckdb)
# Get input


result_directory = "/vast/scratch/users/mangiola.s/pseudobulk_cellNexus_1_0_2"

tar_script({
    
  result_directory = "/vast/scratch/users/mangiola.s/pseudobulk_cellNexus_1_0_2"
  
  
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
        "Seurat", "glue", "qs",  "purrr", "tidybulk", "tidySummarizedExperiment", "edgeR", 
        "crew", "magrittr", "digest", "readr", "forcats", "scuttle", "BiocParallel", "SummarizedExperiment"
      ),
      
      memory = "transient", 
      garbage_collection = 100, 
      storage = "worker", 
      retrieval = "worker", 
      error = "continue", 
      #  debug = "pseudobulk_file_id_quantile_normalised_processed_split_2_HDF5_fe7b6b0f0a71537c", 
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
          verbose = T, 
          seconds_idle = 30
        ),
        crew_controller_slurm(
          name = "slurm_1_10",
          slurm_memory_gigabytes_per_cpu = 10,
          slurm_cpus_per_task = 1,
          workers = 300,
          tasks_max = 20,
          verbose = T, 
          seconds_idle = 30
        ),
        crew_controller_slurm(
          name = "slurm_1_20",
          script_lines = "#SBATCH --mem 20G",
          slurm_cpus_per_task = 1,
          workers = 200,
          tasks_max = 20,
          verbose = T, 
          seconds_idle = 30
        ),
        crew_controller_slurm(
          name = "slurm_1_40",
          script_lines = "#SBATCH --mem 40G",
          slurm_cpus_per_task = 5,
          workers = 200,
          tasks_max = 20,
          verbose = T, 
          seconds_idle = 30
        ),
        crew_controller_slurm(
          name = "slurm_1_80",
          script_lines = "#SBATCH --mem 80G",
          slurm_cpus_per_task = 20,
          workers = 30,
          tasks_max = 1,
          verbose = T, 
          seconds_idle = 30
        ),
        crew_controller_slurm(
          name = "slurm_1_150",
          slurm_memory_gigabytes_per_cpu = 150,
          slurm_cpus_per_task = 1,
          workers = 2,
          tasks_max = 1,
          verbose = T, 
          seconds_idle = 30
        ),
        crew_controller_slurm(
          name = "slurm_2_400",
          slurm_memory_gigabytes_per_cpu = 400,
          slurm_cpus_per_task = 2,
          workers = 1,
          tasks_max = 1,
          verbose = T, 
          seconds_idle = 30
        )
      ),
      # debug = "pseudobulk_file_id_quantile_normalised_2_7240767f602a5810",
      
      resources = tar_resources(crew = tar_resources_crew("slurm_1_20")) 
      #, # Set the target you want to debug.
      #
    )
    
    split_metadata = function(metadata_parquet){
      
      # CAQ_directory = "/vast/projects/cellxgene_curated/cellNexus"
      
      my_metadata =
        tbl(
          dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
          sql(glue("SELECT * FROM read_parquet('{metadata_parquet}')"))
        ) |> 

        # 
      
        # this, because for curated Atlas, query, I don't need of the other meta data, 
        # I can add it later from the full meta data 
        select(sample_id, cell__id, file_id_cellNexus_single_Cell, cell_type_unified_ensemble) |>
        filter(cell_type_unified_ensemble!="other", !cell_type_unified_ensemble |> is.na()) |> 
        filter(!file_id_cellNexus_single_Cell |> is.na()) |> 
        
        # THIS SHOULD BE REMOVED IN THE FUTURE
        mutate(file_id_cellNexus_single_Cell = paste0(file_id_cellNexus_single_Cell, ".h5ad")) |> 
       
        # NEST 
        as_tibble() |> 
        mutate(chunk = file_id_cellNexus_single_Cell) |> 
        nest(data = -c(chunk)) |> 
        
        mutate(number_of_cells = map_int(data, nrow)) 
    
    }
    
    get_sce = 	function(tissue_cell_type_metadata, cache.path) {
      
      dir.create(cache.path, recursive = TRUE, showWarnings = FALSE)
      
      tissue_cell_type_metadata |>
        mutate(data = map2(
          data, chunk,
          ~ {
            
            file_id_cellNexus_single_Cell = .y
            
            # THIS IS TEMPORARY BECAUSE I HAVE BRAIN DATASETS WITH 1M CELL THAT ARE HARD TO SAVE IN THE DB
            if(!file.exists(paste0(cache.path, "/original/", file_id_cellNexus_single_Cell))) {
              warning(paste0(cache.path, "/original/", file_id_cellNexus_single_Cell, " does not exist in the cache") )
              return(NULL)
            } 
            
            .x |>
            CuratedAtlasQueryR:::get_data_container(
              repository = NULL, 
              cache_directory = cache.path, 
              grouping_column = "file_id_cellNexus_single_Cell"
            )
          }
          

        )) |> 
        
        # THIS IS TEMPORARY BECAUSE I HAVE BRAIN DATASETS WITH 1M CELL THAT ARE HARD TO SAVE IN THE DB
        filter(!map_lgl(data, is.null) ) |> 
        filter(map_int(data, ncol)>0)
      
    }
    
    get_pseudobulk = 	function(sce_df) {
      
      sce_df |>
      
        
        mutate(data = map2(
          data, number_of_cells,
          ~ {
            pseudobulk = 
              aggregateAcrossCells(
                .x, 
                colData(.x)[,c("sample_id", "cell_type_unified_ensemble")], 
                BPPARAM = MulticoreParam(workers = if_else(.y>50000, 20, 5))
              )
            colnames(pseudobulk) = paste0(colData(pseudobulk)$sample_id, "___", colData(pseudobulk)$cell_type_unified_ensemble)
            
            pseudobulk = pseudobulk |> select(.cell, sample_id, file_id_cellNexus_single_Cell, cell_type_unified_ensemble)
            
            # Decrease size
            # We will reattach rowames later
            assay(pseudobulk) = assay(pseudobulk) |>  as("sparseMatrix")
            
            pseudobulk |> as("SummarizedExperiment")

          }
        ))
      
    }
    
    aggregate = function(se_df, se_rownames){
      
      print("Start aggregate")
      gc()
      
 
        se_df |>
        
        # Add columns and filter
        mutate(data = map2(
          data, file_id_cellNexus_single_Cell,
          ~ {
            
            # ## TEMPORARY FIX BECAUSE I FORGOT TO ADD ASSAY NAME
            # assays(.x) = .x |> assays() |> as.list() |>  set_names("counts")
            # #####
            
            # Add columns
            se = 
              .x |>
              mutate(file_id_cellNexus_single_Cell = .y) |>
              select(-any_of(c("file_id_cellNexus_single_Cell", ".cell", "original_cell_id")))
            
            
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
        
        nest(data = -file_id_cellNexus_single_Cell) |> 
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
    
    #' Process a SummarizedExperiment Object
    #'
    #' This function processes a single \code{SummarizedExperiment} object by adjusting its assay to include all genes,
    #' filling missing genes with zeros, and recording gene presence per sample as a separate assay.
    #' It returns a \code{SummarizedExperiment} with the adjusted assay and gene presence information as a secondary assay.
    #'
    #' **Why This Function Is Needed:**
    #'
    #' When combining multiple \code{SummarizedExperiment} objects, each may have different gene sets.
    #' This function standardizes each object's assay to include all common genes, facilitating seamless combination.
    #' By returning a \code{SummarizedExperiment}, it simplifies downstream processing and integration.
    #'
    #' @param se A \code{SummarizedExperiment} object to be processed.
    #' @param gene_index A named integer vector mapping gene names to their indices.
    #' @param se_id A unique identifier for the \code{SummarizedExperiment}, used to ensure unique sample names.
    #' @param external_directory A directory path to save the processed \code{SummarizedExperiment}.
    #' @return A \code{SummarizedExperiment} object with:
    #'   \item{assays}{A list containing the adjusted assay with all genes, named using the original assay name, and a "gene_presence" assay.}
    #'   \item{colData}{The column data from the original \code{SummarizedExperiment}, with updated sample names.}
    #'
    #' @examples
    #' # Process a single SummarizedExperiment
    #' processed_se <- process_se(se, gene_index, se_id = 1, external_directory = "output_directory")
    #'
    #' @importFrom SummarizedExperiment assay
    #' @importFrom SummarizedExperiment assayNames
    #' @importFrom SummarizedExperiment colData
    #' @importFrom SummarizedExperiment SummarizedExperiment
    #' @importFrom S4Vectors DataFrame
    #' @importFrom Matrix sparseMatrix
    #' @importFrom HDF5Array saveHDF5SummarizedExperiment
    #' @importFrom stringr str_remove
    #' @export
    process_se <- function(se, gene_index, se_id, external_directory) {
      
      dir.create(external_directory, recursive = TRUE, showWarnings = FALSE)
      
      all_genes <- names(gene_index)
      n_genes <- length(gene_index)
      genes_se <- rownames(se)
      n_cells <- ncol(se)
      
      # Update sample names to be unique by prefixing with se_id
      sample_names = colnames(se)
      
      # Get indices of genes in the full gene list
      idx <- gene_index[genes_se]
      
      # Get the assay data and convert to sparse TsparseMatrix if not already
      assay_se <- assay(se)
      if (!is(assay_se, "TsparseMatrix")) {
        assay_se <- as(as(assay_se, "sparseMatrix"), "TsparseMatrix")
      }
      
      # Adjust indices
      i_indices <- idx[assay_se@i + 1]
      j_indices <- assay_se@j + 1
      
      # Create sparse matrix for the assay with all genes
      assay_full <- sparseMatrix(
        i = i_indices,
        j = j_indices,
        x = assay_se@x,
        dims = c(n_genes, n_cells),
        dimnames = list(all_genes, sample_names)
      )
      
      # Record presence of genes per sample in a sparse binary matrix
      gene_presence_i <- Matrix(FALSE, nrow = n_genes, ncol = n_cells, sparse = TRUE)
      gene_presence_i[idx, ] <- TRUE
      colnames(gene_presence_i) <- sample_names
      rownames(gene_presence_i) <- all_genes
      
      # Get the original assay name(s)
      assay_name <- assayNames(se)[1]  # Assuming there is only one assay
      
      # Create processed SummarizedExperiment with original assay name and "gene_presence"
      processed_se <- SummarizedExperiment(
        assays = setNames(list(assay_full, gene_presence_i), c(assay_name, "gene_presence")),
        colData = colData(se)
      )
      
      # Save the processed SummarizedExperiment to an external file
      HDF5Array::saveHDF5SummarizedExperiment(
        processed_se,
        paste0(external_directory, "/", str_remove(se_id, "\\.h5ad")),
        as.sparse = TRUE,
        replace = TRUE
      )
      
  
    }
    
    #-----------------------#
    # Pipeline
    #-----------------------#
    list(
      tar_target(cache.path, "/vast/scratch/users/mangiola.s/cellNexus/cellxgene/15_11_2024/", deployment = "main"),
      tar_target(
        metadata_parquet, 
        "/vast/projects/cellxgene_curated/cellNexus/cell_metadata_cell_type_consensus_v1_0_2.parquet", 
        format = "file", 
        deployment = "main"
      ),
      
  
      # # Get rownames
      # tar_target(
      #   sce_rownames,
      #   tbl(
      #     dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
      #     sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/cellNexus/cell_metadata_cell_type_consensus_v1_0_2.parquet')")
      #   ) |>
      #     head(1) |> 
      #     mutate(file_id_cellNexus_single_Cell = paste0(file_id_cellNexus_single_Cell, ".h5ad")) |> 
      #     CuratedAtlasQueryR:::get_data_container(
      #       repository = NULL, 
      #       cache_directory = cache.path, 
      #       grouping_column = "file_id_cellNexus_single_Cell"
      #     ) |> 
      #     rownames(),
      #   resources = tar_resources(crew = tar_resources_crew("slurm_1_20")), 
      #   packages = c("arrow", "dplyr",  "duckdb")
      # ),
      
      # Do metadata
      tar_target(
        metadata_split,
        split_metadata(metadata_parquet),
        resources = tar_resources(crew = tar_resources_crew("slurm_1_20")), 
        packages = c("arrow", "dplyr",  "duckdb", "tidyr", "dplyr", "purrr", "glue")
      ),

      # Get SCE SMALL
      tarchetypes::tar_group_by(
        metadata_split_SMALL,
        metadata_split[metadata_split$number_of_cells<50000,],
        chunk,
        resources = tar_resources(crew = tar_resources_crew("slurm_1_20"))
      ),

      tarchetypes::tar_group_by(
        metadata_split_BIG,
        metadata_split[metadata_split$number_of_cells>=50000,],
        chunk,
        resources = tar_resources(crew = tar_resources_crew("slurm_1_80"))
      ),

      # Get SCE SMALL
      tar_target(
        metadata_split_pseudobulk_SMALL,
        metadata_split_SMALL |>  get_sce(cache.path) |> get_pseudobulk(),
        pattern = map(metadata_split_SMALL),
        resources = tar_resources(crew = tar_resources_crew("slurm_1_40"))
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
        chunk,
        resources = tar_resources(crew = tar_resources_crew("slurm_1_150")), 
        deployment = "main"
      )  ,

      # Quantile normalisation taking a "good looking distribution"
      tar_target(
        the_best_target_distribution,
        metadata_grouped_pseudobulk |>
          filter(chunk== "829f7db55e2ec1c73c75d83f2b950c73.h5ad") |>
          pull(data) |>
          _[[1]] |>
          assay( "counts") |>
         as.matrix() |>
          preprocessCore::normalize.quantiles.determine.target() ,
        resources = tar_resources(crew = tar_resources_crew("slurm_1_150"))
      ),

      tar_target(
        pseudobulk_file_id_quantile_normalised,
        metadata_grouped_pseudobulk |> mutate(data = map(
          data,
          
          # I DONT KNOW WHY SOME ROWNAMES ARE NA. IT COULD DEPEND ON CELLXGENE?
          ~ {
            if(rownames(.x) |> is.null() & rownames(assay(.x, withDimnames = FALSE)) |> is.null())
              stop("No rownames present anywhere")
            if(rownames(.x) |> is.null()) rownames(.x) = rownames(assay(.x, withDimnames = FALSE))
            .x |>
            quantile_normalise_abundance(
              method="preprocesscore_normalize_quantiles_use_target",
              target_distribution = the_best_target_distribution
            ) |>
            select(-counts)
          }
        )) ,
        pattern = map(metadata_grouped_pseudobulk),
        resources = tar_resources(crew = tar_resources_crew("slurm_1_40"))
      ),
      
      tar_target(
        sample_genes, 
        pseudobulk_file_id_quantile_normalised |> pull(data) |> _[[1]] |> rownames(), 
        pattern = map(pseudobulk_file_id_quantile_normalised), 
        iteration = "list",
        resources = tar_resources(crew = tar_resources_crew("slurm_1_40"))
      ),
      tar_target(
        all_genes,
        {
          # Create a mapping from gene names to indices
          all_genes <- unique(unlist(sample_genes))
          setNames(seq_len(length(all_genes)), all_genes)
        },
        resources = tar_resources(crew = tar_resources_crew("slurm_1_40"))
      ),
      tar_target(
        pseudobulk_file_id_quantile_normalised_processed,
        process_se(
          pseudobulk_file_id_quantile_normalised |> pull(data) |> _[[1]], 
          all_genes, 
          pseudobulk_file_id_quantile_normalised |> pull(chunk) |> _[[1]], 
          paste0(result_directory, "/external")
        ),
        pattern = map(pseudobulk_file_id_quantile_normalised),
        iteration = "list",
        resources = tar_resources(crew = tar_resources_crew("slurm_1_40"))
      ),
      
      # SPLIT 1
      tar_target(
        pseudobulk_file_id_quantile_normalised_processed_split_1,
        pseudobulk_file_id_quantile_normalised_processed |> 
          split(ceiling(seq_along(pseudobulk_file_id_quantile_normalised_processed) / 10)) ,
        iteration = "list",
        resources = tar_resources(crew = tar_resources_crew("slurm_1_80"))
      ),
      tar_target(
        pseudobulk_file_id_quantile_normalised_processed_split_1_HDF5,
        pseudobulk_file_id_quantile_normalised_processed_split_1 |> 
          do.call(cbind, args = _) |> 
          HDF5Array::saveHDF5SummarizedExperiment(
            paste0(result_directory, "/external", digest::digest(pseudobulk_file_id_quantile_normalised_processed_split_1)),
            as.sparse = TRUE,
            replace = TRUE,
            verbose = TRUE
          ) ,
        pattern = map(pseudobulk_file_id_quantile_normalised_processed_split_1),
        iteration = "list",
        resources = tar_resources(crew = tar_resources_crew("slurm_1_40"))
      ),
      
      # SPLIT 2
      tar_target(
        pseudobulk_file_id_quantile_normalised_processed_split_2,
        pseudobulk_file_id_quantile_normalised_processed_split_1 |> 
          split(ceiling(seq_along(pseudobulk_file_id_quantile_normalised_processed_split_1) / 10)) ,
        iteration = "list",
        resources = tar_resources(crew = tar_resources_crew("slurm_1_80"))
      ),
      tar_target(
        pseudobulk_file_id_quantile_normalised_processed_split_2_HDF5,
        pseudobulk_file_id_quantile_normalised_processed_split_2 |> 
          unlist() |> # Not sure why this is needed
          do.call(cbind, args = _) |> 
          HDF5Array::saveHDF5SummarizedExperiment(
            paste0(result_directory, "/external", digest::digest(pseudobulk_file_id_quantile_normalised_processed_split_2)),
            as.sparse = TRUE,
            replace = TRUE,
            verbose = TRUE
          ) ,
        pattern = map(pseudobulk_file_id_quantile_normalised_processed_split_2),
        iteration = "list",
        resources = tar_resources(crew = tar_resources_crew("slurm_1_40"))
      ),
      
      # SPLIT 3
      tar_target(
        pseudobulk_file_id_quantile_normalised_processed_split_3,
        pseudobulk_file_id_quantile_normalised_processed_split_2 |> 
          split(ceiling(seq_along(pseudobulk_file_id_quantile_normalised_processed_split_2) / 10)) ,
        iteration = "list",
        resources = tar_resources(crew = tar_resources_crew("slurm_1_80"))
      ),
      tar_target(
        pseudobulk_file_id_quantile_normalised_processed_split_3_HDF5,
        pseudobulk_file_id_quantile_normalised_processed_split_3 |> 
          unlist() |> # Not sure why this is needed
          do.call(cbind, args = _) |> 
          HDF5Array::saveHDF5SummarizedExperiment(
            paste0(result_directory, "/external", digest::digest(pseudobulk_file_id_quantile_normalised_processed_split_3)),
            as.sparse = TRUE,
            replace = TRUE,
            verbose = TRUE
          ) ,
        pattern = map(pseudobulk_file_id_quantile_normalised_processed_split_3),
        iteration = "list",
        resources = tar_resources(crew = tar_resources_crew("slurm_1_40"))
      )
      
    )
    
    
  }, 
  ask = FALSE, 
  script = glue("{result_directory}/_targets.R")
)

job::job({
  
  tar_make(
    # callr_function = NULL,
    reporter = "verbose_positives",
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
library(tidySummarizedExperiment)

sccomp_counts = readRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_2/cell_metadata_1_0_2_sccomp_input_counts.rds")

# Save into one SE
job::job({
  tar_read(pseudobulk_file_id_quantile_normalised_processed_split_3_HDF5, store = glue("{result_directory}/_targets"))  |>
    map(~ .x |> inner_join(sccomp_counts)) |> 
    do.call(cbind, args = _) |> 
    HDF5Array::saveHDF5SummarizedExperiment(
      "/vast/projects/cellxgene_curated/cellNexus/pseudobulk_joined",
      as.sparse = TRUE,
      replace = TRUE,
      verbose = TRUE
    )
})

system("~/bin/rclone copy /vast/projects/cellxgene_curated/cellNexus/pseudobulk_joined box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")


se = HDF5Array::loadHDF5SummarizedExperiment("/vast/projects/cellxgene_curated/cellNexus/pseudobulk_joined")
job::job({
  zellkonverter::writeH5AD(
  se, 
  verbose = TRUE, 
  file = "/vast/projects/cellxgene_curated/cellNexus/pseudobulk_sample_cell_type_1_0_2.h5ad", 
  X_name = "counts_scaled"
  )
})

#assays(se) = assays(se)[1] # aggregateAcrossCells does not like booleans

se_sample = 
  se |> 
  aggregateAcrossCells(
  colData(se)[,"sample_id"], 
  use.assay.type = c("counts_scaled", "gene_presence"),  # Specify the correct assay name if there are multiple assays
  BPPARAM =  MulticoreParam(workers = 10)
)
se_sample |> 
  HDF5Array::saveHDF5SummarizedExperiment(
    "/vast/projects/cellxgene_curated/cellNexus/pseudobulk_sample",
    as.sparse = TRUE,
    replace = TRUE,
    verbose = TRUE
  )
system("~/bin/rclone copy /vast/projects/cellxgene_curated/cellNexus/pseudobulk_sample box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")

tar_read_raw("pseudobulk_file_id_quantile_normalised_processed_split_1_HDF5_1b4b19de1079ec95", store = glue("{result_directory}/_targets"))

# tar_invalidate(metadata_split, store = glue("{result_directory}/_targets"))
# tar_invalidate(pseudobulk_file_id_quantile_normalised_processed_split_1, store = glue("{result_directory}/_targets"))

tar_workspace(
  "metadata_split_SMALL", 
  script = glue("{result_directory}/_targets.R"),
  store = glue("{result_directory}/_targets")
)
