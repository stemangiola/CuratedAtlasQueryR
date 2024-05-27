library(tidyverse)
library(targets)
library(glue)

result_directory = "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024"

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
        "zellkonverter", "cellxgenedp", "CuratedAtlasQueryR", "stringr", "tibble", "tidySingleCellExperiment", "dplyr", "Matrix",
        "glue", "qs",  "purrr", "tidybulk", "tidySummarizedExperiment",  "crew", "magrittr", "digest", "readr", "forcats"
      ),
      
      garbage_collection = TRUE,
      #trust_object_timestamps = TRUE,
      memory = "transient",
      storage = "worker",
      retrieval = "worker",
      #error = "continue",
      format = "qs",
      
      # Set the target you want to debug.
      #debug = "metadata_sample",

      cue = tar_cue(mode = "never")		,
      
      #-----------------------#
      # SLURM
      #-----------------------#
      
      controller = crew_controller_group(
        
        
        crew_controller_slurm(
          name = "slurm_1_5",
          slurm_memory_gigabytes_per_cpu = 5,
          slurm_cpus_per_task = 1,
          workers = 200,
          tasks_max = 5,
          verbose = T
        ),
        crew_controller_slurm(
          name = "slurm_1_10",
          slurm_memory_gigabytes_per_cpu = 10,
          slurm_cpus_per_task = 1,
          workers = 100,
          tasks_max = 5,
          verbose = T
        ),
        crew_controller_slurm(
          name = "slurm_1_20",
          slurm_memory_gigabytes_per_cpu = 20,
          slurm_cpus_per_task = 1,
          workers = 100,
          tasks_max = 5,
          verbose = T, 
          script_directory = paste0("/vast/scratch/users/mangiola.s/cellxgenedp/crew_cluster/", basename(tempdir()))
        ),
        crew_controller_slurm(
          name = "slurm_1_40",
          slurm_memory_gigabytes_per_cpu = 40,
          slurm_cpus_per_task = 1,
          workers = 50,
          tasks_max = 5,
          verbose = T
        ),
        crew_controller_slurm(
          name = "slurm_1_80",
          slurm_memory_gigabytes_per_cpu = 80,
          slurm_cpus_per_task = 1,
          workers = 30,
          tasks_max = 5,
          verbose = T
        ),
        crew_controller_slurm(
          name = "slurm_1_200",
          slurm_memory_gigabytes_per_cpu = 200,
          slurm_cpus_per_task = 1,
          workers = 5,
          tasks_max = 5,
          verbose = T
        )
      ),
      resources = tar_resources(crew = tar_resources_crew("slurm_1_10"))  
      
    )
    
    sample_heuristics = function(col_data){
      
      col_data |>
        
        # Sort sample ID
        # Fix some sample id missing
        when(unique(.$dataset_id) %in% c(
          "11b86bc3-6d4d-4e28-903a-0361ea8f6bdf",
          "492b0613-ff5b-4fca-a585-503fc4102e4f",
          "11b86bc3-6d4d-4e28-903a-0361ea8f6bdf",
          "0e8f9ce4-46e5-434e-9ca0-e769d1dd27ea"
        ) ~ mutate(., PatientID = glue("{sample} {replicate} {time_point} {target}") |> as.character()) , ~ (.)) %>%
        when(unique(.$dataset_id) %in% c(
          "0273924c-0387-4f44-98c5-2292dbaab11e",
          "a16bec18-5c9f-40ad-8169-12c5199c7506",
          "556bb449-bbef-43d3-9487-87031fc0decb"
        ) ~ mutate(., PatientID = glue("{Collection.ID} {Genotype} {Location}")|> as.character()) , ~ (.)) %>%
        when(unique(.$dataset_id) %in% c(
          "b83afdc1-baa1-42c0-bd5b-cb607084757d"
        ) ~ mutate(., PatientID = glue("{sex} {development_stage} {disease}")|> as.character()) , ~ (.)) %>%
        when(unique(.$dataset_id) %in% c(
          "3fe53a40-38ff-4f25-b33b-e4d60f2289ef",
          "5c1cc788-2645-45fb-b1d9-2f43d368bba8"
        ) ~ mutate(., PatientID = glue("{Batch} {Fetus_id} {Development_day} {sex} {tissue} {disease}")|> as.character()) , ~ (.)) |>
        
        
        mutate_if(is.factor, as.character) |>
        type_convert(guess_integer = TRUE) |>
        mutate_if(is.integer, as.character) |>
        
        # Convert types
        when("donor_id" %in% colnames(.) ~ mutate(., donor_id = donor_id |> as.character() ), ~(.)) |>
        when("Cluster" %in% colnames(.) ~ mutate(., Cluster = Cluster |> as.character() ), ~(.)) |>
        when("cluster_id" %in% colnames(.) ~ mutate(., cluster_id = cluster_id |> as.character() ), ~(.)) |>
        when("Batch" %in% colnames(.) ~ mutate(., Batch = Batch |> as.character() ), ~(.)) |>
        when("batch" %in% colnames(.) ~ mutate(., batch = batch |> as.character() ), ~(.)) |>
        when("age" %in% colnames(.) ~ mutate(., age = age |> as.character() ), ~(.)) |>
        when("BMI" %in% colnames(.) ~ mutate(., BMI = BMI |> as.character() ), ~(.)) |>
        when("donor_BMI" %in% colnames(.) ~ mutate(., donor_BMI = donor_BMI |> as.character() ), ~(.)) |>
        when("author_cell_type" %in% colnames(.) ~ mutate(., author_cell_type = author_cell_type |> as.character() ), ~(.)) |>
        when("time_point" %in% colnames(.) ~ mutate(., time_point = time_point |> as.character() ), ~(.)) |>
        when("cluster" %in% colnames(.) ~ mutate(., cluster = cluster |> as.character() ), ~(.)) |>
        when("ClusterID" %in% colnames(.) ~ mutate(., ClusterID = ClusterID |> as.character() ), ~(.)) |>
        when("Stage" %in% colnames(.) ~ mutate(., Stage = Stage |> as.character() ), ~(.)) |>
        when("individual" %in% colnames(.) ~ mutate(., individual = individual |> as.character() ), ~(.)) |>
        when("recurrent_cluster" %in% colnames(.) ~ mutate(., recurrent_cluster = recurrent_cluster |> as.character() ), ~(.)) |>
        when("PatientID" %in% colnames(.) ~ mutate(., PatientID = PatientID |> as.character() ), ~(.)) |>
        when("PMI" %in% colnames(.) ~ mutate(., PMI = PMI |> as.character() ), ~(.)) |>
        when("n_genes" %in% colnames(.) ~ mutate(., n_genes = n_genes |> as.numeric() ), ~(.)) |>
        when("n_counts" %in% colnames(.) ~ mutate(., n_counts = n_counts |> as.numeric() ), ~(.)) |>
        when("n_genes_by_counts" %in% colnames(.) ~ mutate(., n_genes_by_counts = n_genes_by_counts |> as.numeric() ), ~(.)) |>
        when("nUMI" %in% colnames(.) ~ mutate(., nUMI = nUMI |> as.numeric() ), ~(.)) |>
        when("percent.cortex" %in% colnames(.) ~ mutate(., percent.cortex = percent.cortex |> as.character() ), ~(.)) |>
        when("percent.medulla" %in% colnames(.) ~ mutate(., percent.medulla = percent.medulla |> as.character() ), ~(.)) |>
        when("Age" %in% colnames(.) ~ mutate(., Age = Age |> as.numeric() ), ~(.)) |>
        when("nCount_RNA" %in% colnames(.) ~ mutate(., nCount_RNA = nCount_RNA |> as.numeric() ), ~(.)) |>
        when("is_primary_data" %in% colnames(.) ~ mutate(., is_primary_data = is_primary_data |> as.character() ), ~(.)) |>
        
        #mutate(across(contains("cluster", ignore.case = TRUE), ~ as.character)) |>
        select(-one_of('PCW')) %>%
        
        # Sort sample ID. It works but not elegant.
        # Based on observation of strangely behaving datasets, where sample ID is not clear
        when("sampleID" %in% colnames(.) & !"PatientID" %in% colnames(.) ~
               mutate(., PatientID = as.character(sampleID )) |>  select(-sampleID), ~(.)) %>%
        when("Patient" %in% colnames(.) ~ mutate(., Sample = NA |> as.character()), ~(.)) |>
        mutate(sample_placeholder = NA |> as.character()) %>%
        when(unique(.$dataset_id)=="e40591e7-0e5a-4bef-9b60-7015abe5b17f" ~ mutate(., sample_placeholder = glue("{batch} {development_stage}") |> as.character()), ~ (.)) %>%
        when(unique(.$dataset_id)=="39b6cc45-8c5c-4f7b-944c-58f66da5efb1" ~ mutate(., sample_placeholder =sample_id), ~ (.))  %>%
        when(unique(.$dataset_id)=="443d6a0e-dbcb-4002-8af0-628e7d4a18fa" ~ mutate(., sample_placeholder =sample_id), ~ (.))  %>%
        when(unique(.$dataset_id)=="a91f075b-52d5-4aa3-8ecc-86c4763a49b3" ~ mutate(., sample_placeholder =sample), ~ (.))  %>%
        when(unique(.$dataset_id)=="0af763e1-0e2f-4de6-9563-5abb0ad2b01e" ~ mutate(., sample_placeholder ="only_one_culture"), ~ (.))  %>%
        when(unique(.$dataset_id)=="5c64f247-5b7c-4842-b290-65c722a65952" ~ mutate(., sample_placeholder ="only_one_culture"), ~ (.))  %>%
        when(unique(.$dataset_id)=="d6f92754-e178-4202-b86f-0f430e965d72" ~ mutate(., sample_placeholder =orig.ident), ~ (.))  %>%
        when(unique(.$dataset_id)=="c790ef7a-1523-4627-8603-d6a02f8f4877" ~ mutate(., sample_placeholder =orig.ident), ~ (.))  %>%
        when(unique(.$dataset_id)=="1e81a742-e457-4fc6-9c39-c55189ec9dc2" ~ mutate(., sample_placeholder =orig.ident), ~ (.))  %>%
        when(unique(.$dataset_id)=="351ef284-b59e-43a5-83ba-0eb907dc282c" ~ mutate(., sample_placeholder =orig.ident), ~ (.))  %>%
        when(unique(.$dataset_id)=="f498030e-246c-4376-87e3-90b28c7efb00" ~ mutate(., sample_placeholder =Name), ~ (.))  %>%
        
        # These are the datasets with too few cells per inferred samples, therefore simplifying
        when(unique(.$dataset_id)=="e3a56e00-8417-4d82-9d35-3fab3aac12f2" ~ mutate(., SpecimenID =NA), ~ (.))  %>%
        when(unique(.$dataset_id)=="17b34e42-bbd2-494b-bf32-b9229344a3f6" ~ mutate(., Sample =NA), ~ (.))  %>%
        
        # Fix huge samples for plate experiments
        tidyr::extract(.cell, "experiment___", "(^expr?[0-9]+)", remove = F) |>
        tidyr::extract(.cell, c("run_from_cell_id"), "(run[a-zA-Z0-9_]+)-.+", remove = FALSE) |> 
        
        mutate(experiment___ = if_else(dataset_id=="3fe53a40-38ff-4f25-b33b-e4d60f2289ef", experiment___, "")) |>
        
        # If run-based embrio study get sample ID from cell ID
        
        
        # Empirically infer samples from many characteristics
        unite("sample_heuristic", one_of(
          "sample_placeholder",
          "Sample",
          "SampleID",
          "sample_uuid",
          "Sample_ID",
          "scRNASeq_sample_ID",
          "Sample_Tag",
          "Sample.ID",
          "sample_names",
          "Short_Sample",
          "Sample.ID.short",
          "Sample.name",
          "patient",
          "Donor.ID",
          "donor_id",
          "donor",
          "PatientID",
          "donor_uuid",
          "library_uuid",
          "suspension_uuid",
          "Patient",
          "tissue_section_uuid",
          "DonorID",
          "specimen",
          "SpecimenID",
          "Fetus_id",
          "individual",
          "tissue",
          "development_stage",
          "assay",
          "experiment___",
          "disease",
          "run_from_cell_id"
        ), na.rm = TRUE, sep = "___", remove = F) |>
        
        
        #parquet does not like . prefix
        rename(cell_ = .cell) |> 
        
        # Add sample hash
        mutate(sample_ = getVDigest(algo="md5")(glue("{sample_heuristic}{dataset_id}"))) |>
        
        # make lighter
        mutate_if(is.character, as.factor)
      
    }
    
    get_metadata = function(.x){
      
      cache.path = "/vast/scratch/users/mangiola.s/cellxgenedp"
      dir.create(cache.path, recursive = TRUE, showWarnings = FALSE)

      h5_path = .x |> files_download(dry.run = FALSE, cache.path = cache.path)
      
      sce = 
        h5_path |> 
        readH5AD(use_hdf5 = TRUE,  obs = FALSE, raw = FALSE, skip_assays = TRUE, layers=FALSE, reader = "R"	)
      
      file.remove(h5_path)
      
      metadata = 
        sce |> 
        as_tibble() 
      
      # join the file metadata
      column_to_omit_becuse_duplicated = 
        colnames(.x) |> 
        intersect(colData(sce) |> colnames()) |> 
        str_subset("donor_id", negate = TRUE) |> 
        c("embedding")
      
      rm(sce)
      gc(verbose = FALSE)
      
      metadata = 
        metadata |> 
        left_join(
          .x |> 
            select(!any_of(column_to_omit_becuse_duplicated)) |> 
            unnest(donor_id) |> 
            unnest(donor_id) |> 
            select_if(negate(is.list)) ,
          by = join_by(donor_id)
        ) 
      
      metadata = 
        metadata |> 
        sample_heuristics() 
      
      # # delete raw data
      # sample_column_to_preserve = 
      #   metadata |> 
      #   slice_sample(n = 500, by = donor_id) |> 
      #   tidybulk::pivot_sample(.sample = sample_) |> 
      #   colnames()
      
      # # Select only sample_ columns
      # metadata = 
      #   metadata |> 
      #   select(sample_, any_of(sample_column_to_preserve)) |> 
      #   distinct()
      

      
      metadata
    }
    
    select_sample_columns = function(metadata){
      
      # delete raw data
      sample_column_to_preserve =
        metadata |>
        slice_sample(n = 500, by = donor_id) |>
        tidybulk::pivot_sample(.sample = sample_) |>
        colnames()
      
      # Select only sample_ columns
        metadata |>
        select(sample_, any_of(sample_column_to_preserve)) |>
        distinct()
      
    }
    #-----------------------#
    # Pipeline
    #-----------------------#
    list(
      
      # Get rownames
      tar_target(
        my_db,
        db(overwrite=TRUE),
        resources = tar_resources(crew = tar_resources_crew("slurm_1_20"))
      ),
      
      # Get SCE SMALL
      tar_target(
        files_dataset_id,
        datasets(my_db) |>
          left_join(
            files(my_db) |> filter(filetype=="H5AD"),
            by = "dataset_id"
          ) |> 
          group_split(dataset_id),
        iteration = "list",
        resources = tar_resources(crew = tar_resources_crew("slurm_1_20"))
      ),
      
      # Get SCE SMALL
      tar_target(
        metadata_dataset_id,
        get_metadata(files_dataset_id),
        pattern = map(files_dataset_id),
        iteration = "list",
        resources = tar_resources(crew = tar_resources_crew("slurm_1_20"))
      ),
      
      tar_target(
        common_columns,
        metadata_dataset_id |>
          map_dfr(~ .x |> colnames() |> as_tibble()) |> 
          dplyr::count(value) |>
          mutate(n_datasets = length(metadata_dataset_id)) |>
          filter(n > (n_datasets / 2)) |>
          pull(value) ,
        resources = tar_resources(crew = tar_resources_crew("slurm_1_200"))
      ),
      
      tar_target(
        metadata_dataset_id_common_sample_columns,
        metadata_dataset_id |> 
          mutate(cell_ =  as.character(cell_)) |> 
          select(any_of(common_columns)) |> 
          select_sample_columns(),
        pattern = map(metadata_dataset_id),
        resources = tar_resources(crew = tar_resources_crew("slurm_1_20") )
      ),
      
      tar_target(
        metadata_dataset_id_cell_to_sample_mapping,
        metadata_dataset_id |> 
          mutate(cell_ =  as.character(cell_)) |> 
          select(cell_, sample_, donor_id),
        pattern = map(metadata_dataset_id),
        resources = tar_resources(crew = tar_resources_crew("slurm_1_20"))
      )
      
    )
    
    
  }, 
  ask = FALSE, 
  script = glue("{result_directory}/_targets.R")
)


tar_make(
  #callr_function = NULL,
  reporter = "verbose_positives",
  script = glue("{result_directory}/_targets.R"),
  store = glue("{result_directory}/_targets")
)

# Sample metadata
tar_read(metadata_dataset_id_common_sample_columns, store = glue("{result_directory}/_targets"))

# Sample to cell link
tar_read(metadata_dataset_id_cell_to_sample_mapping, store = glue("{result_directory}/_targets"))

# 
# # a test to see whether donor ID is present in the new metadata
# test = 
#   files_metadata |> 
#   slice(1:50) |> 
#   nest(data = c(dataset_id, dataset_version_id, filetype, url)) |> 
#   mutate(has_donor_id = map_lgl(
#     data,
#     ~ {
#       browser()
#       h5_path = .x |> files_download(dry.run = FALSE)
#       has_donor_id = 
#         h5_path |> 
#         readH5AD(use_hdf5 = TRUE	) |> 
#         colData() |> 
#         as_tibble() |> 
#         select(any_of("donor_id")) |> 
#         ncol() >
#         0
#       file.remove(h5_path)
#       has_donor_id
#     }
#   )) |> 
#   unnest(data) |> 
#   select(dataset_version_id, has_donor_id)
# 
# 


# files_metadata |>
#   
#   # Get organism list and filter human
#   mutate(organism_name = map_chr(organism, ~ .x |> map(~.x$label) |> paste(collapse=", ") )) |>
#   filter(organism_name |> str_detect("Homo sapiens")) |>
#   
#   # Download
#   files_download(dry.run = FALSE, cache_path = "{root_directory}/raw_data/") |>
#   
#   # Save file list
#   saveRDS("{root_directory}/file_location.rds")
# 
