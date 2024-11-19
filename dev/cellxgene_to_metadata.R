.rs.restartR()
library(tidyverse)
library(targets)
library(glue)

result_directory = "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024"


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
      "zellkonverter", "cellxgenedp", "CuratedAtlasQueryR", "stringr", "tibble", "tidySingleCellExperiment", "dplyr", "Matrix",
      "glue", "qs",  "purrr", "tidybulk", "tidySummarizedExperiment",  "crew", "magrittr", "digest", "readr", "forcats"
    ),
    
    memory = "transient", 
    garbage_collection = 100, 
    storage = "worker", 
    retrieval = "worker", 
    error = "continue", 
    #  debug = "dataset_id_sce_b5312463451d7ee3", 
    cue = tar_cue(mode = "never"),
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
        tasks_max = 5,
        verbose = T, 
        seconds_idle = 30
      ),
      crew_controller_slurm(
        name = "slurm_1_10",
        slurm_memory_gigabytes_per_cpu = 10,
        slurm_cpus_per_task = 1,
        workers = 100,
        tasks_max = 5,
        verbose = T, 
        seconds_idle = 30
      ),
      crew_controller_slurm(
        name = "slurm_1_20",
        slurm_memory_gigabytes_per_cpu = 20,
        slurm_cpus_per_task = 1,
        workers = 100,
        tasks_max = 5,
        verbose = T, , 
        seconds_idle = 30
      ),
      crew_controller_slurm(
        name = "slurm_1_40",
        slurm_memory_gigabytes_per_cpu = 40,
        slurm_cpus_per_task = 1,
        workers = 50,
        tasks_max = 5,
        verbose = T, 
        seconds_idle = 30
      ),
      crew_controller_slurm(
        name = "slurm_1_80",
        slurm_memory_gigabytes_per_cpu = 80,
        slurm_cpus_per_task = 1,
        workers = 30,
        tasks_max = 5,
        verbose = T, 
        seconds_idle = 30
      ),
      crew_controller_slurm(
        name = "slurm_1_200",
        slurm_memory_gigabytes_per_cpu = 200,
        slurm_cpus_per_task = 1,
        workers = 5,
        tasks_max = 5,
        verbose = T, 
        seconds_idle = 30
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
        "run_from_cell_id",
        "is_primary_data"
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
      readH5AD(use_hdf5 = TRUE,  raw = FALSE, skip_assays = TRUE, layers=FALSE, reader = "R"	) 
    
    
    if(is.null(sce) || !"donor_id" %in% colnames(colData(sce)))
      sce = 
      h5_path |> 
      readH5AD(use_hdf5 = TRUE, raw = FALSE, skip_assays = TRUE, layers=FALSE	)
    
    
    
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
      resources = tar_resources(crew = tar_resources_crew("slurm_1_80"))
    ),
    
    # Get SCE SMALL
    tar_target(
      metadata_dataset_id,
      get_metadata(files_dataset_id),
      pattern = map(files_dataset_id),
      iteration = "list",
      resources = tar_resources(crew = tar_resources_crew("slurm_1_20")), 
      deployment = "main"
    ),
    
    # select column that are present in half of the datasets at least, so the common column
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
        
        # Only get primary data
        # filter(is_primary_data=="TRUE") |> 
        
        mutate(cell_ =  as.character(cell_)) |> 
        select(any_of(common_columns)) |> 
        
        # Drop some clearly cell-wise columns
        select(-any_of(c("observation_joinid", "cell_")), -contains("cell_type")) |> 
        
        select_sample_columns(),
      pattern = map(metadata_dataset_id),
      resources = tar_resources(crew = tar_resources_crew("slurm_1_80") )
    ),
    
    tar_target(
      metadata_dataset_id_cell_to_sample_mapping,
      metadata_dataset_id |> 
        
        # Only get primary data
        # filter(is_primary_data=="TRUE") |> 
        
        mutate(
          cell_ =  as.character(cell_), 
          observation_joinid = as.character(observation_joinid)
        ) |> 
        # select(cell_, observation_joinid, sample_, donor_id),
        select(observation_joinid, cell_, sample_, donor_id, dataset_id, is_primary_data, sample_heuristic, cell_type, cell_type_ontology_term_id),
      pattern = map(metadata_dataset_id),
      resources = tar_resources(crew = tar_resources_crew("slurm_1_80")),
      deployment = "main"
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


library(arrow)
library(dplyr)
library(duckdb)

# # Sample metadata
# tar_read(metadata_dataset_id_common_sample_columns, store = glue("{result_directory}/_targets")) |>
#   write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/sample_metadata.parquet")

# # Sample to cell link
# tar_read(metadata_dataset_id_cell_to_sample_mapping, store = glue("{result_directory}/_targets")) |> 
#   write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_ids_for_metadata.parquet")


get_tissue_grouped = function(tissue){
  
  list(
    
    # Respiratory System
    "respiratory system" = c(
      "lung", "lung parenchyma", "alveolus of lung",  "bronchus",
      "respiratory airway", "pleura", "pleural effusion", "middle lobe of right lung",
      "upper lobe of left lung", "lower lobe of left lung", "upper lobe of right lung",
      "lower lobe of right lung", "lingula of left lung", "right lung", "left lung"
    ),
    
    trachea = c( "epithelium of trachea", "trachea"),
    
    # Cardiovascular System
    "cardiovascular system" = c(
      "heart", "heart left ventricle", "heart right ventricle", "cardiac ventricle",
      "cardiac atrium", "right cardiac atrium", "left cardiac atrium", "apex of heart",
      "aorta", "coronary artery", 
      "venous blood", "anterior wall of left ventricle", "myocardium", "interventricular septum", "ventricular tissue", "basal zone of heart"
    ),
    
    vasculature = c("kidney blood vessel", "artery", "vein", "vasculature", "mesenteric artery"),
    # Umbilical Cord Blood
    "umbilical cord blood" = "umbilical cord blood",
    
    # Oesophagus
    "oesophagus" = c(
      "esophagus", "lower esophagus", "esophagus muscularis mucosa",
      "submucosal esophageal gland"
    ),
    
    # Stomach
    "stomach" = c(
      "stomach", "body of stomach", "cardia of stomach"
    ),
    
    # Small Intestine
    "small intestine" = c(
      "small intestine", "duodenum", "jejunum", "ileum"
    ),
    
    # Large Intestine
    "large intestine" = c(
      "large intestine", "colon", "left colon", "right colon",
      "sigmoid colon", "descending colon", "transverse colon",
      "ascending colon", "hepatic flexure of colon", "caecum",
      "rectum", "appendix", "vermiform appendix"
    ),
    
    # Digestive System (General)
    "digestive system (general)" = c(
      "intestine", "hindgut"
    ),
    
    # Nasal, Oral, and Pharyngeal Regions
    "nasal, oral, and pharyngeal regions" = c(
      "nasal cavity", "nasopharynx", "oral mucosa", "tongue", "anterior part of tongue",
      "posterior part of tongue", "gingiva", "nose", "saliva"
    ),
    
    # Cerebral Lobes and Cortical Areas
    "cerebral lobes and cortical areas" = c(
      "frontal lobe", "left frontal lobe", "right frontal lobe", "primary motor cortex",
      "dorsolateral prefrontal cortex", "superior frontal gyrus", "orbitofrontal cortex",
      "medial orbital frontal cortex", "Broca's area", "prefrontal cortex",
      "temporal lobe", "left temporal lobe", "right temporal lobe", 
      "angular gyrus", "entorhinal cortex",
      "parietal lobe", "left parietal lobe", "right parietal lobe", "primary somatosensory cortex",
      "occipital lobe", "right occipital lobe", "primary visual cortex",
      "occipital cortex", "insular cortex", "parietal cortex", "temporal cortex",
      "frontal cortex", "Brodmann (1909) area 4", "temporoparietal junction",
      "middle temporal gyrus", "cingulate cortex", "brain", "brain white matter", "cerebral cortex", "cerebral nuclei"
    ),
    
    # Limbic and Basal Systems
    "limbic and basal systems" = c(
      "anterior cingulate cortex", "anterior cingulate gyrus", "hippocampal formation",
      "hypothalamus", "thalamic complex", "dentate nucleus", "basal ganglion",
      "caudate nucleus", "putamen", "substantia nigra pars compacta",
      "lateral ganglionic eminence", "medial ganglionic eminence",
      "caudal ganglionic eminence", "ganglionic eminence"
    ),
    
    # Brainstem and Cerebellar Structures
    "brainstem and cerebellar structures" = c(
      "pons", "midbrain", "myelencephalon", "telencephalon", "forebrain",
      "cerebellum", "cerebellum vermis lobule", "cerebellar cortex",
      "hemisphere part of cerebellar posterior lobe", "white matter of cerebellum"
    ),
    
    # General Brain and Major Structures
    "general brain and major structures" = c(
      "spinal cord", "neural tube", "cervical spinal cord white matter"
    ),
    
    # Muscular System (Skeletal Muscles)
    "muscular system (skeletal muscles)" = c(
      "rectus abdominis muscle", "gastrocnemius", "muscle of abdomen", "muscle organ",
      "muscle tissue", "pelvic diaphragm muscle", "skeletal muscle tissue", "muscle of pelvic diaphragm"
    ),
    
    # Connective Tissue
    "connective tissue" = c(
      "connective tissue", "tendon of semitendinosus", "vault of skull", "bone spine",
      "rib"
    ),
    
    # Adipose Tissue
    "adipose tissue" = c(
      "adipose tissue", "subcutaneous adipose tissue", "visceral abdominal adipose tissue",
      "perirenal fat", "omental fat pad", "subcutaneous abdominal adipose tissue",
      "abdominal adipose tissue"
    ),
    
    # Endocrine System
    "endocrine system" = c(
      "thyroid gland", "adrenal tissue", "adrenal gland", "islet of Langerhans",
      "endocrine pancreas", "pineal gland"
    ),
    
    # Lymphatic System
    "lymphatic system" = c(
      "lymph node", "mesenteric lymph node", "thoracic lymph node",
      "cervical lymph node", "bronchopulmonary lymph node", "tonsil", "inguinal lymph node"
    ),
    
    # Integumentary System (Skin)
    "integumentary system (skin)" = c(
      "skin of abdomen", "skin of forearm", "skin of scalp", "skin of face", "skin of leg",
      "skin of chest", "skin of back", "skin of hip", "skin of body", "skin of cheek",
      "skin of temple", "skin of shoulder", "skin of external ear", "skin of trunk",
      "skin of prepuce of penis", "skin epidermis", "arm skin", "lower leg skin",
      "hindlimb skin", "zone of skin", "dermis", "skin of nose", "skin of forehead",
      "skin of pes", "axilla"
    ),
    
    # Gastrointestinal Accessory Organs
    "gastrointestinal accessory organs" = c(
      "liver", "gallbladder", "pancreas", "exocrine pancreas", "caudate lobe of liver",
      "hepatic cecum"
    ),
    
    # Spleen
    "spleen" = "spleen",
    
    # Thymus
    "thymus" = "thymus",
    
    # Blood
    "blood" = "blood",
    
    # Bone Marrow
    "bone marrow" = "bone marrow",
    
    # Female Reproductive System
    "female reproductive system" = c(
      "uterus", "myometrium", "fallopian tube", "ampulla of uterine tube",
      "fimbria of uterine tube", "uterine cervix", "endometrium",
      "decidua", "decidua basalis", "placenta", "yolk sac", "isthmus of fallopian tube"
    ),
    "ovary" = "ovary", 
    
    # Male Reproductive System
    "male reproductive system (other)" = c(
      "testis", "gonad"
    ),
    
    # Prostate
    "prostate" = c(
      "prostate gland", "transition zone of prostate", "peripheral zone of prostate"
    ),
    
    # Renal System
    "renal system" = c(
      "kidney", "cortex of kidney", "renal medulla", "renal papilla",
      "renal pelvis", "ureter", "bladder organ"
    ),
    
    # Miscellaneous Glands
    "miscellaneous glands" = c(
      "parotid gland", "lacrimal gland", "sublingual gland", "mammary gland",
      "chorionic villus"
    ),
    
    # Epithelium and Mucosal Tissues
    "epithelium and mucosal tissues" = c(
      "epithelium of small intestine", "epithelium of esophagus", "caecum epithelium",
      "jejunal epithelium", "ileal epithelium", "colonic epithelium",
      "submucosa of ascending colon", "submucosa of ileum", "lamina propria",
      "lamina propria of large intestine", "lamina propria of small intestine",
      "mucosa", "mucosa of colon", "lamina propria of mucosa of colon"
    ),
    
    # Eye and Visual-Related Structures
    "sensory-related structures" = c(
      "retina",
      "retinal neural layer",
      "macula lutea",
      "macula lutea proper",
      "sclera",
      "trabecular meshwork",
      "conjunctiva",
      "pigment epithelium of eye",
      "cornea",
      "iris",
      "ciliary body",
      "peripheral region of retina",
      "eye trabecular meshwork",
      "perifoveal part of retina",
      "choroid plexus",
      "lens of camera-type eye",
      "corneo-scleral junction",
      "fovea centralis",
      "eye",
      "inner ear",
      "vestibular system",
      "primary auditory cortex"
    ),
    
    # Digestive Tract Junctions and Connections
    "digestive tract junctions and connections" = c(
      "esophagogastric junction", "duodeno-jejunal junction", "hepatopancreatic ampulla",
      "hepatopancreatic duct", "pyloric antrum"
    ),
    
    # Peritoneal and Abdominal Cavity Structures
    "peritoneal and abdominal cavity structures" = c(
      "peritoneum", "omentum", "retroperitoneum", "mesentery"
    ),
    
    # Breast
    "breast" = c(
      "breast", "upper outer quadrant of breast"
    )
  ) |> 
    enframe(name ="tissue_groups") |> 
    distinct() |> 
    unnest(value) |> 
    rename(tissue = value) |> 
    mutate()
  
  # #check
  # distinct_tissue = 
  #   tissue |> 
  #   enframe(name = "tissue") |> 
  #   distinct(tissue)
  # 
  # if(nrow(distinct_tissue) != distinct_tissue |> left_join(tissue_grouped_df, copy = TRUE))
  # 
  # 
  # tissue |> 
  #   enframe(name = "tissue") |> 
  #   left_join(tissue_grouped_df)
}

convert_age_labels_to_days <- function(labels) {
  # Initialize vector to store age in days
  age_days <- rep(NA, length(labels))
  
  # Define the mapping for Carnegie stages
  carnegie_stages <- c(
    '9' = 20,
    '10' = 22,
    '11' = 24,
    '12' = 26,
    '13' = 28,
    '14' = 32,
    '16' = 37,
    '17' = 41,
    '18' = 44,
    '19' = 46,
    '20' = 49,
    '21' = 51,
    '22' = 53,
    '23' = 56
  )
  
  # Map words to numbers for decades
  word_to_num <- c(
    'first' = 0,
    'second' = 10,
    'third' = 20,
    'fourth' = 30,
    'fifth' = 40,
    'sixth' = 50,
    'seventh' = 60,
    'eighth' = 70,
    'ninth' = 80,
    'tenth' = 90
  )
  
  # Map ordinal words to numbers
  word_ordinal_to_num <- c(
    'first' = 1,
    'second' = 2,
    'third' = 3,
    'fourth' = 4,
    'fifth' = 5,
    'sixth' = 6,
    'seventh' = 7,
    'eighth' = 8,
    'ninth' = 9,
    'tenth' = 10,
    'eleventh' = 11,
    'twelfth' = 12,
    'thirteenth' = 13,
    'fourteenth' = 14,
    'fifteenth' = 15,
    'sixteenth' = 16,
    'seventeenth' = 17,
    'eighteenth' = 18,
    'nineteenth' = 19,
    'twentieth' = 20,
    'twenty-first' = 21,
    'twenty-second' = 22,
    'twenty-third' = 23
  )
  
  # Map stages to approximate ages in days
  stage_to_age <- list(
    'newborn human' = 0,
    'infant' = 0.5 * 365,
    'child' = 6 * 365,  # Midpoint of 2-12 years
    'adolescent' = 15 * 365,  # Midpoint of 12-18 years
    'young adult' = 25 * 365,  # Approximate age
    'human early adulthood' = 25 * 365,
    'human middle aged' = 50 * 365,
    'human late adulthood' = 70 * 365,
    'human adult' = 40 * 365,
    'human aged' = 75 * 365,
    'mature' = 40 * 365,
    'immature' = 1 * 365,
    'embryonic human' = 28,  # Midpoint of embryonic stage
    'organogenesis' = 28,
    'unknown' = NA
  )
  
  # Loop over labels
  for (i in seq_along(labels)) {
    label <- labels[i]
    
    # Initialize age variable
    age <- NA
    
    # Remove leading and trailing whitespaces
    label <- trimws(label)
    
    # 1. Match "unknown"
    if (grepl("^unknown$", label, ignore.case = TRUE)) {
      age <- NA
    }
    # 2. Match "[number]-month-old human stage"
    else if (grepl("^(\\d+)-month-old human stage$", label)) {
      num <- as.numeric(sub("^(\\d+)-month-old human stage$", "\\1", label))
      age <- num * 30  # Average days in a month
    }
    # 3. Match "[number]-year-old human stage"
    else if (grepl("^(\\d+)-year-old human stage$", label)) {
      num <- as.numeric(sub("^(\\d+)-year-old human stage$", "\\1", label))
      age <- num * 365  # Average days in a year
    }
    # 4. Match "[number]th week post-fertilization human stage"
    else if (grepl("^(\\d+)(?:st|nd|rd|th) week post-fertilization human stage$", label)) {
      num <- as.numeric(sub("^(\\d+)(?:st|nd|rd|th) week post-fertilization human stage$", "\\1", label))
      age <- num * 7
    }
    # 5. Match "Carnegie stage [number]"
    else if (grepl("^Carnegie stage (\\d+)$", label)) {
      num <- sub("^Carnegie stage (\\d+)$", "\\1", label)
      if (num %in% names(carnegie_stages)) {
        age <- carnegie_stages[[num]]
      } else {
        age <- NA
      }
    }
    # 6. Match "[number]-[number] year-old human stage"
    else if (grepl("^(\\d+)-(\\d+) year-old human stage$", label)) {
      num1 <- as.numeric(sub("^(\\d+)-(\\d+) year-old human stage$", "\\1", label))
      num2 <- as.numeric(sub("^(\\d+)-(\\d+) year-old human stage$", "\\2", label))
      avg_years <- (num1 + num2) / 2
      age <- avg_years * 365
    }
    # 7. Match "[number]-[number] year-old child stage"
    else if (grepl("^(\\d+)-(\\d+) year-old child stage$", label)) {
      num1 <- as.numeric(sub("^(\\d+)-(\\d+) year-old child stage$", "\\1", label))
      num2 <- as.numeric(sub("^(\\d+)-(\\d+) year-old child stage$", "\\2", label))
      avg_years <- (num1 + num2) / 2
      age <- avg_years * 365
    }
    # 8. Match "under-1-year-old human stage"
    else if (grepl("^under-1-year-old human stage$", label)) {
      age <- 0.5 * 365  # Assume 0.5 years
    }
    # 9. Match "[number]-month-old human stage" (again)
    else if (grepl("^(\\d+)-month-old human stage$", label)) {
      num <- as.numeric(sub("^(\\d+)-month-old human stage$", "\\1", label))
      age <- num * 30
    }
    # 10. Match "[ordinal] LMP month human stage"
    else if (grepl("^(\\w+) LMP month human stage$", label)) {
      ordinal_word <- tolower(sub("^(\\w+) LMP month human stage$", "\\1", label))
      if (ordinal_word %in% names(word_ordinal_to_num)) {
        num <- word_ordinal_to_num[ordinal_word]
        age <- num * 30
      } else {
        age <- NA
      }
    }
    # 11. Match "[ordinal] decade human stage"
    else if (grepl("^(\\w+) decade human stage$", label)) {
      decade_word <- tolower(sub("^(\\w+) decade human stage$", "\\1", label))
      if (decade_word %in% names(word_to_num)) {
        num1 <- word_to_num[decade_word]
        num2 <- num1 + 9
        avg_years <- (num1 + num2) / 2
        age <- avg_years * 365
      } else {
        age <- NA
      }
    }
    # 12. Match "80 year-old and over human stage"
    else if (grepl("^(\\d+).*year-old and over human stage$", label)) {
      num <- as.numeric(sub("^(\\d+).*year-old and over human stage$", "\\1", label))
      age <- num * 365
    }
    # 13. Match developmental stages
    else if (grepl("^(.*) stage$", label)) {
      stage <- tolower(sub("^(.*) stage$", "\\1", label))
      if (stage %in% names(stage_to_age)) {
        age <- stage_to_age[[stage]]
      } else {
        age <- NA
      }
    }
    # 14. Default case
    else {
      age <- NA
    }
    
    # Assign to age_days vector
    age_days[i] <- age
  }
  
  return(age_days |> as.integer())
}

age_days_tbl = 
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/sample_metadata.parquet')")
  ) |> 
  distinct(development_stage) |> 
  as_tibble() |> 
  mutate(age_days = convert_age_labels_to_days(development_stage)) 

age_days_tbl |>
  write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/age_days.parquet")



# # 
# cell_ids_for_metadata <- tbl(
#   dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
#   sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_ids_for_metadata.parquet')")
# )
# 
# cell_to_refined_sample_from_Mengyuan <- tbl(
#   dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
#   sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/census_samples_to_download.parquet')")
# ) |>
#   select(cell_, observation_joinid, dataset_id, sample_id = sample_2) |> 
#   

# 
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

# Establish a connection to DuckDB in memory
job::job({
  
  con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  # Create views for each of the datasets in DuckDB
  dbExecute(con, "
  CREATE VIEW cell_to_refined_sample_from_Mengyuan AS
  SELECT cell_, observation_joinid, dataset_id, sample_2 AS sample_id, cell_type, cell_type_ontology_term_id
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/census_samples_to_download.parquet')
")
  
  dbExecute(con, "
  CREATE VIEW cell_ids_for_metadata AS
  SELECT cell_, observation_joinid, dataset_id, sample_, donor_id
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_ids_for_metadata.parquet')
")
  
  dbExecute(con, "
  CREATE VIEW sample_metadata AS
  SELECT *
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/sample_metadata.parquet')
")
  
  dbExecute(con, "
  CREATE VIEW age_days_tbl AS
  SELECT development_stage, age_days
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/age_days.parquet')
")
  
  dbExecute(con, "
  CREATE VIEW tissue_grouped AS
  SELECT tissue, tissue_groups
  FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/tissue_grouped.parquet')
")
  
  # Perform optimised joins within DuckDB
  copy_query <- "
COPY (
  SELECT 
    cell_to_refined_sample_from_Mengyuan.cell_,
    cell_to_refined_sample_from_Mengyuan.observation_joinid,
    cell_to_refined_sample_from_Mengyuan.dataset_id,
    cell_to_refined_sample_from_Mengyuan.sample_id,
    cell_to_refined_sample_from_Mengyuan.cell_type,
    cell_to_refined_sample_from_Mengyuan.cell_type_ontology_term_id,
    sample_metadata.*,
    age_days_tbl.age_days,
    tissue_grouped.tissue_groups
  
  FROM cell_to_refined_sample_from_Mengyuan
  
  LEFT JOIN cell_ids_for_metadata
    ON cell_ids_for_metadata.cell_ = cell_to_refined_sample_from_Mengyuan.cell_
    AND cell_ids_for_metadata.observation_joinid = cell_to_refined_sample_from_Mengyuan.observation_joinid
    AND cell_ids_for_metadata.dataset_id = cell_to_refined_sample_from_Mengyuan.dataset_id
    
  LEFT JOIN sample_metadata
    ON cell_ids_for_metadata.sample_ = sample_metadata.sample_
    AND cell_ids_for_metadata.donor_id = sample_metadata.donor_id
    AND cell_ids_for_metadata.dataset_id = sample_metadata.dataset_id
    
  LEFT JOIN age_days_tbl
    ON age_days_tbl.development_stage = sample_metadata.development_stage

  LEFT JOIN tissue_grouped
    ON tissue_grouped.tissue = sample_metadata.tissue
    
) TO '/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata.parquet'
(FORMAT PARQUET, COMPRESSION 'gzip');
"
  
  # Execute the final query to write the result to a Parquet file
  dbExecute(con, copy_query)
  
  # Disconnect from the database
  dbDisconnect(con, shutdown = TRUE)
  
})

# system("~/bin/rclone copy /vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata.parquet box_adelaide:/Mangiola_ImmuneAtlas/reannotation_consensus/")


cell_metadata = tbl(
  dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata.parquet')")
) 

tissues_grouped = get_tissue_grouped() 

tissues_grouped |>
  write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/tissue_grouped.parquet")

##############
#  PLOTS     #
##############

cell_metadata |> 
  left_join(tissues_grouped, copy = TRUE) |> 
  mutate(age_days = convert_age_labels_to_days(development_stage))
  distinct(donor_id, tissue_groups) |> 
  ggplot(aes(fct_infreq(tissue_groups))) +
  geom_bar() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

cell_metadata |> 
  left_join(tissues_grouped, copy = TRUE) |> 
  left_join(days_df, copy = TRUE) |> 
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

cell_metadata |> 
  left_join(tissues_grouped, copy = TRUE) |> 
  left_join(days_df, copy = TRUE) |> 
  write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/cell_metadata_temp.parquet")



# Disconnect from the database
dbDisconnect(con, shutdown = TRUE)


# Get Dharmesh metadata consensus
# system("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/reannotation_consensus/consensus_output_new.parquet /vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/ ")
# system("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/reannotation_consensus/data_driven_consensus_new.parquet /vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/ ")
# system("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/reannotation_consensus/data_driven_consensus.parquet /vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/ ")

# Non immune harmonisation to Dharmesh immune harmonisation
non_immune_harmonisation = 
  read_csv("/vast/projects/mangiola_immune_map/PostDoc/CuratedAtlasQueryR/dev/cell_type_harmonisation_non_immune.csv") 

# system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/CuratedAtlasQueryR/dev/cell_type_harmonisation_non_immune.csv box_adelaide:/Mangiola_ImmuneAtlas/reannotation_consensus/")


tbl(
  dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/consensus_output_new.parquet')")
) |> 
  
  # Add non immune harmonisation to Dharmesh immune harmonisation
  mutate(is_immune = reannotation_consensus == "non immune") |> 
  left_join(non_immune_harmonisation, copy = TRUE
  ) |> 
  mutate(reannotation_consensus = case_when(reannotation_consensus=="non immune" ~ non_immune_harmonised_cell_type, TRUE ~ reannotation_consensus)) |> 
  select(-non_immune_harmonised_cell_type) |> 
  write_parquet_to_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/consensus_output_new_plus_non_immune_harmonisation.parquet")


annotation_with_harmonised <- tbl(
  dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Apr_2024/consensus_output_new_plus_non_immune_harmonisation.parquet')")
)


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


