# how do we want to structure meta and counts container on the cloud?
# how does file_id_db get generated? if load all counts into one container, there would be an issue if file_id_db in new atlas duplicate with that on the cloud 
# Ideally a container should contain meta, and directories for raw counts and cpm

assay_map <- c(counts = "original",
               cpm = "cpm")

testing <- read_parquet("testing.0.2.3.parquet")
testing<- tibble(testing)

upload_cloud <- function(meta,
                         assay = "counts",
                         version = "0.2.3",
                         cache_dir = CuratedAtlasQueryR:::get_default_cache_dir()) {
  #check if the software is in the directory
  check_file_exists(file.path(cache_dir, "harmonised_human_cell_atlas-openrc.sh")) 
  meta_input <- paste(meta, version, sep=".")
  upload_meta <- glue("swift upload {meta_input} {cache_dir}/{meta_input}.parquet)")
  upload_counts <- glue("swift upload {meta_input} {cache_dir}/{assay})")
  # Swift software can be installed when python>=3.11
  system("module load python/3.11", intern = TRUE,)
  check_swift <- system("swift list --help", intern = TRUE, ignore.stderr = TRUE)
  
  if (length(check_swift) == 0){
    system("python - pip install python-swiftclient", intern = TRUE)
  }
  system("source harmonised_human_cell_atlas-openrc.sh", intern = TRUE)
  system(upload_meta)
  system(upload_counts)
}
