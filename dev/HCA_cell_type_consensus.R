library(arrow)
library(dplyr)
library(duckdb)
library(HPCell)

# Read the Parquet file into an R data frame
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
parquet_file = "/vast/projects/cellxgene_curated/census_samples/concensus_input.parquet"

data_tbl <- tbl(con, sql(paste0("SELECT * FROM read_parquet('", parquet_file, "')")))

# annotation_combination = 
#   data_tbl |> 
#   #select(azimuth_predicted.celltype.l2, monaco_first.labels.fine, blueprint_first.labels.fine) |> 
#   select(cell_, dataset_id, cell_type, cell_type_ontology_term_id, azimuth_predicted.celltype.l2, monaco_first.labels.fine, blueprint_first.labels.fine)  
#   #arrange(desc(n)) |> 

annotation_consensus  = 
  data_tbl |>
  distinct(azimuth_predicted.celltype.l2, monaco_first.labels.fine, blueprint_first.labels.fine) |> 
  as_tibble() |> 
  mutate(data_driven_consensus = reference_annotation_to_consensus(azimuth_input = azimuth_predicted.celltype.l2, monaco_input = monaco_first.labels.fine, blueprint_input = blueprint_first.labels.fine )) 


data_tbl = 
  data_tbl |> 
  left_join(annotation_consensus, copy = TRUE)

output_parquet <- "/vast/projects/mangiola_immune_map/PostDoc/CuratedAtlasQueryR/dev/consensus_output.parquet"

con_write <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")

# Use DuckDB's COPY TO command to write the data back to Parquet
# We need to execute a SQL command using dbExecute()
copy_query <- paste0("
  COPY (
    SELECT *
    FROM (
      ", dbplyr::sql_render(data_tbl), "
    )
  ) TO '", output_parquet, "' (FORMAT PARQUET);
")

# Execute the COPY command
dbExecute(con_write, copy_query)

# Disconnect from the database
dbDisconnect(con_write, shutdown = TRUE)

