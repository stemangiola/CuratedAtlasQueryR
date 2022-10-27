library(RSQLite)
library(DBI)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
# Path to the metadata file used as input
metadata_input <- args[[1]]
# Path to the sqlite file which this script will output
output_file <- args[[2]]

metadata <- readRDS(metadata_input)

con <- dbConnect(SQLite(), dbname=output_file)
dbWriteTable(con, "metadata", metadata)
dbDisconnect(con)
