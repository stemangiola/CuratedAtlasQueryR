library(tidyverse)
library(forcats)
library(HCAquery)
library(dittoSeq)
library(sccomp)

source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel")

cell_metadata_with_harmonised_annotation = readRDS("dev/cell_metadata_with_harmonised_annotation.rds")

data_for_plot_1 = 
  cell_metadata_with_harmonised_annotation |> 
  
  left_join(
    get_metadata("dev/metadata.SQLite") |> 
      select(.cell, is_primary_data.y, name, cell_type, file_id, assay) |> 
      as_tibble()
  ) 

# - Number of datasets per tissue
plot_count_dataset = 
  data_for_plot_1 |> 
  distinct(file_id, tissue_harmonised) |> 
  count(tissue_harmonised, name = "Number of datasets") |> 
  ggplot(aes(fct_reorder(tissue_harmonised, desc(`Number of datasets`)), `Number of datasets`)) +
  geom_bar(stat = "identity") +
  xlab("Tissue") +
  ylab("Number of datasets (log10)") +
  scale_y_log10() +
  theme_multipanel +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# - Number of samples per tissue
plot_sample_dataset = 
  data_for_plot_1 |> 
  distinct(.sample, tissue_harmonised) |> 
  count(tissue_harmonised, name = "Number of samples") |> 
  ggplot(aes(fct_reorder(tissue_harmonised, desc(`Number of samples`)), `Number of samples`)) +
  geom_bar(stat = "identity") +
  xlab("Tissue") +
  ylab("Number of samples (log10)") +
  scale_y_log10() +
  theme_multipanel +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# - Histogram of cells per sample
plot_cell_dataset = 
  data_for_plot_1 |> 
  count(.sample, assay) |> 
  ggplot(aes(n)) +
  geom_histogram(aes(fill=assay), bins = 100) +
  scale_fill_manual(values = dittoSeq::dittoColors()) +
  xlab("Number of cells in sample (log10)") +
  ylab("Count instances") +
  scale_x_log10() +
  theme_multipanel 

# - Immune proportion per tissue
data_for_immune_proportion = 
  cell_metadata_with_harmonised_annotation |> 
  
  left_join(
    get_metadata("dev/metadata.SQLite") |> 
      select(.cell, is_primary_data.y, name, cell_type, file_id) |> 
      as_tibble()
  ) |> 
  
  # # Filter only whole tissue
  # filter(
  #   !name |> str_detect(regex('immune', ignore_case = T)) | 
  #     tissue_harmonised %in% c("blood", "lymph node", "bone") |
  #     is_primary_data.y == "PRIMARY"
  # ) |> 
  
  # Filter Immune enriched dataset
  filter(file_id != "e756c34a-abe7-4822-9a35-55ef12270247") |> 
  filter(file_id != "ca4a7d56-739b-4e3c-8ecd-28704914cc14") |> 
  filter(file_id != "59dfc135-19c1-4380-a9e8-958908273756" | tissue_harmonised != "intestine") |> 
  
  # nest(data = -c(.sample, tissue_harmonised)) |> 
  # filter(map_int(data, ~ .x |> filter(cell_type_harmonised == "non_immune") |> nrow()) > 0 | tissue_harmonised %in% c("blood", "lymph node", "bone")) |> 
  # unnest(data) |> 
  
  mutate(is_immune = cell_type_harmonised!="non_immune") |> 
  
  # Fix hematopoietic misclassificsation
  mutate(is_immune = if_else(!is_immune & cell_type |> str_detect("hematopoietic"), TRUE, is_immune)) |> 
  
  # Filter out
  filter(!cell_type |> str_detect("erythrocyte")) |> 
  filter(!cell_type |> str_detect("platelet")) 

data_for_immune_proportion_count = 
  data_for_immune_proportion |> 
  
  # Stats
  count(.sample, tissue_harmonised, is_immune, file_id) |> 
  with_groups(.sample, ~ .x |> mutate(proportion = n/sum(n), sum = sum(n))) |> 
  filter(is_immune) |> 
  with_groups(tissue_harmonised, ~ .x |> mutate( median_proportion = mean(proportion)))

dropLeadingZero <- function(l){  stringr::str_replace(l, '0(?=.)', '') }
S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
IS_sqrt <- function(x){x^2*sign(x)}
S_sqrt_trans <- function() scales::trans_new("S_sqrt",S_sqrt,IS_sqrt)


plot_immune_proportion_dataset = 
  data_for_immune_proportion_count |> 
  ggplot(aes(fct_reorder(tissue_harmonised, desc(median_proportion)), proportion)) +
  geom_point(aes(size = sum, color=file_id)) +
  guides(color="none") +
  scale_size(trans = "log10", range = c(0.1, 2.5), limits = c(1000, 10000)) +
  scale_color_manual(values = dittoSeq::dittoColors()) +
  scale_y_continuous(trans=S_sqrt_trans(), labels = dropLeadingZero) +
  theme_multipanel +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# - scatter plot of abundance vs variability per tissue


# - Confidence class per cell type
# - 

# Study annotation
res =
  data_for_immune_proportion |>
  mutate(is_immune = as.character(is_immune)) |> 
  filter(tissue_harmonised %in% c("blood", "intestine")) |> 
  sccomp_glm(
    formula_composition = ~ 0 + tissue_harmonised,
    formula_variability = ~ 1,
    .sample, is_immune, 
    check_outliers = FALSE,
    approximate_posterior_inference = FALSE,
    #contrasts = c("typecancer - typehealthy", "typehealthy - typecancer"),
    cores = 20,
    mcmc_seed = 42, verbose = T
  )


res |> plot_summary()

cell_metadata_with_harmonised_annotation |> 
  

  
  mutate(is_immune = cell_type_harmonised!="non_immune") |> 
  
  # Stats
  count(.sample, tissue_harmonised, is_immune) |> 
  with_groups(.sample, ~ .x |> mutate(proportion = n/sum(n))) |> 
  filter(tissue_harmonised=="heart" & proportion > 0.75)
  