# This file creates the whole makefile

# Create makeflow file for doublet identification
# makeflow -T slurm -J 200 analysis/empty_droplet_identification/empty_droplet_identification.makeflow


# commands
commands = c()
tab = "\t"

# Read arguments
args = commandArgs(trailingOnly=TRUE)
root_directory = args[[1]]

raw_data_directory = glue("{root_directory}/raw_data")
input_directory = args[[3]]

metadata = glue("{root_directory}/metadata.rds") |> readRDS()

reference_azimuth_path = args[[6]]

# Check modality

# Create dir
result_directory |> dir.create( showWarnings = FALSE, recursive = TRUE)

renv::activate(project = code_directory)

library(tidyverse)
library(glue)
library(here)
library(stringr)
library(Seurat)
library(tidyseurat)

R_code_directory = glue("{code_directory}/R")


# Input demultiplexed THOSE FILES MUST EXIST!
input_directory_demultiplexed = input_directory
input_files_demultiplexed = dir(input_directory_demultiplexed, pattern = ".rds")
input_path_demultiplexed = glue("{input_directory_demultiplexed}/{input_files_demultiplexed}")
samples = input_files_demultiplexed |> str_remove(".rds")




# >>> Parse data
parsed_data_directory = glue("{root_directory}/parsed_data")

# Create input
commands =
  commands |>
  c(
    glue("CATEGORY={suffix}\nMEMORY=30024\nCORES=1\nWALL_TIME=10000"),
    glue("{output_paths_annotation_label_transfer}:{input_path_demultiplexed} {output_path_empty_droplets}\n{tab}Rscript {R_code_directory}/run{suffix}.R {code_directory} {input_path_demultiplexed} {output_path_empty_droplets} {reference_azimuth_path} {output_paths_annotation_label_transfer}")
  )




# >>> ALIVE CELLS
suffix = "__alive_identification"
output_directory = glue("{result_directory}/preprocessing_results/alive_identification")
output_path_alive =   glue("{output_directory}/{samples}{suffix}_output.rds")

# Create input
commands =
  commands |> c(
    glue("CATEGORY={suffix}\nMEMORY=10024\nCORES=1\nWALL_TIME=30000"),
    glue("{output_path_alive}:{input_path_demultiplexed} {output_path_empty_droplets} {output_paths_annotation_label_transfer}\n{tab}Rscript {R_code_directory}/run{suffix}.R {code_directory} {input_path_demultiplexed} {output_path_empty_droplets} {output_paths_annotation_label_transfer} {output_path_alive}")
  )



# >>> DOUBLET IDENTIFICATION
suffix = "__doublet_identification"
output_directory = glue("{result_directory}/preprocessing_results/doublet_identification")
output_path_doublet_identification =   glue("{output_directory}/{samples}{suffix}_output.rds")

# Create input
commands =
  commands |> c(
    glue("CATEGORY={suffix}\nMEMORY=10024\nCORES=2\nWALL_TIME=30000"),
    glue("{output_path_doublet_identification}:{input_path_demultiplexed} {output_path_empty_droplets} {output_path_alive} {output_paths_annotation_label_transfer}\n{tab}Rscript {R_code_directory}/run{suffix}.R {code_directory} {input_path_demultiplexed} {output_path_empty_droplets} {output_path_alive} {output_paths_annotation_label_transfer} {output_path_doublet_identification}")

  )





# >>> VARIABLE GENE SELECTION, BY BATCH AND BROAD CELL TYPE
metadata = readRDS(metadata_path)
suffix = "__variable_gene_identification"
output_directory = glue("{result_directory}/preprocessing_results/variable_gene_identification")


commands_variable_gene =
  tibble(sample = samples, input_path_demultiplexed, output_path_empty_droplets, output_path_alive, output_path_doublet_identification, output_paths_annotation_label_transfer) |>

  # Add batch
  left_join(
    readRDS(metadata_path) |>
      distinct(sample, batch),
    by="sample"
  ) |>

  # Create list of files
  with_groups(
    batch,
    ~ summarise(.x, input_files = paste(
      c(input_path_demultiplexed, output_path_empty_droplets, output_path_alive, output_path_doublet_identification, output_paths_annotation_label_transfer),
      collapse = " "
    ))
  ) |>

  # Output file
  mutate(output_file = glue("{output_directory}/{batch}{suffix}_output.rds") ) |>

  # create command
  mutate(command =  glue("{output_file}:{input_files}\n{tab}Rscript {R_code_directory}/run{suffix}.R {code_directory} {input_files} predicted.celltype.l1 {output_file}"))

output_path_marged_variable_genes =   glue("{output_directory}/merged_{suffix}_output.rds")

commands =
  commands |> c(
    glue("CATEGORY={suffix}\nMEMORY=50024\nCORES=11\nWALL_TIME=30000"),
    commands_variable_gene  |> pull(command) |> unlist(),
    glue("{output_path_marged_variable_genes}:{paste(commands_variable_gene |> pull(output_file), collapse=\" \")}\n{tab}Rscript {R_code_directory}/run{suffix}_merge.R {code_directory} {paste(commands_variable_gene |> pull(output_file), collapse=\" \")} {output_path_marged_variable_genes}")
  )




# >>> NORMALISATION
suffix = "__non_batch_variation_removal"
output_directory = glue("{result_directory}/preprocessing_results/non_batch_variation_removal")
output_path_non_batch_variation_removal =   glue("{output_directory}/{samples}{suffix}_output.rds")

# Create input
commands =
  commands |> c(
    glue("CATEGORY={suffix}\nMEMORY=40024\nCORES=1\nWALL_TIME=30000"),
    glue("{output_path_non_batch_variation_removal}:{input_path_demultiplexed} {output_path_empty_droplets} {output_path_marged_variable_genes}\n{tab}Rscript {R_code_directory}/run{suffix}.R {code_directory} {input_path_demultiplexed} {output_path_empty_droplets} {output_path_marged_variable_genes} {output_path_non_batch_variation_removal}")
  )



# >>> WRITE PREPROCESSING RESULT
suffix = "__preprocessing_output"
output_directory = glue("{result_directory}/preprocessing_results/preprocessing_output")
output_path_preprocessing_results =   glue("{output_directory}/{samples}{suffix}_output.rds")

# Create input
commands =
  commands |> c(
    glue("CATEGORY={suffix}\nMEMORY=10024\nCORES=2\nWALL_TIME=30000"),
    glue("{output_path_preprocessing_results}:{output_path_non_batch_variation_removal} {output_path_alive} {output_paths_annotation_label_transfer} {output_path_doublet_identification}\n{tab}Rscript {R_code_directory}/run{suffix}.R {code_directory} {output_path_non_batch_variation_removal} {output_path_alive} {output_paths_annotation_label_transfer} {output_path_doublet_identification} {output_path_preprocessing_results}")

  )





# >>> PSEUDOBULK PREPROCESSING
suffix = "__pseudobulk_preprocessing"
output_directory = glue("{result_directory}/preprocessing_results/pseudobulk_preprocessing")
output_path_pseudobulk_preprocessing_sample_cell_type =   glue("{output_directory}/pseudobulk_preprocessing_sample_cell_type_output.rds")
output_path_pseudobulk_preprocessing_sample =   glue("{output_directory}/pseudobulk_preprocessing_sample_output.rds")

# Create input
commands =
  commands |> c(
    glue("CATEGORY={suffix}\nMEMORY=50024\nCORES=11\nWALL_TIME=30000"),
    glue("{output_path_pseudobulk_preprocessing_sample_cell_type} {output_path_pseudobulk_preprocessing_sample}:{paste(output_path_preprocessing_results, collapse=\" \")}\n{tab}Rscript {R_code_directory}/run{suffix}.R {code_directory} {paste(output_path_preprocessing_results, collapse=\" \")} {output_path_pseudobulk_preprocessing_sample_cell_type} {output_path_pseudobulk_preprocessing_sample}")

  )








if(modality %in% c("fast_pipeline", "complete_pipeline")){



  # >>> COMMUNICATION FAST

  cores = 4

  # Output
  output_directory = glue("{result_directory}/fast_pipeline_results/communication")
  output_path_ligand_receptor_count =   glue("{output_directory}/{samples}_ligand_receptor_count_output.rds")


  # Create input
  commands =
    commands |>
    c(
      "CATEGORY=communication1\nMEMORY=10024\nCORES=4\nWALL_TIME=9000",
      glue("{output_path_ligand_receptor_count}:{output_path_preprocessing_results}\n{tab}Rscript {R_code_directory}/run__ligand_receptor_count.R {code_directory} {output_path_preprocessing_results} {output_path_ligand_receptor_count}")

    )

  # CellChat results
  #output_path_communication_hypothesis_testing =   glue("{output_directory}/communication_output.rds")
  output_path_plot_overall =  glue("{output_directory}/plot_communication_overall.pdf")
  output_path_plot_heatmap =  glue("{output_directory}/plot_communication_heatmap.pdf")
  output_path_plot_circle =  glue("{output_directory}/plot_communication_circle.pdf")
  output_path_values_communication =  glue("{output_directory}/values_communication.rds")

  commands =
    commands |>
    c(
      "CATEGORY=communication2\nMEMORY=50024\nCORES=2\nWALL_TIME=18000",
      glue("{output_path_plot_overall} {output_path_plot_heatmap} {output_path_plot_circle} {output_path_values_communication}:{paste(output_path_ligand_receptor_count, collapse=\" \")}\n{tab}Rscript {R_code_directory}/run__differential_communication.R {code_directory} {paste(output_path_ligand_receptor_count, collapse=\" \")} {result_directory}/metadata.rds {output_path_plot_overall} {output_path_plot_heatmap} {output_path_plot_circle} {output_path_values_communication}")

    )






  # >>> PSEUDO BULK DIFFERENTIAL GENE TRANSCRIPT ABUNDANCE
  suffix = "__differential_transcript_abundance_fast"
  output_directory = glue("{result_directory}/fast_pipeline_results/differential_transcript_abundance")
  output_path_differential_transcript_abundance_fast=   glue("{output_directory}/differential_transcript_abundance_output.rds")
  output_path_plot_densities =   glue("{output_directory}/plot_densities.rds")
  output_path_plot_significant =   glue("{output_directory}/plot_significant.rds")

  # Create input
  commands =
    commands |> c(
      glue("CATEGORY={suffix}\nMEMORY=70024\nCORES=1\nWALL_TIME=30000"),
      glue("{output_path_differential_transcript_abundance_fast} {output_path_plot_densities} {output_path_plot_significant}:{output_path_pseudobulk_preprocessing_sample_cell_type} {metadata_path}\n{tab}Rscript {R_code_directory}/run__differential_transcript_abundance.R {code_directory} {output_path_pseudobulk_preprocessing_sample_cell_type} {metadata_path} {output_path_differential_transcript_abundance_fast} {output_path_plot_densities} {output_path_plot_significant}")

    )



  # >>> PSEUDO BULK DIFFERENTIAL COMPOSITION
  suffix = "__differential_composition_fast"
  output_directory = glue("{result_directory}/fast_pipeline_results/differential_composition")
  output_path_differential_composition_fast = glue("{output_directory}/differential_composition_output.rds")
  output_path_plot_credible_intervals = glue("{output_directory}/plot_credible_intervals.rds")
  output_path_plot_boxplot = glue("{output_directory}/plot_boxplot.rds")


  # Create input
  commands =
    commands |> c(
      glue("CATEGORY={suffix}\nMEMORY=70024\nCORES=1\nWALL_TIME=30000"),
      glue("{output_path_differential_composition_fast} {output_path_plot_credible_intervals} {output_path_plot_boxplot}:{metadata_path} {paste(output_paths_annotation_label_transfer, collapse=\" \") }\n{tab}Rscript {R_code_directory}/run__differential_composition.R {code_directory} {metadata_path} {paste(output_paths_annotation_label_transfer, collapse=\" \") } {output_path_differential_composition_fast} {output_path_plot_credible_intervals} {output_path_plot_boxplot}")

    )



}





# >>> WRITE TO FILE
commands |>
  write_lines(glue("{result_directory}/pipeline.makeflow"))
