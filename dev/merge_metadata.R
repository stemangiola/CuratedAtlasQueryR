suppressPackageStartupMessages({
	library(zellkonverter)
	library(SingleCellExperiment) # load early to avoid masking dplyr::count()
	library(dplyr)
	library(cellxgenedp)
})

library(tidyverse)
#library(tidySingleCellExperiment)
library(stringr)
library(scMerge)
library(glue)
library(DelayedArray)
library(HDF5Array)
library(openssl)

# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file_paths = args[1:(length(args)-2)]
files_metadata = args[[length(args)-1]]
output_file = args[[length(args)]]

# Read metadata
# metadata =
# 	input_file_paths|>
# 	map(readRDS)

# Get which column don't have too many NAs
common_colnames =
	input_file_paths |>
	imap_dfr(
		~ .x %>%
			{print(.y); (.)} |>
			readRDS() |>
			colnames() |>
			as_tibble()
	) |>
	dplyr::count(value) |>
	mutate(n_datasets = length(input_file_paths)) |>
	filter(n > (n_datasets / 2)) |>
	pull(value)

print(common_colnames)

# Get all metadata
input_file_paths  |>

	# Select core columns
	imap(~ .x %>%
			 	{print(.y); (.)} |>
			 	readRDS() |>
			 	select(.cell, .sample, .sample_name, one_of(common_colnames)) |>
			 	mutate_if(is.factor, as.character)
					) |>
	bind_rows() |>

	unite("file_id_db", c(.sample, cell_type), remove = FALSE) |>
	mutate(file_id_db = file_id_db |> md5()) |>

	# Curate tissue
	left_join(
		read_csv("tissue_label_curated.csv"),
		by="tissue"
	) |>

	# Curate age
	 mutate(time_unit = case_when(
		development_stage |> str_detect("year") ~ "year",
		development_stage |> str_detect("month") ~ "month",
		development_stage |> str_detect("week") ~ "week",
		development_stage |> str_detect("day") ~ "day"
	)) |>
	tidyr::extract(development_stage, "age", "^([0-9]+).*", remove = FALSE, convert = TRUE) |>
	mutate(age = case_when(
		development_stage |> str_detect("third decade") ~ 30,
		development_stage |> str_detect("fifth decade") ~ 50,
		development_stage |> str_detect( "eighth decade") ~ 80,
		development_stage |> str_detect( "fourth decade") ~ 40,
		development_stage |> str_detect( "seventh decade") ~ 70,
		development_stage |> str_detect( "ninth decade") ~ 90,
		development_stage |> str_detect( "tenth decade") ~ 100,
		development_stage |> str_detect( "human aged") ~ 80,
		development_stage |> str_detect( "sixth decade") ~ 60,
		development_stage |> str_detect( "human middle") ~ 45,
		development_stage |> str_detect( "65-79") ~ 72,
		development_stage |> str_detect( "human adult") ~ 30,
		development_stage |> str_detect( "human late adulthood") ~ 50,
		development_stage |> str_detect( "mature stage") ~ 50,
		development_stage |> str_detect( "adolescent stage") ~ 16,
		development_stage |> str_detect( "young adult stage") ~ 25,
		development_stage |> str_detect( "child stage") ~ 10,
		development_stage |> str_detect( "infant stage") ~ 2,
		development_stage |> str_detect( "newborn human stage") ~ 1,
		development_stage |> str_detect( "human early adulthood stage") ~ 25,
		development_stage |> str_detect( "post-juvenile adult stage") ~ 35,
		development_stage |> str_detect( "early adult stage") ~ 25,
		development_stage |> str_detect( "immature stage") ~ 12,

		# Carnegie
		development_stage |> str_detect( "Carnegie stage 17") ~ 41,
		development_stage |> str_detect( "Carnegie stage 23") ~ 53,
		development_stage |> str_detect( "Carnegie stage 18") ~ 44,
		development_stage |> str_detect( "Carnegie stage 19") ~ 46,
		development_stage |> str_detect( "Carnegie stage 22") ~ 56,

		TRUE ~ as.numeric(age)
	)) |>

	mutate(age = case_when(
		development_stage |> str_detect( "organogenesis stage") ~ 60,
		development_stage |> str_detect( "embryonic human stage" )~ 30,
		TRUE ~ age
	)) |>

	mutate(age = case_when(
		development_stage |> str_detect( "fetal stage" )~ 15,
		TRUE ~ age
	)) |>
	mutate(time_unit = case_when(
		!is.na(time_unit) ~ time_unit,
		development_stage |> str_detect( "organogenesis stage") ~ "day",
		development_stage |> str_detect( "embryonic human stage" )~ "day",
		development_stage |> str_detect( "fetal stage" ) ~ "week",
		development_stage |> str_detect( "decade" ) ~ "year",
		development_stage |> str_detect( "aged stage" ) ~ "year",
		development_stage |> str_detect( "adulthood" ) ~ "year",
		development_stage |> str_detect( "adult|mature|adolescent|child|newborn|immature|infant" ) ~ "year",
		development_stage |> str_detect( "Carnegie" ) ~ "day"

	)) |>

	mutate(age_days = case_when(
		time_unit == "day" ~ age,
		time_unit == "week" ~ age * 7,
		time_unit == "month" ~ age * 30,
		time_unit == "year" ~ age * 365
	)) |>
	select(-age, -time_unit) |>

	#mutate_if(is.character, as.factor) |>
	# Add files metadata
	left_join(readRDS(files_metadata) |> select_if(function(x) !is.list(x)), by="file_id") |>


	# Make original cell
	mutate(.cell = .cell |> str_remove(.sample) |> str_remove("_$")) |>

	distinct() |>
	saveRDS(output_file)



