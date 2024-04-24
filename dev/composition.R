library(tidyverse)
library(forcats)
library(CuratedAtlasQueryR)
library(dittoSeq)
library(sccomp)
library(magrittr)
library(patchwork)

source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel")
metadata = readRDS("~/PostDoc/CuratedAtlasQueryR/dev/metadata_annotated.rds")

## from http://tr.im/hH5A


softmax <- function (x) {
	logsumexp <- function (x) {
		y = max(x)
		y + log(sum(exp(x - y)))
	}

	exp(x - logsumexp(x))
}

cell_metadata_with_harmonised_annotation = readRDS("~/PostDoc/CuratedAtlasQueryR/dev/cell_metadata_with_harmonised_annotation.rds")

data_for_plot_1 =
  cell_metadata_with_harmonised_annotation |>

  left_join(
    get_metadata() |>
      dplyr::select(.cell, is_primary_data.y, name, cell_type, file_id, assay) |>
      as_tibble()
  )


# - Number of datasets per tissue
plot_count_dataset =
  data_for_plot_1 |>
  distinct(file_id, tissue_harmonised) |>
  count(tissue_harmonised, name = "Number of datasets") |>
  ggplot(aes(fct_reorder(tissue_harmonised, dplyr::desc(`Number of datasets`)), `Number of datasets`)) +
  geom_bar(stat = "identity") +
  xlab("Tissue") +
  ylab("Number of datasets (log10)") +
  theme_multipanel +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# - Number of samples per tissue
plot_sample_dataset =
  data_for_plot_1 |>
  distinct(.sample, tissue_harmonised) |>
  count(tissue_harmonised, name = "Number of samples") |>
  ggplot(aes(fct_reorder(tissue_harmonised, dplyr::desc(`Number of samples`)), `Number of samples`)) +
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
  xlab("Number of cells in sample (log10, capped at 1e2-1e6)") +
  ylab("Count instances") +
  scale_x_log10(limits = c(10, 1e6)) +
  theme_multipanel

# - Immune proportion per tissue
data_for_immune_proportion =
  cell_metadata_with_harmonised_annotation |>

  left_join(
    # get_metadata() |>
  	metadata |>
    	dplyr::select(.cell, cell_type, file_id, assay, age_days, development_stage, sex, ethnicity) |>
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
  filter(!cell_type |> str_detect("platelet")) |>

	# Frmat sme covatriates
	mutate(assay = if_else(assay |> str_detect("10x|scRNA-seq"), "10x", assay)) |>
	mutate(ethnicity = case_when(
		ethnicity |> str_detect("Chinese|Asian") ~ "Chinese",
		ethnicity |> str_detect("African") ~ "African",
		TRUE ~ ethnicity
	)) |>

	# Fix samples with multiple assays
	unite(".sample", c(.sample , assay), remove = FALSE) |>

	# Fix groups
	unite("group", c(tissue_harmonised , file_id), remove = FALSE)


# - Confidence class per cell type
# -

# # Study annotation
# job::job({
#
# 	res_absolute =
# 		data_for_immune_proportion |>
# 		mutate(is_immune = as.character(is_immune)) |>
#
# 		# Mutate days
# 		mutate(age_days = age_days  |> scale(center = FALSE) |> as.numeric()) |>
# 		filter(development_stage!="unknown") |>
#
#
# 		sccomp_glm(
# 			formula_composition = ~ 0 + tissue_harmonised + sex + ethnicity  + age_days + assay + (tissue_harmonised | group) + (age_days | tissue_harmonised),
# 			formula_variability = ~ 0 + tissue_harmonised + sex + ethnicity ,
# 			.sample, is_immune,
# 			check_outliers = F,
# 			approximate_posterior_inference = FALSE,
# 			cores = 20,
# 			mcmc_seed = 42, verbose = T
# 		)
#
# 	res_absolute |> saveRDS("~/PostDoc/CuratedAtlasQueryR/dev/immune_non_immune_differential_composition.rds")
#
# })

res_absolute = readRDS("~/PostDoc/CuratedAtlasQueryR/dev/immune_non_immune_differential_composition.rds")

res_generated_proportions =
	res_absolute |>
	replicate_data(number_of_draws = 20) |>
	filter(is_immune=="TRUE") |>
	left_join(
		data_for_immune_proportion |>
			unite(".sample", c(.sample , assay), remove = FALSE) |>
			select(.sample, tissue_harmonised)
	) |>
	with_groups(tissue_harmonised, ~ .x |> sample_n(30, replace = T))

# -- Rank immune
dropLeadingZero <- function(l){  stringr::str_replace(l, '0(?=.)', '') }
S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
IS_sqrt <- function(x){x^2*sign(x)}
S_sqrt_trans <- function() scales::trans_new("S_sqrt",S_sqrt,IS_sqrt)

library(ggforestplot)

plot_immune_proportion_dataset =
	data_for_immune_proportion |>

	# Stats
	count(.sample, tissue_harmonised, is_immune, file_id) |>
	with_groups(.sample, ~ .x |> mutate(proportion = n/sum(n), sum = sum(n))) |>
	filter(is_immune) |>
	with_groups(tissue_harmonised, ~ .x |> mutate( median_proportion = mean(proportion))) |>

	# Add multilevel proportion medians
	left_join(
		res_generated_proportions |>
			with_groups(tissue_harmonised, ~ .x |> summarise(median_generated = median(generated_proportions, na.rm = TRUE)))
	) |>

	# Plot
	ggplot(aes( proportion, fct_reorder(tissue_harmonised, median_generated))) +
	ggforestplot::geom_stripes(odd = "#33333333", even = "#00000000") +
	geom_point(aes(size = sum, color=file_id)) +
	geom_boxplot(aes(generated_proportions, fct_reorder(tissue_harmonised, median_generated)), color="red", data =
							 	res_generated_proportions |>
							 	with_groups(tissue_harmonised, ~ .x |> mutate(median_generated = median(generated_proportions, na.rm = TRUE))),
							 fill = NA, outlier.shape = NA,
							 ) +
	guides(color="none") +
	scale_size(trans = "log10", range = c(0.1, 2.5), limits = c(1000, 10000)) +
	scale_color_manual(values = dittoSeq::dittoColors()) +
	scale_x_continuous(trans=S_sqrt_trans(), labels = dropLeadingZero) +
	xlab("Immune proportion (sqrt)") +
	ylab("Tissue") +
	theme_multipanel

coefficients_regression =
	res_absolute |>
	filter(is_immune == "TRUE") |>
	filter(!parameter %in% c("spleen", "bone", "blood", "lymph node") ) %>%
	lm( v_effect ~ c_effect, data = .) %$%
	coefficients

library(ggrepel)

median_composition =
	res_absolute |>
	filter(is_immune == "TRUE") |>
	filter(!parameter %in% c("spleen", "bone", "blood", "lymph node") ) |>
	pull(c_effect) |>
	median()

# - scatter plot of abundance vs variability per tissue
res_for_plot =
	res_absolute |>
	filter(is_immune == "TRUE") |>
	mutate(parameter = parameter |> str_remove("tissue_harmonised")) |>
	mutate(intercept = coefficients_regression[1], slope = coefficients_regression[2]) |>

	# Normalise effects
	mutate(
		v_effect = v_effect - (c_effect * slope + intercept),
		v_lower = v_lower - (c_effect * slope + intercept),
		v_upper = v_upper - (c_effect * slope + intercept)
	) |>
	mutate(
		c_effect = c_effect - median_composition,
		c_lower = c_lower - median_composition,
		c_upper = c_upper - median_composition
	) |>

	# Define significance
	mutate(
		v_significant = v_lower * v_upper > 0,
		c_significant = c_lower * c_upper > 0
	) |>

	# Define quadrants
	mutate(quadrant = case_when(
		c_effect > 0 & v_effect > 0 ~ "Hot and variable",
		c_effect > 0 & v_effect < 0 ~ "Hot and consistent",
		c_effect < 0 & v_effect > 0 ~ "Cold and variable",
		c_effect < 0 & v_effect < 0 ~ "Cold and consistent"

	)) |>

	# Limit the values
	mutate(
		v_effect = pmax(v_effect, -2),
		c_effect = pmin(c_effect, 1),
		v_lower = pmax(v_lower, -2),
		v_upper = pmax(v_upper, -2),
		c_lower = pmin(c_lower, 1),
		c_upper = pmin(c_upper, 1),
	)


res_for_plot |>
	ggplot(aes(c_effect, v_effect, label = parameter)) +
	geom_vline(xintercept = 0, linetype = "dashed") +
	geom_hline(yintercept = 0, linetype = "dashed") +
	geom_errorbar(aes(ymin = v_lower, ymax = v_upper, color = v_significant), alpha = 0.5) +
	geom_errorbar(aes(xmin = c_lower, xmax = c_upper, color = c_significant), alpha = 0.5) +
	geom_point(aes(fill = quadrant), shape = 21) +
	geom_text_repel(data =  res_for_plot |> filter(!v_significant & !c_significant)) +
	geom_label_repel(data = res_for_plot |> filter(v_significant | c_significant)) +
	xlab("Composition") +
	ylab("Variability") +
	scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
		scale_fill_brewer(palette = "Set1") +
	theme_multipanel



# Stats
cell_metadata_with_harmonised_annotation |>

  mutate(is_immune = cell_type_harmonised!="non_immune") |>

  # Stats
  count(.sample, tissue_harmonised, is_immune) |>
  with_groups(.sample, ~ .x |> mutate(proportion = n/sum(n))) |>
  filter(tissue_harmonised=="heart" & proportion > 0.75)


# # Relative
# data_for_immune_proportion_relative =
# 	cell_metadata_with_harmonised_annotation |>
#
# 	left_join(
# 		#get_metadata() |>
# 		metadata |>
# 			dplyr::select(.cell, cell_type, file_id, assay, age_days, development_stage, sex, ethnicity) |>
# 			as_tibble()
# 	) |>
#
# 	# Fix hematopoietic misclassification
# 	mutate(cell_type_harmonised = if_else(cell_type_harmonised=="non_immune" & cell_type |> str_detect("hematopoietic"), "stem", cell_type_harmonised)) |>
#
# 	# Filter out
# 	filter(!cell_type |> str_detect("erythrocyte")) |>
# 	filter(!cell_type |> str_detect("platelet")) |>
#
# 	mutate(is_immune = cell_type_harmonised!="non_immune") |>
#
# 	# Filter only immune
# 	filter(is_immune ) |>
#
# 	# Frmat sme covatriates
# 	mutate(assay = if_else(assay |> str_detect("10x|scRNA-seq"), "10x", assay)) |>
# 	mutate(ethnicity = case_when(
# 		ethnicity |> str_detect("Chinese|Asian") ~ "Chinese",
# 		ethnicity |> str_detect("African") ~ "African",
# 		TRUE ~ ethnicity
# 	)) |>
#
# 	filter(development_stage!="unknown") |>
#
# 	# Fix samples with multiple assays
# 	unite(".sample", c(.sample , assay), remove = FALSE) |>
#
# 	# Fix groups
# 	unite("group", c(tissue_harmonised , file_id), remove = FALSE)
#
# data_for_immune_proportion_relative |> saveRDS("dev/data_for_immune_proportion_relative.rds")

data_for_immune_proportion_relative = readRDS("~/PostDoc/CuratedAtlasQueryR/dev/data_for_immune_proportion_relative.rds")

# Analysis of counts relative
data_for_immune_proportion_relative |>
	distinct(.sample, cell_type_harmonised) |>
	count(.sample) |>
	add_count(n) |>
	arrange(n)

# job::job({
# 	res_relative =
# 		data_for_immune_proportion_relative |>
#
# 		# Scale days
# 		mutate(age_days = age_days  |> scale(center = FALSE) |> as.numeric()) |>
#
# 		# Estimate
# 		sccomp_glm(
# 			formula_composition = ~ 0 + tissue_harmonised + sex + ethnicity  + age_days + assay + (tissue_harmonised | group) + (age_days | tissue_harmonised),
# 			formula_variability = ~ 0 + tissue_harmonised + sex + ethnicity,
# 			.sample, cell_type_harmonised,
# 			check_outliers = F,
# 			approximate_posterior_inference = FALSE,
# 			cores = 20,
# 			mcmc_seed = 42, verbose = T
# 		)
#
# 	res_relative |> saveRDS("dev/immune_non_immune_differential_composition_relative_5.rds")
# })


res_relative = readRDS("~/PostDoc/CuratedAtlasQueryR/dev/immune_non_immune_differential_composition_relative_5.rds")

# job::job({
# 	readRDS("~/PostDoc/CuratedAtlasQueryR/dev/immune_non_immune_differential_composition_relative_5.rds") |>
# 		remove_unwanted_variation(~ 0 + tissue_harmonised) |>
# 		saveRDS("~/PostDoc/CuratedAtlasQueryR/dev/proportions_tissue_adjusted_5.rds")
# })
# job::job({ readRDS("~/PostDoc/CuratedAtlasQueryR/dev/immune_non_immune_differential_composition_relative_4.rds") |> remove_unwanted_variation(~ age_days) })


proportions_tissue_replicate = readRDS("~/PostDoc/CuratedAtlasQueryR/dev/immune_non_immune_differential_composition_relative_5.rds") |> replicate_data (~ 0 + tissue_harmonised)

 # PCA
library(tidybulk)
observed_proportion_PCA =
	data_for_immune_proportion_relative |>
	add_count(tissue_harmonised) |>
	filter(n > 5) |>
	select(-n) |>
	count(.sample, cell_type_harmonised, tissue_harmonised, assay, sex, file_id) |>
	with_groups(.sample, ~ .x |> mutate(observed_proportion = n/sum(n))) |>
	complete(nesting(.sample, tissue_harmonised, assay, sex, file_id), cell_type_harmonised, fill = list(observed_proportion = 0)) |>
	reduce_dimensions(.sample , cell_type_harmonised, observed_proportion, method="tSNE", action="get") |>
	ggplot(aes(tSNE1, tSNE2)) +
	geom_point(aes(color = tissue_harmonised)) +
	guides(color="none")  +
	theme_multipanel


adjusted_proportion_PCA =
	proportions_tissue_adjusted |>
	left_join(
		data_for_immune_proportion_relative |>
			distinct(.sample, tissue_harmonised, assay, file_id, sex, ethnicity)
	) |>
	add_count(tissue_harmonised) |>
	filter(n > 5) |>
	reduce_dimensions(.sample , cell_type_harmonised, adjusted_proportion, method="tSNE", action="get") |>
	ggplot(aes(tSNE1, tSNE2)) +
	geom_point(aes(color = tissue_harmonised))  +
	#ggdensity::geom_hdr_lines(aes(color = tissue_harmonised)) +
	theme_multipanel


library(ggforce)
library(ggpubr)

circle_plot = function(res){

	logsumexp <- function (x) {
		y = max(x)
		y + log(sum(exp(x - y)))
	}
	softmax <- function (x) {

		exp(x - logsumexp(x))
	}

	res_relative_for_plot =
		res |>
		filter(!parameter |> str_detect("group___")) |>

		# Cell type abundance
		with_groups(tissue_harmonised, ~ .x |>  mutate(proportion = softmax(c_effect))) |>
		with_groups(cell_type_harmonised, ~ .x |>  mutate(cell_type_mean_abundance = mean(proportion))) |>

		# Filter for visualisation
		filter(!cell_type_harmonised %in% c("non_immune", "immune_unclassified")) |>

		# Tissue diversity
		with_groups(tissue_harmonised, ~ .x |>  mutate(inter_type_diversity = sd(c_effect))) |>

		# First rank
		with_groups(cell_type_harmonised, ~ .x |> arrange(desc(c_effect)) |>  mutate(rank = 1:n())) |>

		# Cap
		mutate(c_effect = c_effect |> pmax(-5) |> pmin(5))

	inter_type_diversity_plot =
		res_relative_for_plot |>
		distinct(inter_type_diversity, tissue_harmonised) |>
		ggplot(aes(inter_type_diversity, fct_reorder(tissue_harmonised, inter_type_diversity))) +
		geom_bar(stat = "identity") +
		scale_x_reverse() +
		theme_multipanel

	cell_type_mean_abundance_plot =
		res_relative_for_plot |>
		distinct(cell_type_mean_abundance, cell_type_harmonised) |>
		ggplot(aes(fct_reorder(cell_type_harmonised, dplyr::desc(cell_type_mean_abundance)), cell_type_mean_abundance)) +
		geom_bar(stat = "identity") +
		scale_y_continuous(position = "right") +
		theme_multipanel +
		theme(
			axis.text.x = element_blank(),
			axis.title.x = element_blank(),
			axis.ticks.x = element_blank(),
			axis.title.y = element_blank(),
		)

	circle_plot =
		res_relative_for_plot |>
		arrange(rank==1) |>
		ggplot() +
		geom_point(aes(
			fct_reorder(cell_type_harmonised, dplyr::desc(cell_type_mean_abundance)),
			fct_reorder(tissue_harmonised, inter_type_diversity) ,
			fill = rank, size = c_effect, stroke=rank==1), shape=21
		) +
		scale_fill_viridis_c(direction = -1) +
		theme_multipanel +
		theme(
			axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5),
			axis.text.y = element_blank(),
			axis.title.y = element_blank(),
			axis.ticks.y = element_blank()
		)


	plot_spacer() +
		cell_type_mean_abundance_plot +
		inter_type_diversity_plot +
		circle_plot +
		plot_layout(guides = 'collect', height = c(1,5), width = c(1, 5)) &
		theme( plot.margin = margin(0, 0, 0, 0, "pt"), legend.position = "bottom")

}

res_relative |>
	filter(covariate == "tissue_harmonised") |>
	mutate(tissue_harmonised = parameter |> str_remove("tissue_harmonised")) |>
	circle_plot()

res_relative_plot = res_relative |> plot_summary()
res_relative_plot |> saveRDS("dev/immune_non_immune_differential_composition_relative_summary_plot.rds")



# Track of immune system in life

# # Study age absolute
# job::job({
#
# 	data_for_immune_proportion |>
#
# 		# Drop only-immune organs
# 		filter(!tissue_harmonised %in% c("blood", "lymph node", "spleen", "bone")) |>
# 		mutate(is_immune = as.character(is_immune)) |>
#
# 		# Mutate days
# 		mutate(age_days = age_days  |> scale(center = FALSE) |> as.numeric()) |>
# 		filter(development_stage!="unknown") |>
#
# 		sccomp_glm(
# 			formula_composition = ~ age_days + tissue_harmonised + sex + ethnicity  + assay + (1 | group) + (age_days | tissue_harmonised),
# 			formula_variability = ~ age_days + tissue_harmonised ,
# 			.sample, is_immune,
# 			check_outliers = F,
# 			approximate_posterior_inference = FALSE,
# 			cores = 20,
# 			mcmc_seed = 42, verbose = T
# 		) |>
#
# 		saveRDS("~/PostDoc/CuratedAtlasQueryR/dev/immune_non_immune_differential_composition_age2.rds")
#
# })

# Study age relative
job::job({

		data_for_immune_proportion_relative |>

		# Scale days
		mutate(age_days = age_days  |> scale(center = FALSE) |> as.numeric()) |>

		# Estimate
		sccomp_glm(
			formula_composition = ~ age_days + tissue_harmonised + sex + ethnicity  + assay + (1 | group) + (age_days | tissue_harmonised),
			formula_variability = ~ age_days + tissue_harmonised ,
			.sample, cell_type_harmonised,
			check_outliers = F,
			approximate_posterior_inference = FALSE,
			cores = 20,
			mcmc_seed = 42, verbose = T
		) |> saveRDS("~/PostDoc/CuratedAtlasQueryR/dev/immune_non_immune_differential_composition_age_relative.rds")
})

# job::job({
# 	readRDS("~/PostDoc/CuratedAtlasQueryR/dev/immune_non_immune_differential_composition_age_relative.rds") |>
# 		remove_unwanted_variation(~ age_days) |>
# 		saveRDS("~/PostDoc/CuratedAtlasQueryR/dev/proportions_age_adjusted_relative.rds")
# 	})
# job::job({
# 	readRDS("~/PostDoc/CuratedAtlasQueryR/dev/immune_non_immune_differential_composition_age2.rds") |>
# 		remove_unwanted_variation(~ age_days) |>
# 		saveRDS("~/PostDoc/CuratedAtlasQueryR/dev/proportions_age_adjusted_absolute.rds")
# 	})
proportions_age5 = readRDS("~/PostDoc/CuratedAtlasQueryR/dev/proportions_age_adjusted_relative.rds")
proportions_age_absolute = readRDS("~/PostDoc/CuratedAtlasQueryR/dev/proportions_age_adjusted_absolute.rds")

proportions_age_absolute |>
	left_join(
		data_for_immune_proportion |>
			tidybulk::pivot_sample(.sample)
	) |>

	filter(development_stage!="unknown") |>
	filter(is_immune == "TRUE") |>
	#filter(tissue_harmonised != "blood") |>
	# Fix samples with multiple assays
	unite(".sample", c(.sample , assay), remove = FALSE) |>

	# Fix groups
	unite("group", c(tissue_harmonised , file_id), remove = FALSE)  |>
	ggplot(aes(age_days, adjusted_proportion)) +
	geom_point(aes(color = tissue_harmonised), shape = ".") +
	geom_smooth(method="lm") +
	facet_wrap(~ is_immune, ncol=9) +
	scale_y_continuous(trans=S_sqrt_trans(), labels = dropLeadingZero) +
	theme_multipanel

proportions_age5 |>
	left_join(
		data_for_immune_proportion_relative |>
			tidybulk::pivot_sample(.sample)
	) |>

	filter(development_stage!="unknown") |>
	filter(cell_type_harmonised!="immune_unclassified") |>
	#filter(tissue_harmonised != "blood") |>
	# Fix samples with multiple assays
	unite(".sample", c(.sample , assay), remove = FALSE) |>

	# Fix groups
	unite("group", c(tissue_harmonised , file_id), remove = FALSE)  |>
	ggplot(aes(age_days, adjusted_proportion)) +
	geom_point(aes(color = tissue_harmonised), shape = ".") +
	geom_smooth(method="lm") +
	facet_wrap(~ cell_type_harmonised, ncol=9) +
	scale_y_continuous(trans=S_sqrt_trans(), labels = dropLeadingZero) +
	theme_multipanel

# Ethnicicy
res_absolute |>
	filter(covariate == "ethnicity") |>
	filter(!parameter |> str_detect("unknown")) |>
	filter(is_immune == "TRUE")


# Technology bias
res_relative |>
	filter(covariate == "assay") |>
	mutate(tissue_harmonised = parameter |> str_remove("assay")) |>
	circle_plot()

