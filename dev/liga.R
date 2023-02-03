library(tidyverse)
library(glue)
library(HCAquery)
library(brms)
library(MASS)
library(magrittr)
library(stringr)
library(purrr)
library(tidybayes)
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel")



root_directory = "/vast/projects/RCP/human_cell_atlas"
ligand_receptor_data_directory = glue("{root_directory}/ligand_receptor_data_0.2")
# x =
# 	dir(ligand_receptor_data_directory, full.names = T) |>
# 	map_dfr(~ {
# 		readRDS(.x) |>
# 			dplyr::select(-any_of("data")) |>
# 			mutate(file = .x)
#
# 	})
#
#
# x = x |>
# 	bind_cols(x |> pull(sample))
#
#
# x |> saveRDS("~/PostDoc/HCAquery/dev/ligand_receptor_count_all_0.2.rds")

x = readRDS("~/PostDoc/HCAquery/dev/ligand_receptor_count_all_0.2.rds")


metadata_df =
	readRDS("/vast/projects/RCP/human_cell_atlas/metadata_annotated_0.2.rds" ) |>
	count(.sample, sex, age_days, assay, ethnicity, tissue_harmonised, file_id, disease)

plotly::ggplotly(
	x |>
	left_join(
		metadata_df |>
			count(.sample, sex, age_days, assay, ethnicity, tissue_harmonised, file_id),
		by = c("sample" = ".sample")
	) |>
	mutate(lymphoid_organ = tissue_harmonised %in% c("blood", "lymph node", "bone", "spleen", "thymus")) |>

	filter(file_id |> str_detect("973fa", negate = TRUE)) |>
		mutate(how_many_cell_types = map_int(cell_cell_weight, ~ .x |> distinct(cell_from) |> nrow())) |>
		mutate( tot_interactions= tot_interactions/ how_many_cell_types) |>
	ggplot(aes(age_days, tot_interactions, color = tissue_harmonised, label=file_id)) +
	geom_point() +
	facet_grid(DB~lymphoid_organ, scales = "free_x") +
	theme(axis.text.x = element_text(angle = 30)) +
	geom_smooth(method = "glm.nb", se = FALSE) +
	#scale_color_continuous(trans="log") +
	scale_y_sqrt()
)




# DC communication


cdc_weight_data =
	x |>
	dplyr::select(DB, sample, cell_cell_weight) |>
	unnest(cell_cell_weight) |>

	#with_groups(c(sample, DB), ~ .x |> summarise(weight = mean(weight))) |>
	#filter(cell_from |> str_detect("cdc")) |>
	#filter(cell_to |> str_detect("cd8|cd4")) |>

	left_join(
		metadata_df ,
		by = c("sample" = ".sample")
	) |>
	filter(disease=="normal") |>

	mutate(lymphoid_organ = tissue_harmonised %in% c("blood", "lymph node", "bone", "spleen", "thymus")) |>
	filter(file_id |> str_detect("973fa", negate = TRUE)) |>

	#filter(DB=="Cell-Cell Contact") |>
	mutate(weight = weight + 0.0001) |>
	filter(lymphoid_organ) |>
	filter(!is.na(age_days))



library(brms)

# library(furrr)
# plan(multisession, workers = 20)

job::job({


	library(future)
	library("future.batchtools")
	library(furrr)

	slurm <- future::tweak(batchtools_slurm,
												 template = glue("~/third_party_sofware/slurm_batchtools.tmpl"),
												 resources=list(
												 	ncpus = 4,
												 	memory = 5000,
												 	walltime = 6000
												 )
	)
	plan(slurm)

	cdc_weight_data |>


		tidyr::separate( cell_to, c("cell_to", "phenotype", sep=" "), remove = FALSE) |>
		tidyr::separate( cell_from, c("cell_from", "phenotype", sep=" "), remove = FALSE) |>

		with_groups(c(DB, cell_to, cell_from, tissue_harmonised, sex , ethnicity , assay, sample, age_days, file_id), ~ .x |> summarise(weight = mean(weight))) |>

		# Filter
		add_count(DB, cell_to, cell_from, tissue_harmonised, assay, name = "count_partitions") |>
		filter(count_partitions>10) |>

		nest(data = -c(DB, cell_to, cell_from, tissue_harmonised, file_id)) |>
		filter(map_int(data, ~ .x |> distinct(age_days) |> nrow())>1) |>

		# Interval > 10 years
		filter(map_dbl(data, ~ max(.x$age_days, na.rm=T)-min(.x$age_days, na.rm=T)) > 3600*2) |>
		unnest(data) |>

		# Scale age
		mutate(age_days = scale(age_days) |> as.numeric()) |>

		nest(data = -c(DB, cell_to, cell_from, tissue_harmonised)) |>
		filter(map_int(data, ~ .x |> distinct(sample) |> nrow()) > 30) |>
		filter(map_int(data, ~ .x |> distinct(weight) |> nrow()) > 1) |>



		# Filter out if no covariate have more than one value
		filter(!map_lgl(data, ~
											(
												c("sex" , "ethnicity" , "assay" ) %in%
													(unlist(lapply(.x, function(x) length(unique(x)))) |> equals(1) %>% .[which(.)] |> names())
											) |> all()
		)) |>

		#slice(1:2) |>
		mutate(fit = future_imap(
			data,
			~ {
				print(.y)
				one_value = .x  |>
					lapply(function(x) table(x)) |>
					lapply(function(x) x[x>1]) |>
					lapply(function(x) length(unique(x))) |>
					unlist() |>
					equals(1) %>%
					.[which(.)] |>
					names()

				other_covariates =
					#c("sex" , "ethnicity" , "assay" ) |>
					c("ethnicity" , "sex", "assay" ) |>
					setdiff(one_value) |>
					str_c(collapse = " + ")

				f = "weight ~ 1 + age_days"
				if(other_covariates!="") f =	glue("{f} + {other_covariates}")
				if(!"file_id" %in% one_value) f =	glue("{f} + (1 + age_days | file_id)")


				# stan_glm(as.formula(f), family = Gamma(link="log"), data= .x, cores=2, chains = 2) |>
				# 	spread_draws(age_days) %>%
				# 	mean_qi()

				brm(
					as.formula(f),
					family = Gamma(link = "log"),
					data = .x,
					backend = "cmdstanr",
					cores = 4,
					prior = prior(normal(0,5), class = "b"),
					warmup = 500, iter = 800
				) |>
					summary(prob=0.95) %$%
					fixed %>%
					.["age_days",, drop=F]

			}
		)) |>

		# Parse
		mutate(estimate = map(fit, ~ .x )) |>
		dplyr::select(-fit) |>
		unnest(estimate) |>
		filter(!cell_to %in% c("immune", "non")) |>
		filter(!cell_from %in% c("immune", "non")) |>


		saveRDS("cdc_weight_fit_random_effect_file_broad_cells.rds")

})

cdc_weight_fit |> saveRDS("cdc_weight_fit_random_effect_assay.rds")

readRDS("cdc_weight_fit_random_effect_assay.rds")  |>
	filter((`l-95% CI` * `u-95% CI`)>0) |>
	print(n=99)


plot_trends =
	readRDS("cdc_weight_fit_random_effect_file_resolved_cells.rds") |>
	#filter((`l-95% CI` * `u-95% CI`)>0) |>
	filter(cell_from == "cdc") |>
	#filter(cell_to |> str_detect("cd8|cd4|b |ilc")) |>
	filter(!(DB=="ECM-Receptor" & tissue_harmonised=="blood")) |>
	filter(tissue_harmonised=="blood") |>
	filter(cell_to  |> str_detect("cd8|cd4|ilc|mait|tgd")) |>
	filter(DB=="Cell-Cell Contact") |>
	unnest(data) |>

	ggplot(aes(age_days, weight)) +
	geom_point(aes(color = assay), size=0.3) +
	stat_smooth(aes(color = assay, group=file_id), method = "glm", method.args = list(family = "Gamma"), se = F, level = 0.8) +
	facet_wrap(cell_to~DB, scales = "free_y") +
	scale_y_sqrt(lim = c(0, NA)) +
	theme_multipanel +
	theme(axis.text.x = element_text(angle = 30))



# CIRCLE
library(CellChat)


draw_cellchat_circle_plot = function (net, color.use = NULL, title.name = NULL, sources.use = NULL,
																			targets.use = NULL, remove.isolate = FALSE, top = 1, top_absolute = NULL, weight.scale = T,
																			vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = 15,
																			vertex.label.cex = 0.8, vertex.label.color = "black", edge.weight.max = NULL,
																			edge.width.max = 8, alpha.edge = 0.6, label.edge = FALSE,
																			edge.label.color = "black", edge.label.cex = 0.8, edge.curved = 0.2,
																			shape = "circle", layout = in_circle(), margin = 0.2, vertex.size = NULL,
																			arrow.width = 1, arrow.size = 0.2)
{
	if (!is.null(vertex.size)) {
		warning("'vertex.size' is deprecated. Use `vertex.weight`")
	}
	options(warn = -1)

	if(!is.null(top_absolute)) {
		thresh = top_absolute
		net[abs(net) < thresh] <- 0
	}

	thresh <- stats::quantile(as.numeric(net) %>% abs %>% .[.>0], probs = 1 - top)

	net[abs(net) < thresh] <- 0

	if(sum(net)==0) return(NULL)

	if ((!is.null(sources.use)) | (!is.null(targets.use))) {
		if (is.null(rownames(net))) {
			stop("The input weighted matrix should have rownames!")
		}
		cells.level <- rownames(net)
		df.net <- reshape2::melt(net, value.name = "value")
		colnames(df.net)[1:2] <- c("source", "target")
		if (!is.null(sources.use)) {
			if (is.numeric(sources.use)) {
				sources.use <- cells.level[sources.use]
			}
			df.net <- subset(df.net, source %in% sources.use)
		}
		if (!is.null(targets.use)) {
			if (is.numeric(targets.use)) {
				targets.use <- cells.level[targets.use]
			}
			df.net <- subset(df.net, target %in% targets.use)
		}
		df.net$source <- factor(df.net$source, levels = cells.level)
		df.net$target <- factor(df.net$target, levels = cells.level)
		df.net$value[is.na(df.net$value)] <- 0
		net <- tapply(df.net[["value"]], list(df.net[["source"]],
																					df.net[["target"]]), sum)
	}
	net[is.na(net)] <- 0
	if (remove.isolate) {
		idx1 <- which(Matrix::rowSums(net) == 0)
		idx2 <- which(Matrix::colSums(net) == 0)
		idx <- intersect(idx1, idx2)
		if(length(idx)>0){
			net <- net[-idx, ,drop=FALSE]
			net <- net[, -idx, drop=FALSE]
		}
	}
	g <- graph_from_adjacency_matrix(net, mode = "directed",
																	 weighted = T)
	edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
	coords <- layout_(g, layout)
	if (nrow(coords) != 1) {
		coords_scale = scale(coords)
	}
	else {
		coords_scale <- coords
	}
	if (is.null(color.use)) {
		color.use = scPalette(length(igraph::V(g)))
	}
	if (is.null(vertex.weight.max)) {
		vertex.weight.max <- max(vertex.weight)
	}
	vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max +
		5
	loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, -atan(coords_scale[igraph::V(g),
																																						 2]/coords_scale[igraph::V(g), 1]), pi - atan(coords_scale[igraph::V(g),
																																						 																													2]/coords_scale[igraph::V(g), 1]))
	igraph::V(g)$size <- vertex.weight

	if (is.null(color.use)) {
		igraph::V(g)$frame.color <- color.use[igraph::V(g)]
		igraph::V(g)$color <- color.use[igraph::V(g)]
	}
	else{
		igraph::V(g)$frame.color <- color.use[igraph::V(g) |> names()]
		igraph::V(g)$color <- color.use[igraph::V(g) |> names()]
	}
	igraph::V(g)$label.color <- vertex.label.color
	igraph::V(g)$label.cex <- vertex.label.cex
	if (label.edge) {
		igraph::E(g)$label <- igraph::E(g)$weight
		igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
	}
	if (is.null(edge.weight.max)) {
		edge.weight.max <- max(abs(igraph::E(g)$weight))
	}
	if (weight.scale == TRUE) {
		igraph::E(g)$width <- 0.3 + abs(igraph::E(g)$weight)/edge.weight.max *
			edge.width.max
	}
	else {
		igraph::E(g)$width <- 0.3 + edge.width.max * abs(igraph::E(g)$weight)
	}
	igraph::E(g)$arrow.width <- arrow.width
	igraph::E(g)$arrow.size <- arrow.size
	igraph::E(g)$label.color <- edge.label.color
	igraph::E(g)$label.cex <- edge.label.cex

	igraph::E(g)$color =
		circlize::colorRamp2(seq(max(abs(igraph::E(g)$weight)), -max(abs(igraph::E(g)$weight)), length.out =11), RColorBrewer::brewer.pal(11, "RdBu"))(igraph::E(g)$weight) %>%
		grDevices::adjustcolor(alpha.edge)


	if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
		igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[,
																																1])] <- loop.angle[edge.start[which(edge.start[,
																																																							 2] == edge.start[, 1]), 1]]
	}
	radian.rescale <- function(x, start = 0, direction = 1) {
		c.rotate <- function(x) (x + start)%%(2 * pi) * direction
		c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
	}
	label.locs <- radian.rescale(x = 1:length(igraph::V(g)),
															 direction = -1, start = 0)
	label.dist <- vertex.weight/max(vertex.weight) + 2
	plot(g, edge.curved = edge.curved, vertex.shape = shape,
			 layout = coords_scale,
			 margin = margin,
			 #vertex.label.dist = label.dist,
			 vertex.label.degree = label.locs,
			 vertex.label.family = "Helvetica",
			 edge.label.family = "Helvetica")
	if (!is.null(title.name)) {
		text(0, 1.5, title.name, cex = 0.8)
	}
	gg <- recordPlot()
	return(gg)
}

source("https://gist.githubusercontent.com/stemangiola/cfa08c45c28fdf223d4996a6c1256a39/raw/7dc960e253bf6c26284dd713935e11ecf5c0cfc9/color_cell_types.R")
cell_type_color =
	color_array %>%
	enframe(name="cell_from", value="color") |>
	filter(cell_from != "") |>
	#tidyr::separate( cell_from, c("cell_from", "phenotype", sep=" "), remove = FALSE) |>
	#with_groups(cell_from, slice, 1) |>
	dplyr::select(cell_from, color) |>
	deframe()


cdc_weight_fit = readRDS("cdc_weight_fit_random_intercept_file_resolved_cells.rds")

		cdc_weight_fit |>
			filter(DB == "Cell-Cell Contact" & tissue_harmonised=="blood") |>
			mutate(Estimate = if_else((`l-95% CI` * `u-95% CI`)>0, Estimate, 0)) |>
			filter(!cell_to %in% c("immune", "non")) |>
			filter(!cell_from %in% c("immune", "non")) |>

			mutate(Estimate = Estimate/Est.Error/2) |>
			#mutate(Estimate = if_else(cell_from=="cdc" | cell_to =="cdc", Estimate, 0)) |>
			dplyr::select(cell_from, cell_to, Estimate) |>
			tidyr::complete(cell_from, cell_to, fill= list(Estimate = 0)) |>
			pivot_wider(names_from = cell_to, values_from =  Estimate) |>
			tidybulk:::as_matrix(rownames = cell_from) |>
			draw_cellchat_circle_plot(
				edge.width.max = 8,
				#remove.isolate = TRUE,
				#top = 0.2,
				arrow.width = 8,
				arrow.size = 0.3,
				edge.weight.max = 2,
				color.use = cell_type_color
			)

		cdc_weight_fit |>
			filter(DB == "Secreted Signaling" & tissue_harmonised=="blood") |>
			mutate(Estimate = if_else((`l-95% CI` * `u-95% CI`)>0, Estimate, 0)) |>
			filter(!cell_to %in% c("immune", "non")) |>
			filter(!cell_from %in% c("immune", "non")) |>

			mutate(Estimate = Estimate/Est.Error/2) |>
			#mutate(Estimate = if_else(cell_from=="cdc" | cell_to =="cdc", Estimate, 0)) |>
			dplyr::select(cell_from, cell_to, Estimate) |>
			tidyr::complete(cell_from, cell_to, fill= list(Estimate = 0)) |>
			pivot_wider(names_from = cell_to, values_from =  Estimate) |>
			tidybulk:::as_matrix(rownames = cell_from) |>
			draw_cellchat_circle_plot(
				edge.width.max = 4,
				#remove.isolate = TRUE,
				#top = 0.2,
				arrow.width = 8,
				arrow.size = 0.3,
				edge.weight.max = 2
			)



# # BIG model
# 		cdc_weight_data |>
#
#
# 			tidyr::separate( cell_to, c("cell_to", "phenotype", sep=" "), remove = FALSE) |>
# 			tidyr::separate( cell_from, c("cell_from", "phenotype", sep=" "), remove = FALSE) |>
#
# 			with_groups(c(DB, cell_to, cell_from, tissue_harmonised, sex , ethnicity , assay, sample, age_days, file_id), ~ .x |> summarise(weight = mean(weight))) |>
#
# 			# Filter
# 			add_count(DB, cell_to, cell_from, tissue_harmonised, assay, name = "count_partitions") |>
# 			filter(count_partitions>10) |>
#
# 			nest(data = -c(DB, cell_to, cell_from, tissue_harmonised, file_id)) |>
# 			filter(map_int(data, ~ .x |> distinct(age_days) |> nrow())>1) |>
#
# 			# Interval > 10 years
# 			filter(map_dbl(data, ~ max(.x$age_days, na.rm=T)-min(.x$age_days, na.rm=T)) > 3600*2) |>
# 			unnest(data) |>
#
# 			# Scale age
# 			mutate(age_days = scale(age_days) |> as.numeric()) |>
#
# 			nest(data = -c(DB, cell_from, tissue_harmonised)) |>
# 			filter(map_int(data, ~ .x |> distinct(sample) |> nrow()) > 30) |>
# 			filter(map_int(data, ~ .x |> distinct(weight) |> nrow()) > 1) |>
#
#
#
# 			# Filter out if no covariate have more than one value
# 			filter(!map_lgl(data, ~
# 												(
# 													c("sex" , "ethnicity" , "assay" ) %in%
# 														(unlist(lapply(.x, function(x) length(unique(x)))) |> equals(1) %>% .[which(.)] |> names())
# 												) |> all()
# 			)) |>
#
# 			#slice(1:2) |>
# 			mutate(fit = future_imap(
# 				data,
# 				~ {
# 					#browser()
# 					print(.y)
# 					one_value = .x  |>
# 						lapply(function(x) table(x)) |>
# 						lapply(function(x) x[x>1]) |>
# 						lapply(function(x) length(unique(x))) |>
# 						unlist() |>
# 						equals(1) %>%
# 						.[which(.)] |>
# 						names()
#
# 					other_covariates =
# 						#c("sex" , "ethnicity" , "assay" ) |>
# 						c("ethnicity" , "sex", "assay" ) |>
# 						setdiff(one_value) |>
# 						str_c(collapse = " + ")
#
# 					f = "weight ~ 1 + age_days"
# 					if(other_covariates!="") f =	glue("{f} + {other_covariates}")
# 					if(!"file_id" %in% one_value) f =	glue("{f} + (1 + age_days | file_id)")
# 					f = glue("{f} + (1 + age_days | cell_to)")
#
#
#
# 					# stan_glm(as.formula(f), family = Gamma(link="log"), data= .x, cores=2, chains = 2) |>
# 					# 	spread_draws(age_days) %>%
# 					# 	mean_qi()
#
# 					brm(
# 						as.formula(f),
# 						family = Gamma(link = "log"),
# 						data = .x,
# 						backend = "cmdstanr",
# 						cores = 4,
# 						prior = prior(normal(0,5), class = "b"),
# 						warmup = 500, iter = 800
# 					) |>
# 						summary(prob=0.95) %$%
# 						fixed %>%
# 						.["age_days",, drop=F]
#
# 				}
# 			)) |>
#
# 			# Parse
# 			mutate(estimate = map(fit, ~ .x )) |>
# 			dplyr::select(-fit) |>
# 			unnest(estimate) |>
# 			# filter(!cell_to %in% c("immune", "non")) |>
# 			# filter(!cell_from %in% c("immune", "non")) |>
#
#
# 			saveRDS("cdc_weight_fit_random_effect_file_broad_cells_cell_from.rds")
#
#
# 		readRDS("cdc_weight_fit_random_effect_file_broad_cells_cell_from.rds")  |>
# 			filter((`l-95% CI` * `u-95% CI`)>0) |>
# 			print(n=99)








# PER GENE ANALYSIS

		compress_zero_one = function(y){
			# https://stats.stackexchange.com/questions/48028/beta-regression-of-proportion-data-including-1-and-0



			n = length(y)
			(y * (n-1) + 0.5) / n
		}

		cdc_weight_data_per_gene =
			x |>
			dplyr::select(DB, sample, cell_vs_all_cells_per_pathway) |>
			unnest(cell_vs_all_cells_per_pathway) |>

			left_join(
				metadata_df ,
				by = c("sample" = ".sample")
			) |>
			filter(disease=="normal") |>

			mutate(lymphoid_organ = tissue_harmonised %in% c("blood", "lymph node", "bone", "spleen", "thymus")) |>
			filter(file_id |> str_detect("973fa", negate = TRUE)) |>

			filter(lymphoid_organ) |>
			filter(!is.na(age_days)) |>
			filter(cell_type == "cdc") |>
			mutate(value = compress_zero_one(value))



		job::job({


			library(future)
			library("future.batchtools")
			library(furrr)

			slurm <- future::tweak(batchtools_slurm,
														 template = glue("~/third_party_sofware/slurm_batchtools.tmpl"),
														 resources=list(
														 	ncpus = 4,
														 	memory = 5000,
														 	walltime = 6000
														 )
			)
			plan(slurm)

			cdc_weight_data_per_gene |>


				# tidyr::separate( cell_to, c("cell_to", "phenotype", sep=" "), remove = FALSE) |>
				# tidyr::separate( cell_from, c("cell_from", "phenotype", sep=" "), remove = FALSE) |>
				# with_groups(c(DB, cell_to, cell_from, tissue_harmonised, sex , ethnicity , assay, sample, age_days, file_id), ~ .x |> summarise(weight = mean(weight))) |>
				#
				# Filter
				add_count(DB, cell_type, tissue_harmonised, assay, gene, name = "count_partitions") |>
				filter(count_partitions>10) |>

				nest(data = -c(DB, cell_type, tissue_harmonised, file_id, gene)) |>
				filter(map_int(data, ~ .x |> distinct(age_days) |> nrow())>1) |>

				# Interval > 10 years
				filter(map_dbl(data, ~ max(.x$age_days, na.rm=T)-min(.x$age_days, na.rm=T)) > 3600*2) |>
				unnest(data) |>

				# Scale age
				mutate(age_days_scaled = scale(age_days) |> as.numeric()) |>

				nest(data = -c(DB, cell_type, tissue_harmonised, gene)) |>
				filter(map_int(data, ~ .x |> distinct(sample) |> nrow()) > 30) |>
				filter(map_int(data, ~ .x |> distinct(value) |> nrow()) > 1) |>



				# Filter out if no covariate have more than one value
				filter(!map_lgl(data, ~
													(
														c("sex" , "ethnicity" , "assay" ) %in%
															(unlist(lapply(.x, function(x) length(unique(x)))) |> equals(1) %>% .[which(.)] |> names())
													) |> all()
				)) |>

				#slice(1:2) |>
				mutate(fit = future_imap(
					data,
					~ {
						print(.y)
						one_value = .x  |>
							lapply(function(x) table(x)) |>
							lapply(function(x) x[x>1]) |>
							lapply(function(x) length(unique(x))) |>
							unlist() |>
							equals(1) %>%
							.[which(.)] |>
							names()

						other_covariates =
							#c("sex" , "ethnicity" , "assay" ) |>
							c("ethnicity" , "sex", "assay" ) |>
							setdiff(one_value) |>
							str_c(collapse = " + ")

						f = "value ~ 1 + age_days_scaled"
						if(other_covariates!="") f =	glue("{f} + {other_covariates}")
						if(!"file_id" %in% one_value) f =	glue("{f} + (1 | file_id)")


						# stan_glm(as.formula(f), family = Gamma(link="log"), data= .x, cores=2, chains = 2) |>
						# 	spread_draws(age_days_scaled) %>%
						# 	mean_qi()

						brm(
							as.formula(f),
							family = Beta(),
							data = .x,
							backend = "cmdstanr",
							cores = 4,
							prior = prior(normal(0,2.5), class = "b")
						)

					}
				)) |>

				saveRDS("cdc_gene_fit_random_intercept_file_resolved_cells.rds")

		})


	ff = 	readRDS("cdc_gene_fit_random_intercept_file_resolved_cells.rds")  |>
			mutate(estimate = map(fit,
														~.x |>
															summary(prob=0.95) %$%
															fixed %>%
															.["age_days_scaled",, drop=F]
															)) |>
			unnest(estimate)

	# Volcano plot
	volcano_plot =
		ff |>
		filter(tissue_harmonised=="blood") |>
		mutate(posterior = map(fit, ~ .x |> spread_draws(b_age_days_scaled) |> pull(b_age_days_scaled))) |>
		dplyr::select(DB, gene, cell_type, tissue_harmonised, posterior, Estimate) |>
		unnest(posterior) |>
		#tidyr::unite("condition", DB, gene, cell_type, tissue_harmonised, sep = "\n", remove = FALSE)  |>
		count(DB, gene, cell_type, tissue_harmonised, Estimate, posterior>0) |>
		complete(nesting(DB, gene, cell_type, tissue_harmonised), `posterior > 0`, fill = list(n=0)) |>
		with_groups(c(DB, gene, cell_type, tissue_harmonised), ~ .x |> mutate(n_tot = sum(n))) |> mutate(prob = n/n_tot) |>
		with_groups(c(DB, gene, cell_type, tissue_harmonised), ~ .x |> arrange(prob) |> slice(1)) |>
		filter(DB!="ECM-Receptor") |>
		mutate(significant = prob<0.05) |>
		mutate(DB = case_when(significant~DB)) |>
		mutate(gene = case_when(significant~gene)) |>
		ggplot(aes(Estimate, prob, label=gene)) +
		geom_point(aes(color=DB, size=significant)) +
		ggrepel::geom_text_repel(size = 1.5 ) +
		geom_hline(yintercept = 0.05, linetype="dashed", color="darkgrey", size = 0.2) +
		geom_vline(xintercept = 0, linetype="dashed", color="darkgrey", size = 0.2) +
		scale_size_discrete(range = c(0, 0.5)) +
		scale_color_manual(values = c("#F1A340", "#998FC3"), na.value = "black") +
		scale_y_continuous(trans = tidybulk::log10_reverse_trans()) +
		xlab("Association communication strenght/age") +
		ylab("Probability") +
		theme_multipanel


	ggsave(
		"volcano_communication.pdf",
		plot = volcano_plot,
		units = c("mm"),
		width = 40 ,
		height = 53 ,
		limitsize = FALSE
	)

	ff |>
		filter((`l-95% CI` * `u-95% CI`)>0) |>

			print(n=99)


		plot_trends =
			readRDS("cdc_gene_fit_random_intercept_file_resolved_cells.rds") |>
			filter((`l-95% CI` * `u-95% CI`)>0) |>
			#filter(cell_from == "cdc") |>
			#filter(cell_to |> str_detect("cd8|cd4|b |ilc")) |>
			#filter(!(DB=="ECM-Receptor" & tissue_harmonised=="blood")) |>
			filter(tissue_harmonised=="blood") |>
			#filter(cell_to  |> str_detect("cd8|cd4|ilc|mait|tgd")) |>
			#filter(DB=="Cell-Cell Contact") |>
			unnest(data) |>

			ggplot(aes(age_days, value)) +
			geom_point(aes(color = assay), size=0.3) +
			stat_smooth(aes(color = assay), method = "glm", method.args = list(family = "binomial"), se = F, level = 0.8) +
			facet_wrap(gene~DB, scales = "free_y") +
			scale_y_sqrt(lim = c(0, NA)) +
			theme_multipanel +
			theme(axis.text.x = element_text(angle = 30))


	# GENE-CELL ANALYSIS


		cdc_weight_data_per_gene_AND_cell =
			x |>
			dplyr::select(DB, sample, cell_vs_cell_per_pathway) |>
			unnest(cell_vs_cell_per_pathway) |>
			filter(map_int(result, length)>0) |>

			filter(gene %in% c("GALECTIN", "ANNEXIN",  "BAFF" ,    "MHC-I" ,
												 "SELPLG",   "ITGB2",    "ICAM"   ,  "ADGRE5"  , "CD23" ,
												 "ADGRE5"   )
			) |>


			left_join(
				metadata_df ,
				by = c("sample" = ".sample")
			) |>
			filter(disease=="normal") |>

			mutate(lymphoid_organ = tissue_harmonised %in% c("blood", "lymph node", "bone", "spleen", "thymus")) |>
			filter(file_id |> str_detect("973fa", negate = TRUE)) |>

			#filter(DB=="Cell-Cell Contact") |>
			#mutate(value = value + 0.0001) |>
			filter(lymphoid_organ) |>
			filter(!is.na(age_days)) |>

			mutate(result = map(
				result,
				~ .x |>
					as_tibble(rownames = "cell_from") |>
					pivot_longer(-cell_from, names_to = "cell_to", values_to = "value")

			)) |>
			unnest(result) |>

			#with_groups(c(sample, DB), ~ .x |> summarise(weight = mean(weight))) |>
			#filter(cell_from |> str_detect("cdc")) |>
			#filter(cell_to |> str_detect("cd8|cd4")) |>

			filter(cell_from == "cdc" | cell_to == "cdc") |>
			mutate(value = compress_zero_one(value))



		job::job({


			library(future)
			library("future.batchtools")
			library(furrr)

			slurm <- future::tweak(batchtools_slurm,
														 template = glue("~/third_party_sofware/slurm_batchtools.tmpl"),
														 resources=list(
														 	ncpus = 4,
														 	memory = 5000,
														 	walltime = 6000
														 )
			)
			plan(slurm)

			cdc_weight_data_per_gene_AND_cell |>


				# tidyr::separate( cell_to, c("cell_to", "phenotype", sep=" "), remove = FALSE) |>
				# tidyr::separate( cell_from, c("cell_from", "phenotype", sep=" "), remove = FALSE) |>
				# with_groups(c(DB, cell_to, cell_from, tissue_harmonised, sex , ethnicity , assay, sample, age_days, file_id), ~ .x |> summarise(weight = mean(weight))) |>
				#
				# Filter
				mutate(log_n = log(n)) |>

				add_count(DB, cell_from, cell_to, tissue_harmonised, assay, gene, name = "count_partitions") |>
				filter(count_partitions>10) |>

				nest(data = -c(DB, cell_from, cell_to, tissue_harmonised, file_id, gene)) |>
				filter(map_int(data, ~ .x |> distinct(age_days) |> nrow())>1) |>

				# Interval > 10 years
				filter(map_dbl(data, ~ max(.x$age_days, na.rm=T)-min(.x$age_days, na.rm=T)) > 3600*2) |>
				unnest(data) |>

				# Scale age
				mutate(age_days = scale(age_days) |> as.numeric()) |>
				mutate(value = scale(value, center=F) |> as.numeric()) |>

				nest(data = -c(DB, cell_from,  cell_to, tissue_harmonised, gene)) |>
				filter(map_int(data, ~ .x |> distinct(sample) |> nrow()) > 30) |>
				filter(map_int(data, ~ .x |> distinct(value) |> nrow()) > 1) |>



				# Filter out if no covariate have more than one value
				filter(!map_lgl(data, ~
													(
														c("sex" , "ethnicity" , "assay" ) %in%
															(unlist(lapply(.x, function(x) length(unique(x)))) |> equals(1) %>% .[which(.)] |> names())
													) |> all()
				)) |>

				#slice(1:2) |>
				mutate(fit = future_imap(
					data,
					~ {
						print(.y)
						one_value = .x  |>
							lapply(function(x) table(x)) |>
							lapply(function(x) x[x>1]) |>
							lapply(function(x) length(unique(x))) |>
							unlist() |>
							equals(1) %>%
							.[which(.)] |>
							names()

						other_covariates =
							#c("sex" , "ethnicity" , "assay" ) |>
							c("ethnicity" , "sex", "assay" ) |>
							setdiff(one_value) |>
							str_c(collapse = " + ")

						f = "value ~ 1 + age_days + log_n "
						if(other_covariates!="") f =	glue("{f} + {other_covariates}")
						if(!"file_id" %in% one_value) f =	glue("{f} + (1 | file_id)")


						# stan_glm(as.formula(f), family = Gamma(link="log"), data= .x, cores=2, chains = 2) |>
						# 	spread_draws(age_days) %>%
						# 	mean_qi()

						brm(
							as.formula(f),
							family = Gamma(link = "log"),
							data = .x,
							backend = "cmdstanr",
							cores = 4,
							prior = prior(normal(0,2.5), class = "b")
						) |>
							summary(prob=0.95) %$%
							fixed %>%
							.["age_days",, drop=F]

					}
				)) |>

				# Parse
				mutate(estimate = map(fit, ~ .x )) |>
				dplyr::select(-fit) |>
				unnest(estimate) |>
				# filter(!cell_to %in% c("immune", "non")) |>
				# filter(!cell_from %in% c("immune", "non")) |>


				saveRDS("cdc_gene_AND_cell_fit_random_intercept_file_resolved_cells.rds")

		})

		readRDS("cdc_gene_AND_cell_fit_random_intercept_file_resolved_cells.rds")  |>
			filter((`l-95% CI` * `u-95% CI`)>0) |>
			print(n=99)



		source("https://gist.githubusercontent.com/stemangiola/cfa08c45c28fdf223d4996a6c1256a39/raw/7dc960e253bf6c26284dd713935e11ecf5c0cfc9/color_cell_types.R")
		cell_type_color =
			color_array %>%
			enframe(name="cell_from", value="color") |>
			filter(cell_from != "") |>
			tidyr::separate( cell_from, c("cell_from", "phenotype", sep=" "), remove = FALSE) |>
			with_groups(cell_from, slice, 1) |>
			dplyr::select(cell_from, color) |>
			deframe()

pdf("circles_communication.pdf", width = 183*1.6*0.0393701, height = 70*1.6*0.0393701)
par(mfrow=c(2, 3))
plot_circles =
		readRDS("cdc_gene_AND_cell_fit_random_intercept_file_resolved_cells.rds")  |>

			mutate(Estimate = Estimate/Est.Error/2) |>

			filter(tissue_harmonised=="blood") |>
			#filter(cell_from=="cdc") |>
			filter(cell_from!="non_immune" & cell_to != "non_immune") |>
			filter(cell_from!="immune_unclassified" & cell_to != "immune_unclassified") |>
			#filter((`l-95% CI` * `u-95% CI`)>0) |>

			nest(gene_data = -c(gene, DB, tissue_harmonised)) |>

	left_join(
		CellChatDB.human$interaction|>
			filter(pathway_name |> str_detect("CD23|MHC-I$|ICAM|ADGRE5|SELPLG|ITGB2")) |>
			nest(data = -pathway_name) |>
			mutate(ligand_receptor = map_chr(
				data,
				~ glue("({.x$ligand |> unique() |> sort() |>  str_c(collapse = ',')}) -> \n ({.x$receptor |> unique() |> sort() |>  str_c(collapse = ',')})")
			)) |>
			mutate(ligand_receptor = ligand_receptor |>
						 	str_replace("HLA-A,HLA-B,HLA-C,HLA-E,HLA-F,HLA-G", "HLA-A-C/E-G") |>
						 	str_replace("CD8A,CD8B,CD8B2", "CD8A/B/B2") |>
						 	str_replace("ITGAL,ITGAL_ITGB2,ITGAM_ITGB2,ITGAX_ITGB2", "ITGAL/M/X_ITGB2") |>
						 	str_replace("RAET1E,RAET1F,RAET1G", "RAET1E-G") |>
						 	str_replace("CD94:NKG2A,CD94:NKG2C,CD94:NKG2E", "CD94:NKG2A/C/E") |>
						 	str_replace("KIR2DL1,KIR2DL2,KIR2DL3,KIR2DS1,KIR2DS4,KIR3DL1,KIR3DL2,KIR3DL3,KIR3DS1,KLRC1,KLRC2,KLRK1", "KIR2DL/DS/C") |>
						 	str_replace("LILRB1,LILRB2", "LILRB1/2") |>
						 	str_replace("ITGAM_ITGB2,ITGAV_ITGB3,ITGAX_ITGB2", "ITGAM/X_ITGB2,ITGAV_ITGB3")
						),
		by=c("gene" = "pathway_name")
	) |>
			mutate(title = glue("{ligand_receptor}")) |>
			filter(gene %in% c("CD23", "MHC-I", "ICAM", "ADGRE5", "SELPLG", "ITGB2")) |>

			mutate(plot = map2(
				gene_data, title,
				~ {

					cell_types = .x |> select(cell_from, cell_to) |> gather() |> pull(value) |>  unique()

					expand_grid(cell_from = cell_types, cell_to = cell_types) |>
						left_join(
							.x |>
								dplyr::select(cell_from, cell_to, Estimate)
						) |>
					tidyr::complete(cell_from, cell_to, fill= list(Estimate = 0)) |>
					pivot_wider(names_from = cell_to, values_from =  Estimate) |>
					tidybulk:::as_matrix(rownames = cell_from) |>
						multiply_by(2) |>
					draw_cellchat_circle_plot(
						edge.width.max = 4,
						#remove.isolate = TRUE,
						#top = 0.2,
						arrow.width = 3,
						arrow.size = 0.3,
						edge.weight.max = 8,
						color.use = cell_type_color,
						title.name = .y
					)
				}
			))
dev.off()


plot_circles |>
			pull(plot)


		cdc_weight_fit |>
			filter(DB == "Secreted Signaling" & tissue_harmonised=="blood") |>
			mutate(Estimate = if_else((`l-95% CI` * `u-95% CI`)>0, Estimate, 0)) |>
			filter(!cell_to %in% c("immune", "non")) |>
			filter(!cell_from %in% c("immune", "non")) |>

			mutate(Estimate = Estimate/Est.Error/2) |>
			#mutate(Estimate = if_else(cell_from=="cdc" | cell_to =="cdc", Estimate, 0)) |>
			dplyr::select(cell_from, cell_to, Estimate) |>
			tidyr::complete(cell_from, cell_to, fill= list(Estimate = 0)) |>
			pivot_wider(names_from = cell_to, values_from =  Estimate) |>
			tidybulk:::as_matrix(rownames = cell_from) |>
			draw_cellchat_circle_plot(
				edge.width.max = 4,
				#remove.isolate = TRUE,
				#top = 0.2,
				arrow.width = 8,
				arrow.size = 0.3,
				edge.weight.max = 2
			)


















# library(lme4)
# # Analyses cell-cell weight
# fits =  x |>
#  	dplyr::select(DB, sample, cell_cell_weight) |>
#  	unnest(cell_cell_weight) |>
#
#  	#with_groups(c(sample, DB), ~ .x |> summarise(weight = mean(weight))) |>
#  	# filter(cell_from |> str_detect("cdc")) |>
#  	# filter(cell_to |> str_detect("cd8|cd4")) |>
#
#  	left_join(
#  		metadata_df |>
#  			count(.sample, sex, age_days, assay, ethnicity, tissue_harmonised, file_id),
#  		by = c("sample" = ".sample")
#  	) |>
#  	mutate(lymphoid_organ = tissue_harmonised %in% c("blood", "lymph node", "bone", "spleen")) |>
#  	filter(file_id |> str_detect("973fa", negate = TRUE)) |>
#
#  	nest(data = -c(DB, cell_from, cell_to, tissue_harmonised)) |>
#  	filter(map_int(data, ~ .x |> distinct(sample) |> nrow()) > 10) |>
#
# 	# Filter out if no covariate have more than one value
# 	filter(!map_lgl(data, ~
# 										(
# 											c("sex" , "ethnicity" , "assay" ) %in%
# 										 	(unlist(lapply(.x, function(x) length(unique(x)))) |> equals(1) %>% .[which(.)] |> names())
# 										 ) |> all()
# 					)) |>
#  	mutate(fit = imap(
#  		data,
#  		~ {print(.y)
#  			one_value = .x  |>
#  				lapply(function(x) x |> unique() |> length()) |>
#  				unlist() |>
#  				equals(1) %>%
#  				.[which(.)] |>
#  				names()
#
#  			f =
#  				sprintf(
#  					"weight ~ age_days + %s",
#  				c("sex" , "ethnicity" , "assay" ) |> setdiff(one_value) |> str_c(collapse = " + ")
#  			)
#
#  			lm(as.formula(f), data= .x) |> broom::tidy()
#  		}
#  	))

# pvalue histogram
fits |> unnest(fit) |> pull(p.value) |> hist()

# Volcano plot
# Volcano plot
fits |>
	unnest(fit) |>
	filter(term!="(Intercept)") |>
	mutate(p.value.adjusted = p.adjust(p.value)) |>
	mutate(p.value = p.value |> pmax(1e-10)) |>
	mutate(term_parsed = term |> str_sub(1, 3)) |>
	mutate(term = if_else(p.value.adjusted<0.05, term, "")) |>
	ggplot(aes(estimate, p.value)) +
	geom_point(aes(color=tissue_harmonised, shape=p.value.adjusted<0.05)) +
	ggrepel::geom_text_repel(aes(label=term)) +
	facet_wrap(~ term_parsed, scales = "free") +
	ylim(1,1e-9) +
	scale_y_continuous(trans = tidybulk::log10_reverse_trans())



library(lme4)
# Analyses cell pathway
fits =
	x |>
	dplyr::select(DB, sample, cell_vs_all_cells_per_pathway) |>
	unnest(cell_vs_all_cells_per_pathway) |>

	#with_groups(c(sample, DB), ~ .x |> summarise(weight = mean(weight))) |>
	# filter(cell_from |> str_detect("cdc")) |>
	# filter(cell_to |> str_detect("cd8|cd4")) |>

	left_join(
		metadata_df |>
			count(.sample, sex, age_days, assay, ethnicity, tissue_harmonised, file_id),
		by = c("sample" = ".sample")
	) |>
	filter(!is.na(age_days)) |>
	mutate(age_days = scale(age_days) |> as.numeric()) |>
	mutate(lymphoid_organ = tissue_harmonised %in% c("blood", "lymph node", "bone", "spleen")) |>
	filter(file_id |> str_detect("973fa", negate = TRUE)) |>

	nest(data = -c(DB, cell_type, gene, tissue_harmonised)) |>
	filter(map_int(data, ~ .x |> distinct(sample) |> nrow()) > 10) |>

	# Filter out if no covariate have more than one value
	filter(!map_lgl(data, ~
										(
											c("sex" , "ethnicity" , "assay" ) %in%
												(unlist(lapply(.x, function(x) length(unique(x)))) |> equals(1) %>% .[which(.)] |> names())
										) |> all()
	)) |>
	mutate(fit = imap(
		data,
		~ {print(.y)

			one_value = .x  |>
				lapply(function(x) table(x)) |>
				lapply(function(x) x[x>1]) |>
				lapply(function(x) length(unique(x))) |>
					unlist() |>
				equals(1) %>%
				.[which(.)] |>
				names()

		other_covariates = c("sex" , "ethnicity" , "assay" ) |> setdiff(one_value) |> str_c(collapse = " + ")

		f = "value ~ age_days"
		if(other_covariates!="") f =	glue("{f} + {other_covariates}")

			lm(as.formula(f), data= .x) |> broom::tidy()
		}
	))



# pvalue histogram
fits |>  unnest(fit) |> 	filter(term=="age_days") |>
 pull(p.value) |> hist()

# Volcano plot
# Volcano plot
fits |>
	unnest(fit) |>
	filter(term=="age_days") |>
	mutate(p.value.adjusted = p.adjust(p.value)) |>
	mutate(p.value = p.value |> pmax(1e-10)) |>
	mutate(term_parsed = term |> str_sub(1, 3)) |>
	mutate(term = if_else(p.value.adjusted<0.05, term, "")) |>
	ggplot(aes(estimate, p.value)) +
	geom_point(aes(color=tissue_harmonised, shape=p.value.adjusted<0.05)) +
	ggrepel::geom_text_repel(aes(label=term)) +
	facet_wrap(~ term_parsed, scales = "free") +
	ylim(1,1e-9) +
	scale_y_continuous(trans = tidybulk::log10_reverse_trans())



















 	`filter(DB=="Cell-Cell Contact") |>
 	#filter(cell_to=="cd8 tem") %>%
 with_groups(tissue_harmonised, ~ .x |> mutate(weight = scale(weight, center =FALSE) |> as.numeric())) %>%
 	mutate(age_days = scale(age_days) |> as.numeric()) %>%
 		filter(map_int(data, ~ .x |> distinct(sample) |> nrow) > 10) |>
 		lmer(
 	formula = weight ~ 1 + age_days + tissue_harmonised + sex + ethnicity  + assay + (1 | file_id) + (age_days | tissue_harmonised),
 		 data    = .
 		) |>
 	anova()

 netVisual_heatmap(cellchat)


job::job({

	fits_tot_interactions =
		x |>
		left_join(
			metadata |> count(.sample, sex, age_days, assay, ethnicity, tissue_harmonised, file_id)
		) |>
		mutate(age_days = scale(age_days) |> as.numeric()) |>
		nest(data = -c(DB)) |>
		#filter(cell_to == "cd8 tem") |>
		mutate(fit = map(
			data,
			~  brm(tot_interactions ~ 1 + age_days + sex + ethnicity + assay + (age_days | tissue_harmonised) + (1 | file_id),
						 data = .x,
						 warmup = 500, iter = 1000,
						 cores = 2, chains = 2,
						 seed = 123,
						 family = negbinomial()
			) #to run the model
		))

})



xx =
	x |>
	select(DB, sample, cell_vs_cell_per_pathway) |>
	unnest(cell_vs_cell_per_pathway) |>
	filter(!map_lgl(result, is.null)) |>
	mutate(result =  map(
		result,
		~ .x |>
			as_tibble(rownames = "cell_from") |>
			pivot_longer(-cell_from, names_to = "cell_to", values_to = "score")
	)) |>
	unnest(result)

xx = xx |>
	bind_cols(xx |> pull(sample))

xx |> saveRDS("~/temp.rds")

# Plot from to directionality score
xx |>
	filter(cell_from == "cdc" | cell_to == "cdc") |>
	filter(cell_from  |> str_detect("cd4|cd8") | cell_to  |> str_detect("cd4|cd8") ) |>
	unite("cell_types", c(cell_from, cell_to), sep="->", remove = FALSE) |>
	ggplot(aes(cell_types, score)) +
	geom_boxplot() +
	facet_wrap(~cell_from, scales = "free_x") +
	theme(axis.text.x = element_text(angle = 30))

metadata = get_metadata() |> as_tibble()

xx |>
	filter(cell_from == "cdc" ) |>
	filter( cell_to  |> str_detect("cd4|cd8") ) |>
	unite("cell_types", c(cell_from, cell_to), sep="->", remove = FALSE) |>
	with_groups(c(.sample, DB, cell_from, cell_to), ~ .x |> summarise(mean_score = mean(score))) |>
	#nest(data = -.sample) |>
	left_join(
		metadata |> distinct(.sample, sex, age_days, assay, ethnicity, tissue_harmonised)
	) |>
	ggplot(aes(age_days, mean_score, color = DB)) +
	geom_point() +
	facet_grid(tissue_harmonised~cell_to, scales = "free_x") +
	theme(axis.text.x = element_text(angle = 30)) +
	geom_smooth(method = "lm") +
	ylim(-0.6, 0.6)



fits =
	xx |>
	filter(cell_from == "cdc" ) |>
	filter( cell_to  |> str_detect("cd4|cd8") ) |>
	unite("cell_types", c(cell_from, cell_to), sep="->", remove = FALSE) |>
	with_groups(c(.sample, DB, cell_from, cell_to), ~ .x |> summarise(mean_score = mean(score) |> sqrt())) |>
	#nest(data = -.sample) |>
	left_join(
		metadata |>
			count(.sample, sex, age_days, assay, ethnicity, tissue_harmonised, file_id, cell_type_harmonised) |>
			pivot_wider(names_from = cell_type_harmonised, values_from = n)
	) |>
	mutate(age_days = scale(age_days) |> as.numeric()) |>
	nest(data = -c(cell_to, DB)) |>
	#filter(cell_to == "cd8 tem") |>
	mutate(fit = map(
		data,
		~  brm(mean_score ~ 1 + age_days + sex + ethnicity + assay + (1 + age_days | tissue_harmonised),
					 data = .x,
					 warmup = 500, iter = 1000,
					 cores = 2, chains = 2,
					 seed = 123) #to run the model
	))

fits |> pull(fit)
