library(zellkonverter)
library(Seurat)
library(SingleCellExperiment) # load early to avoid masking dplyr::count()
library(dplyr)
library(cellxgenedp)
library(tidyverse)
library(stringr)
library(scMerge)
library(glue)
library(DelayedArray)
library(HDF5Array)
library(tidyseurat)
library(celldex)
library(SingleR)
library(dplyr); library(tidyr); library(ggplot2)
library(Seurat)
library(glue)
library(CellChat)
library(tidyseurat)
library(tidySingleCellExperiment)
library(stringr)
library(HCAquery)
library(tidySingleCellExperiment)

# BUILD PIPELINE

# # # CREATE MAKEFILE
# tab = "\t"
# root_directory = "/vast/projects/RCP/human_cell_atlas"
# ligand_receptor_data_directory = glue("{root_directory}/ligand_receptor_data_0.2")
# light_data_directory = glue("{root_directory}/splitted_light_data")
# metadata_annotated = "/vast/projects/RCP/human_cell_atlas/metadata_annotated.rds" # glue("{root_directory}/metadata.rds")
# cell_type_df = "metadata_cell_type.csv"
#
# light_file_paths = dir(light_data_directory, full.names = TRUE)
# .sample = basename(light_file_paths) |> tools::file_path_sans_ext()
# ligand_receptor_file_paths = glue("{ligand_receptor_data_directory}/{.sample}.rds")
# file_for_annotation_workflow = glue("{root_directory}/cell_sample_cell_type_df_for_annotation_workflow.rds")
#
#
#
# readRDS(metadata_annotated) |>
#   distinct(.sample, file_id) |>
#   # filter(.sample %in% samples_to_use) |>
#   mutate(
#     input_file_path = glue("{light_data_directory}/{.sample}") |> as.character(),
#     output_file_path = glue("{ligand_receptor_data_directory}/{.sample}.rds" |> as.character())
#   ) |>
#
#   mutate(Mb = map_dbl(input_file_path, ~
#                         ( (file.info(glue("{.x}/se.rds"))$size /1e6) |> as.integer() ) +
#                         ( (file.info(glue("{.x}/assays.h5"))$size /1e6) |> as.integer() )
#   )) |>
#   mutate(memory = pmax(Mb * 20 + 10000, 20000)) |>
#
#   # mutate(memory = case_when(
#   # 	output_file_path %in% c(
#   # 		"/vast/projects/RCP/human_cell_atlas/ligand_receptor_data/305a8bd9b8e529a967feb5f73cc8c4df.rds" ,
#   # 		"/vast/projects/RCP/human_cell_atlas/ligand_receptor_data/087c2093be040a404c9685af1ecb3c65.rds",
#   # 		"/vast/projects/RCP/human_cell_atlas/ligand_receptor_data/86a6d20305d912e98318ad4d1d5d1814.rds",
#   # 		"/vast/projects/RCP/human_cell_atlas/ligand_receptor_data/829b99a569ec9ebb5fdd1b0b29208aaf.rds",
#   # 		"/vast/projects/RCP/human_cell_atlas/ligand_receptor_data/176f8892a21bec1bd7bdbc4181af75ed.rds",
#   # 		"/vast/projects/RCP/human_cell_atlas/ligand_receptor_data/23c822334c194bceb576a9ccb1db5929.rds",
#   # 		"/vast/projects/RCP/human_cell_atlas/ligand_receptor_data/3c20ba18525fb5e0b41cb8ea189b5d33.rds",
#   # 		"/vast/projects/RCP/human_cell_atlas/ligand_receptor_data/522dde7ab389d65b265d4cd598576f31.rds",
# # 		"/vast/projects/RCP/human_cell_atlas/ligand_receptor_data/2cf3bb4ffbb2024a9ca04baec073ae14.rds",
# # 		"/vast/projects/RCP/human_cell_atlas/ligand_receptor_data/c8ff7c63b3152a25c338cc279b31ab07.rds",
# # 		"/vast/projects/RCP/human_cell_atlas/ligand_receptor_data/5be263dbc1384b3cec21c5d3c580f838.rds",
# # 		"/vast/projects/RCP/human_cell_atlas/ligand_receptor_data/024d53b702b1846a476cabe5d691f992.rds",
# # 		"/vast/projects/RCP/human_cell_atlas/ligand_receptor_data/9da244f06591fa49e5649c65ed3b0934.rds",
# # 		)~ 160000,
# # 	TRUE ~ memory
# # )) |>
#
# rowid_to_column() |>
#   mutate(commands = pmap(list(output_file_path, input_file_path,  memory, rowid, file_id), ~
#                            c(
#                              glue("CATEGORY=light_data{..4}\nMEMORY={..3}\nCORES=1\nWALL_TIME=30000"),
#                              glue("{..1}:{..2} {metadata_annotated} {cell_type_df}\n{tab}Rscript ligand_receptor_count.R {..2} {metadata_annotated} {..1}")
#                            )
#   ))  |>
#   pull(commands) |>
#   unlist() |>
#   write_lines(glue("~/PostDoc/HCAquery/dev/ligand_receptor_count.makeflow"))
#


# ANALYSIS

# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_path_preprocessing = args[[1]]
file_for_annotation_workflow = args[[2]]
output_path = args[[3]]


#library(furrr)
#plan(multisession, workers=3)

# Create dir
output_path |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

data =
  loadHDF5SummarizedExperiment(input_path_preprocessing	)  |>

  # add lineage 1
  left_join(
    readRDS(file_for_annotation_workflow) |>
      unite(".cell", c(.cell, .sample), sep="_", remove=FALSE	)
  ) |>

	# Filter
	filter(!is.na(cell_type_harmonised)) |>
	tidySingleCellExperiment::add_count(cell_type_harmonised) |>
	tidySingleCellExperiment::filter(n>=30)

names(data@assays@data) = "counts"

# Filter primary data


# If data with less than 2 cell types left save empty output and exit
if((data |> distinct(cell_type_harmonised) |> nrow()) < 2){

	tibble()  |>

		# Save
		saveRDS(output_path)

} else {

# Otherwise continue



	my_sample = data |> distinct(.sample)

	counts_cellchat =
		data |>
		scater::logNormCounts() |>
		createCellChat(group.by = "cell_type_harmonised") |>
		setIdent( ident.use = "cell_type_harmonised")

	# More robust implementation that does not fail if no results
	computeCommunProbPathway = function (object = NULL, net = NULL, pairLR.use = NULL, thresh = 0.05)
	{
		if (is.null(net)) {
			net <- object@net
		}
		if (is.null(pairLR.use)) {
			pairLR.use <- object@LR$LRsig
		}
		prob <- net$prob
		prob[net$pval > thresh] <- 0
		pathways <- unique(pairLR.use$pathway_name)
		group <- factor(pairLR.use$pathway_name, levels = pathways)
		prob.pathways <- aperm(apply(prob, c(1, 2), by, group, sum),
													 c(2, 3, 1))
		pathways.sig <- pathways[apply(prob.pathways, 3, sum) !=     0]
		prob.pathways.sig <- prob.pathways[, , pathways.sig, drop=FALSE]
		idx <- sort(apply(prob.pathways.sig, 3, sum), decreasing = TRUE,
								index.return = TRUE)$ix
		pathways.sig <- pathways.sig[idx]
		prob.pathways.sig <- prob.pathways.sig[, , idx, drop=FALSE]
		if (is.null(object)) {
			netP = list(pathways = pathways.sig, prob = prob.pathways.sig)
			return(netP)
		}
		else {
			object@netP$pathways <- pathways.sig
			object@netP$prob <- prob.pathways.sig
			return(object)
		}
	}

	cellchat_matrix_for_circle = function (object, signaling, signaling.name = NULL, color.use = NULL,
																				 vertex.receiver = NULL, sources.use = NULL, targets.use = NULL,
																				 top = 1, remove.isolate = FALSE, vertex.weight = NULL, vertex.weight.max = NULL,
																				 vertex.size.max = 15, weight.scale = TRUE, edge.weight.max = NULL,
																				 edge.width.max = 8, layout = c("hierarchy", "circle", "chord"),
																				 thresh = 0.05, from = NULL, to = NULL, bidirection = NULL,
																				 vertex.size = NULL, pt.title = 12, title.space = 6, vertex.label.cex = 0.8,
																				 group = NULL, cell.order = NULL, small.gap = 1, big.gap = 10,
																				 scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20,
																				 legend.pos.y = 20, ...) {

		if(object@LR$LRsig %>% filter(pathway_name == signaling) %>% nrow %>% magrittr::equals(0)) return(NULL)

		layout <- match.arg(layout)
		if (!is.null(vertex.size)) {
			warning("'vertex.size' is deprecated. Use `vertex.weight`")
		}
		if (is.null(vertex.weight)) {
			vertex.weight <- as.numeric(table(object@idents))
		}
		pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig,
												 key = "pathway_name", matching.exact = T, pair.only = T)
		if (is.null(signaling.name)) {
			signaling.name <- signaling
		}
		net <- object@net
		pairLR.use.name <- dimnames(net$prob)[[3]]
		pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
		pairLR <- pairLR[pairLR.name, ]
		prob <- net$prob
		pval <- net$pval
		prob[pval > thresh] <- 0
		if (length(pairLR.name) > 1) {
			pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name],
																					 3, sum) != 0]
		}
		else {
			pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name]) !=
																		 	0]
		}
		if (length(pairLR.name.use) == 0) {
			return(NULL)
			#stop(paste0("There is no significant communication of ", 					signaling.name))
		}
		else {
			pairLR <- pairLR[pairLR.name.use, ]
		}
		nRow <- length(pairLR.name.use)
		prob <- prob[, , pairLR.name.use]
		pval <- pval[, , pairLR.name.use]
		if (length(dim(prob)) == 2) {
			prob <- replicate(1, prob, simplify = "array")
			pval <- replicate(1, pval, simplify = "array")
		}
		prob.sum <- apply(prob, c(1, 2), sum)

		prob.sum
	}

	cellchat_process_sample_signal = function (object, signaling = NULL, pattern = c("outgoing",
																																									 "incoming", "all"), slot.name = "netP", color.use = NULL,
																						 color.heatmap = "BuGn", title = NULL, width = 10, height = 8,
																						 font.size = 8, font.size.title = 10, cluster.rows = FALSE,
																						 cluster.cols = FALSE)
	{
		pattern <- match.arg(pattern)
		if (length(slot(object, slot.name)$centr) == 0) {
			stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
		}
		centr <- slot(object, slot.name)$centr
		outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
		incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
		dimnames(outgoing) <- list(levels(object@idents), names(centr))
		dimnames(incoming) <- dimnames(outgoing)
		for (i in 1:length(centr)) {
			outgoing[, i] <- centr[[i]]$outdeg
			incoming[, i] <- centr[[i]]$indeg
		}
		if (pattern == "outgoing") {
			mat <- t(outgoing)
			legend.name <- "Outgoing"
		}
		else if (pattern == "incoming") {
			mat <- t(incoming)
			legend.name <- "Incoming"
		}
		else if (pattern == "all") {
			mat <- t(outgoing + incoming)
			legend.name <- "Overall"
		}
		if (is.null(title)) {
			title <- paste0(legend.name, " signaling patterns")
		}
		else {
			title <- paste0(paste0(legend.name, " signaling patterns"),
											" - ", title)
		}
		if (!is.null(signaling)) {
			mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
			mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
			idx <- match(rownames(mat1), signaling)
			mat[idx[!is.na(idx)], ] <- mat1
			dimnames(mat) <- list(signaling, colnames(mat1))
		}
		mat.ori <- mat
		mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
		mat[mat == 0] <- NA

		mat %>%
			as_tibble(rownames = "gene") %>%
			gather(cell_type, value, -gene) %>%
			mutate(value = if_else(value %in% c(NaN, NA), 0, value))
	}

	computeCommunProbPathway= function (object = NULL, net = NULL, pairLR.use = NULL, thresh = 0.05)
	{
		if (is.null(net)) {
			net <- object@net
		}
		if (is.null(pairLR.use)) {
			pairLR.use <- object@LR$LRsig
		}
		prob <- net$prob
		prob[net$pval > thresh] <- 0
		pathways <- unique(pairLR.use$pathway_name)
		group <- factor(pairLR.use$pathway_name, levels = pathways)

		# STEFANO FIX
		if(length(levels(group))==1){
			xx = apply(prob, c(1, 2), by, group, sum)
			prob.pathways = xx |> array(dim = c(nrow(xx), ncol(xx), 1), dimnames = list(rownames(xx), colnames(xx), levels(group)))
		}
			else
				prob.pathways <- aperm(apply(prob, c(1, 2), by, group, sum),
															 c(2, 3, 1))



		pathways.sig <- pathways[apply(prob.pathways, 3, sum) !=     0]
		prob.pathways.sig <- prob.pathways[, , pathways.sig, drop=FALSE]
		idx <- sort(apply(prob.pathways.sig, 3, sum), decreasing = TRUE,
								index.return = TRUE)$ix
		pathways.sig <- pathways.sig[idx]
		prob.pathways.sig <- prob.pathways.sig[, , idx, drop=FALSE]
		if (is.null(object)) {
			netP = list(pathways = pathways.sig, prob = prob.pathways.sig)
			return(netP)
		}
		else {
			object@netP$pathways <- pathways.sig
			object@netP$prob <- prob.pathways.sig
			return(object)
		}
	}

	computeCommunProb = function (object, type = c("triMean", "truncatedMean", "thresholdedMean",
														 "median"), trim = 0.1, LR.use = NULL, raw.use = TRUE, population.size = FALSE,
						distance.use = TRUE, interaction.length = 200, scale.distance = 0.01,
						k.min = 10, nboot = 100, seed.use = 1L, Kh = 0.5, n = 1)
	{
		type <- match.arg(type)
		cat(type, "is used for calculating the average gene expression per cell group.",
				"\n")
		FunMean <- switch(type, triMean = triMean, truncatedMean = function(x) mean(x,
																																								trim = trim, na.rm = TRUE), median = function(x) median(x,
																																																																				na.rm = TRUE))
		if (raw.use) {
			data <- as.matrix(object@data.signaling)
		}
		else {
			data <- object@data.project
		}
		if (is.null(LR.use)) {
			pairLR.use <- object@LR$LRsig
		}
		else {
			pairLR.use <- LR.use
		}
		complex_input <- object@DB$complex
		cofactor_input <- object@DB$cofactor
		my.sapply <-  sapply
		ptm = Sys.time()
		pairLRsig <- pairLR.use
		group <- object@idents
		geneL <- as.character(pairLRsig$ligand)
		geneR <- as.character(pairLRsig$receptor)
		nLR <- nrow(pairLRsig)
		numCluster <- nlevels(group)
		if (numCluster != length(unique(group))) {
			stop("Please check `unique(object@idents)` and ensure that the factor levels are correct!\n         You may need to drop unused levels using 'droplevels' function. e.g.,\n         `meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))`")
		}
		data.use <- data/max(data)
		nC <- ncol(data.use)
		data.use.avg <- aggregate(t(data.use), list(group), FUN = FunMean)
		data.use.avg <- t(data.use.avg[, -1])
		colnames(data.use.avg) <- levels(group)
		dataLavg <- computeExpr_LR(geneL, data.use.avg, complex_input)
		dataRavg <- computeExpr_LR(geneR, data.use.avg, complex_input)

		# STEFANO ADDED THIS
		dataRavg[dataRavg=="NaN"]=0
		dataLavg[dataLavg=="NaN"]=0

		dataRavg.co.A.receptor <- computeExpr_coreceptor(cofactor_input,
																										 data.use.avg, pairLRsig, type = "A")
		dataRavg.co.I.receptor <- computeExpr_coreceptor(cofactor_input,
																										 data.use.avg, pairLRsig, type = "I")
		dataRavg <- dataRavg * dataRavg.co.A.receptor/dataRavg.co.I.receptor
		dataLavg2 <- t(replicate(nrow(dataLavg), as.numeric(table(group))/nC))
		dataRavg2 <- dataLavg2
		index.agonist <- which(!is.na(pairLRsig$agonist) & pairLRsig$agonist !=
													 	"")
		index.antagonist <- which(!is.na(pairLRsig$antagonist) &
																pairLRsig$antagonist != "")
		if (object@options$datatype != "RNA") {
			data.spatial <- object@images$coordinates
			spot.size.fullres <- object@images$scale.factors$spot
			spot.size <- object@images$scale.factors$spot.diameter
			d.spatial <- computeRegionDistance(coordinates = data.spatial,
																				 group = group, trim = trim, interaction.length = interaction.length,
																				 spot.size = spot.size, spot.size.fullres = spot.size.fullres,
																				 k.min = k.min)
			if (distance.use) {
				print(paste0(">>> Run CellChat on spatial imaging data using distances as constraints <<< [",
										 Sys.time(), "]"))
				d.spatial <- d.spatial * scale.distance
				diag(d.spatial) <- NaN
				cat("The suggested minimum value of scaled distances is in [1,2], and the calculated value here is ",
						min(d.spatial, na.rm = TRUE), "\n")
				if (min(d.spatial, na.rm = TRUE) < 1) {
					stop("Please increase the value of `scale.distance` and check the suggested values in the parameter description (e.g., 1, 0.1, 0.01, 0.001, 0.11, 0.011)")
				}
				P.spatial <- 1/d.spatial
				P.spatial[is.na(d.spatial)] <- 0
				diag(P.spatial) <- max(P.spatial)
				d.spatial <- d.spatial/scale.distance
			}
			else {
				print(paste0(">>> Run CellChat on spatial imaging data without distances as constraints <<< [",
										 Sys.time(), "]"))
				P.spatial <- matrix(1, nrow = numCluster, ncol = numCluster)
				P.spatial[is.na(d.spatial)] <- 0
			}
		}
		else {
			print(paste0(">>> Run CellChat on sc/snRNA-seq data <<< [",
									 Sys.time(), "]"))
			d.spatial <- matrix(NaN, nrow = numCluster, ncol = numCluster)
			P.spatial <- matrix(1, nrow = numCluster, ncol = numCluster)
			distance.use = NULL
			interaction.length = NULL
			spot.size = NULL
			spot.size.fullres = NULL
			k.min = NULL
		}
		Prob <- array(0, dim = c(numCluster, numCluster, nLR))
		Pval <- array(0, dim = c(numCluster, numCluster, nLR))
		set.seed(seed.use)
		permutation <- replicate(nboot, sample.int(nC, size = nC))
		data.use.avg.boot <- my.sapply(X = 1:nboot, FUN = function(nE) {
			groupboot <- group[permutation[, nE]]
			data.use.avgB <- aggregate(t(data.use), list(groupboot),
																 FUN = FunMean)
			data.use.avgB <- t(data.use.avgB[, -1])
			return(data.use.avgB)
		}, simplify = FALSE)
		pb <- txtProgressBar(min = 0, max = nLR, style = 3, file = stderr())
		for (i in 1:nLR) {
			dataLR <- Matrix::crossprod(matrix(dataLavg[i, ], nrow = 1),
																	matrix(dataRavg[i, ], nrow = 1))
			P1 <- dataLR^n/(Kh^n + dataLR^n)
			P1_Pspatial <- P1 * P.spatial
			if (sum(P1_Pspatial) == 0) {
				Pnull = P1_Pspatial
				Prob[, , i] <- Pnull
				p = 1
				Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster,
															byrow = FALSE)
			}
			else {
				if (is.element(i, index.agonist)) {
					data.agonist <- computeExpr_agonist(data.use = data.use.avg,
																							pairLRsig, cofactor_input, index.agonist = i,
																							Kh = Kh, n = n)
					P2 <- Matrix::crossprod(matrix(data.agonist,
																				 nrow = 1))
				}
				else {
					P2 <- matrix(1, nrow = numCluster, ncol = numCluster)
				}
				if (is.element(i, index.antagonist)) {
					data.antagonist <- computeExpr_antagonist(data.use = data.use.avg,
																										pairLRsig, cofactor_input, index.antagonist = i,
																										Kh = Kh, n = n)
					P3 <- Matrix::crossprod(matrix(data.antagonist,
																				 nrow = 1))
				}
				else {
					P3 <- matrix(1, nrow = numCluster, ncol = numCluster)
				}
				if (population.size) {
					P4 <- Matrix::crossprod(matrix(dataLavg2[i,
					], nrow = 1), matrix(dataRavg2[i, ], nrow = 1))
				}
				else {
					P4 <- matrix(1, nrow = numCluster, ncol = numCluster)
				}
				Pnull = P1 * P2 * P3 * P4 * P.spatial
				Prob[, , i] <- Pnull
				Pnull <- as.vector(Pnull)
				Pboot <- sapply(X = 1:nboot, FUN = function(nE) {
					data.use.avgB <- data.use.avg.boot[[nE]]
					dataLavgB <- computeExpr_LR(geneL[i], data.use.avgB,
																			complex_input)
					dataRavgB <- computeExpr_LR(geneR[i], data.use.avgB,
																			complex_input)
					dataRavgB.co.A.receptor <- computeExpr_coreceptor(cofactor_input,
																														data.use.avgB, pairLRsig[i, , drop = FALSE],
																														type = "A")
					dataRavgB.co.I.receptor <- computeExpr_coreceptor(cofactor_input,
																														data.use.avgB, pairLRsig[i, , drop = FALSE],
																														type = "I")
					dataRavgB <- dataRavgB * dataRavgB.co.A.receptor/dataRavgB.co.I.receptor
					dataLRB = Matrix::crossprod(dataLavgB, dataRavgB)
					P1.boot <- dataLRB^n/(Kh^n + dataLRB^n)
					if (is.element(i, index.agonist)) {
						data.agonist <- computeExpr_agonist(data.use = data.use.avgB,
																								pairLRsig, cofactor_input, index.agonist = i,
																								Kh = Kh, n = n)
						P2.boot <- Matrix::crossprod(matrix(data.agonist,
																								nrow = 1))
					}
					else {
						P2.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
					}
					if (is.element(i, index.antagonist)) {
						data.antagonist <- computeExpr_antagonist(data.use = data.use.avgB,
																											pairLRsig, cofactor_input, index.antagonist = i,
																											Kh = Kh, n = n)
						P3.boot <- Matrix::crossprod(matrix(data.antagonist,
																								nrow = 1))
					}
					else {
						P3.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
					}
					if (population.size) {
						groupboot <- group[permutation[, nE]]
						dataLavg2B <- as.numeric(table(groupboot))/nC
						dataLavg2B <- matrix(dataLavg2B, nrow = 1)
						dataRavg2B <- dataLavg2B
						P4.boot = Matrix::crossprod(dataLavg2B, dataRavg2B)
					}
					else {
						P4.boot = matrix(1, nrow = numCluster, ncol = numCluster)
					}
					Pboot = P1.boot * P2.boot * P3.boot * P4.boot *
						P.spatial
					return(as.vector(Pboot))
				})
				Pboot <- matrix(unlist(Pboot), nrow = length(Pnull),
												ncol = nboot, byrow = FALSE)
				nReject <- rowSums(Pboot - Pnull > 0)
				p = nReject/nboot
				Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster,
															byrow = FALSE)
			}
			setTxtProgressBar(pb = pb, value = i)
		}
		close(con = pb)
		Pval[Prob == 0] <- 1
		dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
		dimnames(Pval) <- dimnames(Prob)
		net <- list(prob = Prob, pval = Pval)
		execution.time = Sys.time() - ptm
		object@options$run.time <- as.numeric(execution.time, units = "secs")
		object@options$parameter <- list(type.mean = type, trim = trim,
																		 raw.use = raw.use, population.size = population.size,
																		 nboot = nboot, seed.use = seed.use, Kh = Kh, n = n,
																		 distance.use = distance.use, interaction.length = interaction.length,
																		 spot.size = spot.size, spot.size.fullres = spot.size.fullres,
																		 k.min = k.min)
		if (object@options$datatype != "RNA") {
			object@images$distance <- d.spatial
		}
		object@net <- net
		print(paste0(">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [",
								 Sys.time(), "]"))
		return(object)
	}

	netAnalysis_computeCentrality = function (object = NULL, slot.name = "netP", net = NULL, net.name = NULL,
						thresh = 0.05)
	{
		if (is.null(net)) {
			prob <- methods::slot(object, slot.name)$prob
			pval <- methods::slot(object, slot.name)$pval
			pval[prob == 0] <- 1
			prob[pval >= thresh] <- 0
			net = prob
		}
		if (is.null(net.name)) {
			net.name <- dimnames(net)[[3]]
		}
		if (length(dim(net)) == 3) {
			nrun <- dim(net)[3]
			my.sapply <- ifelse(test = future::nbrOfWorkers() ==
														1, yes = pbapply::pbsapply, no = future.apply::future_sapply)
			centr.all = my.sapply(X = 1:nrun, FUN = function(x) {
				net0 <- net[, , x]

				# ADDED BY STEFANO
				net0[net0<0] = 0

				return(CellChat:::computeCentralityLocal(net0))
			}, simplify = FALSE)
		}
		else {
			centr.all <- as.list(CellChat:::computeCentralityLocal(net))
		}
		names(centr.all) <- net.name
		if (is.null(object)) {
			return(centr.all)
		}
		else {
			slot(object, slot.name)$centr <- centr.all
			return(object)
		}
	}

	communication_results =
		tibble(DB = c("Secreted Signaling", "ECM-Receptor" , "Cell-Cell Contact" )) |>
		mutate(data = list(counts_cellchat)) |>
		mutate(data = map2(
			data, DB,
			~ {
				print(.y)
				.x@DB <- subsetDB(CellChatDB.human, search = .y)

				x = .x |>
					subsetData() |>
					identifyOverExpressedGenes() |>
					identifyOverExpressedInteractions() |>
					projectData(PPI.human)

				if(nrow(x@LR$LRsig)==0) return(NA)

				x |>
					computeCommunProb() |>
					filterCommunication() |>
					computeCommunProbPathway() |>
					aggregateNet()

			}
		)) |>

		# Record sample
		mutate(sample = my_sample) |>

		# Add histogram
		mutate(tot_interactions = map2_dbl(
			data, DB,
			~  .x |> when(
				!is.na(.) ~ sum(.x@net$count),
				~ 0
			)
		)) |>

		# values_df_for_heatmap
		# Scores for each cell types across all others. How communicative is each cell type
		mutate(cell_vs_all_cells_per_pathway = map2(
			data  , sample,
			~ when(
				.x,
				!is.na(.x) && length(.x@netP$pathways) > 0 ~
					netAnalysis_computeCentrality(., slot.name = "netP") |>
					cellchat_process_sample_signal(
						pattern = "all", signaling = .x@netP$pathways,
						title = .y, width = 5, height = 6, color.heatmap = "OrRd"
					),
				~ tibble(gene = character(), cell_type = character(), value = double())
			)
		))

	genes = communication_results |> select(cell_vs_all_cells_per_pathway) |> unnest(cell_vs_all_cells_per_pathway) |> distinct(gene) |> pull(gene)

	# Hugh resolution
	communication_results |>
		mutate(cell_vs_cell_per_pathway = map(
			data,
			~ {
				my_data = .x

				# Return empty if no results
				if(is.na(my_data)) return(tibble(gene = character(),  result = list()))

				tibble(gene = genes) |>
					mutate(result = map(gene, ~ {

						unparsed_result = cellchat_matrix_for_circle(my_data,  layout = "circle", signaling = .x)

						if(!is.null(unparsed_result))
							unparsed_result |>
							as_tibble(rownames = "cell_type_from") |>
							pivot_longer(-cell_type_from, names_to = "cell_type_to", values_to = "score")

						unparsed_result

					}))
			}
		)) |>

		# Save
		saveRDS(output_path)

}

