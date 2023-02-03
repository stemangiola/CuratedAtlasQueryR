
library(zellkonverter)
library(Seurat)
library(SingleCellExperiment) # load early to avoid masking dplyr::count()
library(tidySingleCellExperiment)
library(dplyr)
library(cellxgenedp)
library(tidyverse)
#library(tidySingleCellExperiment)
library(stringr)
library(scMerge)
library(glue)
library(DelayedArray)
library(HDF5Array)
library(tidyseurat)
library(celldex)
library(SingleR)
library(glmGamPoi)
library(magrittr)
# source("utility.R")


# # CREATE MAKEFILE
tab = "\t"
root_directory = "/vast/projects/RCP/human_cell_atlas"
annotated_data_directory = glue("{root_directory}/annotated_data_0.2")
light_data_directory = glue("{root_directory}/splitted_light_data_0.2")
metadata = glue("{root_directory}/metadata_0.2.rds")
cell_type_df = "metadata_cell_type.csv"

light_file_paths = dir(light_data_directory, full.names = TRUE)
.sample = basename(light_file_paths) |> tools::file_path_sans_ext()
annotated_file_paths = glue("{annotated_data_directory}/{.sample}.rds")
file_for_annotation_workflow = glue("{root_directory}/cell_sample_cell_type_df_for_annotation_workflow.rds")

metadata_df = readRDS(metadata)

metadata_df |> distinct(.cell, .sample, file_id, cell_type) |>
	saveRDS(file_for_annotation_workflow)

metadata_df |>
	left_join(read_csv(cell_type_df)) |>
	filter(lineage_1=="immune") |>
	distinct(.sample, file_id) |>

	# Build chunks
	with_groups(file_id, ~ .x |> rowid_to_column("chunk") ) |>
	mutate(chunk = chunk |> divide_by(10) |> ceiling()) |>

	mutate(
		input_files = glue("{light_data_directory}/{.sample}") |> as.character(),
		output_files =  glue("{annotated_data_directory}/{.sample}.rds") |> as.character()
	) |>

	# Get memory
	mutate(Mb = map_dbl(input_files, ~
												( (file.info(glue("{.x}/se.rds"))$size /1e6) |> as.integer() ) +
												( (file.info(glue("{.x}/assays.h5"))$size /1e6) |> as.integer() )
	)) |>
	mutate(memory = pmax(Mb * 20 )) |>

	#
	# mutate(memory = case_when(
	# 	chunk %in%
	#
	# 		c('80507b531a4665ba331ff13f7996f2e0', 'fcbada6e796bb52c7bc3b9cef6532b1e', 'ff612ffdd90b7c4a082ab85e9b54b924', '6d812c89125c54cf915ba8330d74273a', '7a1c50a4cfd567e8559ad46bbd01c1e0', '8bb69081caaf38a3abbb4407589ab5b6', 'b93d5a0d64f34aada9d372ca70f681e8', 'c265cb6e68b8f747bddb10a56fee19c7', '6e9d4cef82b994569d2c98ce75218743', 'a77d2ea6a66417cf8cd63e7e726cf5f0', '078f35a53e9d14b71a4e2b716f557138', '6327cd460c7e5f3b50e8ec07b45389b8', 'a232ee361076b8092ac20118b6dd5aae', 'cb899529c5bf8e6db8cdab2bb24ce0ca', '7cd9204e50d1b12fb1884ab8b9924960', 'ffbfb187f9d7e743a74b03f30d819d34', 'b8eafa88efa0e3c600e36b52a7da3e35', '73ffcdfb290afd536614a159410c6267', 'bb5f7468db4e3110f3399e5c2ab09350', 'b7d54c16cfad21869c8166c410a62e5e', '486eebeda7bdebce5927d7918eb12df0', '25de28cdd79d796cfecb3d3180e1c677', '2bb2ad830713fc3d7d1bfdfddf4ff742', 'b4f6e5212241973692c9cc2f7aeebb3d', 'cb2d45b523a31d3e3a33dab5dc0bd342', '88ca1f6a444a590e695a95e763ce6dc3', '3d2ef9fd4c49b8eef3634cfdd6391ff7',
	# 			'bd72dd64cddad9218397a1cad2f26d67', 'f18fba81d42fe306475c95a239b675c4', '7c9a311a483f8e09e855d093a64f1cc9', '71b9c02e332f0ac479e6e7172e6f5888', 'e22d7bb7b5fdf9dc5a559092767f43aa', '828d4e3d3dcd8de36fcf21fdb0871112', 'ca70d83f2e1df93a98499d18230cb4f5', 'ce3dbb49a1d2d56bc86e67f3b5cae390', '622b67c70869c4f501308efaa92f08a7', '0a6eddf091d95fccb1d5b6660f8e19ce', '70a021e3c1128b46eba2acb77dc02a05', 'f3b104a8f147e678859ae13fb1e92e39', '9750e782e1085296a5a350819f78af97', 'de87abf02112258da54326de0dad1f13', 'edb8ccfa9f2156b02546824c89e8b4ac', 'd060d78733d18e6d4a1e91b562bbb83a', '91d9efed34f1b27d595a479d2e50c886', 'd5f9ab4d42dac4898b7fd2dfdb90f7f0', '6ff2e5eb72d1e258edbef4d7fa000307', '41564528ffbab06b4c91010e95d7c172', 'e1825cb9f73f0cbe40b0b3a85aa350d5', 'b2dca115cc815b7c5a072b7ae06da968', '0cffaa5dbc6dc5d4039b5738cc8452eb', '42759990c5b59cdc1df3a8b92f0d7514', 'f578ec86c5856bad34aaa117524e0815', '54ba38ed501bd4c8af57f1ed587a2431', '14973806ff12d5a58612a4a34ab6e859', 'c9609b8b49020d9072240ac9a4411770', '6b68d3a7e2a6e5f444853f04202e7406', '0ab56037d72ed9139e60a4b9a02aaae0', 'ed2cf4b0d9be99a1b92fb0063909ac66', 'c4b0f90b9e961b9c1ab2006d587c6c6b', '4fcfaf21f34ffc5af7b0ce633ab0ab80', 'c83f157deace6e1931d274c8368648c0', '6e57e21b10cd40ef13fe63c22b58b2b0', 'fbc148d5367f63bb52cc70cd4d0ffd68', '30f109f93bff686b8c1f7bb3f48ff8af', 'f6d1dc92da6cca9960875e392cc8b420', 'fd96b94edeea6ab72ccd7d57a4b503d4', '6464c5c5220928b6cc0b38ada7fb191a', '2d0cfac6f34789c4bf9ca33a983c37a0', 'fdc44855792439d966accdbf7df9cac9',
	# 			'85275069150b6110ad5db5410f1a89e7', '516350a4f6a2d659b2204c1b0b1ac533', '1175a3cd6b2579e263324b639f12b1b1', 'fd6ca7ea36b6b478f905f091c62583c3', '9fd65311516c59f984276101d0b80318', 'c41cb335fc85785d8eca7513c35ebd2d', '2db63e5254393a64cbddec0a6721afee', '5c3d0504c1226f51ae021e6a190d7704', '680d389dbbe7cadea204672c28e9a449', '5021486f9dbc708f8c33f088f26e3758', '7dcceaa2fb865568f73494f5e2b3aa8f', 'e8fe4460a4af317fc2849c4f2077547a', '1afb07303d1ae1d0b52a5ba7a245075d', 'ef6876826a441cbbb2c7c84ebe85d325', 'a657082b863a036686415659a295a66a', '0adc338efd606eb4691aac68bffb97f8', 'd835a01a7f879b62db0a61b8ab574d78', '79c02d7a10ff735fe64dfc5bf88b8f41', '40fba5da714de22ae3aca0d1ca41a80b', '944b2d55e96bd62f205bd006754d3322', '9843e641df41e746a68c3df99a0da789', 'a245885261b384ce0261082eb7a0b229', 'e74a84462d89b7db4a7d07011c42b03c', 'a77e0fb892c34bdcb561da344a13b8ef', '33d9a38728c6c3c0fde2ff58bf86c44a', '7a9ee4c76bf5f27e1eba348b1d191ea6', '3220304cec9621329939086695c360a4', '256686225636c8e7dfb32c69d8dcfa3c', 'd31eb3a29018429c09baeb78c6edea43', '05e6a4d998c123286beb1dca770a741c', '1f86910976e3c2b983a2823fc427834e', 'd6a43e21614625dddede6dc631fd82f0', '783a3b9af01864747d7de20dc55262f6', '88d0a7ca73b95c483a879d67b57cb436', '3b3ed94e2932f27653e814135f55b2c9', 'cbd9ddf771e2f427d2054d491b47ee57', '93d1f9df324981115727aae376278726', '4ff1f53e149edbd0bae956d5d93943e3', '7cfef32e9318fcff73c2bd2884e3cff2', '0583b837fee41bd285c9f29462566c04', 'b4e8336d9e5e4d6bbf30293c07e001b8', '860fdf1db1b377ee854557afbf5c797d', '22d4b286b892ef8b5138072bd9566cd2', '18631218fc1e8c8b8201f5eebe637eeb', '98d25b1a524d63fac90f670dbd2a99a8', 'a5401d0c2f3a6ec0de28a5bf02845fc7', 'd7056438ec8650c2862f0417903863cd',
	# 			'e2830f5e44b5d390e7d7bbf5ba4fc9fe', '0460dc7b199f9265ed173f3ac6ffdd03', '871a9b8c0ec355163d594a4da160bc0c', '1f291c12284a4ae96e154793c811d1d8', 'd2480ab8bce610edee93722824dd7c6c', '6f1aad9a7d84c9ca6a107a6b530175a6', 'bc5f32bc5bf3f46142123f146f1ca136', 'c6906bd160cf94942422974ebf3b57e3', '5572c1c9acaeb29de4764f22cdc7e566', '05aed84debb5078d9aa6c8b20ff62663', '67bb75769570b3a47ee3e83f3998bad2', 'a1d98925e8acd3e7131e677e8718c325', '4ed49a161d65b6dae132877112dbccd9', '6d505bad20d48653f1d095856ef2eb81', 'dd28df166c8bc00f4af59c5687edc5c9', 'a7a375fb168bd68b5eae2a9b69116a4c', 'c7df9c073269959dcea83877c2bea3b7', '8f47d0f5f0a50856b0aebfd0d32c3024', 'a9852fc70bb43641991126f95231051c', '918f1c89b13b74148f21f6b5697c509e', '094b9ba158a3c0728c7d7ace42e0b8f8', '0d5d079fd4fc3743b288fa472e8053d8', 'e525bb6af01cd8659e94df65d6b05cfc', 'a97bf1010beb6c6916316a782cdda57e', '43f3cced1cbbfb710499449703d63d21', 'e98c033fd7a510119bfdb7b48e88f0af', '15e5c907583d76232d7ff970c3a1c7cd', '45160921927eda11f62a7d06129c66ee', '572d4065692d96070e3dcac2868fd206', 'd792a7d0a4d3ce02711f08b5c556913d', '13aed27c8e15d8683a26590275977d14', '896ac37da8f88e33f05ebd5283d4b806', '34b4c632e6befd44d9741d3b5571fa29', '64fa87c0da688c9b4e1722c5ba66b555', '7c5a341956a86aae089954a057c010f9', '1bf6aa8f66fe2383c0fc7e9cb68d328f', '38d2b7349e1f0771f3205725d1f0dbc8', '474fe9670783f8d1a551ff96f0cb5c7a', '4c9f3c2b964ab525ed5b8f68d9ec9d71', 'e9f517c60c53c26fa13b479b5bd4ac43', 'a27097a860f150681d686fc7b8aeb3b4', '8cb6533eeea450f4b2d706c8bb8f3c2c', '58d0f39601731f290778a3b51766ec6d', 'aafb785305edb6bf99bba386dfb7edef', 'b25a5f6f8196daafbe940daf95691e97', '10000a7099938e8483d1cc832aa17b81', '6bfcb65bb93b7690c41bc57ae6f271c0',
	# 			'd461b1af552977ce8755d9d15c15308e',
	# 			'ac8b2e1209f000a1aabb0919ba08ae73', 'd899b5bc64bc0ea8defa003750ab96e7'	) ~ 40000,
	# 	TRUE ~ memory
	# )) |>

	nest(data = -c(file_id, chunk)) |>
	mutate(
		input_file_path = map_chr(data, ~ .x$input_files |> str_c(collapse = " ")),
		output_file_path = map_chr(data, ~ .x$output_files |> str_c(collapse = " "))
	) |>

	# Sum memory
	mutate(total_memory = map_dbl(data, ~ .x$memory |> sum() |> sum(10000) |> max(30000) |>  min(160000))) |>

	left_join(
		tibble::tribble(
			~chunk,                               ~file_id,
			21, "07beec85-51be-4d73-bb80-8f85b7b643d5",
			29, "07beec85-51be-4d73-bb80-8f85b7b643d5",
			32, "07beec85-51be-4d73-bb80-8f85b7b643d5",
			36, "08247324-7ab7-45c5-8bd6-6c22676761ed",
			38, "08247324-7ab7-45c5-8bd6-6c22676761ed",
			40, "08247324-7ab7-45c5-8bd6-6c22676761ed",
			41, "08247324-7ab7-45c5-8bd6-6c22676761ed",
			52, "09132373-0ea7-4d8b-add8-9b0717781109",
			54, "09132373-0ea7-4d8b-add8-9b0717781109",
			57, "09132373-0ea7-4d8b-add8-9b0717781109",
			88, "1042ba0a-98c5-4816-897d-e192eb9303e3",
			91, "12271708-2ba9-4073-9d03-d0416f157f65",
			92, "14b3d1f9-4949-4330-9f9c-e64b68b2197a",
			94, "14b3d1f9-4949-4330-9f9c-e64b68b2197a",
			96, "14b3d1f9-4949-4330-9f9c-e64b68b2197a",
			98, "14b3d1f9-4949-4330-9f9c-e64b68b2197a",
			100, "14b3d1f9-4949-4330-9f9c-e64b68b2197a",
			101, "14b3d1f9-4949-4330-9f9c-e64b68b2197a",
			102, "14b3d1f9-4949-4330-9f9c-e64b68b2197a",
			105, "14b3d1f9-4949-4330-9f9c-e64b68b2197a",
			129, "3431ab62-b11d-445f-a461-1408d2b29f8c",
			175, "5500774a-6ebe-4ddf-adce-90302b7cd007",
			176, "5500774a-6ebe-4ddf-adce-90302b7cd007",
			304, "6661ab3a-792a-4682-b58c-4afb98b2c016",
			582, "c0dca32e-fa64-4448-bc12-2b8f16702c29"
		) |> mutate(memory_up = "yes")
	) |>

	mutate(total_memory = if_else(memory_up=="yes", max(100000, total_memory), total_memory)) |>

	# mutate(
	# 	output_file_path = glue("{annotated_data_directory}/{.sample}.rds" |> as.character())
	# ) |>

	rowid_to_column() |>
	mutate(commands = pmap(list(output_file_path, input_file_path,  total_memory, rowid, file_id), ~
												 	c(
												 		glue("CATEGORY=light_data{..4}\nMEMORY={..3}\nCORES=2"),
												 		glue("{..1}:{..2} {file_for_annotation_workflow} {cell_type_df}\n{tab}Rscript annotate_files.R {..2} {..5} {light_data_directory} {file_for_annotation_workflow} {cell_type_df} {annotated_data_directory} {..3}")
												 	)
	))  |>
	pull(commands) |>
	unlist() |>
	write_lines(glue("~/PostDoc/HCAquery/dev/annotate_files.makeflow"))



