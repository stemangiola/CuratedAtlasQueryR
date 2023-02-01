library(zellkonverter)
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
library(openssl)
library(stringr)
library(HCAquery)
library(purrr)
## source("utility.R")
options(scipen = 999)
#

# # CREATE MAKEFILE
# tab = "\t"
# root_directory = "/vast/projects/RCP/human_cell_atlas"
# raw_data_directory = glue("{root_directory}/raw_data")
# splitted_DB2_data_directory = glue("{root_directory}/splitted_DB2_data_0.2")
# file_cell_types_directory = glue("{root_directory}/file_cell_types")
# input_files_path = dir(file_cell_types_directory, full.names = TRUE)
# gene_names = glue("{root_directory}/gene_names.rds")
#
#
# ## metadata = readRDS(metadata_path)
#
# get_metadata("/vast/projects/RCP/human_cell_atlas/metadata_annotated_0.2.sqlite") |>
# 	distinct(file_id, cell_type) |>
# 	as_tibble() |>
#
# 	unite("file_id_db2", c(file_id, cell_type), remove = FALSE) |>
# 	mutate(file_id_db2 = file_id_db2 |> md5() |> as.character()) |>
#
# 	mutate(
# 		input_file_path = glue("{raw_data_directory}/{file_id}.H5AD") |> as.character(),
# 		output_file_path = glue("{splitted_DB2_data_directory}/{file_id_db2}" |> as.character())
# 	) |>
#
# 	mutate(Mb = map_dbl(input_file_path, ~ (file.info(.x)$size / 1e6) |> as.integer())) |>
# 	mutate(memory = pmax(Mb * 2, 30000)) |>
# 	mutate(
# 		memory = case_when(
# 			file_id_db2 %in% c(
# 		'643d1a28d35ed6d9f087ed8732d0fbf0', 'f7cb9ed27af8bd19b8b05feeaf84ec3c', 'abb6ff154b70c74e395add7dca444bbc', '6690e8e6b89123b1e11eef9d4a153e49', '5734fb9bd18f48c88f1bab6e1ea6e5f7', 'd9f28c596bb5813d7da2e690b0473570', 'fea1eba9397c3763e77dc321ac1c35b8', '6b312db0bbe47249bf6dbb382a0f25aa', 'a252f45beddcd62c5ccc0894461861e2', '2d7cc7e650275abf80b672febbdbb47f', '64b9770f1679dc9a6d8b7b12613391ef', '133f0de6e86d712c079af21344ae4501', '2d6fa434060cd4e16987245c35294fdb', '5754ab342ca2f630f442c6fb11f7977f', 'b7411ae6d2d5025fe8735c9f5635bed3', '5a7c941fcf37216ce181e5af9857a897', 'b48aab951c5ac53a9b9b4250c309cb04', '9fe44e0afee2d48bbe9b9922490c4fdf', 'c4cda7f034432e9f88588905e96b8878', 'cb0fa793cff47ff26de5ebcf9d75d8d8', 'ec2cf3480224c44c2f87875b8a383416', '6c3915e70f391f4cc6d8460b8237e686', 'ecc5c83172b38a4c2d703cef9fee6600', '19bc78159fda9d0795ff82201a493627', '5e0dba8e536570c43c0e79bb4695671b', '680201e2dae91ed66c4cb2da45590607', '94c522ee71123a3b115bdebce647fb88', 'b448edda111f93495eb9eb0488c95613', 'ad33c246db845178ffa1d6d1de15946e', '672f90e815461cf1fde0fed7d83528b3', '95b234bb577816c848af6bdbf4e5f639', '5bf8f63f1b8435f51db88701578c9944', 'e604de2b6508a0056cf187b1790e2a2e', '7df8990ca100ac47df3dc857a739a531', 'a983aca109e8ef8cdfa2befbfec8a85d', 'e6cf0c9d608465fe69fc477e5d24a4da', 'a224c2660000d6fbc79296b31ddb05dd', '321ab8ef8102d8e9817588539216fcb2', 'd8fe285941ee94128c466ffe2fc72214', 'f464369e48799c7225ce613d4fac64e4', '60981c2a4aba42bfb81d4d52eb7d91c4', '2b2ca0d6dce081902657f28e6942a546', 'd48585f19116758ac55715d36af29faf', '2763074a04eafc1590ec8e4195c2d30b', '894c26dfae23520d5b1ec13681b973d4', 'bb270e6b0aea2bba1f8f044402b50d0f', '3837b82ce7110024d6010eb668c4d65f', 'afd38219dfc9fa010b181a6f1dc167df', 'f29f36400dcaa5b550e055dee3a6bfab', '0031d2758f70b4602aa5980a96a96a12', '044fc974c7f546b2f61d15ae5c804ddf', '562594615b5d41f9ebe52c2489325045', '9f7827abf98794f8c37d87e5de01691b', '1994ae237838e3271488f98a7a5afc37', 'aa9b46796d3d0e7cf774c12e34cc872f', '2bc6ef5c0f153cdaba618788f75dff13', '124aa1104d8fed6fa9e98fa2ca5c00eb', 'c032b8738a8d06e9f261f76564de0ec8', 'c54f874f27a35cce97fcb3c67701c8c3', 'bc89346be5f9ffc76d8b920f117f244d', '136bcc485075f101e9d109d0222f7211', '9ed4da7076000ffac93a812b075812da', '90099dad83392f82d7cccaff14ab21a8', '2bcc0d793ed9db6ac939c86f4fa67044', 'e29a49f612b72c770276406500392790', '88eeb885ba02a8f5a30202a65f2314d2', '93cfef1122d69293a9f19c054e031b0a', 'cfe2deeeb735797b871828d9860eff41', '11b27ed03f34435b1ed9b513705a99cf', 'b72533ee50311aded927f2d4acc2d696', '63473d2b4aabefe0385fc0d8a0378115', '6909d30b9f20cbb1c4f22a596918bd69', '4c39fdd1e8eb9a146f6e0766a8fbe4c7', 'ff572a8430b657bfc79fb6c6c3760ac1', 'f01c20412cab27a5e8feacd170cfc900', '0ae84745b3261a021cc0cfe055ca8c46', 'b61dc9d6067dc0aae0db5007e4d7f876', '1b855ddca92afe1bed522bb07dc5f788', '8aff79b4657c0aff54c9016ac1eef46e', '1d85feaf473d55854639514693504c02', 'b2377c149a40b99ec331dbf5ab1cf4ce', '77a8cd32e8188e707ed83575fcbed904', '7e8e47f7eb2537368d06cec14ac3a458', '03ff63330952859dbf86f41a48df6a83', '96c252abaf78877e2f400d838bf69166', '42b08eb49a5390b9e66ddc98663eb269', '148e3d4c60692354e3d1948519e4a0ea', 'f92e17f4efd838108af21129becd7b19', '5b56fc4fb694c62d4de406c24c7ca5e7', '63fbfcf952e803f5c39473fafb599251',
# 		'b531a86d7188be31428a03603d7405c6', '768700aa330a2a5891ad22b67ea49dc1', '2f36daf2a3a8aedd48db06acd6727b0a', '01a6385ee03ac917f245d7e52faa9f66', '02ffa11eb4dbc3106eda5dbdec29e817', '542e6a445f3a4d9dad4e32acc525f2c4', '3513dd06f3363335e0561979a5e5f2d1', '3c88cfffed093f566b06783c7033948f', '45c95a5fb1f303c138b0f60b54a94c2c', '848c8bf283e6c39555ed5085b3eba735', '2b9f5b251111939dd8689ddf86309c7b', '93ce665d647aa0319fc444596eb2bb5a', 'e61d313adfe3604cafc773210f461fc1', 'c6ebf4689bb0b3fe66b0139ccaf83ac2', '0b934d2125cf7ca9564ccdab65797c11', '396a6963d3b402548b38676200b9a544', 'dcc9e9def75bc9bae1b893da13daf378', 'c4aafa1124d87da0606bc67bb125d702', '018d85cd9702db9339920f0cb91d7afc', 'b628ccb3d8721fc019ae2ef279375a08', 'ca1de9c53113407dc7373bd33cad2a57', '92852e108c7612c0d7ff4ceccd59bdca', '133f20d89ba49eb68a7e231e590a0ec7', '90679d805829232802b794baf6a09338', '255400de1eda3bec328dcb932a0ce70e', 'ea03efb1ec2466ff0c3925a2c25255bc', '69bd89af9e84e240de4ae1a9c4730b78', '4b92fdac8484757e1ed805f2ef485aff', 'cea187ba19982b64b651046d0c55d683', '09a1bc602a7ad508f10b1c83785e3f5f', 'e27843805c7e5e9743202ad132f0cb74', 'd8216db7dd67453b948025092631e8f0', '76cb634ba33e4438517fa388f472413f', '0e3b0626806414deb561d4c95edf6cc9', 'd6dde5273fdf62bf1776a24007d59ae2', 'f95b23918aa5bcab2c0cc9dcc85e3f2f', 'ba1017fcabbc3f4552554cac3e1e8d51', 'd7096504ae4361e995f97dbcdcdea94e', 'ee701a2474084b3477ce862f2cef54d0', 'b708a95d6bc07e1d79dab68972fde016', '495f8be32619aa791c71947b66a6e44c', 'e48af5a58cf363b65dc2c027616155c3', 'a837f0376b0af3ee2243a26f96003fb6', 'bf522b5a8a05ef06876b1166cf76cd76', 'f2e97e6eb83c6ae6315735c3a58475dd', 'db3dc0285cc2e4866c211e151559df40', '4d957fd55b43b527b17f82a2e8534337', 'c646e9fb94b3331944050967ea820d11', '051579deb85a5b7eabe7c6df3264d302', 'ebdee2282dfd6ea3cb60844087f4c923', '3795117167beb9768034b7443fe5b754', 'cf1cf8debfc581c217707fb09d1c0cc9', '0c34089dbd4db0d48af97873baa8cb10', '867797d630b0a67438a7b1b6ddb9f46f', '44600d9ff423b2f8a674e4fdf928c36f', 'd8036a66b91006fb87312ab73b6f2647', 'd7ab07844e32d4e609e3fae4f6e624b5', '37decab3aa35040d72af0eeff9e557c2', '98435ec4176e33164ac2dd092911518b', 'b957873f07b73bc9fe2858de7560d354', '52643bbd6767cf5904b0bec3eaafe283', '453b598c0fa37d3df507932d8b742eca', '25be440bc59f97942579d25419517781', '48730f8b9ea4a4eb67c4826cec75e6d7', '6d4ceccbec07d4b3f41b8d73a5cc54a4', '099e4f9788ea5646b01899a4565b449b', '3cd5853986f381dc40e02414f76ae0a4', '63f6d25357ac92d889d2a2e08a6d95f1', 'ccf04cce8134efe5e667e2c0caa57474', 'd0f11393c659cf2f70dc989d21a682fa', '5b845723652627151f4ff121333fede1', 'a3802f7db8829109ac98eb4552231d54', '558dd61fb7a6a47ef76b97384898de27', 'cc556cae149faa24f2d011806d65168e', 'c3d01f08f2eb4f769b0fa3497e8222ff', '64bcbbf198f7e191e06173b054cff02a', 'a50e764e1e49e7877bf91864665de354', '4269e17cd671c56a9865b1ed5e316a61', '5a8d3746fadcec538c3d70490cc888ed', '2e5d52ce7ee10ea5bb37769409449e01', 'd0ad238960be64ad4cd8f67fa58a946b', '40d764a54eae749b6852910655e7c35c', '749821212c25ca95abb1831c48764d2a', '1625d480e0c3160c5e34c53006cd85ee', '354f1780bb909fe141b16a5ce155caf8', '51ae4a30442db8f7fcc3f77899370100', '4bb27edd1f01bf191c77271b580a186c', 'e22f00a4abd2d8d9ebe2ca5a6055a0d5', '4d3eeb032b517f91b8ac75babbd8a13b', 'c5db9e3b4ecc205b7d5d09ca660a0c4b', '314dce2d49ecf8b3d83d97f854b471ef', 'd58447f6ceb1c959538846df532bf82a', '346c5e08ea5bafee68a471f0f04fe397', '2c162ef8a4bdfef1ab5dd2d45389b232', 'f55a91bc0057a859f94245b047087dbf', '24dad6be9f859b50fb4f00e1a7981f82', '59a53c1c0078ad9cb88bdbf161c11800', 'c5a05f23f9784a3be3bfa651198a48eb', '18d13cca5e607b917f3c60bd75f863cd', '521b9b73c5b032f8e2809a5032d18266', '18ecc71ecc0df104a1628bdd4e529f85', '7aec18c648d573f34641b8329a2afe3e'
# 			) ~ 60000,
# 			TRUE ~ memory
# 		)
# 	) |>
# 	mutate(
# 		memory = case_when(
# 			file_id_db2 %in% c(
# 				'f7cb9ed27af8bd19b8b05feeaf84ec3c', 'abb6ff154b70c74e395add7dca444bbc', '2d6fa434060cd4e16987245c35294fdb', '5a7c941fcf37216ce181e5af9857a897', '2c162ef8a4bdfef1ab5dd2d45389b232', 'f55a91bc0057a859f94245b047087dbf', '24dad6be9f859b50fb4f00e1a7981f82', 'c5a05f23f9784a3be3bfa651198a48eb'
# 		) ~ 120000,
# 			TRUE ~ memory
# 		)
# 	) |>
# 	rowid_to_column() |>
# 	mutate(commands = pmap(
# 		list(output_file_path, input_file_path,  memory, rowid),
# 		~
# 			c(
# 				glue("CATEGORY=split_data{..4}\nMEMORY={..3}\nCORES=1\nWALL_TIME=30000"),
# 				glue(
# 					"{..1}:{..2}\n{tab}Rscript DB2_files.R {..2} {gene_names} {..1}"
# 				)
# 			)
# 	))  |>
# 	pull(commands) |>
# 	unlist() |>
# 	write_lines(glue("DB2_files.makeflow"))



# Read arguments
args = commandArgs(trailingOnly = TRUE)
input_file = args[[1]]
gene_names = args[[2]]
output_file = args[[3]]

output_file |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)
file_id = basename(input_file) |> tools::file_path_sans_ext() |> str_split("___") %>% .[[1]] %>% .[1]
file_id_db2 = basename(output_file) |> tools::file_path_sans_ext()

cells_to_keep =
	get_metadata() |>
	filter(file_id == !!file_id) |>
	dplyr::select(.cell, .sample, file_id, cell_type) |>
	as_tibble() |>
	unite("file_id_db2", c(file_id, cell_type), remove = FALSE) |>
	mutate(file_id_db2 = file_id_db2 |> md5() |> as.character()) |>
	filter(file_id_db2 == !!file_id_db2) |>
	mutate(.cell_original = .cell |> str_remove(.sample) |> str_remove("_$")) |>
	pull(.cell_original)

# Read
data = readH5AD(input_file,	use_hdf5 = TRUE, reader = "R")
data_for_rows_and_columns = readH5AD(input_file,	use_hdf5 = TRUE, skip_assays = TRUE)
colnames(data) = data_for_rows_and_columns |> colnames()
rownames(data) = data_for_rows_and_columns|> rownames()
rm(data_for_rows_and_columns)
gc()

# Subset cells
data = data[,cells_to_keep] |> select(.cell)
gc()

col_sums = colSums(data@assays@data$X)

# If I have negative values
if((data@assays@data$X[,1:min(10000, ncol(data@assays@data$X))] |> min()) < 0)
	data@assays@data$X[data@assays@data$X<0] <- 0

# If there are huge value cap
if(max(col_sums) > 1e100) {
	temp = data@assays@data$X[,sample(1:ncol(data), size = 10000, replace = T), drop=T]
	q = quantile(temp[temp>0], 0.9)
	data@assays@data$X[data@assays@data$X>1e100] = q
}

# Drop all 0 cells
data = data[,col_sums>0]

# Set gene names
rownames(data) = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
																			 keys = rownames(data),
																			 keytype = "ENSEMBL",
																			 column = "SYMBOL",
																			 multiVals = "first"
)
data = data[!is.na(rownames(data)),]
rownames(data@assays@data$X)  = rownames(data)

# Merge duplicate genes
split_matrix_by_column <- function(matrix, chunk_size) {
  ncols <- ncol(matrix)
  chunks <- lapply(seq(1, ncols, by = chunk_size), function(i) {
    if(i + chunk_size -1 > ncols){
      matrix[, i:ncols, drop = FALSE]
    }else{
      matrix[, i:(i + chunk_size - 1), drop = FALSE]
    }
  })
  return(chunks)
}

duplicated_genes = data[duplicated(rownames(data)),] |> rownames()
deduplicated_matrix =

	data[rownames(data) %in% duplicated_genes,]@assays@data$X %>%
	split_matrix_by_column(10000) |>
	map(~ rowsum(.x, rownames(.x))) |>
	reduce(cbind)	 %>%
	{
		.x = writeHDF5Array(., as.sparse = T)
		rownames(.x) = rownames(.)
		.x
	}

sce = SingleCellExperiment(list(X = rbind(
	data[!rownames(data) %in% duplicated_genes,]@assays@data$X ,
	deduplicated_matrix
	)[rownames(data),,drop=FALSE]
))
colnames(sce) = colnames(data)

rm(data)
rm(deduplicated_matrix)
gc()

# Complete gene set
missing_genes = readRDS(gene_names)  |> setdiff(rownames(sce))

missing_matrix =
	HDF5RealizationSink(c(length(missing_genes),ncol(sce)), as.sparse = TRUE) |>
	as("DelayedArray")
rownames(missing_matrix) = missing_genes

missing_sce = SingleCellExperiment(list(X=missing_matrix))
colnames(missing_sce) = colnames(sce)
missing_sce = missing_sce[!is.na(rownames(missing_sce)),]

# Make cell name unique
sce = sce |> rbind(missing_sce	)
sce = sce[sce |> rownames() |> sort(),]
rm(missing_sce)
gc()

# Tranform
transformation =
	get_metadata() |>
	filter(file_id == !!file_id) |>
	distinct(x_normalization) |>
	as_tibble() |>
	mutate(transformation = case_when(
		x_normalization |> str_detect("log|ln|Log") ~ "log",
		x_normalization |> str_detect("square root") ~ "square_root",
		TRUE ~ "none"
	)) |>
	pull(transformation) |>
	unique()

if(file_id %in% c(
	"1e81a742-e457-4fc6-9c39-c55189ec9dc2",
	"b51bfa2d-22b2-4c65-9803-d36d4de973fa",
	"4e4bbb2d-f341-4523-a5a0-5407d8b03e0e",
	"c48402e4-e7db-4c82-a9e9-51e285e5165c",
	"82ad3285-e5d4-46d1-89c0-3acf91a9e33f",
	"7addb561-c1bf-4fb5-ad10-16dd65b3643a",
	"575513b2-6e53-41e2-85a9-bc08a6233ce4"
)) transformation = "log"



sce@assays@data$X =
	sce@assays@data$X |>
	when(
		transformation=="log" ~ {
			.x = expm1(.)
			#type(.x) = "integer"
			.x
		},
		transformation=="square_root" ~ {
			.x = (.)^2
			#type(.x) = "integer"
			.x
		},
		~ {
			.x = (.)
			#type(.x) = "integer"
			.x
		}
	)


sce |>	saveHDF5SummarizedExperiment(output_file, replace=TRUE)




