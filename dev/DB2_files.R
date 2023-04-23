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
library(CuratedAtlasQueryR)
library(purrr)

library(dbplyr)
library(DBI)
library(duckdb)

## source("utility.R")
options(scipen = 999)

metadata_DB = "/vast/projects/cellxgene_curated/metadata_annotated_0.2.3.parquet"
root_directory = "/vast/projects/cellxgene_curated"


### GET COMMON GENES
# library(zellkonverter)
# library(SingleCellExperiment)
# library(tidyverse)
# library(purrr)
# library(glue)
# library(CuratedAtlasQueryR)
# library(HDF5Array)
# 

# 
# 
# # Read arguments
# 
# samples =
#   duckdb() |>
#   dbConnect(drv = _, read_only = TRUE) |>
#   tbl(metadata_DB) |>
#   distinct(file_id, sample_) |>
#   as_tibble() |>
#   group_by(file_id) |>
#   slice(1) |>
#   pull(sample_)
# 
# # Read gene names
# glue("{root_directory}/splitted_data_0.2") |> 
#   dir( full.names = TRUE) |>
#   str_subset(samples  |> str_c(collapse = "|")) |>
#   imap(	~ {
#     print(.y)
# 
#     AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
#                           keys = rownames(loadHDF5SummarizedExperiment(.x)),
#                           keytype = "ENSEMBL",
#                           column = "SYMBOL",
#                           multiVals = "first"
#     )
#   }
#   )  |>
#   unlist() |>
#   unique() %>%
#   .[!is.na(.)] |>
#   saveRDS(glue("{root_directory}/gene_names.rds"))


#
# # CREATE MAKEFILE
# tab = "\t"
# raw_data_directory = glue("{root_directory}/raw_data")
# splitted_DB2_data_directory = glue("{root_directory}/splitted_DB2_data_0.2.1")
# file_cell_types_directory = glue("{root_directory}/file_cell_types")
# input_files_path = dir(file_cell_types_directory, full.names = TRUE)
# gene_names = glue("{root_directory}/gene_names.rds")
# 
# 
# duckdb() |>
#   dbConnect(drv = _, read_only = TRUE) |>
#   tbl(metadata_DB) |>
# 	distinct(file_id, cell_type) |>
# 	as_tibble() |>
# 
# 	unite("file_id_db2", c(file_id, cell_type), remove = FALSE) |>
# 	mutate(file_id_db2 = file_id_db2 |> md5() |> as.character()) |>
# 
#   #filter(!file_id_db2 %in% dir(splitted_DB2_data_directory)) |>
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
# 		'b531a86d7188be31428a03603d7405c6', '768700aa330a2a5891ad22b67ea49dc1', '2f36daf2a3a8aedd48db06acd6727b0a', '01a6385ee03ac917f245d7e52faa9f66', '02ffa11eb4dbc3106eda5dbdec29e817', '542e6a445f3a4d9dad4e32acc525f2c4', '3513dd06f3363335e0561979a5e5f2d1', '3c88cfffed093f566b06783c7033948f', '45c95a5fb1f303c138b0f60b54a94c2c', '848c8bf283e6c39555ed5085b3eba735', '2b9f5b251111939dd8689ddf86309c7b', '93ce665d647aa0319fc444596eb2bb5a', 'e61d313adfe3604cafc773210f461fc1', 'c6ebf4689bb0b3fe66b0139ccaf83ac2', '0b934d2125cf7ca9564ccdab65797c11', '396a6963d3b402548b38676200b9a544', 'dcc9e9def75bc9bae1b893da13daf378', 'c4aafa1124d87da0606bc67bb125d702', '018d85cd9702db9339920f0cb91d7afc', 'b628ccb3d8721fc019ae2ef279375a08', 'ca1de9c53113407dc7373bd33cad2a57', '92852e108c7612c0d7ff4ceccd59bdca', '133f20d89ba49eb68a7e231e590a0ec7', '90679d805829232802b794baf6a09338', '255400de1eda3bec328dcb932a0ce70e', 'ea03efb1ec2466ff0c3925a2c25255bc', '69bd89af9e84e240de4ae1a9c4730b78', '4b92fdac8484757e1ed805f2ef485aff', 'cea187ba19982b64b651046d0c55d683', '09a1bc602a7ad508f10b1c83785e3f5f', 'e27843805c7e5e9743202ad132f0cb74', 'd8216db7dd67453b948025092631e8f0', '76cb634ba33e4438517fa388f472413f', '0e3b0626806414deb561d4c95edf6cc9', 'd6dde5273fdf62bf1776a24007d59ae2', 'f95b23918aa5bcab2c0cc9dcc85e3f2f', 'ba1017fcabbc3f4552554cac3e1e8d51', 'd7096504ae4361e995f97dbcdcdea94e', 'ee701a2474084b3477ce862f2cef54d0', 'b708a95d6bc07e1d79dab68972fde016', '495f8be32619aa791c71947b66a6e44c', 'e48af5a58cf363b65dc2c027616155c3', 'a837f0376b0af3ee2243a26f96003fb6', 'bf522b5a8a05ef06876b1166cf76cd76', 'f2e97e6eb83c6ae6315735c3a58475dd', 'db3dc0285cc2e4866c211e151559df40', '4d957fd55b43b527b17f82a2e8534337', 'c646e9fb94b3331944050967ea820d11', '051579deb85a5b7eabe7c6df3264d302', 'ebdee2282dfd6ea3cb60844087f4c923', '3795117167beb9768034b7443fe5b754', 'cf1cf8debfc581c217707fb09d1c0cc9', '0c34089dbd4db0d48af97873baa8cb10', '867797d630b0a67438a7b1b6ddb9f46f', '44600d9ff423b2f8a674e4fdf928c36f', 'd8036a66b91006fb87312ab73b6f2647', 'd7ab07844e32d4e609e3fae4f6e624b5', '37decab3aa35040d72af0eeff9e557c2', '98435ec4176e33164ac2dd092911518b', 'b957873f07b73bc9fe2858de7560d354', '52643bbd6767cf5904b0bec3eaafe283', '453b598c0fa37d3df507932d8b742eca', '25be440bc59f97942579d25419517781', '48730f8b9ea4a4eb67c4826cec75e6d7', '6d4ceccbec07d4b3f41b8d73a5cc54a4', '099e4f9788ea5646b01899a4565b449b', '3cd5853986f381dc40e02414f76ae0a4', '63f6d25357ac92d889d2a2e08a6d95f1', 'ccf04cce8134efe5e667e2c0caa57474', 'd0f11393c659cf2f70dc989d21a682fa', '5b845723652627151f4ff121333fede1', 'a3802f7db8829109ac98eb4552231d54', '558dd61fb7a6a47ef76b97384898de27', 'cc556cae149faa24f2d011806d65168e', 'c3d01f08f2eb4f769b0fa3497e8222ff', '64bcbbf198f7e191e06173b054cff02a', 'a50e764e1e49e7877bf91864665de354', '4269e17cd671c56a9865b1ed5e316a61', '5a8d3746fadcec538c3d70490cc888ed', '2e5d52ce7ee10ea5bb37769409449e01', 'd0ad238960be64ad4cd8f67fa58a946b', '40d764a54eae749b6852910655e7c35c', '749821212c25ca95abb1831c48764d2a', '1625d480e0c3160c5e34c53006cd85ee', '354f1780bb909fe141b16a5ce155caf8', '51ae4a30442db8f7fcc3f77899370100', '4bb27edd1f01bf191c77271b580a186c', 'e22f00a4abd2d8d9ebe2ca5a6055a0d5', '4d3eeb032b517f91b8ac75babbd8a13b', 'c5db9e3b4ecc205b7d5d09ca660a0c4b', '314dce2d49ecf8b3d83d97f854b471ef', 'd58447f6ceb1c959538846df532bf82a', '346c5e08ea5bafee68a471f0f04fe397', '2c162ef8a4bdfef1ab5dd2d45389b232', 'f55a91bc0057a859f94245b047087dbf', '24dad6be9f859b50fb4f00e1a7981f82', '59a53c1c0078ad9cb88bdbf161c11800', 'c5a05f23f9784a3be3bfa651198a48eb', '18d13cca5e607b917f3c60bd75f863cd', '521b9b73c5b032f8e2809a5032d18266', '18ecc71ecc0df104a1628bdd4e529f85', '7aec18c648d573f34641b8329a2afe3e',
# 		'898f0448395da8aafa068c41aa2a02ae', '49e8ffdc3f384aa6c1a3b9ec6a76d4b9', '29db794970769658e1724d4cb9b38e3f', '38123381b3212f2a34c10b49951c2fd0', 'e71b6525fdf629e1e68cfb5d3ff941d8', 'af9105d9560eb7bbdd6ee6cfe1790bb3', '505582b521f1ffaee1bdcb7af70b67af', 'a51622b2d7542e94f2716874ed40a082', 'df65266474af01b30499b17346905968', '48e9ffb0a654b9d9615328eb92c3890d', '11fbe9d7dc93a6b26fa5660bc341313b', '0512fe8ec40a55e460138bffff54f010', '16a931947d00cd7c77349e0ab3dced0d', '478e48cc0ac669c61a4ee6546d25eeed', '6ca623696596a157aea1a522c7edc208', 'f3a0a6b8a6e8b772de8e782626ee317d', 'b49780be0bc03c2692883cfe184270b4', '69ef823b8cd30eaff590d3e42ae2b6a5', '6a904f7160253fb6d8033216ea864233', '1c6326dbae072c12227a9ecd3a003c96', 'a4c2f15dbc90546f3040d5c7beaa7526', 'a4c2c783a641fe4aefe9d55d4fc1de7e', '84e4a08aca08e75f36a478d3e151e727', '7cdd302b4bd837fd924784ebf91ca1c6', '9f6f76264015d2239c1f9080464b593c', '542243559a5b8c3c8c01339aa020498f', '23f4851a4b408ef507a81f31e916edb7', 'd1770c8af69a5f408dbf00931eed91f3', 'c04c71facfac812a4d811451a37be7fe', '1008438a9414a08c4a5c3268b11399b2', '6eedfb5561bb2d4a74ff9d08ab20ff13', '2e46cc8a9bf75e6ec69513e67673a65f', '567a56e7aa824c10589a1d92b545af83', '499ff6bf88889992abc415d35b3e6770', 'b646b6c748e9516e601332dc5bc54ae4', '4deb52e76dac733ed22c3404e270b0bb', '5b629102e1a1799339b21ac92186a176', 'afc8137d04bf1f7a6005a1eebf3773b6', 'a910b1fca3cc6417c729621bc16f1cf7', '6072005ad56ee5aa6151e5d1f83df8ca', 'f16d3043f8dd35450709f567607d8eb7', '9c40939ad680b87257507e1f689b327d', '27c2872532ad6b3e87490a7e84cf539e', '64e8fb3b32e84b488a26bb3d2664a7e9', 'd131d54a9cb6f5607a58144b2efce675', 'beb76a9060ea63899fdc102db4a4000d', '08f422fc394073e8b49eb17aaaafcc9e', '6fd62220291d6788ed2a69a0dbe7c480', 'c0ae19edebbe47eacae23e32d066393a', 'e020054d8b221312877fcf58f4a934cb', '33d7add4d749dcb4794cc67624c0a8dc', 'bcc2976ded8ed2ebd635e4bd2ae01a61', 'c7966e8aea8204f7b71b1ece22e3884a', '950b3360caf8d94bd02893429ecd2701', 'e32ff5f107e7ac9e26de395836d16b5e', 'f82c6018cd5ba244154242d35f6b7807', '97bb24336744339ae5711949d31edd25', '2a6b59bc275d1f799fbc3e6a82df415e', 'fd00ebebf66981f523e45f8bc62b3e08', '3f21b8fcd26af7023a3f2217c814783d', '75072d36160092526a69695ba89df5f5', '1e21bd7daf93b3dd4e35ae705f93c86d', '66e3ab59897912fa03d71e6c18d05830', 'eac92c2848c10feb1b54d94db228200f', '27686cd8736335a4ba5b38db83c99521', 'e7970013603df7808c7637069993442d', 'e2cf48045e1e4354dea02a091c6bbf15', 'afd38219dfc9fa010b181a6f1dc167df', '7df8990ca100ac47df3dc857a739a531', '158aa7e32c3bbcb984cf98cc36dc5d72', '6afcf63d1de62fc052e1aeb8e813e24f', '36ecd245a36bcf622cb83d6ab5d60133', 'fc7e4c35a8b52a75ee8b3ff4bdf7e9fb', 'e427efe71e8e94de5b3e48eb98236323', '1c70c44177d8f60eefb5338444e87549', '6d9190b20ab932a61569ec9ab8494f43', 'f9c08bb38fc10fb516e9466ebbe5f066', '0636d21d01e58665233956672543a0d9', '51fd154d6a251b2a4c0c880159e6a81e', '9300d739ce878a0d954cd0324df6efb3', '1c812bfc08cfaeb7ecfe6c991bc257dd', '2555fa54536c47a75563bb912fe33f55', '05cd25cf665a14365e0307f3b97229e4', '70c64089bd78379a419665ad4329b7f8', 'a064315f06f44c89de4381b38e67840c', 'cecf531d47cc5bb04e25c18ac27ad115', 'bda4adac022f5533f9d886b9c41024a0', 'bc0315c63076810bbd9d6e16d431dc2b', '074209e2b1aada170dc20d152ba5ba8f', '8fb5d8bae0dfb4b65b2f55c40f9c4f14', 'df87afc07497c146894ecc4aefeccccf', '04b329a774082ef0de750aa44312029a', '25fecbe66aa26fe70891b99adf1ffc0c',
# 		'3bb6de378fb690e55d4c67e29451bd62', '9aeb0b7003dde3e1baad11890f9c2dd4', '303db6a8617095217055dfdbaf079b69', '749331f5277515201850e1e6e22acdf4', 'aaf0033decdd384363799bc476f0e55e', '14e60be208219d814895a2625a1b8de7', '360a49eee0407cb33e4215eb5e8b4572', 'b454c0a557d5f3d6e2ffef899f6c0a21', 'c4673d6932c9d37d80a2e4d024d60e7f', '0549e6e9250e147a79ace4b98899b76a', '294f4c1a209eec3ab114d555ea611f49', 'f764c45e4cd13973274847c75dbd080f', '91cab5a0ca6300d36f7e8c27631f1b07', 'ac64baa8d0ad2302c32fcf1b040b6ad1', 'f2d70953e9665473f34ffc9acd242db8', '05a065416fb12f6bd4d6e7dd0b14fc8e', 'ebc02999fde0682fd8f85954ea23406d', '06026a8049be7fe60326d4bd26280626', 'c43b42dcb8743f070ded3e5bf79b1cd0', 'ba25c0fe3395629980ec815dceb2d463', '12d206a55265f90119c797ee00d719bc', '077309d95bb7df2cc11642a63161e8bc', '56ee9e4ab5ebfdf59071e7755236e42d', '2cd6ab8226b2bdd4310200f22d173aac', 'd8c3f606d92524a41e1ca3383fc39e75', '3d73ca7574ad13e1ae53afb25d82938d', '9c7efc497f98cd25a80a1cdd232b5b26', 'c7464bebfc11c6c47cb5142f2d505f23', '8bb1e49d39bc1b69869bee715b8529eb', '317b03f4caa6af8d147fd74d44c04adf', 'beda5a9dd9ab716431b0e2db194ffbe0', '7db83cd468bfcb5cf125996437412b96', '74e3dc029ccd9e5e1c82a64ece3c616e', '5e75192c96af82a0a715998bcb417cb1', '8c0ee4326521ac2daeb411e161cbf6dd', 'e37562a76ad784eca662bf5c7a3ca409', 'aea9854f73ec194517b3eb525b203cf4', 'e6b742a673feaf596c090b1620bfbd92', '41fa6301923bac1d7d51d912f24ed624', '7e4a75099b9bd368db54bd566c54c263', 'ad7fd5ae3295f634d7f37659f66faf7e', 'a0a56250b7b38ab3c30df2ded38d1c1f', 'fef7902e4fbd978114afbc7536b25e91', '11585c1edaa11913c306477bb3d9d552', 'c92f7f7450298c0660c1397f1d3620d1', '2a98c934e57f6d7438e2af42954cc998', '61d10f2a3c2ea845413b7ccb74974e2f', 'd28a9f52f7b6bdf78a08681fb03533a9', '3bf61ed6e2b1bddadbafde3d9f73e2e1',
# 		'6b8aef89e3919a3f9f0012120e3230ac'
# 			) ~ 60000,
# 			TRUE ~ memory
# 		)
# 	) |>
# 	mutate(
# 		memory = case_when(
# 			file_id_db2 %in% c(
# 			  'e2cf48045e1e4354dea02a091c6bbf15', '317b03f4caa6af8d147fd74d44c04adf', 'beda5a9dd9ab716431b0e2db194ffbe0', '7db83cd468bfcb5cf125996437412b96', '74e3dc029ccd9e5e1c82a64ece3c616e', '5e75192c96af82a0a715998bcb417cb1', 'e37562a76ad784eca662bf5c7a3ca409', 'aea9854f73ec194517b3eb525b203cf4', '11585c1edaa11913c306477bb3d9d552', '8c0ee4326521ac2daeb411e161cbf6dd', 'a0a56250b7b38ab3c30df2ded38d1c1f',
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
# 	write_lines(glue("~/PostDoc/CuratedAtlasQueryR/dev/DB2_files.makeflow"))



# Read arguments
args = commandArgs(trailingOnly = TRUE)
input_file = args[[1]]
gene_names = args[[2]]
output_file = args[[3]]

output_file |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)
file_id = basename(input_file) |> tools::file_path_sans_ext() |> str_split("___") %>% .[[1]] %>% .[1]
file_id_db2 = basename(output_file) |> tools::file_path_sans_ext()

cells_to_keep =
  
  # THIS HAS TO CHANGE IN SOMETHING STABLE
  duckdb() |>
  dbConnect(drv = _, read_only = TRUE) |>
  tbl(metadata_DB) |>
  
  
	dplyr::filter(file_id == !!file_id) |>
	dplyr::select(cell_, sample_, file_id, cell_type) |>
	as_tibble() |>
	unite("file_id_db2", c(file_id, cell_type), remove = FALSE) |>
	mutate(file_id_db2 = file_id_db2 |> md5() |> as.character()) |>
	filter(file_id_db2 == !!file_id_db2) |>
	mutate(.cell_original = cell_ |> str_remove(sample_) |> str_remove("_$")) |>
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

# Tranform
transformation =
  
  # THIS HAS TO CHANGE IN SOMETHING STABLE
  duckdb() |>
  dbConnect(drv = _, read_only = TRUE) |>
  tbl(metadata_DB) |>
  
  
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



data@assays@data$X =
  data@assays@data$X |>
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

# If I have negative values
if((data@assays@data$X[,1:min(10000, ncol(data@assays@data$X))] |> min()) < 0)
	data@assays@data$X[data@assays@data$X<0] <- 0

col_sums = colSums(data@assays@data$X)

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
	)
))
colnames(sce) = colnames(data)


rm(data)
rm(deduplicated_matrix)
gc()

# Complete gene set
missing_genes = readRDS(gene_names) |> unique()  |> setdiff(rownames(sce))

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



sce |>	saveHDF5SummarizedExperiment(output_file, replace=TRUE)




