HCA Harmonised
================

Load the package

``` r
library(CuratedAtlasQuery)
library(dplyr)
```

Load the metadata

``` r
get_metadata()
#> # Source:   table<metadata> [?? x 56]
#> # Database: postgres  [public_access@zki3lfhznsa.db.cloud.edu.au:5432/metadata]
#>    .cell   sampl…¹ .sample .samp…² assay assay…³ file_…⁴ cell_…⁵ cell_…⁶ devel…⁷ devel…⁸ disease disea…⁹ ethni…˟ ethni…˟ file_id is_pr…˟ organ…˟ organ…˟ sampl…˟ sex   sex_o…˟ tissue
#>    <chr>   <chr>   <chr>   <chr>   <chr> <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr> <chr>   <chr> 
#>  1 AAACCT… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s… NCBITa… <NA>    male  PATO:0… perip…
#>  2 AAACCT… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s… NCBITa… <NA>    male  PATO:0… perip…
#>  3 AAACCT… 02eb2e… 5f20d7… D17PrP… 10x … EFO:00… 30f754… lumina… CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s… NCBITa… <NA>    male  PATO:0… perip…
#>  4 AAACCT… 02eb2e… 5f20d7… D17PrP… 10x … EFO:00… 30f754… lumina… CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s… NCBITa… <NA>    male  PATO:0… perip…
#>  5 AAACCT… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s… NCBITa… <NA>    male  PATO:0… perip…
#>  6 AAACCT… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s… NCBITa… <NA>    male  PATO:0… perip…
#>  7 AAACCT… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s… NCBITa… <NA>    male  PATO:0… perip…
#>  8 AAACGG… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s… NCBITa… <NA>    male  PATO:0… perip…
#>  9 AAACGG… 02eb2e… 5f20d7… D17PrP… 10x … EFO:00… 30f754… lumina… CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s… NCBITa… <NA>    male  PATO:0… perip…
#> 10 AAACGG… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s… NCBITa… <NA>    male  PATO:0… perip…
#> # … with more rows, 33 more variables: tissue_ontology_term_id <chr>, tissue_harmonised <chr>, age_days <dbl>, dataset_id <chr>, collection_id <chr>, cell_count <int>,
#> #   dataset_deployments <chr>, is_primary_data.y <chr>, is_valid <int>, linked_genesets <int>, mean_genes_per_cell <dbl>, name <chr>, published <int>, revision <int>,
#> #   schema_version <chr>, tombstone <int>, x_normalization <chr>, created_at.x <dbl>, published_at <dbl>, revised_at <dbl>, updated_at.x <dbl>, filename <chr>, filetype <chr>,
#> #   s3_uri <chr>, user_submitted <int>, created_at.y <dbl>, updated_at.y <dbl>, cell_type_harmonised <chr>, confidence_class <dbl>, cell_annotation_azimuth_l2 <chr>,
#> #   cell_annotation_blueprint_singler <chr>, n_cell_type_in_tissue <int>, n_tissue_in_cell_type <int>, and abbreviated variable names ¹​sample_id_db, ²​.sample_name,
#> #   ³​assay_ontology_term_id, ⁴​file_id_db, ⁵​cell_type, ⁶​cell_type_ontology_term_id, ⁷​development_stage, ⁸​development_stage_ontology_term_id, ⁹​disease_ontology_term_id, ˟​ethnicity,
#> #   ˟​ethnicity_ontology_term_id, ˟​is_primary_data.x, ˟​organism, ˟​organism_ontology_term_id, ˟​sample_placeholder, ˟​sex_ontology_term_id
```

Explore the HCA content

``` r
get_metadata() |>
    distinct(tissue, file_id) |>
    count(tissue) |>
    arrange(desc(n))
#> # Source:     SQL [?? x 2]
#> # Database:   postgres  [public_access@zki3lfhznsa.db.cloud.edu.au:5432/metadata]
#> # Ordered by: desc(n)
#>    tissue                      n
#>    <chr>                 <int64>
#>  1 blood                      47
#>  2 heart left ventricle       46
#>  3 cortex of kidney           31
#>  4 renal medulla              29
#>  5 lung                       27
#>  6 liver                      24
#>  7 middle temporal gyrus      24
#>  8 kidney                     19
#>  9 intestine                  18
#> 10 thymus                     17
#> # … with more rows
```

Query raw counts

``` r
library(stringr)
sce <-
    get_metadata() |>
    filter(
        ethnicity == "African" &
        str_like(assay, "%10x%") &
        tissue == "lung parenchyma" &
        str_like(cell_type, "%CD4%")
    ) |>
    get_SingleCellExperiment()
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Single Cell Experiment.

sce
#> class: SingleCellExperiment 
#> dim: 60661 1571 
#> metadata(0):
#> assays(2): counts cpm
#> rownames(60661): TSPAN6 TNMD ... RP11-175I6.6 PRSS43P
#> rowData names(0):
#> colnames(1571): ACAGCCGGTCCGTTAA_F02526_1 GGGAATGAGCCCAGCT_F02526_1 ... TACAACGTCAGCATTG_SC84_1 CATTCGCTCAATACCG_F02526_1
#> colData names(56): sample_id_db .sample ... n_tissue_in_cell_type original_cell_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

Query counts scaled per million. This is helpful if just few genes are
of interest

``` r
sce <-
    get_metadata() |>
    filter(
        ethnicity == "African" &
        str_like(assay, "%10x%") &
        tissue == "lung parenchyma" &
        str_like(cell_type, "%CD4%")
    ) |>
    get_SingleCellExperiment(assays = "cpm")
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Single Cell Experiment.

sce
#> class: SingleCellExperiment 
#> dim: 60661 1571 
#> metadata(0):
#> assays(1): cpm
#> rownames(60661): TSPAN6 TNMD ... RP11-175I6.6 PRSS43P
#> rowData names(0):
#> colnames(1571): ACAGCCGGTCCGTTAA_F02526_1 GGGAATGAGCCCAGCT_F02526_1 ... TACAACGTCAGCATTG_SC84_1 CATTCGCTCAATACCG_F02526_1
#> colData names(56): sample_id_db .sample ... n_tissue_in_cell_type original_cell_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

Extract only a subset of genes:

``` r
get_metadata() |>
    filter(
        ethnicity == "African" &
        str_like(assay, "%10x%") &
        tissue == "lung parenchyma" &
        str_like(cell_type, "%CD4%")
    ) |>
    get_SingleCellExperiment(features = "PUM1")
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Single Cell Experiment.
#> class: SingleCellExperiment 
#> dim: 1 1571 
#> metadata(0):
#> assays(2): counts cpm
#> rownames(1): PUM1
#> rowData names(0):
#> colnames(1571): ACAGCCGGTCCGTTAA_F02526_1 GGGAATGAGCCCAGCT_F02526_1 ... TACAACGTCAGCATTG_SC84_1 CATTCGCTCAATACCG_F02526_1
#> colData names(56): sample_id_db .sample ... n_tissue_in_cell_type original_cell_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

Extract the counts as a Seurat object:

``` r
get_metadata() |>
    filter(
        ethnicity == "African" &
        str_like(assay, "%10x%") &
        tissue == "lung parenchyma" &
        str_like(cell_type, "%CD4%")
    ) |>
    get_seurat()
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Single Cell Experiment.
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#> An object of class Seurat 
#> 60661 features across 1571 samples within 1 assay 
#> Active assay: originalexp (60661 features, 0 variable features)
```

``` r
sessionInfo()
#> R version 4.2.1 (2022-06-23)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: CentOS Linux 7 (Core)
#> 
#> Matrix products: default
#> BLAS:   /stornext/System/data/apps/R/R-4.2.1/lib64/R/lib/libRblas.so
#> LAPACK: /stornext/System/data/apps/R/R-4.2.1/lib64/R/lib/libRlapack.so
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] stringr_1.5.0   dplyr_1.0.10    CuratedAtlasQuery_0.1.0  testthat_3.1.6  RPostgres_1.4.4
#> 
#> loaded via a namespace (and not attached):
#>   [1] plyr_1.8.8                  igraph_1.3.5                lazyeval_0.2.2              sp_1.5-1                    splines_4.2.1               listenv_0.9.0              
#>   [7] scattermore_0.8             usethis_2.1.6               GenomeInfoDb_1.32.4         ggplot2_3.4.0               digest_0.6.31               htmltools_0.5.4            
#>  [13] fansi_1.0.3                 magrittr_2.0.3              memoise_2.0.1               tensor_1.5                  cluster_2.1.3               ROCR_1.0-11                
#>  [19] remotes_2.4.2               globals_0.16.2              matrixStats_0.63.0          timechange_0.2.0            spatstat.sparse_3.0-0       prettyunits_1.1.1          
#>  [25] colorspace_2.0-3            blob_1.2.3                  rappdirs_0.3.3              ggrepel_0.9.2               xfun_0.36                   callr_3.7.3                
#>  [31] crayon_1.5.2                RCurl_1.98-1.9              jsonlite_1.8.4              roxygen2_7.2.3              progressr_0.12.0            spatstat.data_3.0-0        
#>  [37] survival_3.3-1              zoo_1.8-11                  glue_1.6.2                  polyclip_1.10-4             gtable_0.3.1                zlibbioc_1.42.0            
#>  [43] XVector_0.36.0              leiden_0.4.3                DelayedArray_0.22.0         pkgbuild_1.4.0              Rhdf5lib_1.18.2             future.apply_1.10.0        
#>  [49] SingleCellExperiment_1.18.1 BiocGenerics_0.42.0         HDF5Array_1.24.2            abind_1.4-5                 scales_1.2.1                DBI_1.1.3                  
#>  [55] spatstat.random_3.0-1       miniUI_0.1.1.1              Rcpp_1.0.9                  viridisLite_0.4.1           xtable_1.8-4                reticulate_1.26            
#>  [61] bit_4.0.5                   stats4_4.2.1                profvis_0.3.7               htmlwidgets_1.6.0           httr_1.4.4                  RColorBrewer_1.1-3         
#>  [67] ellipsis_0.3.2              Seurat_4.3.0                ica_1.0-3                   urlchecker_1.0.1            pkgconfig_2.0.3             dbplyr_2.3.0               
#>  [73] uwot_0.1.14                 deldir_1.0-6                utf8_1.2.2                  tidyselect_1.2.0            rlang_1.0.6                 reshape2_1.4.4             
#>  [79] later_1.3.0                 munsell_0.5.0               tools_4.2.1                 cachem_1.0.6                cli_3.5.0                   generics_0.1.3             
#>  [85] RSQLite_2.2.20              devtools_2.4.5              ggridges_0.5.4              evaluate_0.19               fastmap_1.1.0               yaml_2.3.6                 
#>  [91] goftest_1.2-3               knitr_1.41                  processx_3.8.0              bit64_4.0.5                 fs_1.5.2                    fitdistrplus_1.1-8         
#>  [97] purrr_1.0.1                 RANN_2.6.1                  nlme_3.1-157                pbapply_1.6-0               future_1.30.0               mime_0.12                  
#> [103] xml2_1.3.3                  brio_1.1.3                  compiler_4.2.1              rstudioapi_0.14             plotly_4.10.1               png_0.1-8                  
#> [109] spatstat.utils_3.0-1        tibble_3.1.8                stringi_1.7.8               ps_1.7.2                    desc_1.4.2                  lattice_0.20-45            
#> [115] Matrix_1.5-3                vctrs_0.5.1                 pillar_1.8.1                lifecycle_1.0.3             rhdf5filters_1.8.0          spatstat.geom_3.0-3        
#> [121] lmtest_0.9-40               RcppAnnoy_0.0.20            data.table_1.14.6           cowplot_1.1.1               bitops_1.0-7                irlba_2.3.5.1              
#> [127] httpuv_1.6.7                patchwork_1.1.2             GenomicRanges_1.48.0        R6_2.5.1                    promises_1.2.0.1            KernSmooth_2.23-20         
#> [133] gridExtra_2.3               IRanges_2.30.1              parallelly_1.33.0           sessioninfo_1.2.2           codetools_0.2-18            MASS_7.3-57                
#> [139] assertthat_0.2.1            pkgload_1.3.2               rhdf5_2.40.0                SummarizedExperiment_1.26.1 rprojroot_2.0.3             withr_2.5.0                
#> [145] SeuratObject_4.1.3          sctransform_0.3.5           S4Vectors_0.34.0            GenomeInfoDbData_1.2.8      parallel_4.2.1              hms_1.1.2                  
#> [151] grid_4.2.1                  tidyr_1.2.1                 rmarkdown_2.19              MatrixGenerics_1.8.1        Rtsne_0.16                  spatstat.explore_3.0-5     
#> [157] lubridate_1.9.0             Biobase_2.56.0              shiny_1.7.4
```
