CuratedAtlasQueryR
================

## Load the package

``` r
library(CuratedAtlasQueryR)
library(dplyr)
library(stringr)
```

## Load and explore the metadata

### Load the metadata

``` r
get_metadata()
#> # Source:   table<metadata> [?? x 56]
#> # Database: sqlite 3.40.0 [/stornext/Home/data/allstaff/m/mangiola.s/.cache/R/CuratedAtlasQueryR/metadata.sqlite]
#>    .cell  sampl…¹ .sample .samp…² assay assay…³ file_…⁴ cell_…⁵ cell_…⁶ devel…⁷ devel…⁸ disease disea…⁹ ethni…˟ ethni…˟ file_id is_pr…˟ organ…˟
#>    <chr>  <chr>   <chr>   <chr>   <chr> <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>  
#>  1 AAACC… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s…
#>  2 AAACC… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s…
#>  3 AAACC… 02eb2e… 5f20d7… D17PrP… 10x … EFO:00… 30f754… lumina… CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s…
#>  4 AAACC… 02eb2e… 5f20d7… D17PrP… 10x … EFO:00… 30f754… lumina… CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s…
#>  5 AAACC… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s…
#>  6 AAACC… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s…
#>  7 AAACC… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s…
#>  8 AAACG… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s…
#>  9 AAACG… 02eb2e… 5f20d7… D17PrP… 10x … EFO:00… 30f754… lumina… CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s…
#> 10 AAACG… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea… HsapDv… normal  PATO:0… Europe… HANCES… 00d626… FALSE   Homo s…
#> # … with more rows, 38 more variables: organism_ontology_term_id <chr>, sample_placeholder <chr>, sex <chr>, sex_ontology_term_id <chr>,
#> #   tissue <chr>, tissue_ontology_term_id <chr>, tissue_harmonised <chr>, age_days <dbl>, dataset_id <chr>, collection_id <chr>,
#> #   cell_count <int>, dataset_deployments <chr>, is_primary_data.y <chr>, is_valid <int>, linked_genesets <int>, mean_genes_per_cell <dbl>,
#> #   name <chr>, published <int>, revision <int>, schema_version <chr>, tombstone <int>, x_normalization <chr>, created_at.x <dbl>,
#> #   published_at <dbl>, revised_at <dbl>, updated_at.x <dbl>, filename <chr>, filetype <chr>, s3_uri <chr>, user_submitted <int>,
#> #   created_at.y <dbl>, updated_at.y <dbl>, cell_type_harmonised <chr>, confidence_class <dbl>, cell_annotation_azimuth_l2 <chr>,
#> #   cell_annotation_blueprint_singler <chr>, n_cell_type_in_tissue <int>, n_tissue_in_cell_type <int>, and abbreviated variable names …
```

### Explore the tissue

``` r
get_metadata() |>
    dplyr::distinct(tissue, file_id) 
```

``` r
#> # Source:     SQL [?? x 2]
#> # Database:   sqlite 3.40.0 [public_access@zki3lfhznsa.db.cloud.edu.au:5432/metadata]
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

## Download single-cell RNA sequencing counts

### Query raw counts

``` r

single_cell_counts = 
    get_metadata() |>
    dplyr::filter(
        ethnicity == "African" &
        stringr::str_like(assay, "%10x%") &
        tissue == "lung parenchyma" &
        stringr::str_like(cell_type, "%CD4%")
    ) |>
    get_SingleCellExperiment()
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Single Cell Experiment.

single_cell_counts
#> # A SingleCellExperiment-tibble abstraction: 1,571 × 57
#> # [90mFeatures=60661 | Cells=1571 | Assays=counts, cpm[0m
#>    .cell  sampl…¹ .sample .samp…² assay assay…³ file_…⁴ cell_…⁵ cell_…⁶ devel…⁷ devel…⁸ disease disea…⁹ ethni…˟ ethni…˟ file_id is_pr…˟ organ…˟
#>    <chr>  <chr>   <chr>   <chr>   <chr> <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>  
#>  1 ACAGC… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  2 GGGAA… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  3 TCTTC… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  4 CCTTA… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  5 ATCTA… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  6 CATCA… c7f14e… 9ef5ea… VUHD69… 10x … EFO:00… bc380d… CD4-po… CL:000… 31-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  7 AGTCT… c7f14e… 9ef5ea… VUHD69… 10x … EFO:00… bc380d… CD4-po… CL:000… 31-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  8 CGCGT… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  9 TGGCT… c7f14e… 9ef5ea… VUHD69… 10x … EFO:00… bc380d… CD4-po… CL:000… 31-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#> 10 TTCCC… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#> # … with 1,561 more rows, 39 more variables: organism_ontology_term_id <chr>, sample_placeholder <chr>, sex <chr>, sex_ontology_term_id <chr>,
#> #   tissue <chr>, tissue_ontology_term_id <chr>, tissue_harmonised <chr>, age_days <dbl>, dataset_id <chr>, collection_id <chr>,
#> #   cell_count <int>, dataset_deployments <chr>, is_primary_data.y <chr>, is_valid <int>, linked_genesets <int>, mean_genes_per_cell <dbl>,
#> #   name <chr>, published <int>, revision <int>, schema_version <chr>, tombstone <int>, x_normalization <chr>, created_at.x <dbl>,
#> #   published_at <dbl>, revised_at <dbl>, updated_at.x <dbl>, filename <chr>, filetype <chr>, s3_uri <chr>, user_submitted <int>,
#> #   created_at.y <dbl>, updated_at.y <dbl>, cell_type_harmonised <chr>, confidence_class <dbl>, cell_annotation_azimuth_l2 <chr>,
#> #   cell_annotation_blueprint_singler <chr>, n_cell_type_in_tissue <int>, n_tissue_in_cell_type <int>, original_cell_id <chr>, and …
```

### Query counts scaled per million

This is helpful if just few genes are of interest, as they can be
compared across samples.

``` r
single_cell_counts = 
    get_metadata() |>
    dplyr::filter(
        ethnicity == "African" &
        stringr::str_like(assay, "%10x%") &
        tissue == "lung parenchyma" &
        stringr::str_like(cell_type, "%CD4%")
    ) |>
    get_SingleCellExperiment(assays = "cpm")
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Single Cell Experiment.

single_cell_counts
#> # A SingleCellExperiment-tibble abstraction: 1,571 × 57
#> # [90mFeatures=60661 | Cells=1571 | Assays=cpm[0m
#>    .cell  sampl…¹ .sample .samp…² assay assay…³ file_…⁴ cell_…⁵ cell_…⁶ devel…⁷ devel…⁸ disease disea…⁹ ethni…˟ ethni…˟ file_id is_pr…˟ organ…˟
#>    <chr>  <chr>   <chr>   <chr>   <chr> <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>  
#>  1 ACAGC… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  2 GGGAA… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  3 TCTTC… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  4 CCTTA… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  5 ATCTA… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  6 CATCA… c7f14e… 9ef5ea… VUHD69… 10x … EFO:00… bc380d… CD4-po… CL:000… 31-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  7 AGTCT… c7f14e… 9ef5ea… VUHD69… 10x … EFO:00… bc380d… CD4-po… CL:000… 31-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  8 CGCGT… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  9 TGGCT… c7f14e… 9ef5ea… VUHD69… 10x … EFO:00… bc380d… CD4-po… CL:000… 31-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#> 10 TTCCC… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#> # … with 1,561 more rows, 39 more variables: organism_ontology_term_id <chr>, sample_placeholder <chr>, sex <chr>, sex_ontology_term_id <chr>,
#> #   tissue <chr>, tissue_ontology_term_id <chr>, tissue_harmonised <chr>, age_days <dbl>, dataset_id <chr>, collection_id <chr>,
#> #   cell_count <int>, dataset_deployments <chr>, is_primary_data.y <chr>, is_valid <int>, linked_genesets <int>, mean_genes_per_cell <dbl>,
#> #   name <chr>, published <int>, revision <int>, schema_version <chr>, tombstone <int>, x_normalization <chr>, created_at.x <dbl>,
#> #   published_at <dbl>, revised_at <dbl>, updated_at.x <dbl>, filename <chr>, filetype <chr>, s3_uri <chr>, user_submitted <int>,
#> #   created_at.y <dbl>, updated_at.y <dbl>, cell_type_harmonised <chr>, confidence_class <dbl>, cell_annotation_azimuth_l2 <chr>,
#> #   cell_annotation_blueprint_singler <chr>, n_cell_type_in_tissue <int>, n_tissue_in_cell_type <int>, original_cell_id <chr>, and …
```

### Extract only a subset of genes

``` r
single_cell_counts = 
    get_metadata() |>
    dplyr::filter(
        ethnicity == "African" &
        stringr::str_like(assay, "%10x%") &
        tissue == "lung parenchyma" &
        stringr::str_like(cell_type, "%CD4%")
    ) |>
    get_SingleCellExperiment(assays = "cpm", features = "PUM1")
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Single Cell Experiment.

single_cell_counts
#> # A SingleCellExperiment-tibble abstraction: 1,571 × 57
#> # [90mFeatures=1 | Cells=1571 | Assays=cpm[0m
#>    .cell  sampl…¹ .sample .samp…² assay assay…³ file_…⁴ cell_…⁵ cell_…⁶ devel…⁷ devel…⁸ disease disea…⁹ ethni…˟ ethni…˟ file_id is_pr…˟ organ…˟
#>    <chr>  <chr>   <chr>   <chr>   <chr> <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>  
#>  1 ACAGC… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  2 GGGAA… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  3 TCTTC… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  4 CCTTA… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  5 ATCTA… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  6 CATCA… c7f14e… 9ef5ea… VUHD69… 10x … EFO:00… bc380d… CD4-po… CL:000… 31-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  7 AGTCT… c7f14e… 9ef5ea… VUHD69… 10x … EFO:00… bc380d… CD4-po… CL:000… 31-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  8 CGCGT… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#>  9 TGGCT… c7f14e… 9ef5ea… VUHD69… 10x … EFO:00… bc380d… CD4-po… CL:000… 31-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#> 10 TTCCC… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s…
#> # … with 1,561 more rows, 39 more variables: organism_ontology_term_id <chr>, sample_placeholder <chr>, sex <chr>, sex_ontology_term_id <chr>,
#> #   tissue <chr>, tissue_ontology_term_id <chr>, tissue_harmonised <chr>, age_days <dbl>, dataset_id <chr>, collection_id <chr>,
#> #   cell_count <int>, dataset_deployments <chr>, is_primary_data.y <chr>, is_valid <int>, linked_genesets <int>, mean_genes_per_cell <dbl>,
#> #   name <chr>, published <int>, revision <int>, schema_version <chr>, tombstone <int>, x_normalization <chr>, created_at.x <dbl>,
#> #   published_at <dbl>, revised_at <dbl>, updated_at.x <dbl>, filename <chr>, filetype <chr>, s3_uri <chr>, user_submitted <int>,
#> #   created_at.y <dbl>, updated_at.y <dbl>, cell_type_harmonised <chr>, confidence_class <dbl>, cell_annotation_azimuth_l2 <chr>,
#> #   cell_annotation_blueprint_singler <chr>, n_cell_type_in_tissue <int>, n_tissue_in_cell_type <int>, original_cell_id <chr>, and …
```

### Extract the counts as a Seurat object

This convert the H5 SingleCellExperiment to Seurat so it might take long
time and occupy a lot of memory dependeing on how many cells you are
requesting.

``` r
single_cell_counts = 
    get_metadata() |>
    dplyr::filter(
        ethnicity == "African" &
        stringr::str_like(assay, "%10x%") &
        tissue == "lung parenchyma" &
        stringr::str_like(cell_type, "%CD4%")
    ) |>
    get_seurat()
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Single Cell Experiment.
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

single_cell_counts
#> An object of class Seurat 
#> 60661 features across 1571 samples within 1 assay 
#> Active assay: originalexp (60661 features, 0 variable features)
```

## Visualise gene transcription

We can gather all natural killer cells and plot the distribution of CD56
(NCAM1) across all tissues

``` r
library(tidySingleCellExperiment)
library(ggplot2)

get_metadata() |> 
    
  # Filter and subset
  filter(cell_type_harmonised=="nk") |> 
  select(.cell, file_id_db, disease, file_id, tissue_harmonised) |> 
  
  # Get counts per million for NCAM1 gene 
  get_SingleCellExperiment(assays = "cpm", features = "NCAM1") |> 

    # Get transcriptional abundance for plotting with `tidySingleCellExperiment`
  join_features("NCAM1", shape = "wide") |> 
    
    # Plot
  ggplot(aes( tissue_harmonised, NCAM1,color = file_id)) +
  geom_jitter(shape=".") +
    
    # Style
  guides(color="none") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
```

<img src="../inst/NCAM1_figure.png" width="629" />

``` r
sessionInfo()
#> R version 4.2.0 (2022-04-22)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: CentOS Linux 7 (Core)
#> 
#> Matrix products: default
#> BLAS:   /stornext/System/data/apps/R/R-4.2.0/lib64/R/lib/libRblas.so
#> LAPACK: /stornext/System/data/apps/R/R-4.2.0/lib64/R/lib/libRlapack.so
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
#>  [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] here_1.0.1                     stringr_1.5.0                  dplyr_1.1.0                    CuratedAtlasQueryR_0.1.0      
#>  [5] dbplyr_2.3.0                   ggplot2_3.4.0                  tidySingleCellExperiment_1.9.1 ttservice_0.2.2               
#>  [9] SingleCellExperiment_1.20.0    SummarizedExperiment_1.28.0    Biobase_2.58.0                 GenomicRanges_1.50.2          
#> [13] GenomeInfoDb_1.34.7            HDF5Array_1.26.0               rhdf5_2.42.0                   DelayedArray_0.24.0           
#> [17] IRanges_2.32.0                 S4Vectors_0.36.1               MatrixGenerics_1.10.0          matrixStats_0.63.0            
#> [21] BiocGenerics_0.44.0            Matrix_1.5-3                  
#> 
#> loaded via a namespace (and not attached):
#>   [1] plyr_1.8.8             igraph_1.3.5           lazyeval_0.2.2         sp_1.5-1               splines_4.2.0          listenv_0.9.0         
#>   [7] scattermore_0.8        inline_0.3.19          digest_0.6.31          htmltools_0.5.4        rsconnect_0.8.29       fansi_1.0.4           
#>  [13] magrittr_2.0.3         memoise_2.0.1          tensor_1.5             cluster_2.1.4          ROCR_1.0-11            globals_0.16.2        
#>  [19] RcppParallel_5.1.6     spatstat.sparse_3.0-0  prettyunits_1.1.1      colorspace_2.1-0       blob_1.2.3             ggrepel_0.9.2         
#>  [25] xfun_0.36              callr_3.7.3            crayon_1.5.2           RCurl_1.98-1.9         jsonlite_1.8.4         spatstat.data_3.0-0   
#>  [31] progressr_0.13.0       survival_3.5-0         zoo_1.8-11             glue_1.6.2             polyclip_1.10-4        gtable_0.3.1          
#>  [37] zlibbioc_1.44.0        XVector_0.38.0         leiden_0.4.3           V8_4.2.2               pkgbuild_1.4.0         Rhdf5lib_1.20.0       
#>  [43] rstan_2.26.6           future.apply_1.10.0    abind_1.4-5            scales_1.2.1           DBI_1.1.3              spatstat.random_3.0-1 
#>  [49] miniUI_0.1.1.1         Rcpp_1.0.10            viridisLite_0.4.1      xtable_1.8-4           reticulate_1.27        bit_4.0.5             
#>  [55] StanHeaders_2.26.6     htmlwidgets_1.6.1      httr_1.4.4             RColorBrewer_1.1-3     ellipsis_0.3.2         Seurat_4.3.0          
#>  [61] ica_1.0-3              pkgconfig_2.0.3        loo_2.5.1              farver_2.1.1           uwot_0.1.14            deldir_1.0-6          
#>  [67] utf8_1.2.2             tidyselect_1.2.0       labeling_0.4.2         rlang_1.0.6            reshape2_1.4.4         later_1.3.0           
#>  [73] munsell_0.5.0          tools_4.2.0            cachem_1.0.6           cli_3.6.0              generics_0.1.3         RSQLite_2.2.20        
#>  [79] ggridges_0.5.4         evaluate_0.20          fastmap_1.1.0          goftest_1.2-3          yaml_2.3.7             processx_3.8.0        
#>  [85] knitr_1.42             bit64_4.0.5            fitdistrplus_1.1-8     purrr_1.0.1            RANN_2.6.1             nlme_3.1-161          
#>  [91] pbapply_1.7-0          future_1.30.0          mime_0.12              compiler_4.2.0         rstudioapi_0.14        plotly_4.10.1         
#>  [97] curl_5.0.0             png_0.1-8              spatstat.utils_3.0-1   tibble_3.1.8           stringi_1.7.12         highr_0.10            
#> [103] ps_1.7.2               lattice_0.20-45        vctrs_0.5.2            pillar_1.8.1           lifecycle_1.0.3        rhdf5filters_1.10.0   
#> [109] spatstat.geom_3.0-3    lmtest_0.9-40          RcppAnnoy_0.0.20       data.table_1.14.6      cowplot_1.1.1          bitops_1.0-7          
#> [115] irlba_2.3.5.1          httpuv_1.6.8           patchwork_1.1.2        R6_2.5.1               promises_1.2.0.1       KernSmooth_2.23-20    
#> [121] gridExtra_2.3          parallelly_1.34.0      codetools_0.2-18       MASS_7.3-58.1          assertthat_0.2.1       rprojroot_2.0.3       
#> [127] withr_2.5.0            SeuratObject_4.1.3     sctransform_0.3.5      GenomeInfoDbData_1.2.9 parallel_4.2.0         grid_4.2.0            
#> [133] tidyr_1.3.0            rmarkdown_2.20         Rtsne_0.16             spatstat.explore_3.0-5 shiny_1.7.4
```
