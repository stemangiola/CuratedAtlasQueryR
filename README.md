CuratedAtlasQueryR
================

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

`CuratedAtlasQuery` is a query interface that allow the programmatic
exploration and retrieval of the harmonised, curated and reannotated
CELLxGENE single-cell human cell atlas. Data can be retrieved at cell,
sample, or dataset levels based on filtering criteria.

<img src="man/figures/logo.png" width="120x" height="139px" />

<img src="man/figures/svcf_logo.jpeg" width="155x" height="58px" /><img src="man/figures/czi_logo.png" width="129px" height="58px" /><img src="man/figures/bioconductor_logo.jpg" width="202px" height="58px" /><img src="man/figures/vca_logo.png" width="219px" height="58px" />

[website](https://stemangiola.github.io/CuratedAtlasQueryR)

# Query interface

## Installation

``` r
devtools::install_github("stemangiola/CuratedAtlasQueryR")
```

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
#> # Source:   table</stornext/Home/data/allstaff/m/mangiola.s/.cache/R/CuratedAtlasQueryR/metadata.parquet> [?? x 56]
#> # Database: DuckDB 0.6.2-dev1166 [unknown@Linux 3.10.0-1160.81.1.el7.x86_64:R 4.2.0/:memory:]
#>    .cell   sampl…¹ .sample .samp…² assay assay…³ file_…⁴ cell_…⁵ cell_…⁶ devel…⁷
#>    <chr>   <chr>   <chr>   <chr>   <chr> <chr>   <chr>   <chr>   <chr>   <chr>  
#>  1 AAACCT… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea…
#>  2 AAACCT… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea…
#>  3 AAACCT… 02eb2e… 5f20d7… D17PrP… 10x … EFO:00… 30f754… lumina… CL:000… 31-yea…
#>  4 AAACCT… 02eb2e… 5f20d7… D17PrP… 10x … EFO:00… 30f754… lumina… CL:000… 31-yea…
#>  5 AAACCT… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea…
#>  6 AAACCT… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea…
#>  7 AAACCT… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea…
#>  8 AAACGG… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea…
#>  9 AAACGG… 02eb2e… 5f20d7… D17PrP… 10x … EFO:00… 30f754… lumina… CL:000… 31-yea…
#> 10 AAACGG… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea…
#> # … with more rows, 46 more variables:
#> #   development_stage_ontology_term_id <chr>, disease <chr>,
#> #   disease_ontology_term_id <chr>, ethnicity <chr>,
#> #   ethnicity_ontology_term_id <chr>, file_id <chr>, is_primary_data.x <chr>,
#> #   organism <chr>, organism_ontology_term_id <chr>, sample_placeholder <chr>,
#> #   sex <chr>, sex_ontology_term_id <chr>, tissue <chr>,
#> #   tissue_ontology_term_id <chr>, tissue_harmonised <chr>, age_days <dbl>, …
```

### Explore the tissue

``` r
get_metadata() |>
    dplyr::distinct(tissue, file_id) 
#> # Source:   SQL [?? x 2]
#> # Database: DuckDB 0.6.2-dev1166 [unknown@Linux 3.10.0-1160.81.1.el7.x86_64:R 4.2.0/:memory:]
#>    tissue                  file_id                             
#>    <chr>                   <chr>                               
#>  1 cortex of kidney        2977b3fa-e4d6-4929-8540-ae12d33a3c53
#>  2 entorhinal cortex       29d5d028-6f90-4943-91f7-fa3f93731de8
#>  3 middle temporal gyrus   2a689fda-d335-4ac0-81b1-a356fdf939db
#>  4 respiratory airway      2c2d5bea-8be7-4227-8a56-f2a85d57fa56
#>  5 thymus                  2ec94470-8171-4825-8346-34d77383438b
#>  6 ileum                   2f893c30-c91b-487b-b1f8-b386a4ef91a3
#>  7 cortex of kidney        2fe625b6-8c85-4818-a9c5-1189303f6e2b
#>  8 fimbria of uterine tube 3044b5dd-a499-456e-86d9-94769bc3b63e
#>  9 ampulla of uterine tube 3044b5dd-a499-456e-86d9-94769bc3b63e
#> 10 temporal cortex         4e4bbb2d-f341-4523-a5a0-5407d8b03e0e
#> # … with more rows
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
#> class: SingleCellExperiment 
#> dim: 60661 1571 
#> metadata(0):
#> assays(2): counts cpm
#> rownames(60661): TSPAN6 TNMD ... RP11-175I6.6 PRSS43P
#> rowData names(0):
#> colnames(1571): ACAGCCGGTCCGTTAA_F02526_1 GGGAATGAGCCCAGCT_F02526_1 ...
#>   TACAACGTCAGCATTG_SC84_1 CATTCGCTCAATACCG_F02526_1
#> colData names(56): sample_id_db .sample ... n_tissue_in_cell_type
#>   original_cell_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
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
#> class: SingleCellExperiment 
#> dim: 60661 1571 
#> metadata(0):
#> assays(1): cpm
#> rownames(60661): TSPAN6 TNMD ... RP11-175I6.6 PRSS43P
#> rowData names(0):
#> colnames(1571): ACAGCCGGTCCGTTAA_F02526_1 GGGAATGAGCCCAGCT_F02526_1 ...
#>   TACAACGTCAGCATTG_SC84_1 CATTCGCTCAATACCG_F02526_1
#> colData names(56): sample_id_db .sample ... n_tissue_in_cell_type
#>   original_cell_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
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
#> class: SingleCellExperiment 
#> dim: 1 1571 
#> metadata(0):
#> assays(1): cpm
#> rownames(1): PUM1
#> rowData names(0):
#> colnames(1571): ACAGCCGGTCCGTTAA_F02526_1 GGGAATGAGCCCAGCT_F02526_1 ...
#>   TACAACGTCAGCATTG_SC84_1 CATTCGCTCAATACCG_F02526_1
#> colData names(56): sample_id_db .sample ... n_tissue_in_cell_type
#>   original_cell_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

### Extract the counts as a Seurat object

This convert the H5 SingleCellExperiment to Seurat so it might take long
time and occupy a lot of memory depending on how many cells you are
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
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes
#> ('-')

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

<img src="man/figures/NCAM1_figure.png" width="629" />

# Cell metadata

Dataset-specific columns (definitions available at
cellxgene.cziscience.com)

`cell_count`, `collection_id`, `created_at.x`, `created_at.y`,
`dataset_deployments`, `dataset_id`, `file_id`, `filename`, `filetype`,
`is_primary_data.y`, `is_valid`, `linked_genesets`,
`mean_genes_per_cell`, `name`, `published`, `published_at`,
`revised_at`, `revision`, `s3_uri`, `schema_version`, `tombstone`,
`updated_at.x`, `updated_at.y`, `user_submitted`, `x_normalization`

Sample-specific columns (definitions available at
cellxgene.cziscience.com)

`.sample`, `.sample_name`, `age_days`, `assay`,
`assay_ontology_term_id`, `development_stage`,
`development_stage_ontology_term_id`, `ethnicity`,
`ethnicity_ontology_term_id`, `experiment___`, `organism`,
`organism_ontology_term_id`, `sample_placeholder`, `sex`,
`sex_ontology_term_id`, `tissue`, `tissue_harmonised`,
`tissue_ontology_term_id`, `disease`, `disease_ontology_term_id`,
`is_primary_data.x`

Cell-specific columns (definitions available at
cellxgene.cziscience.com)

`.cell`, `cell_type`, `cell_type_ontology_term_idm`,
`cell_type_harmonised`, `confidence_class`,
`cell_annotation_azimuth_l2`, `cell_annotation_blueprint_singler`

Through harmonisation and curation we introduced custom column, not
present in the original CELLxGENE metadata

- `tissue_harmonised`: a coarser tissue name for better filtering
- `age_days`: the number of days corresponding to the age
- `cell_type_harmonised`: the consensus call identity (for immune cells)
  using the original and three novel annotations using Seurat Azimuth
  and SingleR
- `confidence_class`: an ordinal class of how confident
  `cell_type_harmonised` is. 1 is complete consensus, 2 is 3 out of four
  and so on.  
- `cell_annotation_azimuth_l2`: Azimuth cell annotation
- `cell_annotation_blueprint_singler`: SingleR cell annotation using
  Blueprint reference
- `cell_annotation_blueprint_monaco`: SingleR cell annotation using
  Monaco reference
- `sample_id_db`: Sample subdivision for internal use
- `file_id_db`: File subdivision for internal use
- `.sample`: Sample ID
- `.sample_name`: How samples were defined

# RNA abundance

The `raw` assay includes RNA abundance in the positive real scale (not
transformed with non-linear functions, e.g. log sqrt). Originally
CELLxGENE include a mix of scales and transformations specified in the
`x_normalization` column.

The `cpm` assay includes counts per million.

------------------------------------------------------------------------

This project has been funded by

- *Silicon Valley Foundation* CZF2019-002443
- *Bioconductor core funding* NIH NHGRI 5U24HG004059-18
- *Victoria Cancer Agency* ECRF21036
- *Australian National Health and Medical Research Council* 1116955
