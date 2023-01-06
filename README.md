readme
================

Load the package

``` r
library(HCAquery)
library(dplyr)
library(dbplyr)
library(SingleCellExperiment)
library(tidySingleCellExperiment)
options("restore_SingleCellExperiment_show" = TRUE)
```

Load the metadata

``` r
get_metadata()
#> # Source:   table<metadata> [?? x 56]
#> # Database: sqlite 3.40.0 [/vast/scratch/users/milton.m/cache/hca_harmonised/metadata.sqlite]
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
#> # Database:   sqlite 3.40.0 [/vast/scratch/users/milton.m/cache/hca_harmonised/metadata.sqlite]
#> # Ordered by: desc(n)
#>    tissue                    n
#>    <chr>                 <int>
#>  1 blood                    47
#>  2 heart left ventricle     46
#>  3 cortex of kidney         31
#>  4 renal medulla            29
#>  5 lung                     27
#>  6 middle temporal gyrus    24
#>  7 liver                    24
#>  8 kidney                   19
#>  9 intestine                18
#> 10 thymus                   17
#> # … with more rows
```

Query raw counts

``` r
sce <-
  get_metadata() |>
  filter(
    ethnicity == "African" &
      assay %LIKE% "%10x%" &
      tissue == "lung parenchyma" &
      cell_type %LIKE% "%CD4%"
  ) |>
  get_SingleCellExperiment()
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Attaching metadata.
#> ℹ Compiling Single Cell Experiment.

sce
#> class: SingleCellExperiment 
#> dim: 60661 1571 
#> metadata(0):
#> assays(2): counts cpm
#> rownames(60661): TSPAN6 TNMD ... RP11-175I6.6 PRSS43P
#> rowData names(0):
#> colnames(1571): ACAGCCGGTCCGTTAA_F02526 GGGAATGAGCCCAGCT_F02526 ... TACAACGTCAGCATTG_SC84 CATTCGCTCAATACCG_F02526
#> colData names(55): sample_id_db .sample ... n_cell_type_in_tissue n_tissue_in_cell_type
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
      assay %LIKE% "%10x%" &
      tissue == "lung parenchyma" &
      cell_type %LIKE% "%CD4%"
  ) |>
  get_SingleCellExperiment(assays = "cpm")
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Attaching metadata.
#> ℹ Compiling Single Cell Experiment.

sce
#> class: SingleCellExperiment 
#> dim: 60661 1571 
#> metadata(0):
#> assays(1): cpm
#> rownames(60661): TSPAN6 TNMD ... RP11-175I6.6 PRSS43P
#> rowData names(0):
#> colnames(1571): ACAGCCGGTCCGTTAA_F02526 GGGAATGAGCCCAGCT_F02526 ... TACAACGTCAGCATTG_SC84 CATTCGCTCAATACCG_F02526
#> colData names(55): sample_id_db .sample ... n_cell_type_in_tissue n_tissue_in_cell_type
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

Extract only a subset of genes:

``` r
get_metadata() |>
  filter(
    ethnicity == "African" &
      assay %LIKE% "%10x%" &
      tissue == "lung parenchyma" &
      cell_type %LIKE% "%CD4%"
  ) |>
  get_SingleCellExperiment(features = "PUM1")
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Attaching metadata.
#> ℹ Compiling Single Cell Experiment.
#> class: SingleCellExperiment 
#> dim: 1 1571 
#> metadata(0):
#> assays(2): counts cpm
#> rownames(1): PUM1
#> rowData names(0):
#> colnames(1571): ACAGCCGGTCCGTTAA_F02526 GGGAATGAGCCCAGCT_F02526 ... TACAACGTCAGCATTG_SC84 CATTCGCTCAATACCG_F02526
#> colData names(55): sample_id_db .sample ... n_cell_type_in_tissue n_tissue_in_cell_type
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

Extract the counts as a Seurat object:

``` r
get_metadata() |>
  filter(
    ethnicity == "African" &
      assay %LIKE% "%10x%" &
      tissue == "lung parenchyma" &
      cell_type %LIKE% "%CD4%"
  ) |>
  get_seurat()
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Attaching metadata.
#> ℹ Compiling Single Cell Experiment.
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#> An object of class Seurat 
#> 60661 features across 1571 samples within 1 assay 
#> Active assay: originalexp (60661 features, 0 variable features)
```
