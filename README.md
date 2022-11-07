HCAquery
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
```

    ## # Source:   table<metadata> [?? x 56]
    ## # Database: sqlite 3.39.3 [/vast/projects/RCP/human_cell_atlas/metadata.sqlite]
    ##    .cell   sampl…¹ .sample .samp…² assay assay…³ file_…⁴ cell_…⁵ cell_…⁶ devel…⁷
    ##    <chr>   <chr>   <chr>   <chr>   <chr> <chr>   <chr>   <chr>   <chr>   <chr>  
    ##  1 AAACCT… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea…
    ##  2 AAACCT… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea…
    ##  3 AAACCT… 02eb2e… 5f20d7… D17PrP… 10x … EFO:00… 30f754… lumina… CL:000… 31-yea…
    ##  4 AAACCT… 02eb2e… 5f20d7… D17PrP… 10x … EFO:00… 30f754… lumina… CL:000… 31-yea…
    ##  5 AAACCT… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea…
    ##  6 AAACCT… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea…
    ##  7 AAACCT… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea…
    ##  8 AAACGG… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea…
    ##  9 AAACGG… 02eb2e… 5f20d7… D17PrP… 10x … EFO:00… 30f754… lumina… CL:000… 31-yea…
    ## 10 AAACGG… 8a0fe0… 5f20d7… D17PrP… 10x … EFO:00… 1e334b… basal … CL:000… 31-yea…
    ## # … with more rows, 46 more variables:
    ## #   development_stage_ontology_term_id <chr>, disease <chr>,
    ## #   disease_ontology_term_id <chr>, ethnicity <chr>,
    ## #   ethnicity_ontology_term_id <chr>, file_id <chr>, is_primary_data.x <chr>,
    ## #   organism <chr>, organism_ontology_term_id <chr>, sample_placeholder <chr>,
    ## #   sex <chr>, sex_ontology_term_id <chr>, tissue <chr>,
    ## #   tissue_ontology_term_id <chr>, tissue_harmonised <chr>, age_days <dbl>, …

Explore the HCA content

``` r
get_metadata() |> 
    distinct(tissue, file_id) |> 
    count(tissue) |> 
    arrange(desc(n))
```

    ## # Source:     SQL [?? x 2]
    ## # Database:   sqlite 3.39.3 [/vast/projects/RCP/human_cell_atlas/metadata.sqlite]
    ## # Ordered by: desc(n)
    ##    tissue                    n
    ##    <chr>                 <int>
    ##  1 blood                    47
    ##  2 heart left ventricle     46
    ##  3 cortex of kidney         31
    ##  4 renal medulla            29
    ##  5 lung                     27
    ##  6 middle temporal gyrus    24
    ##  7 liver                    24
    ##  8 kidney                   19
    ##  9 intestine                18
    ## 10 thymus                   17
    ## # … with more rows

Query raw counts

``` r
sce = 
    get_metadata() |> 
    filter(
         ethnicity == "African" & 
        assay %LIKE% "%10x%" & 
        tissue == "lung parenchyma" & 
        cell_type %LIKE% "%CD4%"
    ) |> 
    
    get_SingleCellExperiment()
```

    ## Reading 1 files.

    ## .

``` r
sce
```

    ## class: SingleCellExperiment 
    ## dim: 60661 1571 
    ## metadata(0):
    ## assays(1): counts
    ## rownames(60661): TSPAN6 TNMD ... RP11-175I6.6 PRSS43P
    ## rowData names(0):
    ## colnames(1571): ACAGCCGGTCCGTTAA_F02526 GGGAATGAGCCCAGCT_F02526 ...
    ##   TACAACGTCAGCATTG_SC84 CATTCGCTCAATACCG_F02526
    ## colData names(55): sample_id_db .sample ... n_cell_type_in_tissue
    ##   n_tissue_in_cell_type
    ## reducedDimNames(0):
    ## mainExpName: NULL
    ## altExpNames(0):

Query counts scaled per million. This is helpful if just few genes are
of interest

``` r
sce = 
    get_metadata() |> 
    filter(
         ethnicity == "African" & 
        assay %LIKE% "%10x%" & 
        tissue == "lung parenchyma" & 
        cell_type %LIKE% "%CD4%"
    ) |> 
    
    get_SingleCellExperiment(assay = "counts_per_million")
```

    ## Reading 1 files.

    ## .

``` r
sce
```

    ## class: SingleCellExperiment 
    ## dim: 60661 1571 
    ## metadata(0):
    ## assays(1): counts_per_million
    ## rownames(60661): TSPAN6 TNMD ... RP11-175I6.6 PRSS43P
    ## rowData names(0):
    ## colnames(1571): ACAGCCGGTCCGTTAA_F02526 GGGAATGAGCCCAGCT_F02526 ...
    ##   TACAACGTCAGCATTG_SC84 CATTCGCTCAATACCG_F02526
    ## colData names(55): sample_id_db .sample ... n_cell_type_in_tissue
    ##   n_tissue_in_cell_type
    ## reducedDimNames(0):
    ## mainExpName: NULL
    ## altExpNames(0):

Extract only a subset of genes:

``` r
get_metadata() |> 
  filter(
     ethnicity == "African" & 
      assay %LIKE% "%10x%" & 
      tissue == "lung parenchyma" & 
      cell_type %LIKE% "%CD4%"
  ) |> 
get_SingleCellExperiment(genes = "PUM1")
```

    ## Reading 1 files.

    ## .

    ## class: SingleCellExperiment 
    ## dim: 1 1571 
    ## metadata(0):
    ## assays(1): counts
    ## rownames(1): PUM1
    ## rowData names(0):
    ## colnames(1571): ACAGCCGGTCCGTTAA_F02526 GGGAATGAGCCCAGCT_F02526 ...
    ##   TACAACGTCAGCATTG_SC84 CATTCGCTCAATACCG_F02526
    ## colData names(55): sample_id_db .sample ... n_cell_type_in_tissue
    ##   n_tissue_in_cell_type
    ## reducedDimNames(0):
    ## mainExpName: NULL
    ## altExpNames(0):

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
```

    ## Reading 1 files.

    ## .

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## An object of class Seurat 
    ## 60661 features across 1571 samples within 1 assay 
    ## Active assay: originalexp (60661 features, 0 variable features)
