HCAquery
================

Load Metadata

``` r
library(HCAquery)
library(dplyr)
library(dbplyr)
library(SingleCellExperiment)
library(tidySingleCellExperiment)
options("restore_SingleCellExperiment_show" = TRUE)
```

Load the data

``` r
get_metadata()
```

    ## # Source:   table<metadata> [?? x 52]
    ## # Database: sqlite 3.39.2 [/vast/projects/RCP/human_cell_atlas/metadata.sqlite]
    ##    .cell   .sample .samp…¹ assay assay…² cell_…³ cell_…⁴ devel…⁵ devel…⁶ disease
    ##    <chr>   <chr>   <chr>   <chr> <chr>   <chr>   <chr>   <chr>   <chr>   <chr>  
    ##  1 AAACCT… 5f20d7… D17PrP… 10x … EFO:00… basal … CL:000… 31-yea… HsapDv… normal 
    ##  2 AAACCT… 5f20d7… D17PrP… 10x … EFO:00… basal … CL:000… 31-yea… HsapDv… normal 
    ##  3 AAACCT… 5f20d7… D17PrP… 10x … EFO:00… lumina… CL:000… 31-yea… HsapDv… normal 
    ##  4 AAACCT… 5f20d7… D17PrP… 10x … EFO:00… lumina… CL:000… 31-yea… HsapDv… normal 
    ##  5 AAACCT… 5f20d7… D17PrP… 10x … EFO:00… basal … CL:000… 31-yea… HsapDv… normal 
    ##  6 AAACCT… 5f20d7… D17PrP… 10x … EFO:00… basal … CL:000… 31-yea… HsapDv… normal 
    ##  7 AAACCT… 5f20d7… D17PrP… 10x … EFO:00… basal … CL:000… 31-yea… HsapDv… normal 
    ##  8 AAACGG… 5f20d7… D17PrP… 10x … EFO:00… basal … CL:000… 31-yea… HsapDv… normal 
    ##  9 AAACGG… 5f20d7… D17PrP… 10x … EFO:00… lumina… CL:000… 31-yea… HsapDv… normal 
    ## 10 AAACGG… 5f20d7… D17PrP… 10x … EFO:00… basal … CL:000… 31-yea… HsapDv… normal 
    ## # … with more rows, 42 more variables: disease_ontology_term_id <chr>,
    ## #   ethnicity <chr>, ethnicity_ontology_term_id <chr>, file_id <chr>,
    ## #   is_primary_data.x <chr>, organism <chr>, organism_ontology_term_id <chr>,
    ## #   sample_placeholder <chr>, sex <chr>, sex_ontology_term_id <chr>,
    ## #   tissue <chr>, tissue_ontology_term_id <chr>, dataset_id <chr>,
    ## #   collection_id <chr>, cell_count <int>, dataset_deployments <chr>,
    ## #   is_primary_data.y <chr>, is_valid <int>, linked_genesets <int>, …

Explore the HCA content

``` r
get_metadata() |> 
    distinct(tissue, file_id) |> 
    count(tissue) |> 
    arrange(desc(n))
```

    ## # Source:     SQL [?? x 2]
    ## # Database:   sqlite 3.39.2 [/vast/projects/RCP/human_cell_atlas/metadata.sqlite]
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
    ## colData names(51): .sample .sample_name ... cell_annotation_azimuth_l2
    ##   cell_annotation_blueprint_singler
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
    ## colData names(51): .sample .sample_name ... cell_annotation_azimuth_l2
    ##   cell_annotation_blueprint_singler
    ## reducedDimNames(0):
    ## mainExpName: NULL
    ## altExpNames(0):
