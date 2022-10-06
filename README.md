HCAquery
================

Load Metadata

``` r
library(HCAquery)
library(dplyr)
library(dbplyr)
library(SingleCellExperiment)
library(tidySingleCellExperiment)
```

Load the data

``` r
get_metadata()
```

    ## # Source:   table<metadata> [?? x 48]
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
    ## # … with more rows, 38 more variables: disease_ontology_term_id <chr>,
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

Query

``` r
options("restore_SingleCellExperiment_show" = TRUE)
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
    ## dim: 28024 1571 
    ## metadata(0):
    ## assays(1): X
    ## rownames(28024): TSPAN6 TNMD ... RP11-107E5.4 RP11-299P2.2
    ## rowData names(0):
    ## colnames(1571): ACAGCCGGTCCGTTAA_F02526 GGGAATGAGCCCAGCT_F02526 ...
    ##   TACAACGTCAGCATTG_SC84 CATTCGCTCAATACCG_F02526
    ## colData names(47): .sample .sample_name ... file_id_db
    ##   stringr..str_remove.stringr..str_remove..cell...sample...._...
    ## reducedDimNames(0):
    ## mainExpName: NULL
    ## altExpNames(0):
