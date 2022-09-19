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
hca_metadata = get_metadata()
```

Explore the HCA content

``` r
hca_metadata |> 
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
    hca_metadata |> 
    filter(
        ethnicity == "African" & 
        assay %LIKE% "%10x%" & 
        tissue == "lung parenchyma" & 
        cell_type %LIKE% "%CD4%"
    ) |> 
    get_SingleCellExperiment()
```

    ## Reading 9 files.

    ## .........

``` r
sce
```

    ## class: SingleCellExperiment 
    ## dim: 60661 1571 
    ## metadata(63): X_normalization batch_condition ... schema_version title
    ## assays(1): X
    ## rownames(60661): TSPAN6 TNMD ... RP11-175I6.6 PRSS43P
    ## rowData names(0):
    ## colnames(1571):
    ##   ACAGCCGGTCCGTTAA_F02526_4fc10a6b85e5fa688b253db4e0db8ba0
    ##   GGGAATGAGCCCAGCT_F02526_4fc10a6b85e5fa688b253db4e0db8ba0 ...
    ##   TCGTAGAGTATATGGA_SC29_948efa032d715a169f9edae65c696836
    ##   GACGTTACACTGTGTA_SC29_948efa032d715a169f9edae65c696836
    ## colData names(45): .sample .sample_name ... created_at.y updated_at.y
