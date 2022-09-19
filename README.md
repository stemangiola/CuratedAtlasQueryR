HCAquery
================

Load Metadata

``` r
hca_metadata |> 
    distinct(.sample)
```

    ## # Source:   SQL [?? x 1]
    ## # Database: sqlite 3.39.2 [/vast/projects/RCP/human_cell_atlas/metadata.sqlite]
    ##    .sample                         
    ##    <chr>                           
    ##  1 5f20d7daf6c42f4f915fefd93d309120
    ##  2 2b449d364b13d4b045111793c78f53e0
    ##  3 a63ebabb27f83bc4e68db8bb9b2814ae
    ##  4 03c6a35e026a17b3a0e2c480bde49a81
    ##  5 c6a711eee7bc60452e57526dba9c814c
    ##  6 0fd1e7dc5d92dde4afe567186d059a8e
    ##  7 2d110c356cfded19edce4e90143eb40a
    ##  8 b560313b5d45f13770d98ae66e882184
    ##  9 b6e757366b7c70b206c3b8cd553ec355
    ## 10 be30712ad6fb8f58d6c918d51b792fc7
    ## # … with more rows

``` r
hca_metadata |> 
      distinct(file_id)
```

    ## # Source:   SQL [?? x 1]
    ## # Database: sqlite 3.39.2 [/vast/projects/RCP/human_cell_atlas/metadata.sqlite]
    ##    file_id                             
    ##    <chr>                               
    ##  1 00d626ec-c97e-4b2d-bf17-04bc09e52460
    ##  2 016fd62b-75f5-4778-a5f0-76bd378a5c0f
    ##  3 0273924c-0387-4f44-98c5-2292dbaab11e
    ##  4 0332cf8b-fec3-4e63-a6bc-4cca49a110ac
    ##  5 03746fce-b248-49d1-9d00-f078e413bb01
    ##  6 03939949-6033-4ee6-828c-4e2b762a413c
    ##  7 048287c8-b2f3-4183-a62b-0df2fb17c5d4
    ##  8 05ca5886-5db1-4ce7-b5cc-6c0df0a068c4
    ##  9 07beec85-51be-4d73-bb80-8f85b7b643d5
    ## 10 08247324-7ab7-45c5-8bd6-6c22676761ed
    ## # … with more rows

``` r
hca_metadata |> 
      nrow()
```

    ## [1] NA

``` r
hca_metadata |> 
      distinct(cell_type)
```

    ## # Source:   SQL [?? x 1]
    ## # Database: sqlite 3.39.2 [/vast/projects/RCP/human_cell_atlas/metadata.sqlite]
    ##    cell_type                          
    ##    <chr>                              
    ##  1 basal cell of prostate epithelium  
    ##  2 luminal cell of prostate epithelium
    ##  3 epithelial cell of urethra         
    ##  4 secretory cell                     
    ##  5 neuroendocrine cell                
    ##  6 classical monocyte                 
    ##  7 enteric smooth muscle cell         
    ##  8 fibroblast                         
    ##  9 neural cell                        
    ## 10 myofibroblast cell                 
    ## # … with more rows

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
    ## reducedDimNames(0):
    ## mainExpName: NULL
    ## altExpNames(0):
