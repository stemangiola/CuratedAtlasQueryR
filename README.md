
Load Metadata

``` r
library(HCAquery)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(dbplyr)
```

    ## 
    ## Attaching package: 'dbplyr'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     ident, sql

``` r
library(SingleCellExperiment)
```

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## Loading required package: GenomeInfoDb

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

``` r
library(tidySingleCellExperiment)
```

    ## Loading required package: ttservice

    ## 
    ## Attaching package: 'ttservice'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     bind_cols, bind_rows

    ## 
    ## Attaching package: 'tidySingleCellExperiment'

    ## The following objects are masked from 'package:ttservice':
    ## 
    ##     bind_cols, bind_rows

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     slice

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     rename

    ## The following object is masked from 'package:matrixStats':
    ## 
    ##     count

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     add_count, bind_cols, bind_rows, count

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
hca_metadata = get_metadata()
```

``` r
hca_metadata |> 
    distinct(.sample)
```

    ## # Source:   SQL [?? x 1]
    ## # Database: sqlite 3.39.3 [/vast/scratch/users/milton.m/metadata.sqlite]
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
    ## # Database: sqlite 3.39.3 [/vast/scratch/users/milton.m/metadata.sqlite]
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
      distinct(cell_type)
```

    ## # Source:   SQL [?? x 1]
    ## # Database: sqlite 3.39.3 [/vast/scratch/users/milton.m/metadata.sqlite]
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
    ## # Database:   sqlite 3.39.3 [/vast/scratch/users/milton.m/metadata.sqlite]
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

    ## # A SingleCellExperiment-tibble abstraction: 1,571 × 46
    ## # [90mFeatures=60661 | Cells=1571 | Assays=X[0m
    ##    .cell   .sample .samp…¹ assay assay…² cell_…³ cell_…⁴ devel…⁵ devel…⁶ disease
    ##    <chr>   <chr>   <chr>   <chr> <chr>   <chr>   <chr>   <chr>   <chr>   <chr>  
    ##  1 ACAGCC… 4fc10a… VUHD92… 10x … EFO:00… CD4-po… CL:000… 55-yea… HsapDv… normal 
    ##  2 GGGAAT… 4fc10a… VUHD92… 10x … EFO:00… CD4-po… CL:000… 55-yea… HsapDv… normal 
    ##  3 TCTTCG… 4fc10a… VUHD92… 10x … EFO:00… CD4-po… CL:000… 55-yea… HsapDv… normal 
    ##  4 CCTTAC… 4fc10a… VUHD92… 10x … EFO:00… CD4-po… CL:000… 55-yea… HsapDv… normal 
    ##  5 ATCTAC… 4fc10a… VUHD92… 10x … EFO:00… CD4-po… CL:000… 55-yea… HsapDv… normal 
    ##  6 CGCGTT… 4fc10a… VUHD92… 10x … EFO:00… CD4-po… CL:000… 55-yea… HsapDv… normal 
    ##  7 TTCCCA… 4fc10a… VUHD92… 10x … EFO:00… CD4-po… CL:000… 55-yea… HsapDv… normal 
    ##  8 TCGGTA… 4fc10a… VUHD92… 10x … EFO:00… CD4-po… CL:000… 55-yea… HsapDv… normal 
    ##  9 CCGGGA… 4fc10a… VUHD92… 10x … EFO:00… CD4-po… CL:000… 55-yea… HsapDv… normal 
    ## 10 TTGCCG… 4fc10a… VUHD92… 10x … EFO:00… CD4-po… CL:000… 55-yea… HsapDv… normal 
    ## # … with 1,561 more rows, 36 more variables: disease_ontology_term_id <chr>,
    ## #   ethnicity <chr>, ethnicity_ontology_term_id <chr>, file_id <chr>,
    ## #   is_primary_data.x <chr>, organism <chr>, organism_ontology_term_id <chr>,
    ## #   sample_placeholder <chr>, sex <chr>, sex_ontology_term_id <chr>,
    ## #   tissue <chr>, tissue_ontology_term_id <chr>, dataset_id <chr>,
    ## #   collection_id <chr>, cell_count <int>, dataset_deployments <chr>,
    ## #   is_primary_data.y <chr>, is_valid <int>, linked_genesets <int>, …
