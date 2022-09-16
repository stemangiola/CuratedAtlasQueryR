HCAquery
================

Load Metadata

``` r
library(HCAquery)
library(tidyverse)
```

    ## ── Attaching packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ ggplot2 3.3.6.9000     ✔ purrr   0.3.4     
    ## ✔ tibble  3.1.8          ✔ dplyr   1.0.9     
    ## ✔ tidyr   1.2.0          ✔ stringr 1.4.1     
    ## ✔ readr   2.1.2          ✔ forcats 0.5.2

``` r
library(SingleCellExperiment)
```

    ## Loading required package: SummarizedExperiment
    ## Loading required package: methods
    ## Loading required package: MatrixGenerics
    ## Loading required package: matrixStats
    ## 
    ## Attaching package: 'matrixStats'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count
    ## 
    ## 
    ## Attaching package: 'MatrixGenerics'
    ## 
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
    ## 
    ## Loading required package: GenomicRanges
    ## Loading required package: stats4
    ## Loading required package: BiocGenerics
    ## Loading required package: utils
    ## Loading required package: graphics
    ## 
    ## Attaching package: 'graphics'
    ## 
    ## The following object is masked from 'package:stats4':
    ## 
    ##     plot
    ## 
    ## Loading required package: stats
    ## 
    ## Attaching package: 'stats'
    ## 
    ## The following objects are masked from 'package:stats4':
    ## 
    ##     AIC, BIC, coef, confint, logLik, nobs, profile, update, vcov
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     filter, lag
    ## 
    ## 
    ## Attaching package: 'BiocGenerics'
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min
    ## 
    ## Loading required package: S4Vectors
    ## 
    ## Attaching package: 'S4Vectors'
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename
    ## 
    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname
    ## 
    ## Loading required package: IRanges
    ## 
    ## Attaching package: 'IRanges'
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce
    ## 
    ## Loading required package: GenomeInfoDb
    ## Loading required package: Biobase
    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.
    ## 
    ## 
    ## Attaching package: 'Biobase'
    ## 
    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians
    ## 
    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

``` r
library(tidySingleCellExperiment)
```

    ## Loading required package: ttservice
    ## 
    ## Attaching package: 'ttservice'
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     bind_cols, bind_rows
    ## 
    ## 
    ## Attaching package: 'tidySingleCellExperiment'
    ## 
    ## The following objects are masked from 'package:ttservice':
    ## 
    ##     bind_cols, bind_rows
    ## 
    ## The following object is masked from 'package:IRanges':
    ## 
    ##     slice
    ## 
    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     rename
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     filter
    ## 
    ## The following object is masked from 'package:matrixStats':
    ## 
    ##     count
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     add_count, bind_cols, bind_rows, count

``` r
hca_metadata = readRDS("/vast/scratch/users/mangiola.s/human_cell_atlas/metadata.rds")
options("restore_SingleCellExperiment_show" = TRUE)
```

``` r
hca_metadata |> 
    distinct(.sample)
```

    ## # A tibble: 12,429 × 1
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
    ## # … with 12,419 more rows

``` r
hca_metadata |> 
      distinct(file_id)
```

    ## # A tibble: 324 × 1
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
    ## # … with 314 more rows

``` r
hca_metadata |> 
      nrow()
```

    ## [1] 28981744

``` r
hca_metadata |> 
      distinct(cell_type)
```

    ## # A tibble: 497 × 1
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
    ## # … with 487 more rows

``` r
hca_metadata |> 
    distinct(tissue, file_id) |> 
    count(tissue) |> 
    arrange(desc(n))
```

    ## # A tibble: 165 × 2
    ##    tissue                    n
    ##    <chr>                 <int>
    ##  1 blood                    47
    ##  2 heart left ventricle     46
    ##  3 cortex of kidney         31
    ##  4 renal medulla            29
    ##  5 lung                     27
    ##  6 liver                    24
    ##  7 middle temporal gyrus    24
    ##  8 kidney                   19
    ##  9 intestine                18
    ## 10 thymus                   17
    ## # … with 155 more rows

Query

``` r
sce = 
    hca_metadata |> 
    filter(
        ethnicity == "African" & 
        assay |> str_detect("10x") & 
        tissue == "lung parenchyma" & 
        cell_type |> str_detect("CD4")
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
