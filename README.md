CuratedAtlasQueryR
================

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

`CuratedAtlasQuery` is a query interface that allow the programmatic
exploration and retrieval of the harmonised, curated and reannotated
CELLxGENE single-cell human cell atlas. Data can be retrieved at cell,
sample, or dataset levels based on filtering criteria.

Harmonised data is stored in the ARDC Nectar Research Cloud, and most
`CuratedAtlasQuery` functions interact with Nectar via web requests, so
a network connection is required for most functionality.

<img src="man/figures/logo.png" width="120x" height="139px" />

<img src="man/figures/svcf_logo.jpeg" width="155x" height="58px" /><img src="man/figures/czi_logo.png" width="129px" height="58px" /><img src="man/figures/bioconductor_logo.jpg" width="202px" height="58px" /><img src="man/figures/vca_logo.png" width="219px" height="58px" /><img src="man/figures/nectar_logo.png" width="180px" height="58px" />

# Query interface

## Installation

``` r
devtools::install_github("stemangiola/CuratedAtlasQueryR")
```

## Load the package

``` r
library(CuratedAtlasQueryR)
```

## Load and explore the metadata

### Load the metadata

``` r
metadata = get_metadata()

metadata
#> # Source:   table</vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/metadata.0.2.3.parquet> [?? x 56]
#> # Database: DuckDB 0.6.2-dev1166 [unknown@Linux 3.10.0-1160.88.1.el7.x86_64:R 4.2.1/:memory:]
#>    cell_ sample_ cell_…¹ cell_…² confi…³ cell_…⁴ cell_…⁵ cell_…⁶ sampl…⁷ _samp…⁸
#>    <chr> <chr>   <chr>   <chr>     <dbl> <chr>   <chr>   <chr>   <chr>   <chr>  
#>  1 AAAC… 689e2f… basal … basal_…       1 <NA>    <NA>    <NA>    f297c7… D17PrP…
#>  2 AAAC… 689e2f… basal … basal_…       1 <NA>    <NA>    <NA>    f297c7… D17PrP…
#>  3 AAAC… 689e2f… lumina… lumina…       1 <NA>    <NA>    <NA>    930938… D17PrP…
#>  4 AAAC… 689e2f… lumina… lumina…       1 <NA>    <NA>    <NA>    930938… D17PrP…
#>  5 AAAC… 689e2f… basal … basal_…       1 <NA>    <NA>    <NA>    f297c7… D17PrP…
#>  6 AAAC… 689e2f… basal … basal_…       1 <NA>    <NA>    <NA>    f297c7… D17PrP…
#>  7 AAAC… 689e2f… basal … basal_…       1 <NA>    <NA>    <NA>    f297c7… D17PrP…
#>  8 AAAC… 689e2f… basal … basal_…       1 <NA>    <NA>    <NA>    f297c7… D17PrP…
#>  9 AAAC… 689e2f… lumina… lumina…       1 <NA>    <NA>    <NA>    930938… D17PrP…
#> 10 AAAC… 689e2f… basal … basal_…       1 <NA>    <NA>    <NA>    f297c7… D17PrP…
#> # … with more rows, 46 more variables: assay <chr>,
#> #   assay_ontology_term_id <chr>, file_id_db <chr>,
#> #   cell_type_ontology_term_id <chr>, development_stage <chr>,
#> #   development_stage_ontology_term_id <chr>, disease <chr>,
#> #   disease_ontology_term_id <chr>, ethnicity <chr>,
#> #   ethnicity_ontology_term_id <chr>, experiment___ <chr>, file_id <chr>,
#> #   is_primary_data_x <chr>, organism <chr>, organism_ontology_term_id <chr>, …
```

### Explore the number of datasets per tissue

``` r
metadata |>
  dplyr::distinct(tissue, dataset_id) |> 
  dplyr::count(tissue)
#> # Source:   SQL [?? x 2]
#> # Database: DuckDB 0.6.2-dev1166 [unknown@Linux 3.10.0-1160.88.1.el7.x86_64:R 4.2.1/:memory:]
#>    tissue                          n
#>    <chr>                       <dbl>
#>  1 peripheral zone of prostate    10
#>  2 transition zone of prostate    10
#>  3 blood                          47
#>  4 intestine                      18
#>  5 middle temporal gyrus          24
#>  6 heart left ventricle           46
#>  7 apex of heart                  16
#>  8 heart right ventricle          16
#>  9 left cardiac atrium             7
#> 10 interventricular septum        16
#> # … with more rows
```

## Download single-cell RNA sequencing counts

### Query raw counts

``` r
single_cell_counts = 
    metadata |>
    dplyr::filter(
        ethnicity == "African" &
        stringr::str_like(assay, "%10x%") &
        tissue == "lung parenchyma" &
        stringr::str_like(cell_type, "%CD4%")
    ) |>
    get_single_cell_experiment()
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Downloading 0 files, totalling 0 GB
#> ℹ Reading files.
#> ℹ Compiling Single Cell Experiment.

single_cell_counts
#> class: SingleCellExperiment 
#> dim: 36229 1571 
#> metadata(0):
#> assays(1): counts
#> rownames(36229): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
#> rowData names(0):
#> colnames(1571): ACAGCCGGTCCGTTAA_F02526_1 GGGAATGAGCCCAGCT_F02526_1 ...
#>   TACAACGTCAGCATTG_SC84_1 CATTCGCTCAATACCG_F02526_1
#> colData names(56): sample_ cell_type ... updated_at_y original_cell_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

### Query counts scaled per million

This is helpful if just few genes are of interest, as they can be
compared across samples.

``` r
single_cell_counts = 
    metadata |>
    dplyr::filter(
        ethnicity == "African" &
        stringr::str_like(assay, "%10x%") &
        tissue == "lung parenchyma" &
        stringr::str_like(cell_type, "%CD4%")
    ) |>
    get_single_cell_experiment(assays = "cpm")
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Downloading 0 files, totalling 0 GB
#> ℹ Reading files.
#> ℹ Compiling Single Cell Experiment.

single_cell_counts
#> class: SingleCellExperiment 
#> dim: 36229 1571 
#> metadata(0):
#> assays(1): cpm
#> rownames(36229): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
#> rowData names(0):
#> colnames(1571): ACAGCCGGTCCGTTAA_F02526_1 GGGAATGAGCCCAGCT_F02526_1 ...
#>   TACAACGTCAGCATTG_SC84_1 CATTCGCTCAATACCG_F02526_1
#> colData names(56): sample_ cell_type ... updated_at_y original_cell_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

### Extract only a subset of genes

``` r
single_cell_counts = 
    metadata |>
    dplyr::filter(
        ethnicity == "African" &
        stringr::str_like(assay, "%10x%") &
        tissue == "lung parenchyma" &
        stringr::str_like(cell_type, "%CD4%")
    ) |>
    get_single_cell_experiment(assays = "cpm", features = "PUM1")
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Downloading 0 files, totalling 0 GB
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
#> colData names(56): sample_ cell_type ... updated_at_y original_cell_id
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
    metadata |>
    dplyr::filter(
        ethnicity == "African" &
        stringr::str_like(assay, "%10x%") &
        tissue == "lung parenchyma" &
        stringr::str_like(cell_type, "%CD4%")
    ) |>
    get_seurat()
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Downloading 0 files, totalling 0 GB
#> ℹ Reading files.
#> ℹ Compiling Single Cell Experiment.

single_cell_counts
#> An object of class Seurat 
#> 36229 features across 1571 samples within 1 assay 
#> Active assay: originalexp (36229 features, 0 variable features)
```

## Save your `SingleCellExperiment`

The returned `SingleCellExperiment` can be saved with two modalities, as
`.rds` or as `HDF5`.

### Saving as RDS (fast saving, slow reading)

Saving as `.rds` has the advantage that is very fast, the `.rds` file
occupies very little disk space as it only stored the links for the
files in your ache.

However it has the disadvantage that for big `SingleCellExperiment`
objects, which merge a lot of HDF5 from your
`get_single_cell_experiment` the display and manipulation is going to be
slow.

``` r
single_cell_counts |> saveRDS("single_cell_counts.rds")
```

### Saving as HDF5 (slow saving, fast reading)

Saving as `.rds` has the advantage that rewrites on disk a monolithic
`HDF5` and so displaying and manipulating large `SingleCellExperiment`
objects, which merge a lot of HDF5 from your
`get_single_cell_experiment`, is going to be fast.

However it has the disadvantage that the files are going to be larger as
they include the count information, and the saving process is going to
be slow for large objects.

``` r
single_cell_counts |> saveHDF5SummarizedExperiment("single_cell_counts")
```

## Visualise gene transcription

We can gather all CD14 monocytes cells and plot the distribution of
HLA-A across all tissues

    #> Loading required package: ttservice
    #> Loading required package: SingleCellExperiment
    #> Loading required package: SummarizedExperiment
    #> Loading required package: MatrixGenerics
    #> Loading required package: matrixStats
    #> 
    #> Attaching package: 'MatrixGenerics'
    #> The following objects are masked from 'package:matrixStats':
    #> 
    #>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    #>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    #>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    #>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    #>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    #>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    #>     colWeightedMeans, colWeightedMedians, colWeightedSds,
    #>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    #>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    #>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    #>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    #>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    #>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    #>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    #>     rowWeightedSds, rowWeightedVars
    #> Loading required package: GenomicRanges
    #> Loading required package: stats4
    #> Loading required package: BiocGenerics
    #> 
    #> Attaching package: 'BiocGenerics'
    #> The following objects are masked from 'package:stats':
    #> 
    #>     IQR, mad, sd, var, xtabs
    #> The following objects are masked from 'package:base':
    #> 
    #>     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    #>     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    #>     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    #>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    #>     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    #>     union, unique, unsplit, which.max, which.min
    #> Loading required package: S4Vectors
    #> 
    #> Attaching package: 'S4Vectors'
    #> The following objects are masked from 'package:base':
    #> 
    #>     expand.grid, I, unname
    #> Loading required package: IRanges
    #> Loading required package: GenomeInfoDb
    #> Loading required package: Biobase
    #> Welcome to Bioconductor
    #> 
    #>     Vignettes contain introductory material; view with
    #>     'browseVignettes()'. To cite Bioconductor, see
    #>     'citation("Biobase")', and for packages 'citation("pkgname")'.
    #> 
    #> Attaching package: 'Biobase'
    #> The following object is masked from 'package:MatrixGenerics':
    #> 
    #>     rowMedians
    #> The following objects are masked from 'package:matrixStats':
    #> 
    #>     anyMissing, rowMedians
    #> 
    #> Attaching package: 'tidySingleCellExperiment'
    #> The following object is masked from 'package:IRanges':
    #> 
    #>     slice
    #> The following object is masked from 'package:S4Vectors':
    #> 
    #>     rename
    #> The following object is masked from 'package:matrixStats':
    #> 
    #>     count
    #> The following objects are masked from 'package:ttservice':
    #> 
    #>     bind_cols, bind_rows
    #> The following object is masked from 'package:stats':
    #> 
    #>     filter
    #> ℹ Realising metadata.
    #> ℹ Synchronising files
    #> ℹ Downloading 0 files, totalling 0 GB
    #> ℹ Reading files.
    #> ℹ Compiling Single Cell Experiment.
    #> Warning: Transformation introduced infinite values in continuous y-axis

![](README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

    #> ℹ Realising metadata.
    #> ℹ Synchronising files
    #> ℹ Downloading 0 files, totalling 0 GB
    #> ℹ Reading files.
    #> ℹ Compiling Single Cell Experiment.
    #> Warning: Transformation introduced infinite values in continuous y-axis

![](README_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

``` r
library(tidySingleCellExperiment)
library(ggplot2)

metadata |>
  # Filter and subset
  dplyr::filter(cell_type_harmonised=="cd14 mono") |>

  # Get counts per million for HCA-A gene
  get_single_cell_experiment(assays = "cpm", features = "HLA-A") |> 
  
  # Plot (styling code have been omitted)
  join_features("HLA-A", shape = "wide") |> 
  ggplot(aes( disease, `HLA.A`,color = file_id)) +
  geom_jitter(shape=".") 
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Downloading 2 files, totalling 1.15 GB
#> ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/c5a05f23f9784a3be3bfa651198a48eb/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/c5a05f23f9784a3be3bfa651198a48eb/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■                  50% |  ETA:  9s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/c5a05f23f9784a3be3bfa651198a48eb/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/c5a05f23f9784a3be3bfa651198a48eb/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■                  50% |  ETA:  9s                                                                    ℹ Reading files.
#> ℹ Compiling Single Cell Experiment.
```

![](README_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

<img src="man/figures/HLA_A_disease_plot.png" width="525" />

``` r

metadata |> 
    
  # Filter and subset
  filter(cell_type_harmonised=="nk") |> 

  # Get counts per million for HCA-A gene 
  get_single_cell_experiment(assays = "cpm", features = "HLA-A") |> 

    # Plot (styling code have been omitted)
  join_features("HLA-A", shape = "wide") |> 
  ggplot(aes( tissue_harmonised, `HLA.A`,color = file_id)) +
  geom_jitter(shape=".") 
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> Warning in .f(.x[[i]], ...): NAs introduced by coercion to integer range

#> Warning in .f(.x[[i]], ...): NAs introduced by coercion to integer range
#> ℹ Downloading 368 files, totalling NA GB
#> ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/b862367ec3f3dc8588021352491edbd9/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/b862367ec3f3dc8588021352491edbd9/assays.h5
#> ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/b862367ec3f3dc8588021352491edbd9/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/b862367ec3f3dc8588021352491edbd9/se.rds
#> ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ec8be743daa60c684a5b0efbb2d8e7af/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ec8be743daa60c684a5b0efbb2d8e7af/assays.h5
#> ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ec8be743daa60c684a5b0efbb2d8e7af/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ec8be743daa60c684a5b0efbb2d8e7af/se.rds
#> ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/191fa376b1cca987246b55a55e7952f9/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/191fa376b1cca987246b55a55e7952f9/assays.h5
#> ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/191fa376b1cca987246b55a55e7952f9/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/191fa376b1cca987246b55a55e7952f9/se.rds
#> ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/abb6ff154b70c74e395add7dca444bbc/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/abb6ff154b70c74e395add7dca444bbc/assays.h5
#> Downloading files ■■                                 2% |  ETA: 36m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/abb6ff154b70c74e395add7dca444bbc/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/abb6ff154b70c74e395add7dca444bbc/se.rds
#> Downloading files ■■                                 2% |  ETA: 36mDownloading files ■■                                 2% |  ETA: 32m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/c6a230d16750eb88a10e9c3b8bacfb38/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/c6a230d16750eb88a10e9c3b8bacfb38/assays.h5
#> Downloading files ■■                                 2% |  ETA: 32mDownloading files ■■                                 2% |  ETA: 32m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/c6a230d16750eb88a10e9c3b8bacfb38/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/c6a230d16750eb88a10e9c3b8bacfb38/se.rds
#> Downloading files ■■                                 2% |  ETA: 32mDownloading files ■■                                 3% |  ETA: 29m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/d8c3f606d92524a41e1ca3383fc39e75/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/d8c3f606d92524a41e1ca3383fc39e75/assays.h5
#> Downloading files ■■                                 3% |  ETA: 29mDownloading files ■■                                 3% |  ETA: 27m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/d8c3f606d92524a41e1ca3383fc39e75/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/d8c3f606d92524a41e1ca3383fc39e75/se.rds
#> Downloading files ■■                                 3% |  ETA: 27m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/06026a8049be7fe60326d4bd26280626/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/06026a8049be7fe60326d4bd26280626/assays.h5
#> Downloading files ■■                                 3% |  ETA: 27mDownloading files ■■                                 4% |  ETA: 23m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/06026a8049be7fe60326d4bd26280626/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/06026a8049be7fe60326d4bd26280626/se.rds
#> Downloading files ■■                                 4% |  ETA: 23m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/6690e8e6b89123b1e11eef9d4a153e49/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/6690e8e6b89123b1e11eef9d4a153e49/assays.h5
#> Downloading files ■■                                 4% |  ETA: 23mDownloading files ■■                                 4% |  ETA: 21m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/6690e8e6b89123b1e11eef9d4a153e49/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/6690e8e6b89123b1e11eef9d4a153e49/se.rds
#> Downloading files ■■                                 4% |  ETA: 21m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/f2d70953e9665473f34ffc9acd242db8/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/f2d70953e9665473f34ffc9acd242db8/assays.h5
#> Downloading files ■■                                 4% |  ETA: 21mDownloading files ■■                                 5% |  ETA: 19m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/f2d70953e9665473f34ffc9acd242db8/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/f2d70953e9665473f34ffc9acd242db8/se.rds
#> Downloading files ■■                                 5% |  ETA: 19m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/05cd25cf665a14365e0307f3b97229e4/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/05cd25cf665a14365e0307f3b97229e4/assays.h5
#> Downloading files ■■                                 5% |  ETA: 19mDownloading files ■■■                                5% |  ETA: 21m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/05cd25cf665a14365e0307f3b97229e4/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/05cd25cf665a14365e0307f3b97229e4/se.rds
#> Downloading files ■■■                                5% |  ETA: 21mDownloading files ■■■                                5% |  ETA: 20m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/133f0de6e86d712c079af21344ae4501/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/133f0de6e86d712c079af21344ae4501/assays.h5
#> Downloading files ■■■                                5% |  ETA: 20mDownloading files ■■■                                6% |  ETA: 19m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/133f0de6e86d712c079af21344ae4501/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/133f0de6e86d712c079af21344ae4501/se.rds
#> Downloading files ■■■                                6% |  ETA: 19mDownloading files ■■■                                6% |  ETA: 18m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/6b312db0bbe47249bf6dbb382a0f25aa/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/6b312db0bbe47249bf6dbb382a0f25aa/assays.h5
#> Downloading files ■■■                                6% |  ETA: 18mDownloading files ■■■                                6% |  ETA: 18m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/6b312db0bbe47249bf6dbb382a0f25aa/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/6b312db0bbe47249bf6dbb382a0f25aa/se.rds
#> Downloading files ■■■                                6% |  ETA: 18mDownloading files ■■■                                7% |  ETA: 17m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ac7aa01423bb9483c28f49c6053d5d55/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ac7aa01423bb9483c28f49c6053d5d55/assays.h5
#> Downloading files ■■■                                7% |  ETA: 17mDownloading files ■■■                                7% |  ETA: 17m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ac7aa01423bb9483c28f49c6053d5d55/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ac7aa01423bb9483c28f49c6053d5d55/se.rds
#> Downloading files ■■■                                7% |  ETA: 17m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/acedf317b3f8af27496c5ede9988f597/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/acedf317b3f8af27496c5ede9988f597/assays.h5
#> Downloading files ■■■                                7% |  ETA: 17mDownloading files ■■■                                7% |  ETA: 16m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/acedf317b3f8af27496c5ede9988f597/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/acedf317b3f8af27496c5ede9988f597/se.rds
#> Downloading files ■■■                                7% |  ETA: 16m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/6202ef9ac861eae35c6af38b2b607fd6/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/6202ef9ac861eae35c6af38b2b607fd6/assays.h5
#> Downloading files ■■■                                7% |  ETA: 16mDownloading files ■■■                                8% |  ETA: 15m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/6202ef9ac861eae35c6af38b2b607fd6/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/6202ef9ac861eae35c6af38b2b607fd6/se.rds
#> Downloading files ■■■                                8% |  ETA: 15m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/5a44a4df123f5f718653d8fcc261fa91/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/5a44a4df123f5f718653d8fcc261fa91/assays.h5
#> Downloading files ■■■                                8% |  ETA: 15mDownloading files ■■■■                               8% |  ETA: 15m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/5a44a4df123f5f718653d8fcc261fa91/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/5a44a4df123f5f718653d8fcc261fa91/se.rds
#> Downloading files ■■■■                               8% |  ETA: 15m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e0f5bf13c94ae04af8b0dca002705bd0/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e0f5bf13c94ae04af8b0dca002705bd0/assays.h5
#> Downloading files ■■■■                               8% |  ETA: 15mDownloading files ■■■■                               9% |  ETA: 14m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e0f5bf13c94ae04af8b0dca002705bd0/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e0f5bf13c94ae04af8b0dca002705bd0/se.rds
#> Downloading files ■■■■                               9% |  ETA: 14m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/aef6b1cbe98851f255fca5d337632685/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/aef6b1cbe98851f255fca5d337632685/assays.h5
#> Downloading files ■■■■                               9% |  ETA: 14m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/aef6b1cbe98851f255fca5d337632685/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/aef6b1cbe98851f255fca5d337632685/se.rds
#> Downloading files ■■■■                               9% |  ETA: 14mDownloading files ■■■■                              10% |  ETA: 12m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/5fd70160f72771cecc39ac8ac5ed9af4/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/5fd70160f72771cecc39ac8ac5ed9af4/assays.h5
#> Downloading files ■■■■                              10% |  ETA: 12mDownloading files ■■■■                              10% |  ETA: 12m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/5fd70160f72771cecc39ac8ac5ed9af4/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/5fd70160f72771cecc39ac8ac5ed9af4/se.rds
#> Downloading files ■■■■                              10% |  ETA: 12m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/29079a713258550032406bd805ece089/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/29079a713258550032406bd805ece089/assays.h5
#> Downloading files ■■■■                              10% |  ETA: 12mDownloading files ■■■■                              11% |  ETA: 12m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/29079a713258550032406bd805ece089/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/29079a713258550032406bd805ece089/se.rds
#> Downloading files ■■■■                              11% |  ETA: 12m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e918153e603d5df595d559a4db4a5d8d/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e918153e603d5df595d559a4db4a5d8d/assays.h5
#> Downloading files ■■■■                              11% |  ETA: 12mDownloading files ■■■■                              11% |  ETA: 11m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e918153e603d5df595d559a4db4a5d8d/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e918153e603d5df595d559a4db4a5d8d/se.rds
#> Downloading files ■■■■                              11% |  ETA: 11m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/893d8537e318769108b4962020ddd846/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/893d8537e318769108b4962020ddd846/assays.h5
#> Downloading files ■■■■                              11% |  ETA: 11m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/893d8537e318769108b4962020ddd846/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/893d8537e318769108b4962020ddd846/se.rds
#> Downloading files ■■■■                              11% |  ETA: 11mDownloading files ■■■■■                             12% |  ETA: 10m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/da596713eaf8ca89717f9a02d1990f04/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/da596713eaf8ca89717f9a02d1990f04/assays.h5
#> Downloading files ■■■■■                             12% |  ETA: 10mDownloading files ■■■■■                             12% |  ETA: 10m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/da596713eaf8ca89717f9a02d1990f04/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/da596713eaf8ca89717f9a02d1990f04/se.rds
#> Downloading files ■■■■■                             12% |  ETA: 10m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/d5a66283b428cafda33d33f31a08af60/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/d5a66283b428cafda33d33f31a08af60/assays.h5
#> Downloading files ■■■■■                             12% |  ETA: 10mDownloading files ■■■■■                             13% |  ETA: 10m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/d5a66283b428cafda33d33f31a08af60/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/d5a66283b428cafda33d33f31a08af60/se.rds
#> Downloading files ■■■■■                             13% |  ETA: 10mDownloading files ■■■■■                             13% |  ETA: 10m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/2f2504999dc8c56695dcdd1a6e4aa5de/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/2f2504999dc8c56695dcdd1a6e4aa5de/assays.h5
#> Downloading files ■■■■■                             13% |  ETA: 10mDownloading files ■■■■■                             13% |  ETA:  9m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/2f2504999dc8c56695dcdd1a6e4aa5de/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/2f2504999dc8c56695dcdd1a6e4aa5de/se.rds
#> Downloading files ■■■■■                             13% |  ETA:  9mDownloading files ■■■■■                             14% |  ETA:  9m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/678238a66b85415246a8f4ad09f8441c/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/678238a66b85415246a8f4ad09f8441c/assays.h5
#> Downloading files ■■■■■                             14% |  ETA:  9mDownloading files ■■■■■                             14% |  ETA:  9m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/678238a66b85415246a8f4ad09f8441c/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/678238a66b85415246a8f4ad09f8441c/se.rds
#> Downloading files ■■■■■                             14% |  ETA:  9m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a74416e8d11b74a611d95d5bcd910460/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a74416e8d11b74a611d95d5bcd910460/assays.h5
#> Downloading files ■■■■■                             14% |  ETA:  9m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a74416e8d11b74a611d95d5bcd910460/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a74416e8d11b74a611d95d5bcd910460/se.rds
#> Downloading files ■■■■■                             14% |  ETA:  9mDownloading files ■■■■■                             15% |  ETA:  9m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/4662daaf870f56c79448936461afe154/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/4662daaf870f56c79448936461afe154/assays.h5
#> Downloading files ■■■■■                             15% |  ETA:  9mDownloading files ■■■■■                             15% |  ETA:  8m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/4662daaf870f56c79448936461afe154/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/4662daaf870f56c79448936461afe154/se.rds
#> Downloading files ■■■■■                             15% |  ETA:  8m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/d6a5579a96f4ecd1b059b6623c4d98a9/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/d6a5579a96f4ecd1b059b6623c4d98a9/assays.h5
#> Downloading files ■■■■■                             15% |  ETA:  8mDownloading files ■■■■■■                            15% |  ETA:  9m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/d6a5579a96f4ecd1b059b6623c4d98a9/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/d6a5579a96f4ecd1b059b6623c4d98a9/se.rds
#> Downloading files ■■■■■■                            15% |  ETA:  9mDownloading files ■■■■■■                            16% |  ETA:  9m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/6291fe62ea33386b7ae347f7baf9b7ab/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/6291fe62ea33386b7ae347f7baf9b7ab/assays.h5
#> Downloading files ■■■■■■                            16% |  ETA:  9mDownloading files ■■■■■■                            16% |  ETA:  8m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/6291fe62ea33386b7ae347f7baf9b7ab/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/6291fe62ea33386b7ae347f7baf9b7ab/se.rds
#> Downloading files ■■■■■■                            16% |  ETA:  8m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/3c71a704f05efd16c117c56b0977df38/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/3c71a704f05efd16c117c56b0977df38/assays.h5
#> Downloading files ■■■■■■                            16% |  ETA:  8mDownloading files ■■■■■■                            17% |  ETA:  9m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/3c71a704f05efd16c117c56b0977df38/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/3c71a704f05efd16c117c56b0977df38/se.rds
#> Downloading files ■■■■■■                            17% |  ETA:  9m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/7b22a2aaaffed0d0f40b2b8a0e2022fe/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/7b22a2aaaffed0d0f40b2b8a0e2022fe/assays.h5
#> Downloading files ■■■■■■                            17% |  ETA:  9mDownloading files ■■■■■■                            17% |  ETA:  9m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/7b22a2aaaffed0d0f40b2b8a0e2022fe/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/7b22a2aaaffed0d0f40b2b8a0e2022fe/se.rds
#> Downloading files ■■■■■■                            17% |  ETA:  9mDownloading files ■■■■■■                            17% |  ETA:  9m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/6667c08f3ef242dd0db7910755bc3e71/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/6667c08f3ef242dd0db7910755bc3e71/assays.h5
#> Downloading files ■■■■■■                            17% |  ETA:  9mDownloading files ■■■■■■                            18% |  ETA:  8m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/6667c08f3ef242dd0db7910755bc3e71/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/6667c08f3ef242dd0db7910755bc3e71/se.rds
#> Downloading files ■■■■■■                            18% |  ETA:  8m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/1e3b62e0491436a77998c8fbdc9f1db1/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/1e3b62e0491436a77998c8fbdc9f1db1/assays.h5
#> Downloading files ■■■■■■                            18% |  ETA:  8mDownloading files ■■■■■■                            18% |  ETA:  8m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/1e3b62e0491436a77998c8fbdc9f1db1/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/1e3b62e0491436a77998c8fbdc9f1db1/se.rds
#> Downloading files ■■■■■■                            18% |  ETA:  8m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/59a27bb4b97d861d23e9b6f3a83c1d80/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/59a27bb4b97d861d23e9b6f3a83c1d80/assays.h5
#> Downloading files ■■■■■■                            18% |  ETA:  8m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/59a27bb4b97d861d23e9b6f3a83c1d80/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/59a27bb4b97d861d23e9b6f3a83c1d80/se.rds
#> Downloading files ■■■■■■                            18% |  ETA:  8m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/98753cbddb01a0c6a1fda9ef7935a11e/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/98753cbddb01a0c6a1fda9ef7935a11e/assays.h5
#> Downloading files ■■■■■■                            18% |  ETA:  8mDownloading files ■■■■■■■                           19% |  ETA:  8m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/98753cbddb01a0c6a1fda9ef7935a11e/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/98753cbddb01a0c6a1fda9ef7935a11e/se.rds
#> Downloading files ■■■■■■■                           19% |  ETA:  8m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/84c8432647d0a781608da4bdeb88e216/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/84c8432647d0a781608da4bdeb88e216/assays.h5
#> Downloading files ■■■■■■■                           19% |  ETA:  8m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/84c8432647d0a781608da4bdeb88e216/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/84c8432647d0a781608da4bdeb88e216/se.rds
#> Downloading files ■■■■■■■                           19% |  ETA:  8mDownloading files ■■■■■■■                           20% |  ETA:  7m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/79ffeada8f92cf0feec712c94543f43f/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/79ffeada8f92cf0feec712c94543f43f/assays.h5
#> Downloading files ■■■■■■■                           20% |  ETA:  7m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/79ffeada8f92cf0feec712c94543f43f/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/79ffeada8f92cf0feec712c94543f43f/se.rds
#> Downloading files ■■■■■■■                           20% |  ETA:  7mDownloading files ■■■■■■■                           21% |  ETA:  7m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/134e0eb2d115cf884b946f8ef47cfe2d/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/134e0eb2d115cf884b946f8ef47cfe2d/assays.h5
#> Downloading files ■■■■■■■                           21% |  ETA:  7m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/134e0eb2d115cf884b946f8ef47cfe2d/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/134e0eb2d115cf884b946f8ef47cfe2d/se.rds
#> Downloading files ■■■■■■■                           21% |  ETA:  7m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/50e3fe44f76e63bff6b19d18d2e4d22e/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/50e3fe44f76e63bff6b19d18d2e4d22e/assays.h5
#> Downloading files ■■■■■■■                           21% |  ETA:  7mDownloading files ■■■■■■■                           21% |  ETA:  7m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/50e3fe44f76e63bff6b19d18d2e4d22e/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/50e3fe44f76e63bff6b19d18d2e4d22e/se.rds
#> Downloading files ■■■■■■■                           21% |  ETA:  7mDownloading files ■■■■■■■■                          22% |  ETA:  7m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/cc07d91d2bd35da21aed37cefae2b0a7/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/cc07d91d2bd35da21aed37cefae2b0a7/assays.h5
#> Downloading files ■■■■■■■■                          22% |  ETA:  7mDownloading files ■■■■■■■■                          22% |  ETA:  7m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/cc07d91d2bd35da21aed37cefae2b0a7/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/cc07d91d2bd35da21aed37cefae2b0a7/se.rds
#> Downloading files ■■■■■■■■                          22% |  ETA:  7m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/b1f49796da75efb68753ba09e9bf4c47/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/b1f49796da75efb68753ba09e9bf4c47/assays.h5
#> Downloading files ■■■■■■■■                          22% |  ETA:  7mDownloading files ■■■■■■■■                          23% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/b1f49796da75efb68753ba09e9bf4c47/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/b1f49796da75efb68753ba09e9bf4c47/se.rds
#> Downloading files ■■■■■■■■                          23% |  ETA:  6mDownloading files ■■■■■■■■                          23% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a5f0308e5f7c43a1ce8ac06a053e7605/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a5f0308e5f7c43a1ce8ac06a053e7605/assays.h5
#> Downloading files ■■■■■■■■                          23% |  ETA:  6mDownloading files ■■■■■■■■                          23% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a5f0308e5f7c43a1ce8ac06a053e7605/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a5f0308e5f7c43a1ce8ac06a053e7605/se.rds
#> Downloading files ■■■■■■■■                          23% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/3e242480ed6e12fbc6cf89f6b9f9f96e/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/3e242480ed6e12fbc6cf89f6b9f9f96e/assays.h5
#> Downloading files ■■■■■■■■                          23% |  ETA:  6mDownloading files ■■■■■■■■                          24% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/3e242480ed6e12fbc6cf89f6b9f9f96e/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/3e242480ed6e12fbc6cf89f6b9f9f96e/se.rds
#> Downloading files ■■■■■■■■                          24% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a6cae194bbcaffcd8e814b4918098415/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a6cae194bbcaffcd8e814b4918098415/assays.h5
#> Downloading files ■■■■■■■■                          24% |  ETA:  6mDownloading files ■■■■■■■■                          24% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a6cae194bbcaffcd8e814b4918098415/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a6cae194bbcaffcd8e814b4918098415/se.rds
#> Downloading files ■■■■■■■■                          24% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/c6997be738caf6c6967129b02553e492/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/c6997be738caf6c6967129b02553e492/assays.h5
#> Downloading files ■■■■■■■■                          24% |  ETA:  6mDownloading files ■■■■■■■■                          25% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/c6997be738caf6c6967129b02553e492/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/c6997be738caf6c6967129b02553e492/se.rds
#> Downloading files ■■■■■■■■                          25% |  ETA:  6mDownloading files ■■■■■■■■■                         25% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/bf32a401275affb54bdd8145779cfb81/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/bf32a401275affb54bdd8145779cfb81/assays.h5
#> Downloading files ■■■■■■■■■                         25% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/bf32a401275affb54bdd8145779cfb81/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/bf32a401275affb54bdd8145779cfb81/se.rds
#> Downloading files ■■■■■■■■■                         25% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/dd40ab1cff3d95c4a66c4fb8cb8956c3/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/dd40ab1cff3d95c4a66c4fb8cb8956c3/assays.h5
#> Downloading files ■■■■■■■■■                         25% |  ETA:  6mDownloading files ■■■■■■■■■                         26% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/dd40ab1cff3d95c4a66c4fb8cb8956c3/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/dd40ab1cff3d95c4a66c4fb8cb8956c3/se.rds
#> Downloading files ■■■■■■■■■                         26% |  ETA:  6mDownloading files ■■■■■■■■■                         26% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ecc5c83172b38a4c2d703cef9fee6600/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ecc5c83172b38a4c2d703cef9fee6600/assays.h5
#> Downloading files ■■■■■■■■■                         26% |  ETA:  6mDownloading files ■■■■■■■■■                         26% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ecc5c83172b38a4c2d703cef9fee6600/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ecc5c83172b38a4c2d703cef9fee6600/se.rds
#> Downloading files ■■■■■■■■■                         26% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/d8fe285941ee94128c466ffe2fc72214/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/d8fe285941ee94128c466ffe2fc72214/assays.h5
#> Downloading files ■■■■■■■■■                         26% |  ETA:  6mDownloading files ■■■■■■■■■                         27% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/d8fe285941ee94128c466ffe2fc72214/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/d8fe285941ee94128c466ffe2fc72214/se.rds
#> Downloading files ■■■■■■■■■                         27% |  ETA:  6mDownloading files ■■■■■■■■■                         27% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/94c522ee71123a3b115bdebce647fb88/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/94c522ee71123a3b115bdebce647fb88/assays.h5
#> Downloading files ■■■■■■■■■                         27% |  ETA:  5mDownloading files ■■■■■■■■■                         27% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/94c522ee71123a3b115bdebce647fb88/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/94c522ee71123a3b115bdebce647fb88/se.rds
#> Downloading files ■■■■■■■■■                         27% |  ETA:  6m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/bb270e6b0aea2bba1f8f044402b50d0f/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/bb270e6b0aea2bba1f8f044402b50d0f/assays.h5
#> Downloading files ■■■■■■■■■                         27% |  ETA:  6mDownloading files ■■■■■■■■■                         28% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/bb270e6b0aea2bba1f8f044402b50d0f/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/bb270e6b0aea2bba1f8f044402b50d0f/se.rds
#> Downloading files ■■■■■■■■■                         28% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/c032b8738a8d06e9f261f76564de0ec8/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/c032b8738a8d06e9f261f76564de0ec8/assays.h5
#> Downloading files ■■■■■■■■■                         28% |  ETA:  5mDownloading files ■■■■■■■■■■                        29% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/c032b8738a8d06e9f261f76564de0ec8/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/c032b8738a8d06e9f261f76564de0ec8/se.rds
#> Downloading files ■■■■■■■■■■                        29% |  ETA:  5mDownloading files ■■■■■■■■■■                        29% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/4c39fdd1e8eb9a146f6e0766a8fbe4c7/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/4c39fdd1e8eb9a146f6e0766a8fbe4c7/assays.h5
#> Downloading files ■■■■■■■■■■                        29% |  ETA:  5mDownloading files ■■■■■■■■■■                        29% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/4c39fdd1e8eb9a146f6e0766a8fbe4c7/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/4c39fdd1e8eb9a146f6e0766a8fbe4c7/se.rds
#> Downloading files ■■■■■■■■■■                        29% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/cb9dccb0d368871cfe6019dfe074b0b5/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/cb9dccb0d368871cfe6019dfe074b0b5/assays.h5
#> Downloading files ■■■■■■■■■■                        29% |  ETA:  5mDownloading files ■■■■■■■■■■                        30% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/cb9dccb0d368871cfe6019dfe074b0b5/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/cb9dccb0d368871cfe6019dfe074b0b5/se.rds
#> Downloading files ■■■■■■■■■■                        30% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/7b435cb2bf03c171248967cea26f40b1/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/7b435cb2bf03c171248967cea26f40b1/assays.h5
#> Downloading files ■■■■■■■■■■                        30% |  ETA:  5mDownloading files ■■■■■■■■■■                        30% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/7b435cb2bf03c171248967cea26f40b1/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/7b435cb2bf03c171248967cea26f40b1/se.rds
#> Downloading files ■■■■■■■■■■                        30% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/8e89e1787049eb0852763d690ebcbb2c/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/8e89e1787049eb0852763d690ebcbb2c/assays.h5
#> Downloading files ■■■■■■■■■■                        30% |  ETA:  5mDownloading files ■■■■■■■■■■                        31% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/8e89e1787049eb0852763d690ebcbb2c/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/8e89e1787049eb0852763d690ebcbb2c/se.rds
#> Downloading files ■■■■■■■■■■                        31% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/685d38877a6724417e79666a7e01c14f/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/685d38877a6724417e79666a7e01c14f/assays.h5
#> Downloading files ■■■■■■■■■■                        31% |  ETA:  5mDownloading files ■■■■■■■■■■                        31% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/685d38877a6724417e79666a7e01c14f/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/685d38877a6724417e79666a7e01c14f/se.rds
#> Downloading files ■■■■■■■■■■                        31% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e407548a991c9811a22756149481b771/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e407548a991c9811a22756149481b771/assays.h5
#> Downloading files ■■■■■■■■■■                        31% |  ETA:  5mDownloading files ■■■■■■■■■■■                       32% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e407548a991c9811a22756149481b771/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e407548a991c9811a22756149481b771/se.rds
#> Downloading files ■■■■■■■■■■■                       32% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/42d73ccd9afad439b978237d14831c4a/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/42d73ccd9afad439b978237d14831c4a/assays.h5
#> Downloading files ■■■■■■■■■■■                       32% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/42d73ccd9afad439b978237d14831c4a/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/42d73ccd9afad439b978237d14831c4a/se.rds
#> Downloading files ■■■■■■■■■■■                       32% |  ETA:  5mDownloading files ■■■■■■■■■■■                       33% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/316eabbc8257f818249861799a2bbaf4/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/316eabbc8257f818249861799a2bbaf4/assays.h5
#> Downloading files ■■■■■■■■■■■                       33% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/316eabbc8257f818249861799a2bbaf4/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/316eabbc8257f818249861799a2bbaf4/se.rds
#> Downloading files ■■■■■■■■■■■                       33% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9cd0ac3bc4caa46c894b30b1aca6d0d4/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9cd0ac3bc4caa46c894b30b1aca6d0d4/assays.h5
#> Downloading files ■■■■■■■■■■■                       33% |  ETA:  5mDownloading files ■■■■■■■■■■■                       33% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9cd0ac3bc4caa46c894b30b1aca6d0d4/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9cd0ac3bc4caa46c894b30b1aca6d0d4/se.rds
#> Downloading files ■■■■■■■■■■■                       33% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/04ba8db05104a89df380837866e08290/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/04ba8db05104a89df380837866e08290/assays.h5
#> Downloading files ■■■■■■■■■■■                       33% |  ETA:  5mDownloading files ■■■■■■■■■■■                       34% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/04ba8db05104a89df380837866e08290/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/04ba8db05104a89df380837866e08290/se.rds
#> Downloading files ■■■■■■■■■■■                       34% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e1f0618820b20b2f5a7f195be80c3469/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e1f0618820b20b2f5a7f195be80c3469/assays.h5
#> Downloading files ■■■■■■■■■■■                       34% |  ETA:  5mDownloading files ■■■■■■■■■■■                       35% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e1f0618820b20b2f5a7f195be80c3469/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e1f0618820b20b2f5a7f195be80c3469/se.rds
#> Downloading files ■■■■■■■■■■■                       35% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/740ca6249b0054a142a525ef53209619/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/740ca6249b0054a142a525ef53209619/assays.h5
#> Downloading files ■■■■■■■■■■■                       35% |  ETA:  5m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/740ca6249b0054a142a525ef53209619/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/740ca6249b0054a142a525ef53209619/se.rds
#> Downloading files ■■■■■■■■■■■                       35% |  ETA:  5mDownloading files ■■■■■■■■■■■■                      35% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ff6e44761924104cffba9d60d45bba91/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ff6e44761924104cffba9d60d45bba91/assays.h5
#> Downloading files ■■■■■■■■■■■■                      35% |  ETA:  4mDownloading files ■■■■■■■■■■■■                      36% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ff6e44761924104cffba9d60d45bba91/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ff6e44761924104cffba9d60d45bba91/se.rds
#> Downloading files ■■■■■■■■■■■■                      36% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/d28a9f52f7b6bdf78a08681fb03533a9/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/d28a9f52f7b6bdf78a08681fb03533a9/assays.h5
#> Downloading files ■■■■■■■■■■■■                      36% |  ETA:  4mDownloading files ■■■■■■■■■■■■                      36% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/d28a9f52f7b6bdf78a08681fb03533a9/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/d28a9f52f7b6bdf78a08681fb03533a9/se.rds
#> Downloading files ■■■■■■■■■■■■                      36% |  ETA:  4mDownloading files ■■■■■■■■■■■■                      36% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/8bb1e49d39bc1b69869bee715b8529eb/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/8bb1e49d39bc1b69869bee715b8529eb/assays.h5
#> Downloading files ■■■■■■■■■■■■                      36% |  ETA:  4mDownloading files ■■■■■■■■■■■■                      37% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/8bb1e49d39bc1b69869bee715b8529eb/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/8bb1e49d39bc1b69869bee715b8529eb/se.rds
#> Downloading files ■■■■■■■■■■■■                      37% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/735fc9a88f17695c616ba7d956fc3845/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/735fc9a88f17695c616ba7d956fc3845/assays.h5
#> Downloading files ■■■■■■■■■■■■                      37% |  ETA:  4mDownloading files ■■■■■■■■■■■■                      37% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/735fc9a88f17695c616ba7d956fc3845/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/735fc9a88f17695c616ba7d956fc3845/se.rds
#> Downloading files ■■■■■■■■■■■■                      37% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ba2031c6401ced0bad9316eb27d99122/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ba2031c6401ced0bad9316eb27d99122/assays.h5
#> Downloading files ■■■■■■■■■■■■                      37% |  ETA:  4mDownloading files ■■■■■■■■■■■■                      38% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ba2031c6401ced0bad9316eb27d99122/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ba2031c6401ced0bad9316eb27d99122/se.rds
#> Downloading files ■■■■■■■■■■■■                      38% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9b8117a97ed183e9385c0e9956717da8/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9b8117a97ed183e9385c0e9956717da8/assays.h5
#> Downloading files ■■■■■■■■■■■■                      38% |  ETA:  4mDownloading files ■■■■■■■■■■■■                      38% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9b8117a97ed183e9385c0e9956717da8/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9b8117a97ed183e9385c0e9956717da8/se.rds
#> Downloading files ■■■■■■■■■■■■                      38% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/25ba1c61a8e295fff3617b094da1230d/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/25ba1c61a8e295fff3617b094da1230d/assays.h5
#> Downloading files ■■■■■■■■■■■■                      38% |  ETA:  4mDownloading files ■■■■■■■■■■■■■                     39% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/25ba1c61a8e295fff3617b094da1230d/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/25ba1c61a8e295fff3617b094da1230d/se.rds
#> Downloading files ■■■■■■■■■■■■■                     39% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ce4ebaa3ce6cec7f8efa2d1ed9a49d76/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ce4ebaa3ce6cec7f8efa2d1ed9a49d76/assays.h5
#> Downloading files ■■■■■■■■■■■■■                     39% |  ETA:  4mDownloading files ■■■■■■■■■■■■■                     39% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ce4ebaa3ce6cec7f8efa2d1ed9a49d76/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ce4ebaa3ce6cec7f8efa2d1ed9a49d76/se.rds
#> Downloading files ■■■■■■■■■■■■■                     39% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a15e91afd99e3fb782c84b881c350c5d/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a15e91afd99e3fb782c84b881c350c5d/assays.h5
#> Downloading files ■■■■■■■■■■■■■                     39% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a15e91afd99e3fb782c84b881c350c5d/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a15e91afd99e3fb782c84b881c350c5d/se.rds
#> Downloading files ■■■■■■■■■■■■■                     39% |  ETA:  4mDownloading files ■■■■■■■■■■■■■                     40% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/1b147255e0d8840ecc3941f4d3e140ef/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/1b147255e0d8840ecc3941f4d3e140ef/assays.h5
#> Downloading files ■■■■■■■■■■■■■                     40% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/1b147255e0d8840ecc3941f4d3e140ef/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/1b147255e0d8840ecc3941f4d3e140ef/se.rds
#> Downloading files ■■■■■■■■■■■■■                     40% |  ETA:  4mDownloading files ■■■■■■■■■■■■■                     41% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/f6dc528d7af79d09c3cdd5f272962db3/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/f6dc528d7af79d09c3cdd5f272962db3/assays.h5
#> Downloading files ■■■■■■■■■■■■■                     41% |  ETA:  4mDownloading files ■■■■■■■■■■■■■                     41% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/f6dc528d7af79d09c3cdd5f272962db3/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/f6dc528d7af79d09c3cdd5f272962db3/se.rds
#> Downloading files ■■■■■■■■■■■■■                     41% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/674dff612b0a2a575a7249ca1ddf8e55/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/674dff612b0a2a575a7249ca1ddf8e55/assays.h5
#> Downloading files ■■■■■■■■■■■■■                     41% |  ETA:  4mDownloading files ■■■■■■■■■■■■■                     42% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/674dff612b0a2a575a7249ca1ddf8e55/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/674dff612b0a2a575a7249ca1ddf8e55/se.rds
#> Downloading files ■■■■■■■■■■■■■                     42% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■                    42% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a4aa1b10f63f711e14c824629175a86d/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a4aa1b10f63f711e14c824629175a86d/assays.h5
#> Downloading files ■■■■■■■■■■■■■■                    42% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■                    42% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a4aa1b10f63f711e14c824629175a86d/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a4aa1b10f63f711e14c824629175a86d/se.rds
#> Downloading files ■■■■■■■■■■■■■■                    42% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9a91665737dd9d894e1295bf73ce4783/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9a91665737dd9d894e1295bf73ce4783/assays.h5
#> Downloading files ■■■■■■■■■■■■■■                    42% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■                    43% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9a91665737dd9d894e1295bf73ce4783/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9a91665737dd9d894e1295bf73ce4783/se.rds
#> Downloading files ■■■■■■■■■■■■■■                    43% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■                    43% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/1239b9115140fa94d499f3014c1b273b/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/1239b9115140fa94d499f3014c1b273b/assays.h5
#> Downloading files ■■■■■■■■■■■■■■                    43% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■                    43% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/1239b9115140fa94d499f3014c1b273b/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/1239b9115140fa94d499f3014c1b273b/se.rds
#> Downloading files ■■■■■■■■■■■■■■                    43% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/179b344a46f65d444a39a584137e644d/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/179b344a46f65d444a39a584137e644d/assays.h5
#> Downloading files ■■■■■■■■■■■■■■                    43% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■                    44% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/179b344a46f65d444a39a584137e644d/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/179b344a46f65d444a39a584137e644d/se.rds
#> Downloading files ■■■■■■■■■■■■■■                    44% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/b43070d411e72ee7ff835f79821e0e07/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/b43070d411e72ee7ff835f79821e0e07/assays.h5
#> Downloading files ■■■■■■■■■■■■■■                    44% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■                    44% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/b43070d411e72ee7ff835f79821e0e07/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/b43070d411e72ee7ff835f79821e0e07/se.rds
#> Downloading files ■■■■■■■■■■■■■■                    44% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/73323f5e6eae774b5dc8dce6d446a017/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/73323f5e6eae774b5dc8dce6d446a017/assays.h5
#> Downloading files ■■■■■■■■■■■■■■                    44% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■                    45% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/73323f5e6eae774b5dc8dce6d446a017/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/73323f5e6eae774b5dc8dce6d446a017/se.rds
#> Downloading files ■■■■■■■■■■■■■■                    45% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a8db324766da5a04e0abfe9557ff9420/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a8db324766da5a04e0abfe9557ff9420/assays.h5
#> Downloading files ■■■■■■■■■■■■■■                    45% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■■                   45% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a8db324766da5a04e0abfe9557ff9420/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a8db324766da5a04e0abfe9557ff9420/se.rds
#> Downloading files ■■■■■■■■■■■■■■■                   45% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■■                   46% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/6644ccf4ab45049b83c4e7f4a918ec76/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/6644ccf4ab45049b83c4e7f4a918ec76/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■                   46% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/6644ccf4ab45049b83c4e7f4a918ec76/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/6644ccf4ab45049b83c4e7f4a918ec76/se.rds
#> Downloading files ■■■■■■■■■■■■■■■                   46% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■■                   46% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/8319ac7610d93fd36b827d89960b645a/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/8319ac7610d93fd36b827d89960b645a/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■                   46% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■■                   46% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/8319ac7610d93fd36b827d89960b645a/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/8319ac7610d93fd36b827d89960b645a/se.rds
#> Downloading files ■■■■■■■■■■■■■■■                   46% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/02136ca55c87ceeffea8241b9e394f20/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/02136ca55c87ceeffea8241b9e394f20/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■                   46% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■■                   47% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/02136ca55c87ceeffea8241b9e394f20/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/02136ca55c87ceeffea8241b9e394f20/se.rds
#> Downloading files ■■■■■■■■■■■■■■■                   47% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e746e1e6c617304549e1b78c1abaef85/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e746e1e6c617304549e1b78c1abaef85/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■                   47% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■■                   48% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e746e1e6c617304549e1b78c1abaef85/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e746e1e6c617304549e1b78c1abaef85/se.rds
#> Downloading files ■■■■■■■■■■■■■■■                   48% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■■                   48% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/8415c560953210dfe19b2aa624d98719/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/8415c560953210dfe19b2aa624d98719/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■                   48% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■■                   48% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/8415c560953210dfe19b2aa624d98719/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/8415c560953210dfe19b2aa624d98719/se.rds
#> Downloading files ■■■■■■■■■■■■■■■                   48% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/f543c85ad36ea704998d1e482902dae0/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/f543c85ad36ea704998d1e482902dae0/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■                   48% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■■■                  49% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/f543c85ad36ea704998d1e482902dae0/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/f543c85ad36ea704998d1e482902dae0/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■                  49% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/055146c01d6189bebd37d53e3b87ab6f/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/055146c01d6189bebd37d53e3b87ab6f/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■                  49% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■■■                  49% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/055146c01d6189bebd37d53e3b87ab6f/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/055146c01d6189bebd37d53e3b87ab6f/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■                  49% |  ETA:  4m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/c3e0ab923f1e25f191f37339d3094ba9/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/c3e0ab923f1e25f191f37339d3094ba9/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■                  49% |  ETA:  4mDownloading files ■■■■■■■■■■■■■■■■                  50% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/c3e0ab923f1e25f191f37339d3094ba9/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/c3e0ab923f1e25f191f37339d3094ba9/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■                  50% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/c6570618e43e0101fc463e9c2ee5e4b7/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/c6570618e43e0101fc463e9c2ee5e4b7/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■                  50% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■                  50% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/c6570618e43e0101fc463e9c2ee5e4b7/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/c6570618e43e0101fc463e9c2ee5e4b7/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■                  50% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/0ddb709253714902f235e1c2b74af824/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/0ddb709253714902f235e1c2b74af824/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■                  50% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■                  51% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/0ddb709253714902f235e1c2b74af824/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/0ddb709253714902f235e1c2b74af824/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■                  51% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/2b1624d83e7e9464270169961b8b7970/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/2b1624d83e7e9464270169961b8b7970/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■                  51% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■                  51% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/2b1624d83e7e9464270169961b8b7970/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/2b1624d83e7e9464270169961b8b7970/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■                  51% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■                  52% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/afef95b754f5a4051363fc97e8d7a5a3/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/afef95b754f5a4051363fc97e8d7a5a3/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■                  52% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■                 52% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/afef95b754f5a4051363fc97e8d7a5a3/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/afef95b754f5a4051363fc97e8d7a5a3/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■                 52% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/df689ee2619c0183e96c2ef50cc06037/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/df689ee2619c0183e96c2ef50cc06037/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■                 52% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■                 52% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/df689ee2619c0183e96c2ef50cc06037/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/df689ee2619c0183e96c2ef50cc06037/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■                 52% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/cd2569c37c318777f5d6a82a38680362/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/cd2569c37c318777f5d6a82a38680362/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■                 52% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■                 53% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/cd2569c37c318777f5d6a82a38680362/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/cd2569c37c318777f5d6a82a38680362/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■                 53% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■                 53% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/8702ff3fb0f02bf4cb84225f438fa2bc/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/8702ff3fb0f02bf4cb84225f438fa2bc/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■                 53% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■                 54% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/8702ff3fb0f02bf4cb84225f438fa2bc/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/8702ff3fb0f02bf4cb84225f438fa2bc/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■                 54% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■                 54% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/28d4b90ef093453cb6c813a8ceda0094/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/28d4b90ef093453cb6c813a8ceda0094/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■                 54% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■                 54% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/28d4b90ef093453cb6c813a8ceda0094/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/28d4b90ef093453cb6c813a8ceda0094/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■                 54% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/dc521a0aa16399812f6331c6d067c14a/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/dc521a0aa16399812f6331c6d067c14a/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■                 54% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■                 55% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/dc521a0aa16399812f6331c6d067c14a/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/dc521a0aa16399812f6331c6d067c14a/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■                 55% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ebd3b8cf48d5fa2a671f3c0cd7ce1b02/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ebd3b8cf48d5fa2a671f3c0cd7ce1b02/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■                 55% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■■                55% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ebd3b8cf48d5fa2a671f3c0cd7ce1b02/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ebd3b8cf48d5fa2a671f3c0cd7ce1b02/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■                55% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■■                55% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/61024c80c683713497df268c6b3b7ff7/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/61024c80c683713497df268c6b3b7ff7/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■                55% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/61024c80c683713497df268c6b3b7ff7/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/61024c80c683713497df268c6b3b7ff7/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■                55% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/50263076f53ef20a957dda489eadd56a/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/50263076f53ef20a957dda489eadd56a/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■                55% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■■                56% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/50263076f53ef20a957dda489eadd56a/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/50263076f53ef20a957dda489eadd56a/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■                56% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■■                57% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/c9c75d7d75926c99eec6b4d505d04d29/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/c9c75d7d75926c99eec6b4d505d04d29/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■                57% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■■                57% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/c9c75d7d75926c99eec6b4d505d04d29/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/c9c75d7d75926c99eec6b4d505d04d29/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■                57% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/3e26ff7794eda8e798fa1308f1c6b388/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/3e26ff7794eda8e798fa1308f1c6b388/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■                57% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■■                57% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/3e26ff7794eda8e798fa1308f1c6b388/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/3e26ff7794eda8e798fa1308f1c6b388/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■                57% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/61d4c23535cb7cf9bddac49bc31934ad/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/61d4c23535cb7cf9bddac49bc31934ad/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■                57% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■■                58% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/61d4c23535cb7cf9bddac49bc31934ad/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/61d4c23535cb7cf9bddac49bc31934ad/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■                58% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/83162c9b53107e2e565112e57ade39bc/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/83162c9b53107e2e565112e57ade39bc/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■                58% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■■■               58% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/83162c9b53107e2e565112e57ade39bc/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/83162c9b53107e2e565112e57ade39bc/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■               58% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/d57c40c40dd8b6259a200c4c09b562de/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/d57c40c40dd8b6259a200c4c09b562de/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■               58% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■■■               59% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/d57c40c40dd8b6259a200c4c09b562de/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/d57c40c40dd8b6259a200c4c09b562de/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■               59% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■■■               59% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e44c837618eae5dcb44bda61dcb68d3e/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e44c837618eae5dcb44bda61dcb68d3e/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■               59% |  ETA:  3m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e44c837618eae5dcb44bda61dcb68d3e/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e44c837618eae5dcb44bda61dcb68d3e/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■               59% |  ETA:  3mDownloading files ■■■■■■■■■■■■■■■■■■■               60% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/96d09aa430de8a97140313504ae84988/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/96d09aa430de8a97140313504ae84988/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■               60% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/96d09aa430de8a97140313504ae84988/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/96d09aa430de8a97140313504ae84988/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■               60% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a85e40c7755e60d4ff2eb64e68388b4d/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a85e40c7755e60d4ff2eb64e68388b4d/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■               60% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■               61% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a85e40c7755e60d4ff2eb64e68388b4d/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a85e40c7755e60d4ff2eb64e68388b4d/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■               61% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/2609cd3d6deabfd24d250e8b6e394493/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/2609cd3d6deabfd24d250e8b6e394493/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■               61% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■               61% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/2609cd3d6deabfd24d250e8b6e394493/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/2609cd3d6deabfd24d250e8b6e394493/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■               61% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/5988c2dcae795773471174d653569d53/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/5988c2dcae795773471174d653569d53/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■               61% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■              62% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/5988c2dcae795773471174d653569d53/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/5988c2dcae795773471174d653569d53/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■              62% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/3c092f8812fc95ff93edaf8c3e27dbc5/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/3c092f8812fc95ff93edaf8c3e27dbc5/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■              62% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■              62% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/3c092f8812fc95ff93edaf8c3e27dbc5/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/3c092f8812fc95ff93edaf8c3e27dbc5/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■              62% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9a88b76aaac24efdf680a715ff16a5e9/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9a88b76aaac24efdf680a715ff16a5e9/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■              62% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■              63% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9a88b76aaac24efdf680a715ff16a5e9/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9a88b76aaac24efdf680a715ff16a5e9/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■              63% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/4a4f73fc70640ee2c916c31612a9f40d/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/4a4f73fc70640ee2c916c31612a9f40d/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■              63% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■              63% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/4a4f73fc70640ee2c916c31612a9f40d/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/4a4f73fc70640ee2c916c31612a9f40d/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■              63% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/b57f38cf5a7eff68842e1680564c4b29/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/b57f38cf5a7eff68842e1680564c4b29/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■              63% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■              64% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/b57f38cf5a7eff68842e1680564c4b29/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/b57f38cf5a7eff68842e1680564c4b29/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■              64% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/bc6484e9b3b6b6e3531775ebd0c05c98/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/bc6484e9b3b6b6e3531775ebd0c05c98/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■              64% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■              64% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/bc6484e9b3b6b6e3531775ebd0c05c98/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/bc6484e9b3b6b6e3531775ebd0c05c98/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■              64% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/b8f881f6d39c1a975f00449dced89283/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/b8f881f6d39c1a975f00449dced89283/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■              64% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■              65% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/b8f881f6d39c1a975f00449dced89283/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/b8f881f6d39c1a975f00449dced89283/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■              65% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/33598e6f0129327cdc696363884d6e58/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/33598e6f0129327cdc696363884d6e58/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■              65% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■■             65% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/33598e6f0129327cdc696363884d6e58/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/33598e6f0129327cdc696363884d6e58/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■             65% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/6b8aef89e3919a3f9f0012120e3230ac/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/6b8aef89e3919a3f9f0012120e3230ac/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■             65% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■■             66% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/6b8aef89e3919a3f9f0012120e3230ac/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/6b8aef89e3919a3f9f0012120e3230ac/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■             66% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/47950ec1bcd99ac93c721761d00d625e/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/47950ec1bcd99ac93c721761d00d625e/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■             66% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■■             67% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/47950ec1bcd99ac93c721761d00d625e/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/47950ec1bcd99ac93c721761d00d625e/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■             67% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/183aa2cec358635423b794ae5602d1b8/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/183aa2cec358635423b794ae5602d1b8/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■             67% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■■             67% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/183aa2cec358635423b794ae5602d1b8/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/183aa2cec358635423b794ae5602d1b8/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■             67% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ea24cb90b425c32ebe0f50e34fd10204/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ea24cb90b425c32ebe0f50e34fd10204/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■             67% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■■             68% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ea24cb90b425c32ebe0f50e34fd10204/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ea24cb90b425c32ebe0f50e34fd10204/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■             68% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9b8223a9ae8b96a374757cba511992cc/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9b8223a9ae8b96a374757cba511992cc/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■             68% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■■             68% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9b8223a9ae8b96a374757cba511992cc/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9b8223a9ae8b96a374757cba511992cc/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■             68% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■■■            68% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a73e4022cab664b71fe8b8414187d07b/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a73e4022cab664b71fe8b8414187d07b/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■            68% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■■■            69% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a73e4022cab664b71fe8b8414187d07b/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a73e4022cab664b71fe8b8414187d07b/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■            69% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ffdb2056e59118e09bafc37205d55ebd/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ffdb2056e59118e09bafc37205d55ebd/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■            69% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■■■            69% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ffdb2056e59118e09bafc37205d55ebd/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ffdb2056e59118e09bafc37205d55ebd/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■            69% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■■■            70% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/209c8d536474811fa0d75df88fe1701d/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/209c8d536474811fa0d75df88fe1701d/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■            70% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■■■            70% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/209c8d536474811fa0d75df88fe1701d/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/209c8d536474811fa0d75df88fe1701d/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■            70% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e6a8fee32a99139daedbf58cea5d16a4/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e6a8fee32a99139daedbf58cea5d16a4/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■            70% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■■■            70% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e6a8fee32a99139daedbf58cea5d16a4/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e6a8fee32a99139daedbf58cea5d16a4/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■            70% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/b10c64c6fe5054819a216ca0a9ae71e2/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/b10c64c6fe5054819a216ca0a9ae71e2/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■            70% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■■■            71% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/b10c64c6fe5054819a216ca0a9ae71e2/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/b10c64c6fe5054819a216ca0a9ae71e2/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■            71% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/89212f7797b56548e49a25d3f47613be/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/89212f7797b56548e49a25d3f47613be/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■            71% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■■■            71% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/89212f7797b56548e49a25d3f47613be/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/89212f7797b56548e49a25d3f47613be/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■            71% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/d2af6e2a113627d9f64b84cca306e79c/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/d2af6e2a113627d9f64b84cca306e79c/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■            71% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/d2af6e2a113627d9f64b84cca306e79c/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/d2af6e2a113627d9f64b84cca306e79c/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■            71% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■           72% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/7f1a897455cc0270f5b17870fc999eb8/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/7f1a897455cc0270f5b17870fc999eb8/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■           72% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/7f1a897455cc0270f5b17870fc999eb8/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/7f1a897455cc0270f5b17870fc999eb8/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■           72% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/98cdc2d51ea9146190ad8b20d8a12be3/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/98cdc2d51ea9146190ad8b20d8a12be3/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■           72% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■           73% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/98cdc2d51ea9146190ad8b20d8a12be3/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/98cdc2d51ea9146190ad8b20d8a12be3/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■           73% |  ETA:  2m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/51fd154d6a251b2a4c0c880159e6a81e/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/51fd154d6a251b2a4c0c880159e6a81e/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■           73% |  ETA:  2mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■           74% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/51fd154d6a251b2a4c0c880159e6a81e/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/51fd154d6a251b2a4c0c880159e6a81e/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■           74% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e427efe71e8e94de5b3e48eb98236323/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e427efe71e8e94de5b3e48eb98236323/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■           74% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■           74% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e427efe71e8e94de5b3e48eb98236323/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e427efe71e8e94de5b3e48eb98236323/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■           74% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■           74% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/0636d21d01e58665233956672543a0d9/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/0636d21d01e58665233956672543a0d9/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■           74% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■           75% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/0636d21d01e58665233956672543a0d9/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/0636d21d01e58665233956672543a0d9/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■           75% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■           75% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a38594e1bb4b2131f31ae1ea8162cd6c/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a38594e1bb4b2131f31ae1ea8162cd6c/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■           75% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■          75% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a38594e1bb4b2131f31ae1ea8162cd6c/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a38594e1bb4b2131f31ae1ea8162cd6c/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■          75% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■          76% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/639778bf6a165eb46f3f725f4c4c911a/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/639778bf6a165eb46f3f725f4c4c911a/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■          76% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■          76% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/639778bf6a165eb46f3f725f4c4c911a/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/639778bf6a165eb46f3f725f4c4c911a/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■          76% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a6f4db2a37e02adbfbad553b6cc62be8/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a6f4db2a37e02adbfbad553b6cc62be8/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■          76% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/a6f4db2a37e02adbfbad553b6cc62be8/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/a6f4db2a37e02adbfbad553b6cc62be8/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■          76% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■          77% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/450d1967561f1ace2695bb7b1bfa1581/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/450d1967561f1ace2695bb7b1bfa1581/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■          77% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■          77% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/450d1967561f1ace2695bb7b1bfa1581/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/450d1967561f1ace2695bb7b1bfa1581/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■          77% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■          77% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/d1d53a6bf2a4cf2cbda959990c017b1c/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/d1d53a6bf2a4cf2cbda959990c017b1c/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■          77% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■          77% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/d1d53a6bf2a4cf2cbda959990c017b1c/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/d1d53a6bf2a4cf2cbda959990c017b1c/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■          77% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■          78% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/101fc095c54909bd45ffb7ff244cf04e/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/101fc095c54909bd45ffb7ff244cf04e/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■          78% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■          78% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/101fc095c54909bd45ffb7ff244cf04e/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/101fc095c54909bd45ffb7ff244cf04e/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■          78% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/36ce88a4452613c5b0f9aee7a529a963/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/36ce88a4452613c5b0f9aee7a529a963/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■          78% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         79% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/36ce88a4452613c5b0f9aee7a529a963/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/36ce88a4452613c5b0f9aee7a529a963/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         79% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         79% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/7e439918352084138c7ce8107e38bddb/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/7e439918352084138c7ce8107e38bddb/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         79% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         79% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/7e439918352084138c7ce8107e38bddb/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/7e439918352084138c7ce8107e38bddb/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         79% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/5093764246714f95d68b26f467bca198/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/5093764246714f95d68b26f467bca198/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         79% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         80% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/5093764246714f95d68b26f467bca198/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/5093764246714f95d68b26f467bca198/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         80% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/da20483e9528153eeaa421dfcc5ef25d/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/da20483e9528153eeaa421dfcc5ef25d/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         80% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         80% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/da20483e9528153eeaa421dfcc5ef25d/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/da20483e9528153eeaa421dfcc5ef25d/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         80% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e02b3d23fe966909c1bea99c583abea8/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e02b3d23fe966909c1bea99c583abea8/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         80% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         81% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e02b3d23fe966909c1bea99c583abea8/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e02b3d23fe966909c1bea99c583abea8/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         81% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/fed91cfdfca482dfddc57060a4417fb7/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/fed91cfdfca482dfddc57060a4417fb7/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         81% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         81% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/fed91cfdfca482dfddc57060a4417fb7/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/fed91cfdfca482dfddc57060a4417fb7/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         81% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/aeb1fe87751f198f669012a1865936ed/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/aeb1fe87751f198f669012a1865936ed/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         81% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/aeb1fe87751f198f669012a1865936ed/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/aeb1fe87751f198f669012a1865936ed/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         81% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/1d971b8c083349e082be4cc40dfdf115/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/1d971b8c083349e082be4cc40dfdf115/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■         81% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■        82% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/1d971b8c083349e082be4cc40dfdf115/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/1d971b8c083349e082be4cc40dfdf115/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■        82% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/66e4da7aba4b59873f34e5fdc6f0002c/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/66e4da7aba4b59873f34e5fdc6f0002c/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■        82% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■        83% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/66e4da7aba4b59873f34e5fdc6f0002c/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/66e4da7aba4b59873f34e5fdc6f0002c/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■        83% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/eaed3c31d1ec0e539f7c24e8addd1edf/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/eaed3c31d1ec0e539f7c24e8addd1edf/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■        83% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■        83% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/eaed3c31d1ec0e539f7c24e8addd1edf/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/eaed3c31d1ec0e539f7c24e8addd1edf/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■        83% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9f7d4ca403d0313c5b4e55b26bc982e6/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9f7d4ca403d0313c5b4e55b26bc982e6/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■        83% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■        84% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9f7d4ca403d0313c5b4e55b26bc982e6/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9f7d4ca403d0313c5b4e55b26bc982e6/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■        84% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/4983467a9a14ebc24176035d99e57732/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/4983467a9a14ebc24176035d99e57732/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■        84% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■        85% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/4983467a9a14ebc24176035d99e57732/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/4983467a9a14ebc24176035d99e57732/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■        85% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/07a2104348139589da4e7611f6707dd9/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/07a2104348139589da4e7611f6707dd9/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■        85% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/07a2104348139589da4e7611f6707dd9/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/07a2104348139589da4e7611f6707dd9/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■        85% |  ETA:  1m                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/1ce3d511f3dabd31c132944aca053956/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/1ce3d511f3dabd31c132944aca053956/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■        85% |  ETA:  1mDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■       86% |  ETA: 47s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/1ce3d511f3dabd31c132944aca053956/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/1ce3d511f3dabd31c132944aca053956/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■       86% |  ETA: 47s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/91451df1d1c45dddba49fdc7b2b06b2b/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/91451df1d1c45dddba49fdc7b2b06b2b/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■       86% |  ETA: 47sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■       86% |  ETA: 45s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/91451df1d1c45dddba49fdc7b2b06b2b/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/91451df1d1c45dddba49fdc7b2b06b2b/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■       86% |  ETA: 45s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/524e25975ab47e06bfa2795fb91055f6/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/524e25975ab47e06bfa2795fb91055f6/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■       86% |  ETA: 45sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■       87% |  ETA: 43s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/524e25975ab47e06bfa2795fb91055f6/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/524e25975ab47e06bfa2795fb91055f6/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■       87% |  ETA: 43s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/84d342e7cdcbaf7d1f5943a1161b3214/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/84d342e7cdcbaf7d1f5943a1161b3214/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■       87% |  ETA: 43sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■       87% |  ETA: 41s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/84d342e7cdcbaf7d1f5943a1161b3214/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/84d342e7cdcbaf7d1f5943a1161b3214/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■       87% |  ETA: 41s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/b8a12e588ab2159d549ab1ef8c63899b/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/b8a12e588ab2159d549ab1ef8c63899b/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■       87% |  ETA: 41s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/b8a12e588ab2159d549ab1ef8c63899b/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/b8a12e588ab2159d549ab1ef8c63899b/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■       87% |  ETA: 41sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■       88% |  ETA: 38s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e43a27de346b8c38b51c1d31587c9db9/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e43a27de346b8c38b51c1d31587c9db9/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■       88% |  ETA: 38s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e43a27de346b8c38b51c1d31587c9db9/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e43a27de346b8c38b51c1d31587c9db9/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■       88% |  ETA: 38s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e6d69bdd53917b8e12fc77995195ed46/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e6d69bdd53917b8e12fc77995195ed46/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■       88% |  ETA: 38sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      89% |  ETA: 35s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/e6d69bdd53917b8e12fc77995195ed46/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/e6d69bdd53917b8e12fc77995195ed46/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      89% |  ETA: 35s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/377fdf24566c23c4134e632e08622647/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/377fdf24566c23c4134e632e08622647/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      89% |  ETA: 35sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      89% |  ETA: 34s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/377fdf24566c23c4134e632e08622647/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/377fdf24566c23c4134e632e08622647/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      89% |  ETA: 34s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9e2c405b49db89a05a56654ed5ef0c54/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9e2c405b49db89a05a56654ed5ef0c54/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      89% |  ETA: 34sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      90% |  ETA: 32s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9e2c405b49db89a05a56654ed5ef0c54/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9e2c405b49db89a05a56654ed5ef0c54/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      90% |  ETA: 32sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      90% |  ETA: 31s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/f54cc3f57867cde8dad85b3943caf40c/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/f54cc3f57867cde8dad85b3943caf40c/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      90% |  ETA: 31sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      90% |  ETA: 30s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/f54cc3f57867cde8dad85b3943caf40c/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/f54cc3f57867cde8dad85b3943caf40c/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      90% |  ETA: 30s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/652decf843229baad3756757d725b8ce/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/652decf843229baad3756757d725b8ce/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      90% |  ETA: 30s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/652decf843229baad3756757d725b8ce/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/652decf843229baad3756757d725b8ce/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      90% |  ETA: 30sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      91% |  ETA: 27s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/2f4c65a10efceb9c740f2530d2dac22d/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/2f4c65a10efceb9c740f2530d2dac22d/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      91% |  ETA: 27sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      92% |  ETA: 26s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/2f4c65a10efceb9c740f2530d2dac22d/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/2f4c65a10efceb9c740f2530d2dac22d/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      92% |  ETA: 26s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ce990c823897f594c2f67020ce4528cf/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ce990c823897f594c2f67020ce4528cf/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      92% |  ETA: 26sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     92% |  ETA: 25s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ce990c823897f594c2f67020ce4528cf/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ce990c823897f594c2f67020ce4528cf/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     92% |  ETA: 25sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     92% |  ETA: 24s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/476a0495f54da0364300288f83faf945/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/476a0495f54da0364300288f83faf945/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     92% |  ETA: 24sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     93% |  ETA: 23s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/476a0495f54da0364300288f83faf945/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/476a0495f54da0364300288f83faf945/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     93% |  ETA: 23s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/f32e31925a088673e71bfd625c260ad0/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/f32e31925a088673e71bfd625c260ad0/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     93% |  ETA: 23sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     93% |  ETA: 21s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/f32e31925a088673e71bfd625c260ad0/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/f32e31925a088673e71bfd625c260ad0/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     93% |  ETA: 21sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     93% |  ETA: 20s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/5ab6a562a16d683284ee060cc4d13b08/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/5ab6a562a16d683284ee060cc4d13b08/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     93% |  ETA: 20sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     94% |  ETA: 20s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/5ab6a562a16d683284ee060cc4d13b08/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/5ab6a562a16d683284ee060cc4d13b08/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     94% |  ETA: 20s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/aa53a0672ea59137a5b3e9d70fc481e2/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/aa53a0672ea59137a5b3e9d70fc481e2/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     94% |  ETA: 20sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     94% |  ETA: 18s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/aa53a0672ea59137a5b3e9d70fc481e2/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/aa53a0672ea59137a5b3e9d70fc481e2/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     94% |  ETA: 18s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ebe20a6b052d077401bb1cb8bc5ad33f/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ebe20a6b052d077401bb1cb8bc5ad33f/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     94% |  ETA: 18sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     95% |  ETA: 16s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ebe20a6b052d077401bb1cb8bc5ad33f/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ebe20a6b052d077401bb1cb8bc5ad33f/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     95% |  ETA: 16sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    95% |  ETA: 15s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/efcffe370e99dcbdfa56b951916e7b1c/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/efcffe370e99dcbdfa56b951916e7b1c/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    95% |  ETA: 15s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/efcffe370e99dcbdfa56b951916e7b1c/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/efcffe370e99dcbdfa56b951916e7b1c/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    95% |  ETA: 15sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    96% |  ETA: 13s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9389f1e414a052523882287d2b3a9a44/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9389f1e414a052523882287d2b3a9a44/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    96% |  ETA: 13sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    96% |  ETA: 13s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9389f1e414a052523882287d2b3a9a44/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9389f1e414a052523882287d2b3a9a44/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    96% |  ETA: 13s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9de6fe7873eba1cb6d43ad64c6c8b9d5/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9de6fe7873eba1cb6d43ad64c6c8b9d5/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    96% |  ETA: 13sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    96% |  ETA: 11s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9de6fe7873eba1cb6d43ad64c6c8b9d5/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9de6fe7873eba1cb6d43ad64c6c8b9d5/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    96% |  ETA: 11s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/feca086e987d6a6a039afce899d659f4/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/feca086e987d6a6a039afce899d659f4/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    96% |  ETA: 11sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    97% |  ETA:  9s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/feca086e987d6a6a039afce899d659f4/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/feca086e987d6a6a039afce899d659f4/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    97% |  ETA:  9s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9eabc7393d6c2531dbb16717dbb63348/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9eabc7393d6c2531dbb16717dbb63348/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    97% |  ETA:  9sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    98% |  ETA:  7s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/9eabc7393d6c2531dbb16717dbb63348/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/9eabc7393d6c2531dbb16717dbb63348/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    98% |  ETA:  7s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/8e4ca69983ff6bd0c41f23f343093832/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/8e4ca69983ff6bd0c41f23f343093832/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    98% |  ETA:  7sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    98% |  ETA:  6s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/8e4ca69983ff6bd0c41f23f343093832/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/8e4ca69983ff6bd0c41f23f343093832/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    98% |  ETA:  6s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ee04286fca57bfb1351183683353bcc4/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ee04286fca57bfb1351183683353bcc4/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    98% |  ETA:  6sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   99% |  ETA:  4s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/ee04286fca57bfb1351183683353bcc4/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/ee04286fca57bfb1351183683353bcc4/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   99% |  ETA:  4s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/cafae7395ca3e28058321f628f3fd188/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/cafae7395ca3e28058321f628f3fd188/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   99% |  ETA:  4sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   99% |  ETA:  2s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/cafae7395ca3e28058321f628f3fd188/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/cafae7395ca3e28058321f628f3fd188/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   99% |  ETA:  2s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/210419160b0e0cdcd61b162e9500e78b/assays.h5 to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/210419160b0e0cdcd61b162e9500e78b/assays.h5
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   99% |  ETA:  2sDownloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% |  ETA:  1s                                                                    ℹ Downloading https://swift.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/harmonised-human-atlas/cpm/210419160b0e0cdcd61b162e9500e78b/se.rds to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/cpm/210419160b0e0cdcd61b162e9500e78b/se.rds
#> Downloading files ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% |  ETA:  1s                                                                    ℹ Reading files.
#> ℹ Compiling Single Cell Experiment.
```

![](README_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

<img src="man/figures/HLA_A_tissue_plot.png" width="525" />

## Obtain Unharmonised Metadata

Various metadata fields are *not* common between datasets, so it does
not make sense for these to live in the main metadata table. However, we
can obtain it using the `get_unharmonised_metadata()` function. This
function returns a data frame with one row per dataset, including the
`unharmonised` column which contains unharmnised metadata as a nested
data frame.

``` r
harmonised <- get_metadata() |> dplyr::filter(tissue == "kidney blood vessel")
unharmonised <- get_unharmonised_metadata(harmonised)
#> ℹ Downloading https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/unharmonised_metadata/63523aa3-0d04-4fc6-ac59-5cadd3e73a14.parquet to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/unharmonised/63523aa3-0d04-4fc6-ac59-5cadd3e73a14.parquet
#> ℹ Downloading https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/unharmonised_metadata/8fee7b82-178b-4c04-bf23-04689415690d.parquet to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/unharmonised/8fee7b82-178b-4c04-bf23-04689415690d.parquet
#> ℹ Downloading https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/unharmonised_metadata/dc9d8cdd-29ee-4c44-830c-6559cb3d0af6.parquet to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/unharmonised/dc9d8cdd-29ee-4c44-830c-6559cb3d0af6.parquet
#> ℹ Downloading https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/unharmonised_metadata/f7e94dbb-8638-4616-aaf9-16e2212c369f.parquet to /vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/0.2/unharmonised/f7e94dbb-8638-4616-aaf9-16e2212c369f.parquet
unharmonised
#> # A tibble: 4 × 2
#>   file_id                              unharmonised   
#>   <chr>                                <list>         
#> 1 63523aa3-0d04-4fc6-ac59-5cadd3e73a14 <tbl_dck_[,17]>
#> 2 8fee7b82-178b-4c04-bf23-04689415690d <tbl_dck_[,12]>
#> 3 dc9d8cdd-29ee-4c44-830c-6559cb3d0af6 <tbl_dck_[,14]>
#> 4 f7e94dbb-8638-4616-aaf9-16e2212c369f <tbl_dck_[,14]>
```

Notice that the columns differ between each dataset’s data frame:

``` r
dplyr::pull(unharmonised, unharmonised) |> head(2)
#> [[1]]
#> # Source:   SQL [?? x 17]
#> # Database: DuckDB 0.6.2-dev1166 [unknown@Linux 3.10.0-1160.88.1.el7.x86_64:R 4.2.1/:memory:]
#>    cell_ file_id donor…¹ donor…² libra…³ mappe…⁴ sampl…⁵ suspe…⁶ suspe…⁷ autho…⁸
#>    <chr> <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>  
#>  1 4602… 63523a… 27 mon… a8536b… 5ddaea… GENCOD… 61bf84… cell    d8a44f… Pelvic…
#>  2 4602… 63523a… 27 mon… a8536b… 5ddaea… GENCOD… 61bf84… cell    d8a44f… Pelvic…
#>  3 4602… 63523a… 27 mon… a8536b… 5ddaea… GENCOD… 61bf84… cell    d8a44f… Pelvic…
#>  4 4602… 63523a… 27 mon… a8536b… 5ddaea… GENCOD… 61bf84… cell    d8a44f… Pelvic…
#>  5 4602… 63523a… 27 mon… a8536b… 5ddaea… GENCOD… 61bf84… cell    d8a44f… Pelvic…
#>  6 4602… 63523a… 27 mon… a8536b… 5ddaea… GENCOD… 61bf84… cell    d8a44f… Pelvic…
#>  7 4602… 63523a… 27 mon… a8536b… 5ddaea… GENCOD… 61bf84… cell    d8a44f… Pelvic…
#>  8 4602… 63523a… 27 mon… a8536b… 5ddaea… GENCOD… 61bf84… cell    d8a44f… Pelvic…
#>  9 4602… 63523a… 19 mon… 463181… 671785… GENCOD… 125234… cell    c7485e… CD4 T …
#> 10 4602… 63523a… 27 mon… a8536b… 5ddaea… GENCOD… 61bf84… cell    d8a44f… Pelvic…
#> # … with more rows, 7 more variables: cell_state <chr>,
#> #   reported_diseases <chr>, Short_Sample <chr>, Project <chr>,
#> #   Experiment <chr>, compartment <chr>, broad_celltype <chr>, and abbreviated
#> #   variable names ¹​donor_age, ²​donor_uuid, ³​library_uuid,
#> #   ⁴​mapped_reference_annotation, ⁵​sample_uuid, ⁶​suspension_type,
#> #   ⁷​suspension_uuid, ⁸​author_cell_type
#> 
#> [[2]]
#> # Source:   SQL [?? x 12]
#> # Database: DuckDB 0.6.2-dev1166 [unknown@Linux 3.10.0-1160.88.1.el7.x86_64:R 4.2.1/:memory:]
#>    cell_ file_id orig.…¹ nCoun…² nFeat…³ seura…⁴ Project donor…⁵ compa…⁶ broad…⁷
#>    <chr> <chr>   <chr>     <dbl> <chr>   <chr>   <chr>   <chr>   <chr>   <chr>  
#>  1 1069  8fee7b… 4602ST…   16082 3997    25      Experi… Wilms3  non_PT  Pelvic…
#>  2 1214  8fee7b… 4602ST…    1037 606     25      Experi… Wilms3  non_PT  Pelvic…
#>  3 2583  8fee7b… 4602ST…    3028 1361    25      Experi… Wilms3  non_PT  Pelvic…
#>  4 2655  8fee7b… 4602ST…    1605 859     25      Experi… Wilms3  non_PT  Pelvic…
#>  5 3609  8fee7b… 4602ST…    1144 682     25      Experi… Wilms3  non_PT  Pelvic…
#>  6 3624  8fee7b… 4602ST…    1874 963     25      Experi… Wilms3  non_PT  Pelvic…
#>  7 3946  8fee7b… 4602ST…    1296 755     25      Experi… Wilms3  non_PT  Pelvic…
#>  8 5163  8fee7b… 4602ST…   11417 3255    25      Experi… Wilms3  non_PT  Pelvic…
#>  9 5446  8fee7b… 4602ST…    1769 946     19      Experi… Wilms2  lympho… CD4 T …
#> 10 6275  8fee7b… 4602ST…    3750 1559    25      Experi… Wilms3  non_PT  Pelvic…
#> # … with more rows, 2 more variables: author_cell_type <chr>, Sample <chr>, and
#> #   abbreviated variable names ¹​orig.ident, ²​nCount_RNA, ³​nFeature_RNA,
#> #   ⁴​seurat_clusters, ⁵​donor_id, ⁶​compartment, ⁷​broad_celltype
```

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

`sample_`, `sample_name`, `age_days`, `assay`, `assay_ontology_term_id`,
`development_stage`, `development_stage_ontology_term_id`, `ethnicity`,
`ethnicity_ontology_term_id`, `experiment___`, `organism`,
`organism_ontology_term_id`, `sample_placeholder`, `sex`,
`sex_ontology_term_id`, `tissue`, `tissue_harmonised`,
`tissue_ontology_term_id`, `disease`, `disease_ontology_term_id`,
`is_primary_data.x`

Cell-specific columns (definitions available at
cellxgene.cziscience.com)

`cell_`, `cell_type`, `cell_type_ontology_term_idm`,
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
- `sample_`: Sample ID
- `sample_name`: How samples were defined

# RNA abundance

The `raw` assay includes RNA abundance in the positive real scale (not
transformed with non-linear functions, e.g. log sqrt). Originally
CELLxGENE include a mix of scales and transformations specified in the
`x_normalization` column.

The `cpm` assay includes counts per million.

# Installation and getting-started problems

**Problem:** Default R cache path including non-standard characters
(e.g. dash)

``` r
get_metadata()

# Error in `db_query_fields.DBIConnection()`:
# ! Can't query fields.
# Caused by error:
# ! Parser Error: syntax error at or near "/"
# LINE 2: FROM /Users/bob/Library/Cach...
```

**Solution:** Setup custom cache path (e.g. user home directory)

``` r
get_metadata(cache_directory = path.expand('~'))
```

**Problem:** namespace ‘dbplyr’ 2.2.1 is being loaded, but \>= 2.3.0 is
required

**Solution:** Install new dbplyr

``` r
install.packages("dbplyr")
```

------------------------------------------------------------------------

This project has been funded by

- *Silicon Valley Foundation* CZF2019-002443
- *Bioconductor core funding* NIH NHGRI 5U24HG004059-18
- *Victoria Cancer Agency* ECRF21036
- *Australian National Health and Medical Research Council* 1116955
- *The Lorenzo and Pamela Galli Medical Research Trust*
