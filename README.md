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

<img src="man/figures/svcf_logo.jpeg" width="155x" height="58px" /><img src="man/figures/czi_logo.png" width="129px" height="58px" /><img src="man/figures/bioconductor_logo.jpg" width="202px" height="58px" /><img src="man/figures/vca_logo.png" width="219px" height="58px" /><img src="man/figures/nectar_logo.png" width="180px" height="58px" />

[website](https://stemangiola.github.io/CuratedAtlasQueryR)

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
metadata  = get_metadata()

metadata
#> # Source:   table</vast/scratch/users/milton.m/cache/R/CuratedAtlasQueryR/metadata.0.2.3.parquet> [?? x 56]
#> # Database: DuckDB 0.6.2-dev1166 [unknown@Linux 3.10.0-1160.81.1.el7.x86_64:R 4.2.1/:memory:]
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
#> # Database: DuckDB 0.6.2-dev1166 [unknown@Linux 3.10.0-1160.81.1.el7.x86_64:R 4.2.1/:memory:]
#>    tissue                                           n
#>    <chr>                                        <dbl>
#>  1 blood                                           47
#>  2 respiratory airway                              16
#>  3 mammary gland epithelial cell (cell culture)     1
#>  4 colon                                            3
#>  5 intestine                                       18
#>  6 pleural effusion                                11
#>  7 lymph node                                      15
#>  8 lung                                            27
#>  9 liver                                           24
#> 10 axilla                                          10
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
    get_SingleCellExperiment()
#> ℹ Realising metadata.
#> ℹ Synchronising files
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
    get_SingleCellExperiment(assays = "cpm")
#> ℹ Realising metadata.
#> ℹ Synchronising files
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
objects, which merge a lot of HDF5 from your `get_SingleCellExperiment`
the display and manipulation is going to be slow.

``` r
single_cell_counts |> saveRDS("single_cell_counts.rds")
```

### Saving as HDF5 (slow saving, fast reading)

Saving as `.rds` has the advantage that rewrites on disk a monolithic
`HDF5` and so displaying and manipulating large `SingleCellExperiment`
objects, which merge a lot of HDF5 from your `get_SingleCellExperiment`,
is going to be fast.

However it has the disadvantage that the files are going to be larger as
they include the count information, and the saving process is going to
be slow for large objects.

``` r
single_cell_counts |> saveHDF5SummarizedExperiment("single_cell_counts")
```

## Visualise gene transcription

We can gather all CD14 monocytes cells and plot the distribution of
HLA-A across all tissues

``` r
library(tidySingleCellExperiment)
library(ggplot2)

metadata |>
  # Filter and subset
  filter(cell_type_harmonised=="cd14 mono") |>

  # Get counts per million for HCA-A gene
  get_SingleCellExperiment(assays = "cpm", features = "HLA-A") |> 
  
  # Plot (styling code have been omitted)
  join_features("HLA-A", shape = "wide") |> 
  ggplot(aes( disease, `HLA.A`,color = file_id)) +
  geom_jitter(shape=".") 
```

<img src="man/figures/HLA_A_disease_plot.png" width="525" />

``` r

metadata |> 
    
  # Filter and subset
  filter(cell_type_harmonised=="nk") |> 

  # Get counts per million for HCA-A gene 
  get_SingleCellExperiment(assays = "cpm", features = "HLA-A") |> 

    # Plot (styling code have been omitted)
  join_features("HLA-A", shape = "wide") |> 
  ggplot(aes( tissue_harmonised, `HLA.A`,color = file_id)) +
  geom_jitter(shape=".") 
```

<img src="man/figures/HLA_A_tissue_plot.png" width="525" />

## Obtain Unharmonised Metadata

Various metadata fields are *not* common between datasets, so it does
not make sense for these to live in the main metadata table. However, we
can obtain it using the `get_unharmonised_metadata()` function.

Note how this table has additional columns that are not in the normal
metadata:

``` r
dataset = "838ea006-2369-4e2c-b426-b2a744a2b02b"
unharmonised_meta = get_unharmonised_metadata(dataset)
unharmonised_tbl = dplyr::collect(unharmonised_meta[[dataset]])
unharmonised_tbl
#> # A tibble: 168,860 × 23
#>    cell_     file_id Neuro…¹ Class Subcl…² Super…³ Age.a…⁴ Years…⁵ Cogni…⁶ ADNC 
#>    <chr>     <chr>   <lgl>   <chr> <chr>   <chr>   <chr>   <chr>   <chr>   <chr>
#>  1 GGACGAAG… 838ea0… FALSE   Neur… L4 IT   L4 IT_2 90+ ye… 16 to … Dement… High 
#>  2 TCACGGGA… 838ea0… FALSE   Neur… L4 IT   L4 IT_1 90+ ye… 12 to … Dement… Inte…
#>  3 TCAGTTTT… 838ea0… FALSE   Neur… L4 IT   L4 IT_2 78 to … 16 to … No dem… Low  
#>  4 TCAGTCCT… 838ea0… FALSE   Neur… L4 IT   L4 IT_4 78 to … 16 to … Dement… Inte…
#>  5 AGCCACGC… 838ea0… FALSE   Neur… L4 IT   L4 IT_2 78 to … 19 to … No dem… Inte…
#>  6 CCTCAACC… 838ea0… TRUE    Neur… L4 IT   L4 IT_2 Less t… Refere… Refere… Refe…
#>  7 CTCGACAA… 838ea0… FALSE   Neur… L4 IT   L4 IT_2 78 to … 12 to … No dem… Inte…
#>  8 AGCTACAG… 838ea0… FALSE   Neur… L4 IT   L4 IT_4 90+ ye… 16 to … Dement… High 
#>  9 CTCGAGGG… 838ea0… FALSE   Neur… L4 IT   L4 IT_2 65 to … 16 to … Dement… High 
#> 10 AGTGCCGT… 838ea0… FALSE   Neur… L4 IT   L4 IT_4 90+ ye… 16 to … Dement… High 
#> # … with 168,850 more rows, 13 more variables: Braak.stage <chr>,
#> #   Thal.phase <chr>, CERAD.score <chr>, APOE4.status <chr>,
#> #   Lewy.body.disease.pathology <chr>, LATE.NC.stage <chr>,
#> #   Microinfarct.pathology <chr>, Specimen.ID <chr>, Donor.ID <chr>, PMI <chr>,
#> #   Number.of.UMIs <dbl>, Genes.detected <dbl>,
#> #   Fraction.mitochrondrial.UMIs <dbl>, and abbreviated variable names
#> #   ¹​Neurotypical.reference, ²​Subclass, ³​Supertype, ⁴​Age.at.death, …
```

If we have metadata from the normal metadata table that is from a single
dataset, we can even join this additional metadata into one big data
frame:

``` r
harmonised_meta = get_metadata() |> dplyr::filter(file_id == dataset) |> dplyr::collect()
dplyr::left_join(harmonised_meta, unharmonised_tbl, by=c("file_id", "cell_"))
#> # A tibble: 168,860 × 77
#>    cell_ sample_ cell_…¹ cell_…² confi…³ cell_…⁴ cell_…⁵ cell_…⁶ sampl…⁷ _samp…⁸
#>    <chr> <chr>   <chr>   <chr>     <dbl> <chr>   <chr>   <chr>   <chr>   <chr>  
#>  1 GGAC… f63cb4… L2/3-6… neuron        1 <NA>    <NA>    <NA>    168593… H21.33…
#>  2 TCAC… 0d4d1f… L2/3-6… neuron        1 <NA>    <NA>    <NA>    f7d747… H21.33…
#>  3 TCAG… 3e5a3b… L2/3-6… neuron        1 <NA>    <NA>    <NA>    3417a9… H20.33…
#>  4 TCAG… 7010a3… L2/3-6… neuron        1 <NA>    <NA>    <NA>    246a59… H20.33…
#>  5 AGCC… 82bb9a… L2/3-6… neuron        1 <NA>    <NA>    <NA>    7a8f35… H21.33…
#>  6 CCTC… a233eb… L2/3-6… neuron        1 <NA>    <NA>    <NA>    188243… H18.30…
#>  7 CTCG… 27f104… L2/3-6… neuron        1 <NA>    <NA>    <NA>    a62943… H20.33…
#>  8 AGCT… 0190a2… L2/3-6… neuron        1 <NA>    <NA>    <NA>    c508a8… H20.33…
#>  9 CTCG… 95d846… L2/3-6… neuron        1 <NA>    <NA>    <NA>    29285d… H21.33…
#> 10 AGTG… b0e1c5… L2/3-6… neuron        1 <NA>    <NA>    <NA>    cd7823… H21.33…
#> # … with 168,850 more rows, 67 more variables: assay <chr>,
#> #   assay_ontology_term_id <chr>, file_id_db <chr>,
#> #   cell_type_ontology_term_id <chr>, development_stage <chr>,
#> #   development_stage_ontology_term_id <chr>, disease <chr>,
#> #   disease_ontology_term_id <chr>, ethnicity <chr>,
#> #   ethnicity_ontology_term_id <chr>, experiment___ <chr>, file_id <chr>,
#> #   is_primary_data_x <chr>, organism <chr>, organism_ontology_term_id <chr>, …
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
