HCA Harmonised
================

Load the package

``` r
library(HCAquery)
library(dplyr)
library(dbplyr)
library(SingleCellExperiment)
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
#> ℹ Compiling Single Cell Experiment.
#> ℹ Attaching metadata.

sce
#> # A SingleCellExperiment-tibble abstraction: 1,571 × 56
#> # [90mFeatures=60661 | Cells=1571 | Assays=counts, cpm[0m
#>    .cell   sampl…¹ .sample .samp…² assay assay…³ file_…⁴ cell_…⁵ cell_…⁶ devel…⁷ devel…⁸ disease disea…⁹ ethni…˟ ethni…˟ file_id is_pr…˟ organ…˟ organ…˟ sampl…˟ sex   sex_o…˟ tissue
#>    <chr>   <chr>   <chr>   <chr>   <chr> <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr> <chr>   <chr> 
#>  1 ACAGCC… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#>  2 GGGAAT… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#>  3 TCTTCG… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#>  4 CCTTAC… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#>  5 ATCTAC… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#>  6 CATCAG… c7f14e… 9ef5ea… VUHD69… 10x … EFO:00… bc380d… CD4-po… CL:000… 31-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    fema… PATO:0… lung …
#>  7 AGTCTT… c7f14e… 9ef5ea… VUHD69… 10x … EFO:00… bc380d… CD4-po… CL:000… 31-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    fema… PATO:0… lung …
#>  8 CGCGTT… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#>  9 TGGCTG… c7f14e… 9ef5ea… VUHD69… 10x … EFO:00… bc380d… CD4-po… CL:000… 31-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    fema… PATO:0… lung …
#> 10 TTCCCA… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#> # … with 1,561 more rows, 33 more variables: tissue_ontology_term_id <chr>, tissue_harmonised <chr>, age_days <dbl>, dataset_id <chr>, collection_id <chr>, cell_count <int>,
#> #   dataset_deployments <chr>, is_primary_data.y <chr>, is_valid <int>, linked_genesets <int>, mean_genes_per_cell <dbl>, name <chr>, published <int>, revision <int>,
#> #   schema_version <chr>, tombstone <int>, x_normalization <chr>, created_at.x <dbl>, published_at <dbl>, revised_at <dbl>, updated_at.x <dbl>, filename <chr>, filetype <chr>,
#> #   s3_uri <chr>, user_submitted <int>, created_at.y <dbl>, updated_at.y <dbl>, cell_type_harmonised <chr>, confidence_class <dbl>, cell_annotation_azimuth_l2 <chr>,
#> #   cell_annotation_blueprint_singler <chr>, n_cell_type_in_tissue <int>, n_tissue_in_cell_type <int>, and abbreviated variable names ¹​sample_id_db, ²​.sample_name,
#> #   ³​assay_ontology_term_id, ⁴​file_id_db, ⁵​cell_type, ⁶​cell_type_ontology_term_id, ⁷​development_stage, ⁸​development_stage_ontology_term_id, ⁹​disease_ontology_term_id, ˟​ethnicity,
#> #   ˟​ethnicity_ontology_term_id, ˟​is_primary_data.x, ˟​organism, ˟​organism_ontology_term_id, ˟​sample_placeholder, ˟​sex_ontology_term_id
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
#> ℹ Compiling Single Cell Experiment.
#> ℹ Attaching metadata.

sce
#> # A SingleCellExperiment-tibble abstraction: 1,571 × 56
#> # [90mFeatures=60661 | Cells=1571 | Assays=cpm[0m
#>    .cell   sampl…¹ .sample .samp…² assay assay…³ file_…⁴ cell_…⁵ cell_…⁶ devel…⁷ devel…⁸ disease disea…⁹ ethni…˟ ethni…˟ file_id is_pr…˟ organ…˟ organ…˟ sampl…˟ sex   sex_o…˟ tissue
#>    <chr>   <chr>   <chr>   <chr>   <chr> <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr> <chr>   <chr> 
#>  1 ACAGCC… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#>  2 GGGAAT… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#>  3 TCTTCG… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#>  4 CCTTAC… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#>  5 ATCTAC… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#>  6 CATCAG… c7f14e… 9ef5ea… VUHD69… 10x … EFO:00… bc380d… CD4-po… CL:000… 31-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    fema… PATO:0… lung …
#>  7 AGTCTT… c7f14e… 9ef5ea… VUHD69… 10x … EFO:00… bc380d… CD4-po… CL:000… 31-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    fema… PATO:0… lung …
#>  8 CGCGTT… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#>  9 TGGCTG… c7f14e… 9ef5ea… VUHD69… 10x … EFO:00… bc380d… CD4-po… CL:000… 31-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    fema… PATO:0… lung …
#> 10 TTCCCA… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#> # … with 1,561 more rows, 33 more variables: tissue_ontology_term_id <chr>, tissue_harmonised <chr>, age_days <dbl>, dataset_id <chr>, collection_id <chr>, cell_count <int>,
#> #   dataset_deployments <chr>, is_primary_data.y <chr>, is_valid <int>, linked_genesets <int>, mean_genes_per_cell <dbl>, name <chr>, published <int>, revision <int>,
#> #   schema_version <chr>, tombstone <int>, x_normalization <chr>, created_at.x <dbl>, published_at <dbl>, revised_at <dbl>, updated_at.x <dbl>, filename <chr>, filetype <chr>,
#> #   s3_uri <chr>, user_submitted <int>, created_at.y <dbl>, updated_at.y <dbl>, cell_type_harmonised <chr>, confidence_class <dbl>, cell_annotation_azimuth_l2 <chr>,
#> #   cell_annotation_blueprint_singler <chr>, n_cell_type_in_tissue <int>, n_tissue_in_cell_type <int>, and abbreviated variable names ¹​sample_id_db, ²​.sample_name,
#> #   ³​assay_ontology_term_id, ⁴​file_id_db, ⁵​cell_type, ⁶​cell_type_ontology_term_id, ⁷​development_stage, ⁸​development_stage_ontology_term_id, ⁹​disease_ontology_term_id, ˟​ethnicity,
#> #   ˟​ethnicity_ontology_term_id, ˟​is_primary_data.x, ˟​organism, ˟​organism_ontology_term_id, ˟​sample_placeholder, ˟​sex_ontology_term_id
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
#> ℹ Compiling Single Cell Experiment.
#> ℹ Attaching metadata.
#> # A SingleCellExperiment-tibble abstraction: 1,571 × 56
#> # [90mFeatures=1 | Cells=1571 | Assays=counts, cpm[0m
#>    .cell   sampl…¹ .sample .samp…² assay assay…³ file_…⁴ cell_…⁵ cell_…⁶ devel…⁷ devel…⁸ disease disea…⁹ ethni…˟ ethni…˟ file_id is_pr…˟ organ…˟ organ…˟ sampl…˟ sex   sex_o…˟ tissue
#>    <chr>   <chr>   <chr>   <chr>   <chr> <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr> <chr>   <chr> 
#>  1 ACAGCC… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#>  2 GGGAAT… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#>  3 TCTTCG… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#>  4 CCTTAC… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#>  5 ATCTAC… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#>  6 CATCAG… c7f14e… 9ef5ea… VUHD69… 10x … EFO:00… bc380d… CD4-po… CL:000… 31-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    fema… PATO:0… lung …
#>  7 AGTCTT… c7f14e… 9ef5ea… VUHD69… 10x … EFO:00… bc380d… CD4-po… CL:000… 31-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    fema… PATO:0… lung …
#>  8 CGCGTT… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#>  9 TGGCTG… c7f14e… 9ef5ea… VUHD69… 10x … EFO:00… bc380d… CD4-po… CL:000… 31-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    fema… PATO:0… lung …
#> 10 TTCCCA… 33cdeb… 4fc10a… VUHD92… 10x … EFO:00… bc380d… CD4-po… CL:000… 55-yea… HsapDv… normal  PATO:0… African HANCES… 6661ab… TRUE    Homo s… NCBITa… <NA>    male  PATO:0… lung …
#> # … with 1,561 more rows, 33 more variables: tissue_ontology_term_id <chr>, tissue_harmonised <chr>, age_days <dbl>, dataset_id <chr>, collection_id <chr>, cell_count <int>,
#> #   dataset_deployments <chr>, is_primary_data.y <chr>, is_valid <int>, linked_genesets <int>, mean_genes_per_cell <dbl>, name <chr>, published <int>, revision <int>,
#> #   schema_version <chr>, tombstone <int>, x_normalization <chr>, created_at.x <dbl>, published_at <dbl>, revised_at <dbl>, updated_at.x <dbl>, filename <chr>, filetype <chr>,
#> #   s3_uri <chr>, user_submitted <int>, created_at.y <dbl>, updated_at.y <dbl>, cell_type_harmonised <chr>, confidence_class <dbl>, cell_annotation_azimuth_l2 <chr>,
#> #   cell_annotation_blueprint_singler <chr>, n_cell_type_in_tissue <int>, n_tissue_in_cell_type <int>, and abbreviated variable names ¹​sample_id_db, ²​.sample_name,
#> #   ³​assay_ontology_term_id, ⁴​file_id_db, ⁵​cell_type, ⁶​cell_type_ontology_term_id, ⁷​development_stage, ⁸​development_stage_ontology_term_id, ⁹​disease_ontology_term_id, ˟​ethnicity,
#> #   ˟​ethnicity_ontology_term_id, ˟​is_primary_data.x, ˟​organism, ˟​organism_ontology_term_id, ˟​sample_placeholder, ˟​sex_ontology_term_id
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
#> ℹ Compiling Single Cell Experiment.
#> ℹ Attaching metadata.
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#> An object of class Seurat 
#> 60661 features across 1571 samples within 1 assay 
#> Active assay: originalexp (60661 features, 0 variable features)
```

``` r
sessionInfo()
#> R version 4.2.1 (2022-06-23)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: CentOS Linux 7 (Core)
#> 
#> Matrix products: default
#> BLAS:   /stornext/System/data/apps/R/R-4.2.1/lib64/R/lib/libRblas.so
#> LAPACK: /stornext/System/data/apps/R/R-4.2.1/lib64/R/lib/libRlapack.so
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] SingleCellExperiment_1.18.1 dbplyr_2.2.1                dplyr_1.0.10                HCAquery_0.1.0              SummarizedExperiment_1.26.1 Biobase_2.56.0             
#>  [7] GenomicRanges_1.48.0        GenomeInfoDb_1.32.4         IRanges_2.30.1              S4Vectors_0.34.0            BiocGenerics_0.42.0         MatrixGenerics_1.8.1       
#> [13] matrixStats_0.63.0          testthat_3.1.6             
#> 
#> loaded via a namespace (and not attached):
#>   [1] utf8_1.2.2                     spatstat.explore_3.0-5         reticulate_1.26                tidyselect_1.2.0               RSQLite_2.2.20                
#>   [6] htmlwidgets_1.6.0              grid_4.2.1                     Rtsne_0.16                     devtools_2.4.5                 munsell_0.5.0                 
#>  [11] xmlparsedata_1.0.5             codetools_0.2-18               ica_1.0-3                      future_1.30.0                  miniUI_0.1.1.1                
#>  [16] withr_2.5.0                    spatstat.random_3.0-1          colorspace_2.0-3               progressr_0.12.0               knitr_1.41                    
#>  [21] rstudioapi_0.14                Seurat_4.3.0                   ROCR_1.0-11                    tensor_1.5                     listenv_0.9.0                 
#>  [26] GenomeInfoDbData_1.2.8         polyclip_1.10-4                bit64_4.0.5                    rhdf5_2.40.0                   rprojroot_2.0.3               
#>  [31] parallelly_1.33.0              vctrs_0.5.1                    generics_0.1.3                 xfun_0.36                      R6_2.5.1                      
#>  [36] bitops_1.0-7                   rhdf5filters_1.8.0             spatstat.utils_3.0-1           cachem_1.0.6                   DelayedArray_0.22.0           
#>  [41] assertthat_0.2.1               promises_1.2.0.1               scales_1.2.1                   gtable_0.3.1                   globals_0.16.2                
#>  [46] processx_3.8.0                 goftest_1.2-3                  rlang_1.0.6                    codegrip_0.0.0.9000            splines_4.2.1                 
#>  [51] lazyeval_0.2.2                 spatstat.geom_3.0-3            yaml_2.3.6                     reshape2_1.4.4                 abind_1.4-5                   
#>  [56] httpuv_1.6.7                   tools_4.2.1                    usethis_2.1.6                  ggplot2_3.4.0                  ellipsis_0.3.2                
#>  [61] RColorBrewer_1.1-3             sessioninfo_1.2.2              ggridges_0.5.4                 Rcpp_1.0.9                     plyr_1.8.8                    
#>  [66] zlibbioc_1.42.0                purrr_1.0.0                    RCurl_1.98-1.9                 ps_1.7.2                       prettyunits_1.1.1             
#>  [71] deldir_1.0-6                   pbapply_1.6-0                  cowplot_1.1.1                  urlchecker_1.0.1               zoo_1.8-11                    
#>  [76] SeuratObject_4.1.3             ggrepel_0.9.2                  cluster_2.1.3                  fs_1.5.2                       magrittr_2.0.3                
#>  [81] data.table_1.14.6              scattermore_0.8                lmtest_0.9-40                  RANN_2.6.1                     fitdistrplus_1.1-8            
#>  [86] pkgload_1.3.2                  patchwork_1.1.2                mime_0.12                      evaluate_0.19                  xtable_1.8-4                  
#>  [91] gridExtra_2.3                  compiler_4.2.1                 tibble_3.1.8                   KernSmooth_2.23-20             crayon_1.5.2                  
#>  [96] htmltools_0.5.4                later_1.3.0                    tidyr_1.2.1                    DBI_1.1.3                      MASS_7.3-57                   
#> [101] rappdirs_0.3.3                 Matrix_1.5-3                   brio_1.1.3                     cli_3.5.0                      parallel_4.2.1                
#> [106] igraph_1.3.5                   pkgconfig_2.0.3                ttservice_0.2.2                sp_1.5-1                       plotly_4.10.1                 
#> [111] spatstat.sparse_3.0-0          xml2_1.3.3                     roxygen2_7.2.3                 XVector_0.36.0                 stringr_1.5.0                 
#> [116] callr_3.7.3                    digest_0.6.31                  sctransform_0.3.5              RcppAnnoy_0.0.20               spatstat.data_3.0-0           
#> [121] rmarkdown_2.19                 leiden_0.4.3                   uwot_0.1.14                    curl_4.3.3                     shiny_1.7.4                   
#> [126] lifecycle_1.0.3                nlme_3.1-157                   jsonlite_1.8.4                 Rhdf5lib_1.18.2                desc_1.4.2                    
#> [131] viridisLite_0.4.1              fansi_1.0.3                    pillar_1.8.1                   lattice_0.20-45                fastmap_1.1.0                 
#> [136] httr_1.4.4                     pkgbuild_1.4.0                 survival_3.3-1                 glue_1.6.2                     remotes_2.4.2                 
#> [141] png_0.1-8                      bit_4.0.5                      stringi_1.7.8                  profvis_0.3.7                  HDF5Array_1.24.2              
#> [146] blob_1.2.3                     tidySingleCellExperiment_1.6.3 memoise_2.0.1                  irlba_2.3.5.1                  future.apply_1.10.0
```
