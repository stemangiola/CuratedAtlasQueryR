library(dplyr)

test_that("get_SingleCellExperiment() correctly handles duplicate cell IDs", {
    meta <- get_metadata() |>
        dplyr::filter(cell_ == "868417_1") |>
        dplyr::collect()
    sce <- get_SingleCellExperiment(meta)
    # This query should return multiple cells, despite querying only 1 cell ID
    nrow(meta) |> expect_gt(1)
    # Each of the two ambiguous cell IDs should now be unique
    colnames(sce) |> expect_equal(c("868417_1_1", "868417_1_2"))
    # We should have lots of column data, derived from the metadata
    SummarizedExperiment::colData(sce) |>
        dim() |>
        expect_equal(c(2, 56))
})

test_that("get_default_cache_dir() returns the correct directory on Linux", {
    grepl("linux", version$platform, fixed = TRUE) |>
        skip_if_not()

    "~/.cache/R/CuratedAtlasQueryR" |>
        normalizePath() |>
        expect_equal(
            get_default_cache_dir(),
        )
})

test_that("sync_assay_files() syncs appropriate files", {
    temp <- tempfile()
    test_file <- "00095cb0de0dc9528316b636fc9b3446"

    sync_assay_files(
        subdirs = "cpm",
        cache_dir = temp,
        files = test_file
    )

    temp_subdir <- file.path(temp, "cpm")

    test_file %in% list.files(temp_subdir) |>
        expect(failure_message = "The correct subdirectory was not created")

    downloaded <- file.path(temp_subdir, test_file) |>
        list.files()

    expect_equal(downloaded, c("assays.h5", "se.rds"))
})

test_that("get_SingleCellExperiment() syncs appropriate files", {
    temp <- tempfile()
    test_file <- "00095cb0de0dc9528316b636fc9b3446"

    meta <- get_metadata() |> head(2)

    # The remote dataset should have many genes
    sce <- get_SingleCellExperiment(meta, cache_directory = temp)
    sce |>
        row.names() |>
        length() |>
        expect_gt(1)
})

test_that(
    "The assays argument to get_SingleCellExperiment controls the number
  of returned assays",
    {
        # We need this for the assays() function
        library(SummarizedExperiment)

        meta <- get_metadata() |> head(2)

        # If we request both assays, we get both assays
        get_SingleCellExperiment(meta, assays = c("counts", "cpm")) |>
            assays() |>
            names() |>
            expect_setequal(c("counts", "cpm"))

        # If we request one assay, we get one assays
        get_SingleCellExperiment(meta, assays = "counts") |>
            assays() |>
            names() |>
            expect_setequal("counts")
        get_SingleCellExperiment(meta, assays = "cpm") |>
            assays() |>
            names() |>
            expect_setequal("cpm")
    }
)

test_that("The features argument to get_SingleCellExperiment subsets genes", {
    meta <- get_metadata() |> head(2)

    # The un-subset dataset should have many genes
    sce_full <- get_SingleCellExperiment(meta) |>
        row.names() |>
        length()
    expect_gt(sce_full, 1)

    # The subset dataset should only have one gene
    sce_subset <- get_SingleCellExperiment(meta, features = "PUM1") |>
        row.names() |>
        length()
    expect_equal(sce_subset, 1)

    expect_gt(sce_full, sce_subset)
})

test_that("get_seurat() returns the appropriate data in Seurat format", {
    meta <- get_metadata() |> head(2)

    sce <- get_SingleCellExperiment(meta, features = "PUM1")
    seurat <- get_seurat(meta, features = "PUM1")

    # The output should be a Seurat object
    expect_s4_class(seurat, "Seurat")
    # Both methods should have appropriately subset genes
    expect_equal(
        rownames(sce),
        rownames(seurat)
    )
})

test_that("get_SingleCellExperiment() assigns the right cell ID to each cell", {
    id = "3214d8f8986c1e33a85be5322f2db4a9"
    
    # Force the file to be cached
    get_metadata() |>
        filter(file_id_db == id) |>
        get_SingleCellExperiment()
    
    # Load the SCE from cache directly
    assay_1 = CuratedAtlasQueryR:::get_default_cache_dir() |>
        file.path("0.2/original", id) |>
        HDF5Array::loadHDF5SummarizedExperiment() |>
        assay("X") |>
        as.matrix()
    
    # Make a SCE that has the right column names, but reversed
    assay_2 = 
        assay_1 |>
        colnames() |>
        tibble::tibble(
            file_id_db = id,
            cell_ = _
        ) |>
        arrange(-row_number()) |>
        get_SingleCellExperiment(assays = "counts") |>
        assay("counts") |>
        as.matrix()
        
    colnames(assay_2) = sub("_1", "", x=colnames(assay_2))

    expect_equal(
        assay_1,
        assay_2[, colnames(assay_1)]
    )
})

test_that("get_unharmonised_metadata works with one ID", {
    dataset_id = "838ea006-2369-4e2c-b426-b2a744a2b02b"
    unharmonised_meta = get_unharmonised_metadata(dataset_id)
    unharmonised_tbl = unharmonised_meta[[dataset_id]]
    
    expect_s3_class(unharmonised_meta, "list")
    expect_s3_class(unharmonised_tbl, "tbl")
})

test_that("get_unharmonised_metadata works with multiple IDs", {
    dataset_ids = c("838ea006-2369-4e2c-b426-b2a744a2b02b", "83b9cb97-9ee4-404d-8cdf-ccede8235356.parquet")
    unharmonised_meta = get_unharmonised_metadata(dataset_ids)
    
    expect_equal(names(unharmonised_meta), dataset_ids)
})