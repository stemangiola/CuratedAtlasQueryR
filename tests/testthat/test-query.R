library(HCAquery)

test_that("get_SingleCellExperiment() correctly handles duplicate cell IDs", {
    meta <- get_metadata() |> dplyr::filter(.cell == "868417_1")
    sce <- get_SingleCellExperiment(meta)
    # This query should return multiple cells, despite querying only 1 cell ID
    meta |> dplyr::tally() |> dplyr::pull(n) |> expect_gt(1)
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

    expect_equal(
        get_default_cache_dir(),
        "~/.cache/hca_harmonised"
    )
})

test_that("sync_assay_files() syncs appropriate files", {
    temp <- tempfile()
    test_file <- "ffd2891329f66921dbe361af85b07c51"

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
