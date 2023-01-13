test_that("get_SingleCellExperiment() correctly handles duplicate cell IDs", {
    meta = get_metadata() |> filter(.cell == "868417_1")
    # This query should return multiple cells, despite querying only 1 cell ID
    dplyr::collect(meta) |> nrow() |> expect_gt(1)
    sce <- get_SingleCellExperiment(meta)
    colnames(sce) |> expect_equal(meta$.cell)
})
