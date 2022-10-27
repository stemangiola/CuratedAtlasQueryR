library(HCAquery)

test_that("The genes argument to get_SingleCellExperiment subsets genes", {
  meta = get_metadata() |> head(2)
  
  # The un-subset dataset should have many genes
  sce_full = get_SingleCellExperiment(meta) |> row.names() |> length()
  expect_gt(sce_full, 1)
  
  # The subset dataset should only have one gene
  sce_subset = get_SingleCellExperiment(meta, genes = "PUM1") |> row.names() |> length()
  expect_equal(sce_subset, 1)
  
  expect_gt(sce_full, sce_subset)
})
