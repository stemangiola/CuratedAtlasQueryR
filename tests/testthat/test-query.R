library(HCAquery)

LOCAL_HCA_RAW = Sys.getenv("LOCAL_HCA_RAW")
LOCAL_HCA_SCALED = Sys.getenv("LOCAL_HCA_SCALED")
REMOTE_HCA_RAW = Sys.getenv("REMOTE_HCA_RAW")
REMOTE_HCA_SCALED = Sys.getenv("REMOTE_HCA_SCALED")

test_that("get_default_cache_dir() returns the correct directory on Linux", {
  grepl("linux", version$platform, fixed = TRUE) |>
    skip_if_not()
  
  expect_equal(
    get_default_cache_dir(),
    "~/.cache/hca_harmonised"
  )
})

test_that("sync_remote_files() syncs appropriate files", {
  # We need a remote URL to run this test
  skip_if(REMOTE_HCA_RAW == "")
  
  temp = tempfile()
  test_file = "00095cb0de0dc9528316b636fc9b3446"
  
  sync_remote_files(
    url = httr::parse_url(REMOTE_HCA_RAW),
    cache_dir = temp,
    files = test_file
  )
  
  test_file %in% list.files(temp) |>
    expect(failure_message = "The correct subdirectory was not created")
  
  downloaded = file.path(temp, test_file) |>
    list.files() 
  
  expect_equal(downloaded, c("assays.h5", "se.rds"))
})

test_that("get_SingleCellExperiment() syncs appropriate files", {
  # We need a remote URL to run this test
  skip_if(REMOTE_HCA_RAW == "")
  
  temp = tempfile()
  test_file = "00095cb0de0dc9528316b636fc9b3446"
  
  meta = get_metadata() |> head(2)
  
  # The remote dataset should have many genes
  sce = get_SingleCellExperiment(meta, repository = REMOTE_HCA_RAW, cache_dir = temp)
  sce |> row.names() |> length() |> expect_gt(1)
})

test_that("The genes argument to get_SingleCellExperiment subsets genes", {
  # We need a local copy to run this test
  skip_if(LOCAL_HCA_RAW == "")
  
  meta = get_metadata() |> head(2)
  
  # The un-subset dataset should have many genes
  sce_full = get_SingleCellExperiment(meta, cache_dir = LOCAL_HCA_RAW) |> row.names() |> length()
  expect_gt(sce_full, 1)
  
  # The subset dataset should only have one gene
  sce_subset = get_SingleCellExperiment(meta, cache_dir = LOCAL_HCA_RAW, genes = "PUM1") |> row.names() |> length()
  expect_equal(sce_subset, 1)
  
  expect_gt(sce_full, sce_subset)
})

test_that("get_seurat() returns the appropriate data in Seurat format", {
  # We need a local copy to run this test
  skip_if(LOCAL_HCA_RAW == "")
  
  meta = get_metadata() |> head(2)
  
  sce = get_SingleCellExperiment(meta, cache_dir = LOCAL_HCA_RAW, genes = "PUM1")
  seurat = get_seurat(meta, cache_dir = LOCAL_HCA_RAW, genes = "PUM1")
  
  # The output should be a Seurat object
  expect_s4_class(seurat, "Seurat")
  # Both methods should have appropriately subset genes
  expect_equal(
    rownames(sce),
    rownames(seurat)
  )
})
