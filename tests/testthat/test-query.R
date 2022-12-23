library(HCAquery)

LOCAL_HCA = Sys.getenv("LOCAL_HCA")
REMOTE_HCA = Sys.getenv("REMOTE_HCA")

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
  skip_if(REMOTE_HCA == "")
  
  temp = tempfile()
  test_file = "00095cb0de0dc9528316b636fc9b3446"
  
  sync_remote_files(
    url = httr::parse_url(REMOTE_HCA),
    subdirs = "cpm",
    cache_dir = temp,
    files = test_file
  )
  
  temp_subdir = file.path(temp, "cpm") 
  
  test_file %in% list.files(temp_subdir) |>
    expect(failure_message = "The correct subdirectory was not created")
  
  downloaded = file.path(temp_subdir, test_file) |>
    list.files() 
  
  expect_equal(downloaded, c("assays.h5", "se.rds"))
})

test_that("get_SingleCellExperiment() syncs appropriate files", {
  # We need a remote URL to run this test
  skip_if(REMOTE_HCA == "")
  
  temp = tempfile()
  test_file = "00095cb0de0dc9528316b636fc9b3446"
  
  meta = get_metadata() |> head(2)
  
  # The remote dataset should have many genes
  sce = get_SingleCellExperiment(meta, repository = REMOTE_HCA, cache_directory = temp)
  sce |> row.names() |> length() |> expect_gt(1)
})

test_that("The assays argument to get_SingleCellExperiment controls the number of returned assays", {
  # We need this for the assays() function
  library(SummarizedExperiment)
  
  # We need a local copy to run this test
  skip_if(LOCAL_HCA == "")
  
  meta = get_metadata() |> head(2)

  # If we request both assays, we get both assays  
  get_SingleCellExperiment(meta, cache_directory = LOCAL_HCA, assays = "counts", "cpm") |> assays() |> names() |> expect_mapequal(c("counts", "cpm"))
  
  # If we request one assay, we get one assays  
  get_SingleCellExperiment(meta, cache_directory = LOCAL_HCA, assays = "counts") |> assays() |> names() |> expect_mapequal("counts")
  get_SingleCellExperiment(meta, cache_directory = LOCAL_HCA, assays = "cpm") |> assays() |> names() |> expect_mapequal("cpm")
})

test_that("The features argument to get_SingleCellExperiment subsets genes", {
  # We need a local copy to run this test
  skip_if(LOCAL_HCA == "")
  
  meta = get_metadata() |> head(2)
  
  # The un-subset dataset should have many genes
  sce_full = get_SingleCellExperiment(meta, cache_directory = LOCAL_HCA) |> row.names() |> length()
  expect_gt(sce_full, 1)
  
  # The subset dataset should only have one gene
  sce_subset = get_SingleCellExperiment(meta, cache_directory = LOCAL_HCA, features = "PUM1") |> row.names() |> length()
  expect_equal(sce_subset, 1)
  
  expect_gt(sce_full, sce_subset)
})

test_that("get_seurat() returns the appropriate data in Seurat format", {
  # We need a local copy to run this test
  skip_if(LOCAL_HCA == "")
  
  meta = get_metadata() |> head(2)
  
  sce = get_SingleCellExperiment(meta, cache_directory = LOCAL_HCA, features = "PUM1")
  seurat = get_seurat(meta, cache_directory = LOCAL_HCA, features = "PUM1")
  
  # The output should be a Seurat object
  expect_s4_class(seurat, "Seurat")
  # Both methods should have appropriately subset genes
  expect_equal(
    rownames(sce),
    rownames(seurat)
  )
})
