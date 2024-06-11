# Functions that relate to unharmonised metadata

#' @include utils.R
NULL

#' Base URL for all the unharmonised data
#' @noRd
UNHARMONISED_URL <- single_line_str(
    "https://object-store.rc.nectar.org.au/v1/
    AUTH_06d6e008e3e642da99d806ba3ea629c5/unharmonised_metadata"
)

#' Returns unharmonised metadata for selected datasets.
#'
#' Various metadata fields are *not* common between datasets, so it does not
#' make sense for these to live in the main metadata table. This function is a
#' utility that allows easy fetching of this data if necessary.
#'
#' @param dataset_id A character vector, where each entry is a dataset ID
#'   obtained from the `$file_id` column of the table returned from
#'   [get_metadata()]
#' @param cells An optional character vector of cell IDs. If provided, only
#'   metadata for those cells will be returned.
#' @param conn An optional DuckDB connection object. If provided, it will re-use
#'   the existing connection instead of opening a new one.
#' @param remote_url Optional character vector of length 1. An HTTP URL pointing
#'   to the root URL under which all the unharmonised dataset files are located.
#' @param cache_directory Optional character vector of length 1. A file path on
#'   your local system to a directory (not a file) that will be used to store
#'   the unharmonised metadata files.
#' @importFrom purrr map set_names
#' @importFrom glue glue
#' @importFrom DBI dbConnect
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl filter
#' @importFrom rlang .data
#' @return A named list, where each name is a dataset file ID, and each value is
#'   a "lazy data frame", ie a `tbl`.
#' @examples
#' \dontrun{
#' dataset <- "838ea006-2369-4e2c-b426-b2a744a2b02b"
#' harmonised_meta <- get_metadata() |> 
#'     dplyr::filter(file_id == dataset) |> dplyr::collect()
#' unharmonised_meta <- get_unharmonised_dataset(dataset)
#' unharmonised_tbl <- dplyr::collect(unharmonised_meta[[dataset]])
#' dplyr::left_join(harmonised_meta, unharmonised_tbl, by=c("file_id", "cell_"))
#' }
#' @references Mangiola, S., M. Milton, N. Ranathunga, C. S. N. Li-Wai-Suen, 
#'   A. Odainic, E. Yang, W. Hutchison et al. "A multi-organ map of the human 
#'   immune system across age, sex and ethnicity." bioRxiv (2023): 2023-06.
#'   doi:10.1101/2023.06.08.542671.
#' @source [Mangiola et al.,2023](https://www.biorxiv.org/content/10.1101/2023.06.08.542671v3)
get_unharmonised_dataset <- function(
    dataset_id,
    cells = NULL,
    conn = duckdb() |> dbConnect(drv = _, read_only = TRUE),
    remote_url = UNHARMONISED_URL,
    cache_directory = get_default_cache_dir()
){
    unharmonised_root <- file.path(
      cache_directory,
      "unharmonised"
    )
    file_name <- glue::glue("{dataset_id}.parquet")
    local_path <- file.path(unharmonised_root, file_name)
    glue("{remote_url}/{file_name}") |>
        sync_remote_file(
            local_path,
            progress(type = "down", con = stderr())
        )
    
    read_parquet(conn, local_path) |>
        filter(.data$cell_ %in% cells)
}

#' Returns unharmonised metadata for a metadata query
#' @inherit get_unharmonised_dataset description
#' @param metadata A lazy data frame obtained from [get_metadata()], filtered
#'   down to some cells of interest
#' @inheritDotParams get_unharmonised_dataset
#' @return A tibble with two columns:
#'  * `file_id`: the same `file_id` as the main metadata table obtained from
#'    [get_metadata()]
#'  * `unharmonised`: a nested tibble, with one row per cell in the input
#'    `metadata`, containing unharmonised metadata
#' @export
#' @importFrom dplyr group_by summarise filter collect
#' @importFrom rlang .data
#' @importFrom dbplyr remote_con
#' @examples
#' harmonised <- dplyr::filter(get_metadata(), tissue == "kidney blood vessel")
#' unharmonised <- get_unharmonised_metadata(harmonised)
#' @references Mangiola, S., M. Milton, N. Ranathunga, C. S. N. Li-Wai-Suen, 
#'   A. Odainic, E. Yang, W. Hutchison et al. "A multi-organ map of the human 
#'   immune system across age, sex and ethnicity." bioRxiv (2023): 2023-06.
#'   doi:10.1101/2023.06.08.542671.
#' @source [Mangiola et al.,2023](https://www.biorxiv.org/content/10.1101/2023.06.08.542671v3)
get_unharmonised_metadata <- function(metadata, ...){
    args <- list(...)
    metadata |>
        collect() |>
        group_by(.data$file_id) |>
        summarise(
            unharmonised = list(
              dataset_id=.data$file_id[[1]],
              cells=.data$cell_,
              conn=remote_con(metadata)
            ) |>
                c(args) |> 
                do.call(get_unharmonised_dataset, args=_) |> 
                list()
        )
}
