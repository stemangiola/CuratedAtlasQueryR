# Utility scripts that are used internally by the package at runtime

#' Gets the file size of a number of remote files
#' @param urls A character vector containing URLs
#' @return The file size of each of the files pointed to by the provided URL,
#' in gigabytes, as double vector
#' @importFrom purrr map_dbl
#' @importFrom httr HEAD
#' @keywords internal
url_file_size <- function(urls){
    map_dbl(urls, function(url){
        as.integer(
            HEAD(url)$headers$`content-length` 
        ) / 10^9
    })
}

#' Prints a message indicating the size of a download
#' @inheritParams url_file_size
#' @importFrom cli cli_alert_info
#' @keywords internal
#' @return `NULL`, invisibly
report_file_sizes <- function(urls){
    total_size <- url_file_size(urls) |> 
        sum() |>
        round(digits=2)
    
    "Downloading {length(urls)} file{?s}, totalling {total_size} GB" |>
        cli_alert_info()
    
    invisible(NULL)
}

#' Formats a multi-line string as it it were on one line
#' @param text Any character vector
#' @return The same character vector, with newlines and subsequent whitespace
#'   removed
#' @keywords internal
#' @importFrom stringr str_remove_all
single_line_str <- function(text){
    str_remove_all(text, r"(\n\s*)")
}

#' Returns the default cache directory with a version number
#' @export
#' @return A length one character vector.
#' @importFrom tools R_user_dir
#' @importFrom utils packageName
#' @keywords internal
get_default_cache_dir <- function() {
    packageName() |>
        R_user_dir(
            "cache"
        ) |>
        file.path(COUNTS_VERSION) |>
        normalizePath() |>
        suppressWarnings()
}

#' Clear the default cache directory
#' @export
#' @return A length one character vector.
#' @importFrom tools R_user_dir
#' @importFrom utils packageName
#' @keywords internal
clear_cache <- function() {
  get_default_cache_dir() |> unlink(TRUE, TRUE)
}


#' Synchronises a single remote file with a local path
#' @importFrom httr write_disk GET stop_for_status
#' @importFrom cli cli_abort cli_alert_info
#' @return `NULL`, invisibly
#' @keywords internal
sync_remote_file <- function(full_url, output_file, ...) {
    if (!file.exists(output_file)) {
        output_dir <- dirname(output_file)
        dir.create(output_dir,
                   recursive = TRUE,
                   showWarnings = FALSE
        )
        cli_alert_info("Downloading {full_url} to {output_file}")
        
        tryCatch(
            GET(full_url, write_disk(output_file), ...) |> stop_for_status(),
            error = function(e) {
                # Clean up if we had an error
                file.remove(output_file)
                cli_abort("File {full_url} could not be downloaded. {e}")
            }
        )
    }
    invisible(NULL)
}

#' Returns a tibble from a parquet file path
#' 
#' Since dbplyr 2.4.0, raw file paths aren't handled very well
#' See: https://github.com/duckdb/duckdb-r/issues/38
#' Hence the need for this method
#' @importFrom glue glue
#' @importFrom dplyr tbl
#' @importFrom dbplyr sql
#' @importFrom glue glue_sql
#' @return An SQL data frame
#' @keywords internal
read_parquet <- function(conn, path){
    from_clause <- glue_sql("FROM read_parquet([{`path`*}], union_by_name=true)", .con=conn) |> sql()
    tbl(conn, from_clause)
}


