#' Gets the file size of a number of remote files
#' @param urls A character vector containing URLs
#' @return The file size of each of the files pointed to by the provided URL,
#' in gigabytes, as double vector
#' @importFrom purrr map_dbl
#' @keywords internal
url_file_size = function(urls){
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
report_file_sizes = function(urls){
    total_size = url_file_size(urls) |> 
        sum() |>
        round(digits=2)
    
    cli_alert_info("Downloading {length(urls)} file{?s}, totalling {total_size} GB")
}