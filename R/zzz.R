#' @importFrom purrr walk
#' @importFrom dbplyr remote_con
#' @importFrom DBI dbDisconnect
.onUnload <- function(libname, pkgname){
    # Close connections to all cached tables. This should avoid most of the
    # "Connection is garbage-collected" messages
    cache$metadata_table |>
        as.list() |>
        walk(function(table){
            table |>
                remote_con() |>
                dbDisconnect()
        })
}
