# Utility scripts for development purposes, that are not exported to users

#' Upload a file to the Nectar object store
#' @param source A character scalar indicating the local path to the file to
#'   upload
#' @param container A character scalar indicating the name of the container to
#'   upload to
#' @param name An optional character scalar indicating the name the file should
#'   have after being uploaded. Defaults to being the basename of the source
#'   file.
#' @param credential_id The OpenStack application credential ID as a character
#'   scalar. This is optional because you can alternatively source a
#'   `-openrc.sh` file instead of providing it here.
#' @param credential_id The OpenStack application credential secret as a
#'   character scalar
#' @return NULL
#' @keywords internal
upload_swift = function(source, container, name = basename(source), credential_id = NULL, credential_secret = NULL){
    # Create the basilisk environment
    swift_env <- basilisk::BasiliskEnvironment(
        envname="swift-nectar-upload",
        pkgname=packageName(),
        packages=c("python-swiftclient==4.2.0", "python-keystoneclient==5.1.0", "python==3.10.9")
    )
    proc <- basilisk::basiliskStart(swift_env)
    
    # Build the CLI args
    if (!is.null(credential_id) && !is.null(credential_secret)){
        auth <- c(
            "--os-auth-type",
            "v3applicationcredential",
            "--os-application-credential-id",
            credential_id,
            "--os-application-credential-secret",
            credential_secret
        )
    }
    else {
        auth <- character()
    }
    args = c(
        "-m", 
        "swiftclient.shell",
        "--os-auth-url",
        "https://keystone.rc.nectar.org.au:5000/v3/",
        "--os-project-id",
        "06d6e008e3e642da99d806ba3ea629c5",
        auth,
        "upload",
        container,
        source,
        "--object-name",
        name
    )
    
    # Perform the upload
    system2(reticulate::py_exe(), args=args)
    basilisk::basiliskStop(proc)
    
    invisible(NULL)
}

#' Update the metadata database in nectar using a newly created data frame
#' @param metadata The data frame to upload
#' @param version The version for the new metadata as a character scalar, e.g.
#'   "0.2.3"
#' @inheritDotParams upload_swift
#' @examples
#' \dontrun{
#'  metadata = CuratedAtlasQueryR::get_metadata() |> head(10) |> dplyr::collect()
#'  update_database(metadata, "0.2.3", credential_id = "ABCDEFGHIJK", credential_secret = "ABCD1234EFGH-5678IJK")
#'  # Prints "metadata.0.2.3.parquet" if successful
#' }
#' @keywords internal
update_database = function(metadata, version, ...){
    # These are optional dev packages
    rlang::check_installed(c("arrow", "glue", "basilisk"))
    
    dir <- tempdir()
    parquet_name <- glue::glue("metadata.{version}.parquet")
    parquet_path <- file.path(dir, parquet_name)
    arrow::write_parquet(metadata, sink=parquet_path)
    
    upload_swift(parquet_path, container="metadata", name=parquet_name, ...)
}

#' Update the unharmonised parquet files
#' @param unharmonised_parquet_dir The path to a directory containing parquet
#'   files, one for each dataset, e.g.
#'   /vast/projects/cellxgene_curated/metadata_non_harmonised_parquet_0.2
#' @inheritDotParams upload_swift
#' @keywords internal
#' @examples
#' \dontrun{
#' update_unharmonised("/vast/projects/cellxgene_curated/metadata_non_harmonised_parquet_0.2", credential_id = "ABCDEFGHIJK", credential_secret = "ABCD1234EFGH-5678IJK")
#' }
update_unharmonised = function(unharmonised_parquet_dir, ...){
    # name="/" forces it have no prefix, ie be at the top level in the bucket
    upload_swift(unharmonised_parquet_dir, container="unharmonised_metadata", name="/", ...)
}
