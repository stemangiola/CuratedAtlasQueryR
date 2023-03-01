# Utility scripts for development purposes, that are not exported to users

#' Update the metadata database in nectar using a newly created data frame
#' @param metadata The data frame to upload
#' @param version The version for the new metadata as a character scalar, e.g.
#'   "0.2.3"
#' @param credential_id The OpenStack application credential ID as a character
#'   scalar. This is optional because you can alternatively source a
#'   `-openrc.sh` file instead of providing it here.
#' @param credential_id The OpenStack application credential secret as a
#'   character scalar
#' @noRd
update_database = function(metadata, version, credential_id = NULL, credential_secret = NULL){
    # These are optional dev packages
    rlang::check_installed(c("arrow", "glue", "basilisk"))
    
    # Create parquet
    dir <- tempdir()
    parquet_name <- glue::glue("metadata.{version}.parquet")
    parquet_path <- file.path(dir, parquet_name)
    arrow::write_parquet(metadata, sink=parquet_path)
    
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
        "metadata",
        parquet_path,
        "--object-name",
        parquet_name
    )
    
    # Perform the upload
    system2(reticulate::py_exe(), args=args)
    basilisk::basiliskStop(proc)
}