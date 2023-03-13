test_that("url_file_size() returns the correct sizes", {
    c(
        "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz",
        "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz"
    ) |>
        url_file_size() |>
        expect_equal(c(
            0.973,
            0.944
        ), tolerance = 0.001)
})
