if(requireNamespace('spelling', quietly = TRUE))
  spelling::spell_check_test(
    # We hit https://github.com/ropensci/spelling/issues/62 if the 
    # vignettes are checked
    vignettes = FALSE,
    error = FALSE,
    skip_on_cran = TRUE
 )
