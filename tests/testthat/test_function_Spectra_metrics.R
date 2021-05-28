## create toy example data set
fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
sps_sciex <- Spectra(fls, backend = MsBackendMzR())
set.seed(1)

## add some fake charges
sps_sciex@backend$precursorCharge <- as.integer(
    sample(x = c(1, 2, 3, 4), size = 1862, 
                        replace = TRUE, prob = c(0.7, 0.15, 0.1, 0.05)))

## START unit test RatioCharge1over2 ##
test_that("RatioCharge1over2", {
    expect_error(RatioCharge1over2(NULL), "unable to find an inherited method")
    expect_error(RatioCharge1over2(1:10), "unable to find an inherited method")
    expect_equal(RatioCharge1over2(sps_sciex), 0.7 / 0.15, tolerance = 3e-02)
    
    ## do not include charges 1 or 2
    sps_sciex_foo <- sps_sciex
    sps_sciex_foo@backend$precursorCharge <- as.integer(
        sample(x = c(1, 2, 3, 4), size = 1862, 
               replace = TRUE, prob = c(0.0, 0.85, 0.1, 0.05)))
    expect_equal(RatioCharge1over2(sps_sciex_foo), NA)
    sps_sciex_foo <- sps_sciex
    sps_sciex_foo@backend$precursorCharge <- as.integer(
        sample(x = c(1, 2, 3, 4), size = 1862, 
               replace = TRUE, prob = c(0.85, 0.0, 0.1, 0.05)))
    expect_equal(RatioCharge1over2(sps_sciex_foo), NA)
})
## END unit test RatioCharge1over2 ##

## START unit test RatioCharge3over2 ##
test_that("RatioCharge3over2", {
    expect_error(RatioCharge3over2(NULL), "unable to find an inherited method")
    expect_error(RatioCharge3over2(1:10), "unable to find an inherited method")
    expect_equal(RatioCharge3over2(sps_sciex), 0.1 / 0.15, tolerance = 2e-02)
    
    ## do not include charges 2 or 3
    sps_sciex_foo <- sps_sciex
    sps_sciex_foo@backend$precursorCharge <- as.integer(
        sample(x = c(1, 2, 3, 4), size = 1862, 
               replace = TRUE, prob = c(0.7, 0.0, 0.25, 0.05)))
    expect_equal(RatioCharge3over2(sps_sciex_foo), NA)
    sps_sciex_foo <- sps_sciex
    sps_sciex_foo@backend$precursorCharge <- as.integer(
        sample(x = c(1, 2, 3, 4), size = 1862, 
               replace = TRUE, prob = c(0.7, 0.25, 0.0, 0.05)))
    expect_equal(RatioCharge3over2(sps_sciex_foo), NA)
})
## END unit test RatioCharge3over2 ##

## START unit test RatioCharge4over2 ##
test_that("RatioCharge4over2", {
    expect_error(RatioCharge4over2(NULL), "unable to find an inherited method")
    expect_error(RatioCharge4over2(1:10), "unable to find an inherited method")
    expect_equal(RatioCharge4over2(sps_sciex), 0.05 / 0.15, tolerance = 3e-02)
    
    ## do not include charges 2 or 4
    sps_sciex_foo <- sps_sciex
    sps_sciex_foo@backend$precursorCharge <- as.integer(
        sample(x = c(1, 2, 3, 4), size = 1862, 
               replace = TRUE, prob = c(0.7, 0.0, 0.1, 0.20)))
    expect_equal(RatioCharge4over2(sps_sciex_foo), NA)
    sps_sciex_foo <- sps_sciex
    sps_sciex_foo@backend$precursorCharge <- as.integer(
        sample(x = c(1, 2, 3, 4), size = 1862, 
               replace = TRUE, prob = c(0.7, 0.25, 0.15, 0.0)))
    expect_equal(RatioCharge4over2(sps_sciex_foo), NA)
})
## END unit test RatioCharge4over2 ##
