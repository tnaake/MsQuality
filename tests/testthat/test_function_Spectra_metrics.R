## create toy example data set
fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
sps_sciex <- Spectra(fls, backend = MsBackendMzR())
set.seed(1)

## add some fake charges
sps_sciex@backend$precursorCharge <- as.integer(
    sample(x = c(1, 2, 3, 4), size = 1862, 
                        replace = TRUE, prob = c(0.7, 0.15, 0.1, 0.05)))

## add some fake precursorMZ
sps_sciex@backend$precursorMz <- rnorm(n = 1862, mean = 500, sd = 100)
    

## START unit test rtimeDuration ##
test_that("rtimeDuration", {
    expect_error(rtimeDuration(NULL), "unable to find an inherited method")
    expect_error(rtimeDuration(NULL), "unable to find an inherited method")
    expect_equal(rtimeDuration(sps_sciex), 259.482)
})
## END unit test rtimeDuration ##

## START unit test RToverTICquantile ##
test_that("RToverTICquantile", {
    expect_error(RToverTICquantile(NULL), "unable to find an inherited method")
    expect_error(RToverTICquantile(NULL), "unable to find an inherited method")
    tmp <- RToverTICquantile(sps_sciex)
    expect_equal(as.numeric(tmp),
        c(0.001058682, 0.250272370, 0.500517792, 0.749731480, 1.000000000))
    expect_equal(names(tmp), c("0%", "25%", "50%", "75%", "100%"))
})
## END unit test RToverTICquantile ##

## START unit test RToverMSQuarters ##
test_that("RToverMSQuarters", {
    expect_error(RToverMSQuarters(NULL), "unable to find an inherited method")
    expect_error(RToverMSQuarters(NULL), "unable to find an inherited method")
    expect_equal(RToverMSQuarters(sps_sciex, MSLevel = 1L), 
        c(0.2502724, 0.5005370, 0.7507863, 1.0000000))
    expect_error(RToverMSQuarters(sps_sciex, MSLevel = 2L), 
        "Spectra object does not contain any spectra")
})
## END unit test RToverMSQuarters ##

## START unit test TICquantileToQuantileLogRatio ##
test_that("TICquantileToQuantileLogRatio", {
    expect_error(TICquantileToQuantileLogRatio(NULL), 
        "unable to find an inherited method")
    expect_error(TICquantileToQuantileLogRatio(NULL), 
        "unable to find an inherited method")
    expect_equal(
        TICquantileToQuantileLogRatio(sps_sciex, relativeTo = "previous", 
                                                                MSLevel = 1L), 
        c(-0.05853477, -0.58202994, 0.15361253))
    expect_equal(
        TICquantileToQuantileLogRatio(sps_sciex, relativeTo = "Q1", 
                                                                MSLevel = 1L), 
        c(-0.05853477, -0.64056471, -0.48695218))
    expect_error(TICquantiletoQuantileLogRatio(sps_sciex, MSLevel = 2L), 
                 "Spectra object does not contain any spectra")
})
## END unit test TICquantiletoQuantileLogRatio ##

## START unit test TICquantileToQuantileLogRatio ##
test_that("TICquantileToQuantileLogRatio", {
    expect_error(TICquantileToQuantileLogRatio(NULL), 
                 "unable to find an inherited method")
    expect_error(TICquantileToQuantileLogRatio(1:10), 
                 "unable to find an inherited method")
    expect_equal(
        TICquantileToQuantileLogRatio(sps_sciex, relativeTo = "previous", 
                                      MSLevel = 1L), 
        c(-0.05853477, -0.58202994, 0.15361253))
    expect_equal(
        TICquantileToQuantileLogRatio(sps_sciex, relativeTo = "Q1", 
                                      MSLevel = 1L), 
        c(-0.05853477, -0.64056471, -0.48695218))
    expect_error(TICquantiletoQuantileLogRatio(sps_sciex, MSLevel = 2L), 
                 "Spectra object does not contain any spectra")
})
## END unit test TICquantiletoQuantileLogRatio ##

## START unit test numberSpectra ##
test_that("numberSpectra", {
    expect_error(numberSpectra(NULL), "unable to find an inherited method")
    expect_error(numberSpectra(1:10), "unable to find an inherited method")
    expect_equal(numberSpectra(sps_sciex, MSLevel = 1), 1862)
    expect_equal(numberSpectra(sps_sciex, MSLevel = 2), 0)
})
## END unit test numberSpectra ##

## START unit test medianPrecursorMZ ##
test_that("medianPrecursorMZ", {
    expect_error(medianPrecursorMZ(NULL), "unable to find an inherited method")
    expect_error(medianPrecursorMZ(1:10), "unable to find an inherited method")
    expect_equal(medianPrecursorMZ(sps_sciex, MSLevel = 1), 501.6874)
    expect_error(medianPrecursorMZ(sps_sciex, MSLevel = 2), 
        "Spectra object does not contain any spectra")
})
## END unit test medianPrecursorMZ ##

## START unit test rtimeIQR ##
test_that("rtimeIQR", {
    expect_error(rtimeIQR(NULL), "unable to find an inherited method")
    expect_error(rtimeIQR(1:10), "unable to find an inherited method")
    expect_equal(rtimeIQR(sps_sciex, MSLevel = 1), 129.875)
    expect_error(rtimeIQR(sps_sciex, MSLevel = 2), 
        "Spectra object does not contain any spectra")
})
## END unit test rtimeIQR ##

## START unit test rtimeIQRrate ##
test_that("rtimeIQRrate", {
    expect_error(rtimeIQRrate(NULL), "unable to find an inherited method")
    expect_error(rtimeIQRrate(1:10), "unable to find an inherited method")
    expect_equal(rtimeIQRrate(sps_sciex, MSLevel = 1), 7.160731)
    expect_error(rtimeIQRrate(sps_sciex, MSLevel = 2), 
        "Spectra object does not contain any spectra")
})
## END unit test rtimeIQRrate ##










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
