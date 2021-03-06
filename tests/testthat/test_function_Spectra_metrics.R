## create test sets for test_function_Spectra_metrics.R
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

## add some fake precursorIntensity
sps_sciex@backend$precursorIntensity <- rpois(n = 1862, lambda = 10000)

## START unit test rtDuration ##
test_that("rtDuration", {
    expect_error(rtDuration(NULL), "unable to find an inherited method")
    expect_error(rtDuration(NULL), "unable to find an inherited method")
    expect_equal(rtDuration(sps_sciex), 259.482)
})
## END unit test rtDuration ##

## START unit test rtOverTicQuantile ##
test_that("rtOverTicQuantile", {
    expect_error(rtOverTicQuantile(NULL), "unable to find an inherited method")
    expect_error(rtOverTicQuantile(NULL), "unable to find an inherited method")
    suppressWarnings(tmp <- rtOverTicQuantile(sps_sciex))
    expect_equal(as.numeric(tmp),
        c(0.0, 0.25, 0.5, 0.75, 1), tolerance = 1e-02)
    expect_equal(names(tmp), c("0%", "25%", "50%", "75%", "100%"))
})
## END unit test rtOverTicQuantile ##

## START unit test rtOverMsQuarters ##
test_that("rtOverMsQuarters", {
    expect_error(rtOverMsQuarters(NULL), "unable to find an inherited method")
    expect_error(rtOverMsQuarters(NULL), "unable to find an inherited method")
    expect_equal(as.numeric(rtOverMsQuarters(sps_sciex, msLevel = 1L)), 
        c(0.25, 0.5, 0.75, 1), tolerance = 1e-02)
    expect_equal(names(rtOverMsQuarters(sps_sciex, msLevel = 1L)), 
        c("Quarter1", "Quarter2", "Quarter3", "Quarter4"), tolerance = 1e-06)
    expect_error(rtOverMsQuarters(sps_sciex, msLevel = 2L), 
        "'spectra' does contain less than four spectra")
})
## END unit test rtOverMsQuarters ##

## START unit test ticQuantileToQuantileLogRatio ##
test_that("ticQuantileToQuantileLogRatio", {
    expect_error(ticQuantileToQuantileLogRatio(NULL), 
        "unable to find an inherited method")
    expect_error(ticQuantileToQuantileLogRatio(NULL), 
        "unable to find an inherited method")
    tmp <- suppressWarnings(ticQuantileToQuantileLogRatio(sps_sciex, 
                relativeTo = "previous", msLevel = 1L))
    expect_equal(as.numeric(tmp), c(-0.05853477, -0.58202994, 0.15361253))
    expect_equal(names(tmp), c("Q2/Q1", "Q3/Q2", "Q4/Q3"))
    tmp <- suppressWarnings(ticQuantileToQuantileLogRatio(sps_sciex, 
                relativeTo = "Q1", msLevel = 1L))
    expect_equal(as.numeric(tmp), c(-0.05853477, -0.64056471, -0.48695218))
    expect_equal(names(tmp), c("Q2/Q1", "Q3/Q1", "Q4/Q1"))
    expect_error(ticQuantileToQuantileLogRatio(sps_sciex, msLevel = 2L), 
                 "'spectra' does not contain any spectra")
})
## END unit test ticQuantiletoQuantileLogRatio ##

## START unit test numberSpectra ##
test_that("numberSpectra", {
    expect_error(numberSpectra(NULL), "unable to find an inherited method")
    expect_error(numberSpectra(1:10), "unable to find an inherited method")
    expect_equal(numberSpectra(sps_sciex, msLevel = 1L), 1862)
    expect_equal(numberSpectra(sps_sciex, msLevel = 2L), 0)
})
## END unit test numberSpectra ##

## START unit test medianPrecursorMz ##
test_that("medianPrecursorMz", {
    expect_error(medianPrecursorMz(NULL), "unable to find an inherited method")
    expect_error(medianPrecursorMz(1:10), "unable to find an inherited method")
    expect_equal(medianPrecursorMz(sps_sciex, msLevel = 1L), 496.4041, 
        tolerance = 1e-06)
    expect_error(medianPrecursorMz(sps_sciex, msLevel = 2L), 
        "'spectra' does not contain any spectra")
})
## END unit test medianPrecursorMz ##

## START unit test rtIqr ##
test_that("rtIqr", {
    expect_error(rtIqr(NULL), "unable to find an inherited method")
    expect_error(rtIqr(1:10), "unable to find an inherited method")
    expect_equal(rtIqr(sps_sciex, msLevel = 1L), 129.875)
    expect_error(rtIqr(sps_sciex, msLevel = 2L), 
        "'spectra' does not contain any spectra")
})
## END unit test rtIqr ##

## START unit test rtIqrRate ##
test_that("rtIqrRate", {
    expect_error(rtIqrRate(NULL), "unable to find an inherited method")
    expect_error(rtIqrRate(1:10), "unable to find an inherited method")
    expect_equal(rtIqrRate(sps_sciex, msLevel = 1L), 7.160731, 
        tolerance = 1e-06)
    expect_error(rtIqrRate(sps_sciex, msLevel = 2L), 
        "'spectra' does not contain any spectra")
})
## END unit test rtIqrRate ##

## START unit test areaUnderTic ##
test_that("areaUnderTic", {
    expect_error(areaUnderTic(NULL), "unable to find an inherited method")
    expect_error(areaUnderTic(1:10), "unable to find an inherited method")
    expect_equal(suppressWarnings(areaUnderTic(sps_sciex, msLevel = 1L)), 
        1273927561)
    expect_error(areaUnderTic(sps_sciex, msLevel = 2L), 
                 "'spectra' does not contain any spectra")
})
## END unit test areaUnderTic ##

## START unit test areaUnderTicRtQuantiles ##
test_that("areaUnderTicRtQuantiles", {
    expect_error(areaUnderTicRtQuantiles(NULL), "unable to find an inherited method")
    expect_error(areaUnderTicRtQuantiles(1:10), "unable to find an inherited method")
    suppressWarnings(tmp <- areaUnderTicRtQuantiles(sps_sciex, msLevel = 1L))
    expect_equal(as.numeric(tmp), c(383935723, 368643879, 245834029, 274788917))
    expect_equal(names(tmp), c("25%", "50%", "75%", "100%"))
    expect_error(areaUnderTicRtQuantiles(sps_sciex, msLevel = 2L), 
        "'spectra' does not contain any spectra")
})
## END unit test areaUnderTicRtQuantiles ##

## START unit test extentIdentifiedPrecursorIntensity ##
test_that("extentIdentifiedPrecursorIntensity", {
    expect_error(extentIdentifiedPrecursorIntensity(NULL), 
        "unable to find an inherited method")
    expect_error(extentIdentifiedPrecursorIntensity(1:10), 
        "unable to find an inherited method")
    expect_equal(extentIdentifiedPrecursorIntensity(sps_sciex, msLevel = 1L), 
        1.034276, tolerance = 1e-06)
    expect_error(extentIdentifiedPrecursorIntensity(sps_sciex, msLevel = 2), 
        "'spectra' does not contain any spectra")
})
## END unit test extentIdentifiedPrecursorIntensity ##

## START unit test medianTicRtIqr ##
test_that("medianTicRtIqr", {
    expect_error(medianTicRtIqr(NULL), "unable to find an inherited method")
    expect_error(medianTicRtIqr(1:10), "unable to find an inherited method")
    expect_equal(suppressWarnings(medianTicRtIqr(sps_sciex, msLevel = 1L)), 
        718615)
    expect_error(medianTicRtIqr(sps_sciex, msLevel = 2L), 
        "'spectra' does not contain any spectra")
})
## END unit test medianTicRtIqr ##

## START unit test medianTicOfRtRange ##
test_that("medianTicOfRtRange", {
    expect_error(medianTicOfRtRange(NULL), "unable to find an inherited method")
    expect_error(medianTicOfRtRange(1:10), "unable to find an inherited method")
    expect_equal(suppressWarnings(medianTicOfRtRange(sps_sciex, msLevel = 1L)), 
        804944)
    expect_error(medianTicOfRtRange(sps_sciex, msLevel = 2L), 
                 "'spectra' does not contain any spectra")
})
## END unit test medianTicOfRtRange ##

## START unit test mzAcquisitionRange ##
test_that("mzAcquisitionRange", {
    expect_error(mzAcquisitionRange(NULL), "unable to find an inherited method")
    expect_error(mzAcquisitionRange(1:10), "unable to find an inherited method")
    tmp <- mzAcquisitionRange(sps_sciex, msLevel = 1L)
    expect_equal(as.numeric(tmp), c(105, 134), tolerance = 1e-06)
    expect_equal(names(tmp), c("min", "max"))
    expect_error(mzAcquisitionRange(sps_sciex, msLevel = 2L), 
        "'spectra' does not contain any spectra")
})
## END unit test mzAcquisitionRange ##

## START unit test rtAcquisitionRange ##
test_that("rtAcquisitionRange", {
    expect_error(rtAcquisitionRange(NULL), "unable to find an inherited method")
    expect_error(rtAcquisitionRange(1:10), "unable to find an inherited method")
    tmp <- rtAcquisitionRange(sps_sciex, msLevel = 1L)
    expect_equal(as.numeric(tmp), c(0.275, 259.757))
    expect_equal(names(tmp), c("min", "max"))
    expect_error(rtAcquisitionRange(sps_sciex, msLevel = 2L), 
        "'spectra' does not contain any spectra")
})
## END unit test rtAcquisitionRange ##

## START unit test precursorIntensityRange ##
test_that("precursorIntensityRange", {
    expect_error(precursorIntensityRange(NULL), "unable to find an inherited method")
    expect_error(precursorIntensityRange(1:10), "unable to find an inherited method")
    tmp <- precursorIntensityRange(sps_sciex, msLevel = 1L)
    expect_equal(as.numeric(tmp), c(9679, 10286))
    expect_equal(names(tmp), c("min", "max"))
    expect_error(precursorIntensityRange(sps_sciex, msLevel = 2L), 
        "'spectra' does not contain any spectra")
})
## END unit test precursorIntensityRange ##

## START unit test precursorIntensityQuartiles ##
test_that("precursorIntensityQuartiles", {
    expect_error(precursorIntensityQuartiles(NULL), "unable to find an inherited method")
    expect_error(precursorIntensityQuartiles(1:10), "unable to find an inherited method")
    tmp <- precursorIntensityQuartiles(sps_sciex, msLevel = 1L)
    expect_equal(as.numeric(tmp), c(9934, 9999, 10067))
    expect_equal(names(tmp), c("25%", "50%", "75%"))
    expect_error(precursorIntensityQuartiles(sps_sciex, msLevel = 2L), 
                 "'spectra' does not contain any spectra")
})
## END unit test precursorIntensityRange ##

## START unit test precursorIntensityMean ##
test_that("precursorIntensityMean", {
    expect_error(precursorIntensityMean(NULL), "unable to find an inherited method")
    expect_error(precursorIntensityMean(1:10), "unable to find an inherited method")
    expect_equal(precursorIntensityMean(sps_sciex, msLevel = 1L), 9999.646,
        tolerance = 1e-06)
    expect_error(precursorIntensityMean(sps_sciex, msLevel = 2L), 
                 "'spectra' does not contain any spectra")
})
## END unit test precursorIntensityMean ##

## START unit test precursorIntensitySd ##
test_that("precursorIntensitySd", {
    expect_error(precursorIntensitySd(NULL), "unable to find an inherited method")
    expect_error(precursorIntensitySd(1:10), "unable to find an inherited method")
    expect_equal(precursorIntensitySd(sps_sciex, msLevel = 1L), 101.0341, 
        tolerance = 1e-06)
    expect_error(precursorIntensitySd(sps_sciex, msLevel = 2L), 
        "'spectra' does not contain any spectra")
})
## END unit test precursorIntensitySd ##

## START unit test msSignal10xChange ##
test_that("msSignal10xChange", {
    expect_error(msSignal10xChange(NULL), "unable to find an inherited method")
    expect_error(msSignal10xChange(1:10), "unable to find an inherited method")
    expect_equal(suppressWarnings(
        msSignal10xChange(sps_sciex, change = "jump", msLevel = 1L)), 0)
    expect_equal(suppressWarnings(
        msSignal10xChange(sps_sciex, change = "fall", msLevel = 1L)), 0)
    expect_error(msSignal10xChange(sps_sciex, change = "jump", msLevel = 2L), 
        "'spectra' does not contain any spectra")
    expect_error(msSignal10xChange(sps_sciex, change = "fall", msLevel = 2L), 
        "'spectra' does not contain any spectra")
})
## END unit test msSignal10xChange ##

## START unit test ratioCharge1over2 ##
test_that("ratioCharge1over2", {
    expect_error(ratioCharge1over2(NULL), "unable to find an inherited method")
    expect_error(ratioCharge1over2(1:10), "unable to find an inherited method")
    expect_equal(ratioCharge1over2(sps_sciex), 0.7 / 0.15, tolerance = 3e-02)
    
    ## do not include charges 1 or 2
    sps_sciex_foo <- sps_sciex
    sps_sciex_foo@backend$precursorCharge <- as.integer(
        sample(x = c(1, 2, 3, 4), size = 1862, 
               replace = TRUE, prob = c(0.0, 0.85, 0.1, 0.05)))
    expect_equal(ratioCharge1over2(sps_sciex_foo), NA)
    sps_sciex_foo <- sps_sciex
    sps_sciex_foo@backend$precursorCharge <- as.integer(
        sample(x = c(1, 2, 3, 4), size = 1862, 
               replace = TRUE, prob = c(0.85, 0.0, 0.1, 0.05)))
    expect_equal(ratioCharge1over2(sps_sciex_foo), NA)
})
## END unit test ratioCharge1over2 ##

## START unit test ratioCharge3over2 ##
test_that("ratioCharge3over2", {
    expect_error(ratioCharge3over2(NULL), "unable to find an inherited method")
    expect_error(ratioCharge3over2(1:10), "unable to find an inherited method")
    expect_equal(ratioCharge3over2(sps_sciex), 0.1 / 0.15, tolerance = 2e-02)
    
    ## do not include charges 2 or 3
    sps_sciex_foo <- sps_sciex
    sps_sciex_foo@backend$precursorCharge <- as.integer(
        sample(x = c(1, 2, 3, 4), size = 1862, 
               replace = TRUE, prob = c(0.7, 0.0, 0.25, 0.05)))
    expect_equal(ratioCharge3over2(sps_sciex_foo), NA)
    sps_sciex_foo <- sps_sciex
    sps_sciex_foo@backend$precursorCharge <- as.integer(
        sample(x = c(1, 2, 3, 4), size = 1862, 
               replace = TRUE, prob = c(0.7, 0.25, 0.0, 0.05)))
    expect_equal(ratioCharge3over2(sps_sciex_foo), NA)
})
## END unit test ratioCharge3over2 ##

## START unit test ratioCharge4over2 ##
test_that("ratioCharge4over2", {
    expect_error(ratioCharge4over2(NULL), "unable to find an inherited method")
    expect_error(ratioCharge4over2(1:10), "unable to find an inherited method")
    expect_equal(ratioCharge4over2(sps_sciex), 0.05 / 0.15, tolerance = 3e-02)
    
    ## do not include charges 2 or 4
    sps_sciex_foo <- sps_sciex
    sps_sciex_foo@backend$precursorCharge <- as.integer(
        sample(x = c(1, 2, 3, 4), size = 1862, 
               replace = TRUE, prob = c(0.7, 0.0, 0.1, 0.20)))
    expect_equal(ratioCharge4over2(sps_sciex_foo), NA)
    sps_sciex_foo <- sps_sciex
    sps_sciex_foo@backend$precursorCharge <- as.integer(
        sample(x = c(1, 2, 3, 4), size = 1862, 
               replace = TRUE, prob = c(0.7, 0.25, 0.15, 0.0)))
    expect_equal(ratioCharge4over2(sps_sciex_foo), NA)
})
## END unit test ratioCharge4over2 ##

test_that(".rtOrderSpectra works", {
    tmp <- sps_sciex[sample(seq_along(sps_sciex), 10)]
    res <- .rtOrderSpectra(tmp)
    expect_true(!is.unsorted(rtime(res)))

    tmp$rtime[4] <- NA
    expect_warning(res <- .rtOrderSpectra(tmp))
    expect_equal(rtime(res), rtime(tmp))
})
