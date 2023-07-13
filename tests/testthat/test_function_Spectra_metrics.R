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

test_that(".rtOrderSpectra works properly.", {
    tmp <- sps_sciex[sample(seq_along(sps_sciex), 10)]
    res <- .rtOrderSpectra(tmp)
    expect_true(!is.unsorted(rtime(res)))
    
    tmp$rtime[4] <- NA
    expect_warning(res <- .rtOrderSpectra(tmp))
    expect_equal(rtime(res), rtime(tmp))
})

## START unit test chromatographyDuration ##
test_that("chromatographyDuration works properly.", {
    expect_error(chromatographyDuration(NULL), "unable to find an inherited method")
    expect_error(chromatographyDuration(NULL), "unable to find an inherited method")
    tmp <- chromatographyDuration(sps_sciex)
    expect_equal(as.numeric(tmp), 259.482)
    
    ## test attributes
    expect_equal(names(attributes(tmp)), "chromatographyDuration")
    expect_equal(attr(tmp, "chromatographyDuration"), "MS:4000053")
})
## END unit test chromatographyDuration ##

## START unit test ticQuartersRtFraction ##
test_that("ticQuartersRtFraction works properly.", {
    expect_error(ticQuartersRtFraction(NULL), "unable to find an inherited method")
    expect_error(ticQuartersRtFraction(NULL), "unable to find an inherited method")
    suppressWarnings(tmp <- ticQuartersRtFraction(sps_sciex))
    expect_equal(as.numeric(tmp),
        c(0.0, 0.2010891, 0.4257983, 0.7247362, 1), tolerance = 1e-02)
    
    ## test attributes 
    expect_equal(names(attributes(tmp)), c("names", "ticQuartersRtFraction"))
    expect_equal(names(tmp), c("0%", "25%", "50%", "75%", "100%"))
    expect_equal(attr(tmp, "ticQuartersRtFraction"), "MS:4000054")
})
## END unit test ticQuartersRtFraction ##

## START unit test rtOverMsQuarters ##
test_that("rtOverMsQuarters works properly.", {
    expect_error(rtOverMsQuarters(NULL), "unable to find an inherited method")
    expect_error(rtOverMsQuarters(NULL), "unable to find an inherited method")
    expect_equal(as.numeric(rtOverMsQuarters(sps_sciex, msLevel = 1L)), 
        c(0.2494778, 0.5000077, 0.7505222, 1), tolerance = 1e-06)
    expect_equal(names(rtOverMsQuarters(sps_sciex, msLevel = 1L)), 
        c("Quarter1", "Quarter2", "Quarter3", "Quarter4"))
    expect_equal(as.numeric(rtOverMsQuarters(sps_sciex[1:3,], msLevel = 1L)), 
        c(NaN, NaN, NaN, NaN), tolerance = 1e-06)
    expect_equal(names(rtOverMsQuarters(sps_sciex[1:3, ], msLevel = 1L)), 
        c("Quarter1", "Quarter2", "Quarter3", "Quarter4"))
    expect_equal(as.numeric(rtOverMsQuarters(sps_sciex[1:4,], msLevel = 1L)), 
        c(0, 0.333333, 0.6666667, 1), tolerance = 1e-06)
    expect_equal(names(rtOverMsQuarters(sps_sciex[1:4, ], msLevel = 1L)), 
        c("Quarter1", "Quarter2", "Quarter3", "Quarter4"))
    
    ## test attributes
    tmp <- rtOverMsQuarters(sps_sciex[1:4, ], msLevel = 1L)
    expect_equal(names(attributes(tmp)), 
        c("names", "rtOverMsQuarters"))
    expect_equal(attr(tmp, "rtOverMsQuarters"), "MS:4000055")
    tmp <- rtOverMsQuarters(sps_sciex[1:4, ], msLevel = 2L)
    expect_equal(names(attributes(tmp)), 
                 c("names", "rtOverMsQuarters"))
    expect_equal(attr(tmp, "rtOverMsQuarters"), "MS:4000056")
    tmp <- rtOverMsQuarters(sps_sciex[1:4, ], msLevel = 3L)
    expect_equal(names(attributes(tmp)), c("names"))
    
})
## END unit test rtOverMsQuarters ##

## START unit test ticQuartileToQuartileLogRatio ##
test_that("ticQuartileToQuartileLogRatio works properly.", {
    expect_error(ticQuartileToQuartileLogRatio(NULL, relativeTo = "Q1"), 
        "unable to find an inherited method")
    expect_error(ticQuartileToQuartileLogRatio(NULL, relativeTo = "previous"), 
        "unable to find an inherited method")
    
    ## mode = "TIC_change"
    tmp <- suppressWarnings(ticQuartileToQuartileLogRatio(sps_sciex, 
        relativeTo = "previous", mode = "TIC_change", msLevel = 1L))
    expect_equal(as.numeric(tmp), c(-6.047274, NaN, 2.505220), 
        tolerance = 1e-06)
    expect_equal(names(tmp), c("Q2/Q1", "Q3/Q2", "Q4/Q3"))
    tmp <- suppressWarnings(ticQuartileToQuartileLogRatio(sps_sciex, 
        relativeTo = "Q1", mode = "TIC_change", msLevel = 1L))
    expect_equal(as.numeric(tmp), c(-6.047274, NaN, NaN), tolerance = 1e-06)
    expect_equal(names(tmp), c("Q2/Q1", "Q3/Q1", "Q4/Q1"))
    
    ## mode = "TIC"
    tmp <- suppressWarnings(ticQuartileToQuartileLogRatio(sps_sciex, 
        relativeTo = "previous", mode = "TIC", msLevel = 1L))
    expect_equal(as.numeric(tmp), c(0.4759104, 0.2566362, 0.6382157),
        tolerance = 1e-06)
    expect_equal(names(tmp), c("Q2/Q1", "Q3/Q2", "Q4/Q3"))
    tmp <- suppressWarnings(ticQuartileToQuartileLogRatio(sps_sciex, 
        relativeTo = "Q1", mode = "TIC", msLevel = 1L))
    expect_equal(as.numeric(tmp), c(0.4759104, 0.7325466, 1.3707623), tolerance = 1e-06)
    expect_equal(names(tmp), c("Q2/Q1", "Q3/Q1", "Q4/Q1"))
    
    ## tests when Spectra object has length 0
    tmp <- suppressWarnings(ticQuartileToQuartileLogRatio(sps_sciex, 
        relativeTo = "previous", mode = "TIC", msLevel = 2L))
    expect_equal(as.numeric(tmp), c(NaN, NaN, NaN))
    expect_equal(names(tmp), c("Q2/Q1", "Q3/Q2", "Q4/Q3"))
    tmp <- suppressWarnings(ticQuartileToQuartileLogRatio(sps_sciex, 
        relativeTo = "Q1", mode = "TIC", msLevel = 2L))
    expect_equal(as.numeric(tmp), c(NaN, NaN, NaN))
    
    ## test attributes
    tmp <- suppressWarnings(ticQuartileToQuartileLogRatio(sps_sciex[1:4, ], 
        relativeTo = "Q1", mode = "TIC", msLevel = 1L))
    expect_equal(names(attributes(tmp)), c("names"))
    expect_equal(names(tmp), c("Q2/Q1", "Q3/Q1", "Q4/Q1"))
    tmp <- suppressWarnings(ticQuartileToQuartileLogRatio(sps_sciex[1:4, ], 
        relativeTo = "previous", mode = "TIC", msLevel = 1L))
    expect_equal(names(attributes(tmp)), 
        c("names", "ticQuartileToQuartileLogRatio"))
    expect_equal(names(tmp), c("Q2/Q1", "Q3/Q2", "Q4/Q3"))
    expect_equal(attr(tmp, "ticQuartileToQuartileLogRatio"), "MS:4000058")
    tmp <- suppressWarnings(ticQuartileToQuartileLogRatio(sps_sciex[1:4, ], 
        relativeTo = "Q1", mode = "TIC", msLevel = 2L))
    expect_equal(names(attributes(tmp)), c("names"))
    tmp <- suppressWarnings(ticQuartileToQuartileLogRatio(sps_sciex[1:4, ], 
        relativeTo = "previous", mode = "TIC", msLevel = 2L))
    expect_equal(names(attributes(tmp)), c("names"))
    
})
## END unit test ticQuartileToQuartileLogRatio ##

## START unit test numberSpectra ##
test_that("numberSpectra works properly.", {
    expect_error(numberSpectra(NULL), "unable to find an inherited method")
    expect_error(numberSpectra(1:10), "unable to find an inherited method")
    expect_equal(numberSpectra(sps_sciex, msLevel = 1L), 1862)
    expect_equal(numberSpectra(sps_sciex, msLevel = 2L), 0)
    
    ## test attributes
    tmp <- numberSpectra(sps_sciex, msLevel = 1L)
    expect_equal(names(attributes(tmp)), c("numberSpectra"))
    expect_equal(attr(tmp, "numberSpectra"), "MS:4000059")
    tmp <- numberSpectra(sps_sciex, msLevel = 2L)
    expect_equal(names(attributes(tmp)), c("numberSpectra"))
    expect_equal(attr(tmp, "numberSpectra"), "MS:4000060")
    tmp <- numberSpectra(sps_sciex, msLevel = 3L)
    expect_equal(names(attributes(tmp)), NULL)
    
})
## END unit test numberSpectra ##

## START unit test medianPrecursorMz ##
test_that("medianPrecursorMz works properly.", {
    expect_error(medianPrecursorMz(NULL), "unable to find an inherited method")
    expect_error(medianPrecursorMz(1:10), "unable to find an inherited method")
    expect_equal(
        medianPrecursorMz(sps_sciex, msLevel = 1L, 
            identificationLevel = "all"), 
        496.4041, tolerance = 1e-06)
    expect_equal(
        medianPrecursorMz(sps_sciex, msLevel = 1L, 
            identificationLevel = "identified"), 
        496.4041, tolerance = 1e-06)
    expect_equal(
        medianPrecursorMz(sps_sciex, msLevel = 1L, 
            identificationLevel = "unidentified"), 
        496.4041, tolerance = 1e-06)
    expect_equal(medianPrecursorMz(sps_sciex, msLevel = 2L), NaN)
    
    ## test attributes
    tmp <- medianPrecursorMz(sps_sciex, msLevel = 1L, 
        identificationLevel = "all")
    expect_equal(names(attributes(tmp)), NULL)
    tmp <- medianPrecursorMz(sps_sciex, msLevel = 1L, 
        identificationLevel = "identified")
    expect_equal(names(attributes(tmp)), "medianPrecursorMz")
    expect_equal(attr(tmp, "medianPrecursorMz"), "MS:4000152")
    tmp <- medianPrecursorMz(sps_sciex, msLevel = 1L, 
        identificationLevel = "unidentified")
    expect_equal(names(attributes(tmp)), NULL)
    tmp <- medianPrecursorMz(sps_sciex, msLevel = 2L, 
        identificationLevel = "all")
    expect_equal(names(attributes(tmp)), NULL)
    tmp <- medianPrecursorMz(sps_sciex, msLevel = 2L, 
        identificationLevel = "identified")
    expect_equal(names(attributes(tmp)), NULL)
    tmp <- medianPrecursorMz(sps_sciex, msLevel = 2L, 
        identificationLevel = "unidentified")
    expect_equal(names(attributes(tmp)), NULL)
})
## END unit test medianPrecursorMz ##

## START unit test rtIqr ##
test_that("rtIqr works properly.", {
    expect_error(rtIqr(NULL), "unable to find an inherited method")
    expect_error(rtIqr(1:10), "unable to find an inherited method")
    expect_equal(rtIqr(sps_sciex, msLevel = 1L), 129.875)
    expect_equal(rtIqr(sps_sciex, msLevel = 2L), NaN)
    
    ## test attributes
    tmp <- rtIqr(sps_sciex, msLevel = 1L, identificationLevel = "all")
    expect_equal(names(attributes(tmp)), NULL)
    tmp <- rtIqr(sps_sciex, msLevel = 1L, identificationLevel = "identified")
    expect_equal(names(attributes(tmp)), "rtIqr")
    expect_equal(attr(tmp, "rtIqr"), "MS:4000153")
    tmp <- rtIqr(sps_sciex, msLevel = 1L, identificationLevel = "unidentified")
    expect_equal(names(attributes(tmp)), NULL)
})
## END unit test rtIqr ##

## START unit test rtIqrRate ##
test_that("rtIqrRate works properly.", {
    expect_error(rtIqrRate(NULL), "unable to find an inherited method")
    expect_error(rtIqrRate(1:10), "unable to find an inherited method")
    expect_equal(
        as.numeric(rtIqrRate(sps_sciex, msLevel = 1L, 
            identificationLevel = "all")), 
        7.160731, tolerance = 1e-06)
    expect_equal(
        as.numeric(rtIqrRate(sps_sciex, msLevel = 1L, 
                  identificationLevel = "identified")), 
        7.160731, tolerance = 1e-06)
    expect_equal(
        as.numeric(rtIqrRate(sps_sciex, msLevel = 1L, 
                  identificationLevel = "unidentified")), 
        7.160731, tolerance = 1e-06)
    expect_equal(rtIqrRate(sps_sciex, msLevel = 2L), NaN)
    
    ## test attributes
    tmp <- rtIqrRate(sps_sciex, msLevel = 1L, identificationLevel = "all")
    expect_equal(names(attributes(tmp)), NULL)
    tmp <- rtIqrRate(sps_sciex, msLevel = 1L, identificationLevel = "identified")
    expect_equal(names(attributes(tmp)), "rtIqrRate")
    expect_equal(attr(tmp, "rtIqrRate"), "MS:4000154")
    tmp <- rtIqrRate(sps_sciex, msLevel = 1L, identificationLevel = "unidentified")
    expect_equal(names(attributes(tmp)), NULL)
})
## END unit test rtIqrRate ##

## START unit test areaUnderTic ##
test_that("areaUnderTic works properly.", {
    expect_error(areaUnderTic(NULL), "unable to find an inherited method")
    expect_error(areaUnderTic(1:10), "unable to find an inherited method")
    expect_equal(suppressWarnings(areaUnderTic(sps_sciex, msLevel = 1L)), 
        1273927561)
    expect_equal(areaUnderTic(sps_sciex, msLevel = 2L), NaN)
    
    ## test attributes
    tmp <- areaUnderTic(sps_sciex, msLevel = 1L)
    expect_equal(names(attributes(tmp)), "areaUnderTic")
    expect_equal(attr(tmp, "areaUnderTic"), "MS:4000155")
    tmp <- areaUnderTic(sps_sciex, msLevel = 2L)
    expect_equal(names(attributes(tmp)), "areaUnderTic")
    expect_equal(attr(tmp, "areaUnderTic"), "MS:4000155")
})
## END unit test areaUnderTic ##

## START unit test areaUnderTicRtQuantiles ##
test_that("areaUnderTicRtQuantiles works properly.", {
    expect_error(areaUnderTicRtQuantiles(NULL), "unable to find an inherited method")
    expect_error(areaUnderTicRtQuantiles(1:10), "unable to find an inherited method")
    suppressWarnings(tmp <- areaUnderTicRtQuantiles(sps_sciex, msLevel = 1L))
    expect_equal(as.numeric(tmp), c(383935723, 368643879, 245834029, 274788917))
    expect_equal(names(tmp), c("25%", "50%", "75%", "100%"))
    
    tmp <- areaUnderTicRtQuantiles(sps_sciex, msLevel = 2L)
    expect_equal(as.numeric(tmp), c(NaN, NaN, NaN, NaN))
    expect_equal(names(tmp), c("25%", "50%", "75%", "100%"))
    
    ## test attributes
    tmp <- suppressWarnings(areaUnderTicRtQuantiles(sps_sciex, msLevel = 1L))
    expect_equal(names(attributes(tmp)), c("names", "areaUnderTicRtQuantiles"))
    expect_equal(names(tmp), c("25%", "50%", "75%", "100%"))
    expect_equal(attr(tmp, "areaUnderTicRtQuantiles"), "MS:4000156")
    tmp <- suppressWarnings(areaUnderTicRtQuantiles(sps_sciex, msLevel = 2L))
    expect_equal(names(attributes(tmp)), c("names", "areaUnderTicRtQuantiles"))
    expect_equal(names(tmp), c("25%", "50%", "75%", "100%"))
    expect_equal(attr(tmp, "areaUnderTicRtQuantiles"), "MS:4000156")
})
## END unit test areaUnderTicRtQuantiles ##

## START unit test extentIdentifiedPrecursorIntensity ##
test_that("extentIdentifiedPrecursorIntensity works properly.", {
    expect_error(extentIdentifiedPrecursorIntensity(NULL), 
        "unable to find an inherited method")
    expect_error(extentIdentifiedPrecursorIntensity(1:10), 
        "unable to find an inherited method")
    expect_equal(
        extentIdentifiedPrecursorIntensity(sps_sciex, msLevel = 1L,
            identificationLevel = "all"), 
        1.034276, tolerance = 1e-06)
    expect_equal(
        extentIdentifiedPrecursorIntensity(sps_sciex, msLevel = 1L,
            identificationLevel = "identified"), 
        1.034276, tolerance = 1e-06)
    expect_equal(
        extentIdentifiedPrecursorIntensity(sps_sciex, msLevel = 1L,
            identificationLevel = "unidentified"), 
        1.034276, tolerance = 1e-06)
    expect_equal(
        extentIdentifiedPrecursorIntensity(sps_sciex, msLevel = 2L,
            identificationLevel = "all"), 
        NaN)
    expect_equal(
        extentIdentifiedPrecursorIntensity(sps_sciex, msLevel = 2L,
            identificationLevel = "identified"), 
        NaN)
    expect_equal(
        extentIdentifiedPrecursorIntensity(sps_sciex, msLevel = 2L,
            identificationLevel = "unidentified"), 
        NaN)
    
    ## test attributes
    tmp <- extentIdentifiedPrecursorIntensity(sps_sciex, msLevel = 1L, 
        identificationLevel = "all")
    expect_equal(names(attributes(tmp)), NULL)
    tmp <- extentIdentifiedPrecursorIntensity(sps_sciex, msLevel = 1L, 
        identificationLevel = "identified")
    expect_equal(names(attributes(tmp)), c("extentIdentifiedPrecursorIntensity"))
    expect_equal(attr(tmp, "extentIdentifiedPrecursorIntensity"), "MS:4000157")
    tmp <- extentIdentifiedPrecursorIntensity(sps_sciex, msLevel = 1L, 
        identificationLevel = "unidentified")
    expect_equal(names(attributes(tmp)), NULL)
    tmp <- extentIdentifiedPrecursorIntensity(sps_sciex, msLevel = 2L, 
        identificationLevel = "all")
    expect_equal(names(attributes(tmp)), NULL)
    tmp <- extentIdentifiedPrecursorIntensity(sps_sciex, msLevel = 2L, 
        identificationLevel = "identified")
    expect_equal(names(attributes(tmp)), c("extentIdentifiedPrecursorIntensity"))
    expect_equal(attr(tmp, "extentIdentifiedPrecursorIntensity"), "MS:4000157")
    tmp <- extentIdentifiedPrecursorIntensity(sps_sciex, msLevel = 2L, 
        identificationLevel = "unidentified")
    expect_equal(names(attributes(tmp)), NULL)
})
## END unit test extentIdentifiedPrecursorIntensity ##

## START unit test medianTicRtIqr ##
test_that("medianTicRtIqr works properly.", {
    expect_error(medianTicRtIqr(NULL), "unable to find an inherited method")
    expect_error(medianTicRtIqr(1:10), "unable to find an inherited method")
    expect_equal(suppressWarnings(medianTicRtIqr(sps_sciex, msLevel = 1L)), 
        718615)
    expect_equal(medianTicRtIqr(sps_sciex, msLevel = 2L), NaN)
    
    ## test attributes
    tmp <- suppressWarnings(medianTicRtIqr(sps_sciex[1:10,], msLevel = 1L, 
        identificationLevel = "all"))
    expect_equal(names(attributes(tmp)), NULL)
    tmp <- suppressWarnings(medianTicRtIqr(sps_sciex[1:10,], msLevel = 1L, 
        identificationLevel = "identified"))
    expect_equal(names(attributes(tmp)), c("medianTicRtIqr"))
    expect_equal(attr(tmp, "medianTicRtIqr"), "MS:4000158")
    tmp <- suppressWarnings(medianTicRtIqr(sps_sciex[1:10,], msLevel = 1L, 
        identificationLevel = "unidentified"))
    expect_equal(names(attributes(tmp)), NULL)
    tmp <- suppressWarnings(medianTicRtIqr(sps_sciex[1:10,], msLevel = 2L, 
        identificationLevel = "all"))
    expect_equal(names(attributes(tmp)), NULL)
    tmp <- suppressWarnings(medianTicRtIqr(sps_sciex[1:10,], msLevel = 2L, 
        identificationLevel = "identified"))
    expect_equal(names(attributes(tmp)), c("medianTicRtIqr"))
    expect_equal(attr(tmp, "medianTicRtIqr"), "MS:4000158")
    tmp <- suppressWarnings(medianTicRtIqr(sps_sciex[1:10,], msLevel = 2L, 
        identificationLevel = "unidentified"))
    expect_equal(names(attributes(tmp)), NULL)
})
## END unit test medianTicRtIqr ##

## START unit test medianTicOfRtRange ##
test_that("medianTicOfRtRange works properly.", {
    expect_error(medianTicOfRtRange(NULL), "unable to find an inherited method")
    expect_error(medianTicOfRtRange(1:10), "unable to find an inherited method")
    expect_equal(suppressWarnings(medianTicOfRtRange(sps_sciex, msLevel = 1L)), 
        804944)
    expect_equal(medianTicOfRtRange(sps_sciex, msLevel = 2L), NaN)
    
    ## test attributes
    tmp <- suppressWarnings(medianTicOfRtRange(sps_sciex[1:10,], msLevel = 1L, 
        identificationLevel = "all"))
    expect_equal(names(attributes(tmp)), NULL)
    tmp <- suppressWarnings(medianTicOfRtRange(sps_sciex[1:10,], msLevel = 1L, 
        identificationLevel = "identified"))
    expect_equal(names(attributes(tmp)), c("medianTicOfRtRange"))
    expect_equal(attr(tmp, "medianTicOfRtRange"), "MS:4000159")
    tmp <- suppressWarnings(medianTicOfRtRange(sps_sciex[1:10,], msLevel = 1L, 
        identificationLevel = "unidentified"))
    expect_equal(names(attributes(tmp)), NULL)
    tmp <- suppressWarnings(medianTicOfRtRange(sps_sciex[1:10,], msLevel = 2L, 
        identificationLevel = "all"))
    expect_equal(names(attributes(tmp)), NULL)
    tmp <- suppressWarnings(medianTicOfRtRange(sps_sciex[1:10,], msLevel = 2L, 
        identificationLevel = "identified"))
    expect_equal(names(attributes(tmp)), c("medianTicOfRtRange"))
    expect_equal(attr(tmp, "medianTicOfRtRange"), "MS:4000159")
    tmp <- suppressWarnings(medianTicOfRtRange(sps_sciex[1:10,], msLevel = 2L, 
        identificationLevel = "unidentified"))
    expect_equal(names(attributes(tmp)), NULL)
})
## END unit test medianTicOfRtRange ##

## START unit test mzAcquisitionRange ##
test_that("mzAcquisitionRange works properly.", {
    expect_error(mzAcquisitionRange(NULL), "unable to find an inherited method")
    expect_error(mzAcquisitionRange(1:10), "unable to find an inherited method")
    tmp <- mzAcquisitionRange(sps_sciex, msLevel = 1L)
    expect_equal(as.numeric(tmp), c(105, 134), tolerance = 1e-06)
    expect_equal(names(tmp), c("min", "max"))
    tmp <- mzAcquisitionRange(sps_sciex, msLevel = 2L)
    expect_equal(as.numeric(tmp), c(NaN, NaN))
    expect_equal(names(tmp), c("min", "max"))
    
    ## test attributes
    tmp <- suppressWarnings(mzAcquisitionRange(sps_sciex[1:10,], msLevel = 1L))
    expect_equal(names(attributes(tmp)), c("names", "mzAcquisitionRange"))
    expect_equal(names(tmp), c("min", "max"))
    expect_equal(attr(tmp, "mzAcquisitionRange"), "MS:4000069")
    tmp <- suppressWarnings(mzAcquisitionRange(sps_sciex[1:10,], msLevel = 2L))
    expect_equal(names(attributes(tmp)), c("names", "mzAcquisitionRange"))
    expect_equal(names(tmp), c("min", "max"))
    expect_equal(attr(tmp, "mzAcquisitionRange"), "MS:4000069")
})
## END unit test mzAcquisitionRange ##

## START unit test rtAcquisitionRange ##
test_that("rtAcquisitionRange works properly.", {
    expect_error(rtAcquisitionRange(NULL), "unable to find an inherited method")
    expect_error(rtAcquisitionRange(1:10), "unable to find an inherited method")
    tmp <- rtAcquisitionRange(sps_sciex, msLevel = 1L)
    expect_equal(as.numeric(tmp), c(0.275, 259.757))
    expect_equal(names(tmp), c("min", "max"))
    tmp <- rtAcquisitionRange(sps_sciex, msLevel = 2L)
    expect_equal(as.numeric(tmp), c(NaN, NaN))
    expect_equal(names(tmp), c("min", "max"))
    
    ## test attributes
    tmp <- suppressWarnings(rtAcquisitionRange(sps_sciex[1:10,], msLevel = 1L))
    expect_equal(names(attributes(tmp)), c("names", "rtAcquisitionRange"))
    expect_equal(names(tmp), c("min", "max"))
    expect_equal(attr(tmp, "rtAcquisitionRange"), "MS:4000070")
    tmp <- suppressWarnings(rtAcquisitionRange(sps_sciex[1:10,], msLevel = 2L))
    expect_equal(names(attributes(tmp)), c("names", "rtAcquisitionRange"))
    expect_equal(names(tmp), c("min", "max"))
    expect_equal(attr(tmp, "rtAcquisitionRange"), "MS:4000070")
})
## END unit test rtAcquisitionRange ##

## START unit test precursorIntensityRange ##
test_that("precursorIntensityRange works properly.", {
    expect_error(precursorIntensityRange(NULL), "unable to find an inherited method")
    expect_error(precursorIntensityRange(1:10), "unable to find an inherited method")
    tmp <- precursorIntensityRange(sps_sciex, msLevel = 1L)
    expect_equal(as.numeric(tmp), c(9679, 10286))
    expect_equal(names(tmp), c("min", "max"))
    tmp <- precursorIntensityRange(sps_sciex, msLevel = 2L)
    expect_equal(as.numeric(tmp), c(NaN, NaN))
    expect_equal(names(tmp), c("min", "max"))
    
    ## test attributes
    tmp <- suppressWarnings(precursorIntensityRange(sps_sciex[1:10,], msLevel = 1L))
    expect_equal(names(attributes(tmp)), c("names", "precursorIntensityRange"))
    expect_equal(names(tmp), c("min", "max"))
    expect_equal(attr(tmp, "precursorIntensityRange"), "MS:4000070")
    tmp <- suppressWarnings(precursorIntensityRange(sps_sciex[1:10,], msLevel = 2L))
    expect_equal(names(attributes(tmp)), c("names", "precursorIntensityRange"))
    expect_equal(names(tmp), c("min", "max"))
    expect_equal(attr(tmp, "precursorIntensityRange"), "MS:4000070")
})
## END unit test precursorIntensityRange ##

## START unit test precursorIntensityQuartiles ##
test_that("precursorIntensityQuartiles works properly.", {
    expect_error(precursorIntensityQuartiles(NULL), "unable to find an inherited method")
    expect_error(precursorIntensityQuartiles(1:10), "unable to find an inherited method")
    tmp <- precursorIntensityQuartiles(sps_sciex, msLevel = 1L, 
        identificationLevel = "all")
    expect_equal(as.numeric(tmp), c(9934, 9999, 10067))
    expect_equal(names(tmp), c("Q1", "Q2", "Q3"))
    tmp <- precursorIntensityQuartiles(sps_sciex, msLevel = 1L, 
        identificationLevel = "identified")
    expect_equal(as.numeric(tmp), c(9934, 9999, 10067))
    expect_equal(names(tmp), c("Q1", "Q2", "Q3"))
    tmp <- precursorIntensityQuartiles(sps_sciex, msLevel = 1L, 
        identificationLevel = "unidentified")
    expect_equal(as.numeric(tmp), c(9934, 9999, 10067))
    expect_equal(names(tmp), c("Q1", "Q2", "Q3"))
    tmp <- precursorIntensityQuartiles(sps_sciex, msLevel = 2L,
        identificationLevel = "all")
    expect_equal(as.numeric(tmp), c(NaN, NaN, NaN))
    expect_equal(names(tmp), c("Q1", "Q2", "Q3"))
    tmp <- precursorIntensityQuartiles(sps_sciex, msLevel = 2L,
        identificationLevel = "identified")
    expect_equal(as.numeric(tmp), c(NaN, NaN, NaN))
    expect_equal(names(tmp), c("Q1", "Q2", "Q3"))
    tmp <- precursorIntensityQuartiles(sps_sciex, msLevel = 2L,
        identificationLevel = "unidentified")
    expect_equal(as.numeric(tmp), c(NaN, NaN, NaN))
    expect_equal(names(tmp), c("Q1", "Q2", "Q3"))
    
    ## test attributes
    tmp <- suppressWarnings(precursorIntensityQuartiles(sps_sciex[1:10,], 
        msLevel = 1L, identificationLevel = "all"))
    expect_equal(names(attributes(tmp)), 
        c("names", "precursorIntensityQuartiles"))
    expect_equal(names(tmp), c("Q1", "Q2", "Q3"))
    expect_equal(attr(tmp, "precursorIntensityQuartiles"), "MS:4000116")
    tmp <- suppressWarnings(precursorIntensityQuartiles(sps_sciex[1:10,], 
        msLevel = 1L, identificationLevel = "identified"))
    expect_equal(names(attributes(tmp)), 
        c("names", "precursorIntensityQuartiles"))
    expect_equal(names(tmp), c("Q1", "Q2", "Q3"))
    expect_equal(attr(tmp, "precursorIntensityQuartiles"), "MS:4000161")
    tmp <- suppressWarnings(precursorIntensityQuartiles(sps_sciex[1:10,], 
        msLevel = 1L, identificationLevel = "unidentified"))
    expect_equal(names(attributes(tmp)), 
        c("names", "precursorIntensityQuartiles"))
    expect_equal(names(tmp), c("Q1", "Q2", "Q3"))
    expect_equal(attr(tmp, "precursorIntensityQuartiles"), "MS:4000162")
    tmp <- suppressWarnings(precursorIntensityQuartiles(sps_sciex[1:10,], 
        msLevel = 2L, identificationLevel = "all"))
    expect_equal(names(attributes(tmp)), 
        c("names", "precursorIntensityQuartiles"))
    expect_equal(names(tmp), c("Q1", "Q2", "Q3"))
    expect_equal(attr(tmp, "precursorIntensityQuartiles"), "MS:4000116")
    tmp <- suppressWarnings(precursorIntensityQuartiles(sps_sciex[1:10,], 
        msLevel = 2L, identificationLevel = "identified"))
    expect_equal(names(attributes(tmp)), 
        c("names", "precursorIntensityQuartiles"))
    expect_equal(names(tmp), c("Q1", "Q2", "Q3"))
    expect_equal(attr(tmp, "precursorIntensityQuartiles"), "MS:4000161")
    tmp <- suppressWarnings(precursorIntensityQuartiles(sps_sciex[1:10,], 
        msLevel = 2L, identificationLevel = "unidentified"))
    expect_equal(names(attributes(tmp)), 
        c("names", "precursorIntensityQuartiles"))
    expect_equal(names(tmp), c("Q1", "Q2", "Q3"))
    expect_equal(attr(tmp, "precursorIntensityQuartiles"), "MS:4000162")
})
## END unit test precursorIntensityQuartiles ##

## START unit test precursorIntensityMean ##
test_that("precursorIntensityMean works properly.", {
    expect_error(precursorIntensityMean(NULL), "unable to find an inherited method")
    expect_error(precursorIntensityMean(1:10), "unable to find an inherited method")
    expect_equal(precursorIntensityMean(sps_sciex, msLevel = 1L, 
            identificationLevel = "all"), 
        9999.646, tolerance = 1e-06)
    expect_equal(precursorIntensityMean(sps_sciex, msLevel = 1L, 
            identificationLevel = "identified"), 
        9999.646, tolerance = 1e-06)
    expect_equal(precursorIntensityMean(sps_sciex, msLevel = 1L, 
            identificationLevel = "unidentified"), 
        9999.646, tolerance = 1e-06)
    expect_equal(precursorIntensityMean(sps_sciex, msLevel = 2L, 
            identificationLevel = "all"), NaN)
    expect_equal(precursorIntensityMean(sps_sciex, msLevel = 2L, 
            identificationLevel = "identified"), NaN)
    expect_equal(precursorIntensityMean(sps_sciex, msLevel = 2L, 
            identificationLevel = "unidentified"), NaN)
    
    ## test attributes
    tmp <- suppressWarnings(precursorIntensityMean(sps_sciex[1:10,], 
        msLevel = 1L, identificationLevel = "all"))
    expect_equal(names(attributes(tmp)), c("precursorIntensityMean"))
    expect_equal(attr(tmp, "precursorIntensityMean"), "MS:4000117")
    tmp <- suppressWarnings(precursorIntensityMean(sps_sciex[1:10,], 
        msLevel = 1L, identificationLevel = "identified"))
    expect_equal(names(attributes(tmp)), c("precursorIntensityMean"))
    expect_equal(attr(tmp, "precursorIntensityMean"), "MS:4000163")
    tmp <- suppressWarnings(precursorIntensityMean(sps_sciex[1:10,], 
        msLevel = 1L, identificationLevel = "unidentified"))
    expect_equal(names(attributes(tmp)), c("precursorIntensityMean"))
    expect_equal(attr(tmp, "precursorIntensityMean"), "MS:4000164")
    tmp <- suppressWarnings(precursorIntensityMean(sps_sciex[1:10,], 
        msLevel = 2L, identificationLevel = "all"))
    expect_equal(names(attributes(tmp)), c("precursorIntensityMean"))
    expect_equal(attr(tmp, "precursorIntensityMean"), "MS:4000117")
    tmp <- suppressWarnings(precursorIntensityMean(sps_sciex[1:10,], 
        msLevel = 2L, identificationLevel = "identified"))
    expect_equal(names(attributes(tmp)), c("precursorIntensityMean"))
    expect_equal(attr(tmp, "precursorIntensityMean"), "MS:4000163")
    tmp <- suppressWarnings(precursorIntensityMean(sps_sciex[1:10,], 
        msLevel = 2L, identificationLevel = "unidentified"))
    expect_equal(names(attributes(tmp)), c("precursorIntensityMean"))
    expect_equal(attr(tmp, "precursorIntensityMean"), "MS:4000164")
})
## END unit test precursorIntensityMean ##

## START unit test precursorIntensitySd ##
test_that("precursorIntensitySd works properly.", {
    expect_error(precursorIntensitySd(NULL), "unable to find an inherited method")
    expect_error(precursorIntensitySd(1:10), "unable to find an inherited method")
    expect_equal(as.numeric(precursorIntensitySd(sps_sciex, msLevel = 1L, 
        identificationLevel = "all")), 101.0341, tolerance = 1e-06)
    expect_equal(as.numeric(precursorIntensitySd(sps_sciex, msLevel = 1L, 
        identificationLevel = "identified")), 101.0341, tolerance = 1e-06)
    expect_equal(as.numeric(precursorIntensitySd(sps_sciex, msLevel = 1L, 
        identificationLevel = "unidentified")), 101.0341, tolerance = 1e-06)
    expect_equal(as.numeric(precursorIntensitySd(sps_sciex, msLevel = 2L, 
        identificationLevel = "all")), NaN)
    expect_equal(as.numeric(precursorIntensitySd(sps_sciex, msLevel = 2L, 
        identificationLevel = "identified")), NaN)
    expect_equal(as.numeric(precursorIntensitySd(sps_sciex, msLevel = 2L, 
        identificationLevel = "unidentified")), NaN)
    
    ## test attributes
    tmp <- suppressWarnings(precursorIntensitySd(sps_sciex[1:10,], 
        msLevel = 1L, identificationLevel = "all"))
    expect_equal(names(attributes(tmp)), c("precursorIntensitySd"))
    expect_equal(attr(tmp, "precursorIntensitySd"), "MS:4000118")
    tmp <- suppressWarnings(precursorIntensitySd(sps_sciex[1:10,], 
        msLevel = 1L, identificationLevel = "identified"))
    expect_equal(names(attributes(tmp)), c("precursorIntensitySd"))
    expect_equal(attr(tmp, "precursorIntensitySd"), "MS:4000165")
    tmp <- suppressWarnings(precursorIntensitySd(sps_sciex[1:10,], 
        msLevel = 1L, identificationLevel = "unidentified"))
    expect_equal(names(attributes(tmp)), c("precursorIntensitySd"))
    expect_equal(attr(tmp, "precursorIntensitySd"), "MS:4000166")
    tmp <- suppressWarnings(precursorIntensitySd(sps_sciex[1:10,], 
        msLevel = 2L, identificationLevel = "all"))
    expect_equal(names(attributes(tmp)), c("precursorIntensitySd"))
    expect_equal(attr(tmp, "precursorIntensitySd"), "MS:4000118")
    tmp <- suppressWarnings(precursorIntensitySd(sps_sciex[1:10,], 
        msLevel = 2L, identificationLevel = "identified"))
    expect_equal(names(attributes(tmp)), c("precursorIntensitySd"))
    expect_equal(attr(tmp, "precursorIntensitySd"), "MS:4000165")
    tmp <- suppressWarnings(precursorIntensitySd(sps_sciex[1:10,], 
        msLevel = 2L, identificationLevel = "unidentified"))
    expect_equal(names(attributes(tmp)), c("precursorIntensitySd"))
    expect_equal(attr(tmp, "precursorIntensitySd"), "MS:4000166")
})
## END unit test precursorIntensitySd ##

## START unit test msSignal10xChange ##
test_that("msSignal10xChange works properly.", {
    expect_error(msSignal10xChange(NULL), "unable to find an inherited method")
    expect_error(msSignal10xChange(1:10), "unable to find an inherited method")
    expect_equal(suppressWarnings(
        msSignal10xChange(sps_sciex, change = "jump", msLevel = 1L)), 0)
    expect_equal(suppressWarnings(
        msSignal10xChange(sps_sciex, change = "fall", msLevel = 1L)), 0)
    expect_equal(msSignal10xChange(sps_sciex, change = "jump", msLevel = 2L), 
        NaN)
    expect_equal(msSignal10xChange(sps_sciex, change = "fall", msLevel = 2L), 
        NaN)
    
    ## test attributes
    tmp <- suppressWarnings(
        msSignal10xChange(sps_sciex[1:2, ], change = "jump", msLevel = 1L))
    expect_equal(names(attributes(tmp)), c("msSignal10xChange"))
    expect_equal(attr(tmp, "msSignal10xChange"), "MS:4000097")
    tmp <- suppressWarnings(
        msSignal10xChange(sps_sciex[1:2, ], change = "fall", msLevel = 1L))
    expect_equal(names(attributes(tmp)), c("msSignal10xChange"))
    expect_equal(attr(tmp, "msSignal10xChange"), "MS:4000098")
    tmp <- suppressWarnings(
        msSignal10xChange(sps_sciex[1:2, ], change = "jump", msLevel = 2L))
    expect_equal(names(attributes(tmp)), NULL)
    tmp <- suppressWarnings(
        msSignal10xChange(sps_sciex[1:2, ], change = "fall", msLevel = 2L))
    expect_equal(names(attributes(tmp)), NULL)
})
## END unit test msSignal10xChange ##

## START unit test numberEmptyScans ##
test_that("numberEmptyScans works properly.", {
    expect_error(numberEmptyScans(NULL), "unable to find an inherited method")
    expect_error(numberEmptyScans(1:10), "unable to find an inherited method")
    expect_equal(as.numeric(numberEmptyScans(sps_sciex, msLevel = 1L)), 0)
    expect_equal(as.numeric(numberEmptyScans(sps_sciex, msLevel = 2L)), 0)
    expect_equal(as.numeric(numberEmptyScans(sps_sciex, msLevel = 3L)), 0)
    
    ## create one Spectra object with one missing entry
    spd <- DataFrame(
        msLevel = c(2L, 2L), polarity = c(1L, 1L),
        id = c("unknown", "HMDB0000001"), 
        name = c("unknown", "1-Methylhistidine"))
    ## Assign m/z and intensity values
    spd$mz <- list(
        c(NaN), c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16))
    spd$intensity <- list(
        c(NaN), c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643))
    sps_tmp <- Spectra(spd)
    expect_equal(as.numeric(numberEmptyScans(sps_tmp, msLevel = 2L)), 1)
                 
    ## test attributes
    tmp <- numberEmptyScans(sps_sciex, msLevel = 1L)
    expect_equal(names(attributes(tmp)), "numberEmptyScans")
    expect_equal(attr(tmp, "numberEmptyScans"), "MS:4000099")
    tmp <- numberEmptyScans(sps_sciex, msLevel = 2L)
    expect_equal(names(attributes(tmp)), "numberEmptyScans")
    expect_equal(attr(tmp, "numberEmptyScans"), "MS:4000100")
    tmp <- numberEmptyScans(sps_sciex, msLevel = 3L)
    expect_equal(names(attributes(tmp)), "numberEmptyScans")
    expect_equal(attr(tmp, "numberEmptyScans"), "MS:4000101")
})
## END unit test numberEmptyScans ##

## START unit test ratioCharge1over2 ##
test_that("ratioCharge1over2 works properly.", {
    expect_error(ratioCharge1over2(NULL), "unable to find an inherited method")
    expect_error(ratioCharge1over2(1:10), "unable to find an inherited method")
    expect_equal(ratioCharge1over2(sps_sciex, 
        identificationLevel = "all"), 
        0.7 / 0.15, tolerance = 3e-02)
    expect_equal(ratioCharge1over2(sps_sciex, 
        identificationLevel = "identified"), 
        0.7 / 0.15, tolerance = 3e-02)
    expect_equal(ratioCharge1over2(sps_sciex, 
        identificationLevel = "unidentified"), 
        0.7 / 0.15, tolerance = 3e-02)
    
    ## do not include charges 1 or 2
    sps_sciex_foo <- sps_sciex
    sps_sciex_foo@backend$precursorCharge <- as.integer(
        sample(x = c(1, 2, 3, 4), size = 1862, 
               replace = TRUE, prob = c(0.0, 0.85, 0.1, 0.05)))
    expect_equal(ratioCharge1over2(sps_sciex_foo), NaN)
    sps_sciex_foo <- sps_sciex
    sps_sciex_foo@backend$precursorCharge <- as.integer(
        sample(x = c(1, 2, 3, 4), size = 1862, 
               replace = TRUE, prob = c(0.85, 0.0, 0.1, 0.05)))
    expect_equal(ratioCharge1over2(sps_sciex_foo), NaN)
    
    ## test attributes
    tmp <- ratioCharge1over2(sps_sciex, identificationLevel = "all")
    expect_equal(names(attributes(tmp)), "ratioCharge1over2")
    expect_equal(attr(tmp, "ratioCharge1over2"), "MS:4000168")
    tmp <- ratioCharge1over2(sps_sciex, identificationLevel = "identified")
    expect_equal(names(attributes(tmp)), "ratioCharge1over2")
    expect_equal(attr(tmp, "ratioCharge1over2"), "MS:4000167")
    tmp <- ratioCharge1over2(sps_sciex, identificationLevel = "unidentified")
    expect_equal(names(attributes(tmp)), NULL)
})
## END unit test ratioCharge1over2 ##

## START unit test ratioCharge3over2 ##
test_that("ratioCharge3over2 works properly.", {
    expect_error(ratioCharge3over2(NULL), "unable to find an inherited method")
    expect_error(ratioCharge3over2(1:10), "unable to find an inherited method")
    expect_equal(ratioCharge3over2(sps_sciex, 
            identificationLevel = "all"), 
        0.1 / 0.15, tolerance = 2e-02)
    expect_equal(ratioCharge3over2(sps_sciex, 
            identificationLevel = "identified"), 
        0.1 / 0.15, tolerance = 2e-02)
    expect_equal(ratioCharge3over2(sps_sciex, 
            identificationLevel = "unidentified"), 
        0.1 / 0.15, tolerance = 2e-02)
    
    ## do not include charges 2 or 3
    sps_sciex_foo <- sps_sciex
    sps_sciex_foo@backend$precursorCharge <- as.integer(
        sample(x = c(1, 2, 3, 4), size = 1862, 
               replace = TRUE, prob = c(0.7, 0.0, 0.25, 0.05)))
    expect_equal(ratioCharge3over2(sps_sciex_foo), NaN)
    sps_sciex_foo <- sps_sciex
    sps_sciex_foo@backend$precursorCharge <- as.integer(
        sample(x = c(1, 2, 3, 4), size = 1862, 
               replace = TRUE, prob = c(0.7, 0.25, 0.0, 0.05)))
    expect_equal(ratioCharge3over2(sps_sciex_foo), NaN)
    
    ## test attributes
    tmp <- ratioCharge3over2(sps_sciex, identificationLevel = "all")
    expect_equal(names(attributes(tmp)), "ratioCharge3over2")
    expect_equal(attr(tmp, "ratioCharge3over2"), "MS:4000170")
    tmp <- ratioCharge3over2(sps_sciex, identificationLevel = "identified")
    expect_equal(names(attributes(tmp)), "ratioCharge3over2")
    expect_equal(attr(tmp, "ratioCharge3over2"), "MS:4000169")
    tmp <- ratioCharge3over2(sps_sciex, identificationLevel = "unidentified")
    expect_equal(names(attributes(tmp)), NULL)
})
## END unit test ratioCharge3over2 ##

## START unit test ratioCharge4over2 ##
test_that("ratioCharge4over2 works properly.", {
    expect_error(ratioCharge4over2(NULL), "unable to find an inherited method")
    expect_error(ratioCharge4over2(1:10), "unable to find an inherited method")
    expect_equal(ratioCharge4over2(sps_sciex, 
            identificationLevel = "all"), 
        0.05 / 0.15, tolerance = 3e-02)
    expect_equal(ratioCharge4over2(sps_sciex, 
            identificationLevel = "identified"), 
        0.05 / 0.15, tolerance = 3e-02)
    expect_equal(ratioCharge4over2(sps_sciex, 
            identificationLevel = "unidentified"), 
        0.05 / 0.15, tolerance = 3e-02)
    
    ## do not include charges 2 or 4
    sps_sciex_foo <- sps_sciex
    sps_sciex_foo@backend$precursorCharge <- as.integer(
        sample(x = c(1, 2, 3, 4), size = 1862, 
               replace = TRUE, prob = c(0.7, 0.0, 0.1, 0.20)))
    expect_equal(ratioCharge4over2(sps_sciex_foo), NaN)
    sps_sciex_foo <- sps_sciex
    sps_sciex_foo@backend$precursorCharge <- as.integer(
        sample(x = c(1, 2, 3, 4), size = 1862, 
               replace = TRUE, prob = c(0.7, 0.25, 0.15, 0.0)))
    expect_equal(ratioCharge4over2(sps_sciex_foo), NaN)
    
    ## test attributes
    tmp <- ratioCharge4over2(sps_sciex, identificationLevel = "all")
    expect_equal(names(attributes(tmp)), "ratioCharge4over2")
    expect_equal(attr(tmp, "ratioCharge4over2"), "MS:4000172")
    tmp <- ratioCharge4over2(sps_sciex, identificationLevel = "identified")
    expect_equal(names(attributes(tmp)), "ratioCharge4over2")
    expect_equal(attr(tmp, "ratioCharge4over2"), "MS:4000171")
    tmp <- ratioCharge4over2(sps_sciex, identificationLevel = "unidentified")
    expect_equal(names(attributes(tmp)), NULL)
})
## END unit test ratioCharge4over2 ##


## START unit test meanCharge ##
test_that("meanCharge works properly.", {
    expect_error(meanCharge(NULL), "unable to find an inherited method")
    expect_error(meanCharge(1:10), "unable to find an inherited method")
    expect_equal(as.numeric(meanCharge(sps_sciex, msLevel = 1, 
            identificationLevel = "all")), 
        1.520408, tolerance = 2e-02)
    expect_equal(as.numeric(meanCharge(sps_sciex, msLevel = 1, 
            identificationLevel = "identified")), 
        1.520408, tolerance = 2e-02)
    expect_equal(as.numeric(meanCharge(sps_sciex, msLevel = 1, 
            identificationLevel = "unidentified")), 
        1.520408, tolerance = 2e-02)
    expect_equal(as.numeric(meanCharge(sps_sciex, msLevel = 2, 
            identificationLevel = "all")), 
        NaN)
    expect_equal(as.numeric(meanCharge(sps_sciex, msLevel = 2, 
            identificationLevel = "identified")), 
        NaN)
    expect_equal(as.numeric(meanCharge(sps_sciex, msLevel = 2, 
            identificationLevel = "unidentified")), 
        NaN)
    
    ## test attributes
    tmp <- meanCharge(sps_sciex, msLevel = 1, identificationLevel = "all")
    expect_equal(names(attributes(tmp)), "meanCharge")
    expect_equal(attr(tmp, "meanCharge"), "MS:4000174")
    tmp <- meanCharge(sps_sciex, msLevel = 1, identificationLevel = "identified")
    expect_equal(names(attributes(tmp)), "meanCharge")
    expect_equal(attr(tmp, "meanCharge"), "MS:4000173")
    tmp <- meanCharge(sps_sciex, msLevel = 1, identificationLevel = "unidentified")
    expect_equal(names(attributes(tmp)), NULL)
    tmp <- meanCharge(sps_sciex, msLevel = 2, identificationLevel = "all")
    expect_equal(names(attributes(tmp)), "meanCharge")
    expect_equal(attr(tmp, "meanCharge"), "MS:4000174")
    tmp <- meanCharge(sps_sciex, msLevel = 2, identificationLevel = "identified")
    expect_equal(names(attributes(tmp)), "meanCharge")
    expect_equal(attr(tmp, "meanCharge"), "MS:4000173")
    tmp <- meanCharge(sps_sciex, msLevel = 2, identificationLevel = "unidentified")
    expect_equal(names(attributes(tmp)), NULL)
})
## END unit test meanCharge ##


## START unit test medianCharge ##
test_that("medianCharge works properly.", {
    expect_error(medianCharge(NULL), "unable to find an inherited method")
    expect_error(medianCharge(1:10), "unable to find an inherited method")
    expect_equal(as.numeric(medianCharge(sps_sciex, msLevel = 1, 
            identificationLevel = "all")), 
        1.520408, tolerance = 2e-02)
    expect_equal(as.numeric(medianCharge(sps_sciex, msLevel = 1, 
            identificationLevel = "identified")), 
        1.520408, tolerance = 2e-02)
    expect_equal(as.numeric(medianCharge(sps_sciex, msLevel = 1, 
            identificationLevel = "unidentified")), 
        1.520408, tolerance = 2e-02)
    expect_equal(as.numeric(medianCharge(sps_sciex, msLevel = 2, 
            identificationLevel = "all")), 
        NaN)
    expect_equal(as.numeric(medianCharge(sps_sciex, msLevel = 2, 
            identificationLevel = "identified")), 
        NaN)
    expect_equal(as.numeric(medianCharge(sps_sciex, msLevel = 2, 
            identificationLevel = "unidentified")), 
        NaN)
    
    ## test attributes
    tmp <- medianCharge(sps_sciex, msLevel = 1, identificationLevel = "all")
    expect_equal(names(attributes(tmp)), "medianCharge")
    expect_equal(attr(tmp, "medianCharge"), "MS:4000174")
    tmp <- medianCharge(sps_sciex, msLevel = 1, identificationLevel = "identified")
    expect_equal(names(attributes(tmp)), "medianCharge")
    expect_equal(attr(tmp, "medianCharge"), "MS:4000173")
    tmp <- medianCharge(sps_sciex, msLevel = 1, identificationLevel = "unidentified")
    expect_equal(names(attributes(tmp)), NULL)
    tmp <- medianCharge(sps_sciex, msLevel = 2, identificationLevel = "all")
    expect_equal(names(attributes(tmp)), "medianCharge")
    expect_equal(attr(tmp, "medianCharge"), "MS:4000174")
    tmp <- medianCharge(sps_sciex, msLevel = 2, identificationLevel = "identified")
    expect_equal(names(attributes(tmp)), "medianCharge")
    expect_equal(attr(tmp, "medianCharge"), "MS:4000173")
    tmp <- medianCharge(sps_sciex, msLevel = 2, identificationLevel = "unidentified")
    expect_equal(names(attributes(tmp)), NULL)
})
## END unit test medianCharge ##