################################################################################
################## PART 1: FULL CHROMATOGRAMS METRICS ##########################
################################################################################

#' @name chromatographyDuration
#'
#' @rdname chromatographyDuration
#'
#' @aliases chromatographyDuration,Chromatograms-method chromatographyDuration,Spectra-method
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L, 1L),
#'     mz = c(112.2, 123.3, 134.4),
#'     dataOrigin = c("mem1", "mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = numeric(), intensity = numeric()),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#'
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' chromatographyDuration(object = chr)
NULL

#' @noRd
.chromatographyDuration_chromatograms <- function(chromatograms, ...) {
    rts <- unlist(rtime(chromatograms), use.names = FALSE)
    rts <- as.numeric(rts)

    if (length(rts) == 0 || all(is.na(rts))) {
        res <- NA_real_
    } else {
        res <- max(rts, na.rm = TRUE) - min(rts, na.rm = TRUE)
    }

    attr(res, "chromatographyDuration") <- "MS:4000053"
    res
}

#' @rdname chromatographyDuration
setMethod("chromatographyDuration", "Chromatograms", function(object, ...) {
    .chromatographyDuration_chromatograms(chromatograms = object, ...)
})

#' @title Number of Chromatograms
#'
#' @description
#' The function `chromatogramCount` computes the number of chromatograms in the
#' `Chromatograms` object.
#'
#' @details
#' This metric corresponds to the PSI:MS term:
#' MS:4000071 (number of chromatograms).
#'
#' @param chromatograms `Chromatograms` object
#' @param ... further arguments
#'
#' @return `integer(1)`
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L),
#'     mz = c(112.2, 123.3),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Returns number of chromatograms
#' chromatogramCount(chr)
chromatogramCount <- function(chromatograms, ...) {
    res <- length(chromatograms)
    names(res) <- "chromatogramCount"
    attr(res, "chromatogramCount") <- "MS:4000071"
    res
}

#' @title Number of empty chromatograms
#'
#' @description
#' The function `numberEmptyChrom` returns the number of chromatograms that are empty.
#'
#' @details
#' The function returns the number of empty chromatograms.
#'
#' No specific PSI:MS term exists. It could be created.
#'
#' @param chromatograms `Chromatograms` object
#' @param ... further arguments (currently ignored)
#'
#' @return `integer(1)`
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L, 1L),
#'     mz = c(112.2, 123.3, 134.4),
#'     dataOrigin = c("mem1", "mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = numeric(), intensity = numeric()),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Returns 1 (the second chromatogram is empty)
#' numberEmptyChrom(chr)
numberEmptyChrom <- function(chromatograms, ...) {
    res <- sum(lengths(chromatograms) == 0)
    res
}

#' @name rtAcquisitionRange
#'
#' @rdname rtAcquisitionRange
#'
#' @aliases rtAcquisitionRange,Chromatograms-method
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L),
#'     mz = c(112.2, 123.3),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' rtAcquisitionRange(chr)
NULL

#' @noRd
.rtAcquisitionRange_chromatograms <- function(chromatograms, na.rm = TRUE, ...) {
    rts <- unlist(rtime(chromatograms), use.names = FALSE)
    if (length(rts) == 0) {
        res <- c(min = NA_real_, max = NA_real_)
    } else {
        res <- range(rts, na.rm = na.rm)
        names(res) <- c("min", "max")
    }
    attr(res, "rtAcquisitionRange") <- "MS:4000070"
    res
}

#' @rdname rtAcquisitionRange
setMethod("rtAcquisitionRange", "Chromatograms", function(object, ...) {
    .rtAcquisitionRange_chromatograms(object, ...)
})

#' @title Maximum intensity across chromatograms
#'
#' @description
#' The function `maxIntensity` returns the maximum intensity value observed
#' across all chromatograms.
#'
#' @details
#' The function returns the maximum intensity across all chromatograms.
#'
#' No specific PSI:MS term exists for chromatogram maximum intensity.
#' Note: MS:4000063 is "MS2 known precursor charges fractions", not max intensity.
#'
#' @param chromatograms `Chromatograms` object
#' @param na.rm `logical(1)` whether to remove `NA` values (default `TRUE`)
#' @param ... further arguments passed to `max`
#'
#' @return `numeric(1)`
#'
#' @author Philippine Louail
#'
#' @importFrom Chromatograms intensity
#'
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L),
#'     mz = c(112.2, 123.3),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Returns max intensity across all chromatograms: 1200
#' maxIntensity(chr)
maxIntensity <- function(chromatograms, na.rm = TRUE, ...) {
    ints <- intensity(chromatograms)
    # Ensure we properly extract numeric values even in parallel context
    ints <- as.numeric(unlist(ints, use.names = FALSE))
    res <- max(ints, na.rm = na.rm)
    attr(res, "maxIntensity") <- "custom_metric:max_intensity"
    res
}

#' @title Intensity statistics (Mean) across chromatograms
#'
#' @description
#' The function `intensityMean` calculates the mean of intensity values
#' across all chromatograms.
#'
#' @details
#' The function returns the mean intensity across all chromatograms.
#'
#' No specific PSI:MS term exists. It should be created.
#'
#' @param chromatograms `Chromatograms` object
#' @param na.rm `logical(1)` whether to remove `NA` values (default `TRUE`)
#' @param ... further arguments passed to `mean`
#'
#' @return `numeric(1)`
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L),
#'     mz = c(112.2, 123.3),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Returns mean intensity across all chromatograms
#' intensityMean(chr)
intensityMean <- function(chromatograms, na.rm = TRUE, ...) {
    ints <- intensity(chromatograms)
    # Ensure we properly extract numeric values even in parallel context
    ints <- as.numeric(unlist(ints, use.names = FALSE))
    res <- mean(ints, na.rm = na.rm)
    attr(res, "intensityMean") <- "custom_metric:intensity_mean"
    res
}

#' @title Intensity statistics (Standard Deviation) across chromatograms
#'
#' @description
#' The function `intensitySd` calculates the standard deviation of intensity
#' values across all chromatograms.
#'
#' @details
#' The function returns the standard deviation of the intensity across all chromatograms.
#'
#' No specific PSI:MS term exists. It should be created.
#'
#' @param chromatograms `Chromatograms` object
#' @param na.rm `logical(1)` whether to remove `NA` values (default `TRUE`)
#' @param ... further arguments passed to `sd`
#'
#' @return `numeric(1)`
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L),
#'     mz = c(112.2, 123.3),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Returns standard deviation of intensities across all chromatograms
#' intensitySd(chr)
intensitySd <- function(chromatograms, na.rm = TRUE, ...) {
    ints <- intensity(chromatograms)
    # Ensure we properly extract numeric values even in parallel context
    ints <- as.numeric(unlist(ints, use.names = FALSE))
    res <- sd(ints, na.rm = na.rm)
    attr(res, "intensitySd") <- "custom_metric:intensity_sd"
    res
}

#' @title Intensity Quartiles per chromatogram
#'
#' @description
#' The function `intensityQuartiles` calculates the minimum, 1st quartile,
#' median, 3rd quartile, and maximum of intensity values for each chromatogram.
#'
#' @details
#' The function returns a list of vectors, where each vector contains the
#' summary statistics for a chromatogram.
#'
#' No specific PSI:MS term exists. It should be created.
#'
#' @param chromatograms `Chromatograms` object
#' @param ... further arguments passed to `summary`
#'
#' @return `numeric(6)` named vector with Min, 1st Qu., Median, Mean, 3rd Qu., Max
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L),
#'     mz = c(112.2, 123.3),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Returns Min, 1st Qu., Median, Mean, 3rd Qu., Max for the overall intensities
#' intensityQuartiles(chr)
intensityQuartiles <- function(chromatograms, ...) {
    res <- unlist(intensity(chromatograms), use.names = FALSE)
    ## order overall based on RT
    rts <- unlist(rtime(chromatograms), use.names = FALSE)
    ord <- order(rts)
    res <- res[ord]
    qt <- summary(res)
    ## name the output, add attributes
    names(qt) <- c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max")
    attr(qt, "intensityQuartiles") <- "custom_metric:intensity_quartiles"
    qt
}

#' @title Intensity Range per chromatogram
#'
#' @description
#' The function `intensityRange` calculates the range (min, max) of intensity
#' values for each chromatogram.
#'
#' @details
#' The function returns a list of numeric vectors of length 2.
#'
#' No specific PSI:MS term exists. It should be created.
#'
#' @param chromatograms `Chromatograms` object
#' @param na.rm `logical(1)` whether to remove `NA` values (default `TRUE`)
#' @param ... further arguments passed to `range`
#'
#' @return `numeric(2)` named vector with min and max
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L),
#'     mz = c(112.2, 123.3),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Returns (min, max) intensity across all chromatograms
#' intensityRange(chr)
intensityRange <- function(chromatograms, na.rm = TRUE, ...) {
    ints <- unlist(intensity(chromatograms), use.names = FALSE)
    min_val <- min(ints, na.rm = na.rm)
    max_val <- max(ints, na.rm = na.rm)
    res <- c(min = min_val, max = max_val)
    attr(res, "intensityRange") <- "custom_metric:intensity_range"
    res
}

#' @title Peak Count (Number of data points) per chromatogram
#'
#' @description
#' The function `peakCount` returns the number of data points (retention
#' time - intensity pairs) for each chromatogram.
#'
#' @details
#' The function returns the count of data points per chromatogram.
#'
#' No specific PSI:MS term exists for chromatogram peak count.
#'
#' @param chromatograms `Chromatograms` object
#' @param na.rm `logical(1)` indicating whether `NA` values should be
#' removed before counting (default `FALSE`)
#' @param ... further arguments (currently ignored)
#'
#' @return `integer(1)` total count of data points across all chromatograms
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L, 1L),
#'     mz = c(112.2, 123.3, 134.4),
#'     dataOrigin = c("mem1", "mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = numeric(), intensity = numeric()),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Returns total number of data points: 10 (5 + 0 + 5)
#' peakCount(chr)
peakCount <- function(chromatograms, na.rm = FALSE, ...) {
    int <- unlist(intensity(chromatograms), use.names = FALSE)
    if (na.rm) {
        res <- sum(!is.na(int))
    } else {
        res <- length(int)
    }
    attr(res, "peakCount") <- "custom_metric:peak_count"
    res
}

#' @name rtIqr
#'
#' @rdname rtIqr
#'
#' @aliases rtIqr,Chromatograms-method
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L),
#'     mz = c(112.2, 123.3),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Returns IQR of retention times across all chromatograms
#' rtIqr(chr)
NULL

#' @noRd
.rtIqr_chromatograms <- function(chromatograms, na.rm = TRUE, ...) {
    res <- IQR(unlist(rtime(chromatograms), use.names = FALSE), na.rm = na.rm)
    attr(res, "rtIqr") <- "custom_metric:rt_iqr"
    res
}

#' @rdname rtIqr
setMethod("rtIqr", "Chromatograms", function(object, ...) {
    .rtIqr_chromatograms(object, ...)
})

#' @title Baseline Intensity across chromatograms
#'
#' @description
#' The function `baselineIntensity` estimates the baseline intensity as the
#' 5th percentile of intensity values across all chromatograms.
#'
#' @details
#' The function returns the 5th percentile of intensities across all chromatograms.
#'
#' No specific PSI:MS term exists. It should be created.
#'
#' @param chromatograms `Chromatograms` object
#' @param probs `numeric(1)` quantile probability (default 0.05)
#' @param na.rm `logical(1)` whether to remove `NA` values (default `TRUE`)
#' @param ... further arguments passed to `quantile`
#'
#' @return `numeric(1)`
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L),
#'     mz = c(112.2, 123.3),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Returns 5th percentile of intensity (baseline estimate)
#' baselineIntensity(chr)
#' ## Use different percentile
#' baselineIntensity(chr, probs = 0.10)
baselineIntensity <- function(chromatograms, probs = 0.05, na.rm = TRUE, ...) {
    res <- quantile(
        unlist(intensity(chromatograms), use.names = FALSE),
        probs = probs,
        na.rm = na.rm
    )
    attr(res, "baselineIntensity") <- "custom_metric:baseline_intensity"
    res
}

#' @title Signal-to-Noise Ratio across chromatograms
#'
#' @description
#' The function `signalToNoiseRatio` calculates a robust signal-to-noise ratio
#' approximated as the Maximum Intensity divided by the Median of Non-Zero
#' Intensities.
#'
#' @details
#' The noise is estimated as the median of all intensity values greater than 0.
#' This avoids underestimating noise in sparse chromatograms where the
#' 5th percentile might be 0.
#'
#' @param chromatograms `Chromatograms` object
#' @param ... further arguments (currently ignored)
#'
#' @return `numeric(1)`
#'
#' @author Philippine Louail
#'
#' @export
#'
signalToNoiseRatio <- function(chromatograms, ...) {
    all_ints <- unlist(intensity(chromatograms), use.names = FALSE)

    ## Filter for valid signals
    valid_ints <- all_ints[!is.na(all_ints) & all_ints > 0]

    if (length(valid_ints) == 0) {
        res <- NA_real_
    } else {
        ## Signal: Maximum observed intensity
        signal <- max(valid_ints)

        ## Noise: Median of non-zero intensities (Robust estimator)
        noise <- median(valid_ints)

        if (noise == 0) {
            res <- NA_real_
        } else {
            res <- signal / noise
        }
    }

    attr(res, "signalToNoiseRatio") <- "custom_metric:signal_to_noise_ratio"
    res
}

#' @title Number of 10x Intensity Changes across chromatograms
#'
#' @description
#' The function `intensity10xChange` calculates the total number of times the
#' intensity changes by a factor of 10 or more between consecutive data points
#' across all chromatograms.
#'
#' @details
#' The function returns the total count of such changes across all chromatograms.
#'
#' No specific PSI:MS term exists. It should be created.
#'
#' @param chromatograms `Chromatograms` object
#' @param change `character(1)` "jump" or "fall"
#' @param ... further arguments
#'
#' @return `integer(1)`
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L),
#'     mz = c(112.2, 123.3),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(1.0, 2.0, 3.0, 4.0, 5.0),
#'                intensity = c(50, 100, 2000, 100, 50)),
#'     data.frame(rtime = c(1.0, 2.0, 3.0, 4.0, 5.0),
#'                intensity = c(200, 20, 30, 20, 200))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Count 10x intensity jumps (increases) across all chromatograms
#' intensity10xChange(chr, change = "jump")
#' ## Count 10x intensity falls (decreases) across all chromatograms
#' intensity10xChange(chr, change = "fall")
intensity10xChange <- function(chromatograms, change = c("jump", "fall"), ...) {
change <- match.arg(change)
    ## Extract list once - much faster
    all_ints <- intensity(chromatograms)

    count_changes <- function(ints) {
        if (length(ints) < 2) return(0L)
        ## Assuming RTs are already sorted in Chromatograms objects (usually true)
        ## If not, you must sort ints here.
        valid <- ints > 0 & !is.na(ints)
        ints <- ints[valid]
        if (length(ints) < 2) return(0L)

        ratios <- ints[-1] / ints[-length(ints)]
        if (change == "jump") sum(ratios >= 10) else sum(ratios <= 0.1)
    }

    res <- sum(vapply(all_ints, count_changes, integer(1)))
    if (change == "jump") {
        attr(res, "intensity10xJump") <- "custom_metric:intensity_10x_jump"
    } else {
        attr(res, "intensity10xFall") <- "custom_metric:intensity_10x_fall"
    }
    res

}#' @title Median Intensity in RT IQR across chromatograms
#'
#' @description
#' The function `medianIntensityRtIqr` calculates the median intensity of data
#' points that fall within the interquartile range (IQR) of retention times
#' across all chromatograms.
#'
#' @details
#' The function returns the median intensity within the RT IQR across all
#' chromatograms.
#'
#' No specific PSI:MS term exists. It should be created.
#'
#' @param chromatograms `Chromatograms` object
#' @param na.rm `logical(1)` whether to remove `NA` values (default `TRUE`)
#' @param ... further arguments passed to `median` and `quantile`
#'
#' @return `numeric(1)`
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L),
#'     mz = c(112.2, 123.3),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Returns median intensity of points within RT IQR
#' medianIntensityRtIqr(chr)
medianIntensityRtIqr <- function(chromatograms, na.rm = TRUE, ...) {
    all_rts <- unlist(rtime(chromatograms), use.names = FALSE)
    all_ints <- unlist(intensity(chromatograms), use.names = FALSE)

    if (length(all_rts) == 0 || all(is.na(all_rts))) {
        res <- NA_real_
        attr(res, "medianIntensityRtIqr") <- "custom_metric:median_intensity_rt_iqr"
        return(res)
    }

    q_rt <- quantile(all_rts, probs = c(0.25, 0.75), na.rm = na.rm)
    mask <- all_rts >= q_rt[1] & all_rts <= q_rt[2] & !is.na(all_rts)

    if (!any(mask)) {
        res <- NA_real_
        attr(res, "medianIntensityRtIqr") <- "custom_metric:median_intensity_rt_iqr"
        return(res)
    }
    res <- median(all_ints[mask], na.rm = na.rm)
    attr(res, "medianIntensityRtIqr") <- "custom_metric:median_intensity_rt_iqr"
    res
}

################################################################################
#################### PART 2: XIC/EIC METRICS FROM PSI-MS #######################
################################################################################

#' @title Full Width at Half Maximum (FWHM) for a Chromatogram
#'
#' @description
#' The function `xicFwhm` calculates the Full Width at Half Maximum (FWHM)
#' for a single chromatogram.
#'
#' @details
#' The FWHM is calculated by finding the maximum intensity peak, determining
#' 50\% of that intensity, and linearly interpolating the time difference between
#' the left and right crossing points.
#'
#' If `peakBoundary` is provided, the calculation is restricted to the peak
#' region defined by those boundaries, which can improve accuracy for noisy
#' chromatograms and efficiency when calculating multiple metrics.
#'
#' This metric is analogous to MS:4000051 (XIC-FWHM quantiles).
#'
#' @param chromatograms `Chromatograms` object containing a single chromatogram
#' @param peakBoundary optional `numeric(2)` named vector with `left_boundary`
#'   and `right_boundary` values from a previous call to `peakBoundary()`.
#'   If not provided, the full chromatogram is used.
#' @param ... further arguments
#'
#' @return `numeric(1)` FWHM value. Returns `NA` if FWHM cannot be calculated.
#'
#' @author Philippine Louail
#'
#' @importFrom stats approx
#' @importFrom utils tail
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(msLevel = 1L, mz = 112.2, dataOrigin = "mem1")
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' xicFwhm(chr)
#'
#' ## Use pre-computed peak boundaries for efficiency
#' pb <- peakBoundary(chr)
#' xicFwhm(chr, peakBoundary = pb)
xicFwhm <- function(chromatograms, peakBoundary = NULL, ...) {
    rts <- rtime(chromatograms)[[1L]]
    ints <- intensity(chromatograms)[[1L]]

    ## If peakBoundary provided, subset to peak region
    if (!is.null(peakBoundary) && !any(is.na(peakBoundary))) {
        left_rt <- peakBoundary["left_boundary"]
        right_rt <- peakBoundary["right_boundary"]
        mask <- rts >= left_rt & rts <= right_rt
        rts <- rts[mask]
        ints <- ints[mask]
    }

    if (length(ints) < 3 || all(is.na(ints))) {
        res <- NA_real_
        attr(res, "xicFwhm") <- "custom_metric:xic_fwhm"
        return(res)
    }

    max_int <- max(ints, na.rm = TRUE)
    if (max_int == 0) {
        res <- NA_real_
        attr(res, "xicFwhm") <- "custom_metric:xic_fwhm"
        return(res)
    }

    max_idx <- which.max(ints)
    half_max <- max_int / 2

    left_candidates <- which(ints[1:max_idx] < half_max)
    if (length(left_candidates) == 0) {
        res <- NA_real_
        attr(res, "xicFwhm") <- "custom_metric:xic_fwhm"
        return(res)
    }
    left_idx <- tail(left_candidates, 1)

    right_candidates <- which(ints[max_idx:length(ints)] < half_max)
    if (length(right_candidates) == 0) {
        res <- NA_real_
        attr(res, "xicFwhm") <- "custom_metric:xic_fwhm"
        return(res)
    }
    right_idx <- max_idx + right_candidates[1] - 1
    left_ints <- ints[c(left_idx, left_idx + 1)]
    left_rts <- rts[c(left_idx, left_idx + 1)]
    right_ints <- ints[c(right_idx - 1, right_idx)]
    right_rts <- rts[c(right_idx - 1, right_idx)]
    
    if (any(is.na(left_ints)) || any(is.na(left_rts)) ||
        any(is.na(right_ints)) || any(is.na(right_rts))) {
        res <- NA_real_
        attr(res, "xicFwhm") <- "custom_metric:xic_fwhm"
        return(res)
    }
    rt_left <- approx(x = left_ints, y = left_rts, xout = half_max)$y
    rt_right <- approx(x = right_ints, y = right_rts, xout = half_max)$y
    res <- rt_right - rt_left
    attr(res, "xicFwhm") <- "custom_metric:xic_fwhm"
    res
}

#' @title Peak Boundary for a Chromatogram
#'
#' @description
#' The function `peakBoundary` finds the left and right retention time boundaries
#' of the main peak in a single chromatogram (EIC).
#'
#' @details
#' The function uses an adaptive approach to find peak boundaries:
#'
#' 1. First attempts to find boundaries using local minima (valleys) on each
#'    side of the peak apex via `MsCoreUtils::valleys()`.
#'
#' 2. Validates the boundaries by checking if intensities at boundaries are
#'    near baseline level.
#'
#' 3. If boundaries are not at baseline (e.g., for peaks with elevated baseline
#'    or monotonic slopes), falls back to a relative threshold method that
#'    finds where intensity drops below a fraction of peak height above baseline.
#'
#' The baseline is estimated from a lower quantile of intensities (default 10th
#' percentile), and the threshold is calculated relative to peak height above
#' this baseline.
#'
#' @param chromatograms `Chromatograms` object containing a single chromatogram
#' @param threshold `numeric(1)` fraction of peak height above baseline used
#'   as fallback threshold (default 0.1 = 10\%). The boundary is placed where
#'   intensity drops below: `baseline + (max - baseline) * threshold`.
#' @param baselineThreshold `numeric(1)` maximum acceptable intensity at 
#'   boundaries as a fraction of peak height above baseline. If boundary 
#'   intensity exceeds this, falls back to threshold method. Default is 0.1.
#' @param baselineQuantile `numeric(1)` quantile used to estimate the baseline
#'   intensity (default 0.1 = 10th percentile of intensities).
#' @param ... further arguments (currently unused)
#'
#' @return `numeric(2)` named vector with `left_boundary` and `right_boundary` 
#'   values. Returns `NA` values if boundaries cannot be determined.
#'
#' @author Philippine Louail
#'
#' @importFrom MsCoreUtils valleys
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#'
#' ## Find peak boundaries
#' peakBoundary(chr)
#'
#' ## Adjust threshold for different sensitivity
#' peakBoundary(chr, threshold = 0.05)
peakBoundary <- function(chromatograms, 
                         threshold = 0.1,
                         baselineThreshold = 0.1,
                         baselineQuantile = 0.1,
                         ...) {
    rts <- rtime(chromatograms)[[1L]]
    ints <- intensity(chromatograms)[[1L]]
    n <- length(ints)
    na_result <- c(left_boundary = NA_real_, right_boundary = NA_real_)
    attr(na_result, "peakBoundary") <- "custom_metric:peak_boundary"

    if (n < 3 || all(is.na(ints))) return(na_result)
    
    max_int <- max(ints, na.rm = TRUE)
    if (max_int == 0) return(na_result)
    
    max_idx <- which.max(ints)
    baseline_int <- quantile(ints, probs = baselineQuantile, na.rm = TRUE)
    peak_height <- max_int - baseline_int
    baseline_thresh <- baseline_int + peak_height * baselineThreshold

    ## Try valley-based boundaries first
    v <- MsCoreUtils::valleys(ints, max_idx)
    left_idx <- if ("left" %in% colnames(v)) v[1L, "left"] else 1L
    right_idx <- if ("right" %in% colnames(v)) v[1L, "right"] else n

    ## Check if valleys are valid (at baseline level, not NA, not adjacent to NA)
    left_ok <- !is.na(ints[left_idx]) && ints[left_idx] <= baseline_thresh &&
               !(left_idx > 1 && is.na(ints[left_idx - 1]))
    right_ok <- !is.na(ints[right_idx]) && ints[right_idx] <= baseline_thresh &&
                !(right_idx < n && is.na(ints[right_idx + 1]))

    ## Fallback to threshold method if needed
    if (!left_ok || !right_ok) {
        thresh_val <- baseline_int + peak_height * threshold
        left_cand <- which(ints[seq_len(max_idx)] <= thresh_val)
        right_cand <- which(ints[max_idx:n] <= thresh_val)
        left_idx <- if (length(left_cand)) max(left_cand) else 1L
        right_idx <- if (length(right_cand)) max_idx + min(right_cand) - 1L else n
    }

    res <- c(left_boundary = rts[left_idx], right_boundary = rts[right_idx])
    attr(res, "peakBoundary") <- "custom_metric:peak_boundary"
    res
}

#' @title Peak Width for a Chromatogram
#'
#' @description
#' The function `peakWidth` calculates the width of the main peak in a single
#' chromatogram (EIC).
#'
#' @details
#' The peak width is calculated as the difference between the right and left
#' peak boundaries determined by `peakBoundary()`, which uses
#' `MsCoreUtils::valleys()` to find the local minima (valleys) on each side
#' of the peak apex.
#'
#' This is different from FWHM which measures width at 50\% of max intensity.
#'
#' The user should provide a single chromatogram representing an extracted ion
#' chromatogram (EIC) for a specific compound/feature.
#'
#' @param chromatograms `Chromatograms` object containing a single chromatogram
#' @param peakBoundary optional `numeric(2)` named vector with `left_boundary`
#'   and `right_boundary` values from a previous call to `peakBoundary()`.
#'   If not provided, boundaries are calculated automatically.
#' @param ... further arguments passed to `peakBoundary()`
#'
#' @return `numeric(1)` peak width value. Returns `NA` if width cannot be
#'   calculated.
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Returns peak width for the chromatogram
#' peakWidth(chr)
#'
#' ## Use pre-computed peak boundaries for efficiency
#' pb <- peakBoundary(chr)
#' peakWidth(chr, peakBoundary = pb)
peakWidth <- function(chromatograms, peakBoundary = NULL, ...) {
    if (is.null(peakBoundary)) {
        peakBoundary <- peakBoundary(chromatograms, ...)
    }
    res <- unname(peakBoundary["right_boundary"] - peakBoundary["left_boundary"])
    attr(res, "peakWidth") <- "custom_metric:peak_width"
    res
}

#' @title Peak Beta Values (Peak Shape Quality)
#'
#' @description
#' The function `peakBeta` calculates beta parameters for a chromatographic peak,
#' assessing its similarity to a bell curve (beta distribution) and the quality
#' of the peak shape.
#'
#' @details
#' This function wraps `MetaboCoreUtils::betaValues()` which compares the
#' chromatographic peak to Beta distribution curves. It returns two values:
#' - `beta_cor`: correlation/similarity to the best-fit beta distribution curve
#' - `beta_snr`: signal-to-noise ratio of the peak relative to residuals
#'
#' Higher `beta_cor` values (close to 1.0) indicate more symmetric, bell-shaped peaks.
#' Higher `beta_snr` values indicate stronger peak signal relative to noise/residuals,
#' meaning cleaner, better-defined peaks.
#'
#' If `peakBoundary` is provided, only the peak region is analyzed. Otherwise,
#' boundaries are calculated automatically using `peakBoundary()`.
#'
#' Requires at least 5 data points within the peak region.
#'
#' @references
#' Kumler W, Hazelton B J and Ingalls A E (2023) "Picky with peakpicking:
#' assessing chromatographic peak quality with simple metrics in metabolomics"
#' BMC Bioinformatics 24(1):404. doi: 10.1186/s12859-023-05533-4
#'
#' @param chromatograms `Chromatograms` object containing a single chromatogram
#' @param peakBoundary optional `numeric(2)` named vector with `left_boundary`
#'   and `right_boundary` values from a previous call to `peakBoundary()`.
#'   If not provided, boundaries are calculated automatically.
#' @param ... further arguments (currently unused)
#'
#' @return `numeric(2)` named vector with `beta_cor` and `beta_snr` values.
#'   Returns `NA` for both if values cannot be calculated.
#'
#' @author  William Kumler, Philippine Louail
#'
#' @importFrom MetaboCoreUtils betaValues
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
#' pdata <- list(
#'     data.frame(rtime = seq(1, 20, by = 0.5),
#'                intensity = c(10, 15, 25, 50, 100, 200, 350, 500, 650,
#'                              700, 650, 500, 350, 200, 100, 50, 25, 15,
#'                              10, 8, 6, 5, 4, 3, 2, 2, 1, 1, 1, 1, 1,
#'                              1, 1, 1, 1, 1, 1, 1, 1))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' peakBeta(chr)
#'
#' ## Use pre-computed peak boundaries for efficiency
#' pb <- peakBoundary(chr)
#' peakBeta(chr, peakBoundary = pb)
peakBeta <- function(chromatograms, peakBoundary = NULL, ...) {
    rts <- rtime(chromatograms)[[1L]]
    ints <- intensity(chromatograms)[[1L]]

    ## Calculate peak boundaries if not provided
    if (is.null(peakBoundary)) {
        peakBoundary <- peakBoundary(chromatograms)
    }

    ## Check for valid boundaries
    if (any(is.na(peakBoundary))) {
        res <- c(beta_cor = NA_real_, beta_snr = NA_real_)
        attr(res, "peakBeta") <- "custom_metric:peak_beta"
        return(res)
    }

    ## Subset to peak region
    left_rt <- peakBoundary["left_boundary"]
    right_rt <- peakBoundary["right_boundary"]
    mask <- rts >= left_rt & rts <= right_rt
    peak_rts <- rts[mask]
    peak_ints <- ints[mask]

    ## Need at least 5 points for betaValues
    if (length(peak_ints) < 5) {
        res <- c(beta_cor = NA_real_, beta_snr = NA_real_)
        attr(res, "peakBeta") <- "custom_metric:peak_beta"
        return(res)
    }

    beta_vals <- MetaboCoreUtils::betaValues(
        intensity = peak_ints,
        rtime = peak_rts
    )

    res <- c(beta_cor = unname(beta_vals[1]), beta_snr = unname(beta_vals[2]))
    attr(res, "peakBeta") <- "custom_metric:peak_beta"
    res
}

#' @title Peak Prominence (Peak-to-Baseline Ratio)
#'
#' @description
#' The function `peakProminence` calculates the prominence of a chromatographic
#' peak relative to its baseline, useful for filtering out flat/noisy signals.
#'
#' @details
#' Peak prominence is calculated as the ratio of peak height above baseline
#' to the baseline level:
#' 
#' \deqn{prominence = \frac{max - baseline}{baseline}}{prominence = (max - baseline) / baseline}
#' 
#' Where baseline is estimated from the lower quantile of intensities 
#' (controlled by `baselineQuantile`).
#'
#' Higher values indicate more prominent peaks that stand out clearly from
#' the baseline. Typical good peaks have prominence > 5-10, while noisy
#' plateaus or flat signals have prominence < 3-5.
#'
#' This metric is particularly useful for filtering out:
#' - Noisy chromatograms with no clear peak
#' - Flat plateau signals
#' - Weak peaks barely above baseline
#'
#' @param chromatograms `Chromatograms` object containing a single chromatogram
#' @param peakBoundary optional `numeric(2)` named vector with `left_boundary`
#'   and `right_boundary` values from a previous call to `peakBoundary()`.
#'   If provided, only the peak region is analyzed.
#' @param baselineQuantile `numeric(1)` quantile used to estimate the baseline
#'   intensity (default 0.1 = 10th percentile of intensities).
#' @param ... further arguments (currently unused)
#'
#' @return `numeric(1)` peak prominence value. Returns `NA` if prominence 
#'   cannot be calculated (e.g., baseline is zero or NA).
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @examples
#' library(Chromatograms)
#' 
#' ## Good peak with high prominence
#' cdata <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
#' pdata_good <- list(data.frame(
#'     rtime = 1:20,
#'     intensity = c(100, 100, 100, 200, 500, 1000, 2000, 5000, 10000, 15000,
#'                   10000, 5000, 2000, 1000, 500, 200, 100, 100, 100, 100)
#' ))
#' chr_good <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata_good)
#' peakProminence(chr_good)  # High value (~150)
#' 
#' ## Bad peak (noisy plateau) with low prominence
#' pdata_bad <- list(data.frame(
#'     rtime = 1:20,
#'     intensity = c(3000, 3500, 4000, 5000, 8000, 10000, 12000, 11000, 10000,
#'                   11000, 12000, 10000, 9000, 8000, 7000, 6000, 5000, 4000, 3500, 3000)
#' ))
#' chr_bad <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata_bad)
#' peakProminence(chr_bad)  # Low value (~3)
peakProminence <- function(chromatograms, peakBoundary = NULL, 
                           baselineQuantile = 0.1, ...) {
    rts <- rtime(chromatograms)[[1L]]
    ints <- intensity(chromatograms)[[1L]]
    
    ## If peakBoundary provided, subset to peak region
    if (!is.null(peakBoundary) && !any(is.na(peakBoundary))) {
        left_rt <- peakBoundary["left_boundary"]
        right_rt <- peakBoundary["right_boundary"]
        mask <- rts >= left_rt & rts <= right_rt
        ints <- ints[mask]
    }
    
    if (length(ints) < 3 || all(is.na(ints))) {
        res <- NA_real_
        attr(res, "peakProminence") <- "custom_metric:peak_prominence"
        return(res)
    }
    
    max_int <- max(ints, na.rm = TRUE)
    baseline_int <- quantile(ints, probs = baselineQuantile, na.rm = TRUE)
    
    ## Avoid division by zero or negative baseline
    if (is.na(baseline_int) || baseline_int <= 0) {
        res <- NA_real_
        attr(res, "peakProminence") <- "custom_metric:peak_prominence"
        return(res)
    }
    
    res <- (max_int - baseline_int) / baseline_int
    attr(res, "peakProminence") <- "custom_metric:peak_prominence"
    res
}

#' @title TIC quantile RT fraction
#'
#' @rdname ticQuantileRtFraction
#'
#' @description
#' The function `ticQuantileRtFraction` calculates the intervals when
#' respective quantiles of the TIC accumulate, divided by retention time
#' duration.
#'
#' @details
#' Returns the RT fraction when Q1, Q2 (median), Q3, and Q4 (max) of the
#' cumulative TIC are reached. Provides information on sample flow along the
#' chromatographic run.
#'
#' id: MS:4000183
#' name: TIC quantile RT fraction
#' def: "The interval when the respective quantile of the TIC accumulates
#' divided by retention time duration. The number of values in the tuple
#' implies the quantile mode." [PSI:MS]
#' comment: The metric informs about the dynamic range of the acquisition along
#' the chromatographic separation. The metric provides information on the sample
#' (compound) flow along the chromatographic run, potentially revealing poor
#' chromatographic performance, such as the absence of a signal for a
#' significant portion of the run.
#' is_a: MS:4000004 ! n-tuple
#' relationship: has_metric_category MS:4000009 ! ID free metric
#' relationship: has_metric_category MS:4000012 ! single run based metric
#' relationship: has_metric_category MS:4000016 ! retention time metric
#' relationship: has_metric_category MS:4000017 ! chromatogram metric
#' relationship: has_units UO:0000191 ! fraction
#'
#'
#' @return `numeric(4)` named vector with Q1, Q2, Q3, Q4 RT fractions
#'
#' @author Philippine Louail
#'
#' @importFrom stats setNames
#'
#' @aliases ticQuantileRtFraction,Chromatograms-method ticQuantileRtFraction,Spectra-method
#'
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L),
#'     mz = c(112.2, 123.3),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Returns RT fractions when Q1, Q2, Q3, Q4 of TIC accumulate
#' ticQuantileRtFraction(chr)
setMethod("ticQuantileRtFraction", "Chromatograms", function(
    object,
    probs = seq(0, 1, 0.25),
    ...
) {
    all_rts <- unlist(rtime(object), use.names = FALSE)
    all_ints <- unlist(intensity(object), use.names = FALSE)

    n_probs <- length(probs)
    if (length(all_rts) < 2 || all(is.na(all_ints))) {
        res <- setNames(rep(NA_real_, n_probs), paste0(probs * 100, "%"))
    } else {
        ## Sort by retention time
        ord <- order(all_rts)
        sorted_rts <- all_rts[ord]
        sorted_ints <- all_ints[ord]

        ## Calculate cumulative TIC
        cumsum_tic <- cumsum(sorted_ints)
        total_tic <- sum(sorted_ints)

        ## Find RT when each quantile is reached
        rt_duration <- max(sorted_rts) - min(sorted_rts)

        if (rt_duration == 0 || total_tic == 0) {
            res <- setNames(rep(NA_real_, n_probs), paste0(probs * 100, "%"))
        } else {
            quantile_rts <- vapply(
                probs,
                function(p) {
                    target <- p * total_tic
                    idx <- which(cumsum_tic >= target)[1]
                    if (is.na(idx)) {
                        return(max(sorted_rts))
                    }
                    sorted_rts[idx]
                },
                numeric(1)
            )

            ## Calculate fraction of RT duration
            rt_start <- min(sorted_rts)
            res <- (quantile_rts - rt_start) / rt_duration
            names(res) <- paste0(probs * 100, "%")
        }
    }
    attr(res, "ticQuantileRtFraction") <- "MS:4000183"
    res
})

#' @title Area under TIC in MS1
#'
#' @description
#' The function `areaUnderTicMs1` calculates the area under the total ion
#' current chromatogram for all MS1 spectra.
#'
#' @details
#' The function filters the chromatograms for `msLevel == 1` and calculates
#' the sum of intensities.
#'
#' id: MS:4000029
#' name: area under TIC in MS1
#' def: "The area under the total ion current chromatogram (MS:1000235) of all MS1 spectra." [PSI:MS]
#'
#' @param chromatograms `Chromatograms` object
#' @param na.rm `logical(1)` whether to remove `NA` values (default `TRUE`)
#' @param ... further arguments passed to `sum`
#'
#' @return `numeric(1)`
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @importFrom Chromatograms filterChromData
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 2L),
#'     mz = c(112.2, 123.3),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0), intensity = c(100, 250, 400)),
#'     data.frame(rtime = c(5.1, 5.8), intensity = c(80, 500))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Returns sum of intensities for MS1 only
#' areaUnderTicMs1(chr)
areaUnderTicMs1 <- function(chromatograms, na.rm = TRUE, ...) {
    # Check for empty input
    if (length(chromatograms) == 0) {
        res <- NA_real_
    } else {
        # Filter for MS1 - use tryCatch in case no chromatograms match
        chr_ms1 <- tryCatch(
            filterChromData(chromatograms, variables = c("msLevel"), ranges = c(1L, 1L)),
            error = function(e) NULL
        )
        if (is.null(chr_ms1) || length(chr_ms1) == 0) {
            res <- NA_real_
        } else {
            # Subset and calculate sum
            ints <- unlist(intensity(chr_ms1), use.names = FALSE)
            res <- sum(ints, na.rm = na.rm)
        }
    }

    attr(res, "areaUnderTicMs1") <- "MS:4000029"
    res
}

#' @title Area under TIC in MS2
#'
#' @description
#' The function `areaUnderTicMs2` calculates the area under the total ion
#' current chromatogram for all MS2 spectra.
#'
#' @details
#' The function filters the chromatograms for `msLevel == 2` and calculates
#' the sum of intensities.
#'
#' id: MS:4000030
#' name: area under TIC in MS2
#' def: "The area under the total ion current chromatogram (MS:1000235) of all MS2 spectra." [PSI:MS]
#'
#' @param chromatograms `Chromatograms` object
#' @param na.rm `logical(1)` whether to remove `NA` values (default `TRUE`)
#' @param ... further arguments passed to `sum`
#'
#' @return `numeric(1)`
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 2L),
#'     mz = c(112.2, 123.3),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0), intensity = c(100, 250, 400)),
#'     data.frame(rtime = c(5.1, 5.8), intensity = c(80, 500))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Returns sum of intensities for MS2 only
#' areaUnderTicMs2(chr)
areaUnderTicMs2 <- function(chromatograms, na.rm = TRUE, ...) {
    # Check for empty input
    if (length(chromatograms) == 0) {
        res <- NA_real_
    } else {
        # Filter for MS2 - use tryCatch in case no chromatograms match
        chr_ms2 <- tryCatch(
            filterChromData(chromatograms, variables = c("msLevel"), ranges = c(2L, 2L)),
            error = function(e) NULL
        )
        if (is.null(chr_ms2) || length(chr_ms2) == 0) {
            res <- NA_real_
        } else {
            # Subset and calculate sum
            ints <- unlist(intensity(chr_ms2), use.names = FALSE)
            res <- sum(ints, na.rm = na.rm)
        }
    }

    attr(res, "areaUnderTicMs2") <- "MS:4000030"
    res
}

#' @name areaUnderTicRtQuantiles
#'
#' @rdname areaUnderTicRtQuantiles
#'
#' @aliases areaUnderTicRtQuantiles,Chromatograms-method
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L),
#'     mz = c(112.2, 123.3),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' areaUnderTicRtQuantiles(chr)
NULL

#' @noRd
.areaUnderTicRtQuantiles_chromatograms <- function(chromatograms, msLevel = NULL, ...) {
    if (length(chromatograms) == 0) {
        res <- setNames(rep(NA_real_, 4), c("25%", "50%", "75%", "100%"))
        attr(res, "areaUnderTicRtQuantiles") <- "MS:4000156"
        return(res)
    }

    if (!is.null(msLevel)) {
        chromatograms <- tryCatch(
            filterChromData(chromatograms, variables = c("msLevel"), ranges = c(msLevel, msLevel)),
            error = function(e) NULL
        )

        if (is.null(chromatograms) || length(chromatograms) == 0) {
            res <- setNames(rep(NA_real_, 4), c("25%", "50%", "75%", "100%"))
            attr(res, "areaUnderTicRtQuantiles") <- "MS:4000156"
            return(res)
        }
    }

    ints <- unlist(intensity(chromatograms), use.names = FALSE)
    rts <- unlist(rtime(chromatograms), use.names = FALSE)

    if (length(rts) < 2 || all(is.na(rts))) {
        res <- setNames(rep(NA_real_, 4), c("25%", "50%", "75%", "100%"))
        attr(res, "areaUnderTicRtQuantiles") <- "MS:4000156"
        return(res)
    }

    ord <- order(rts)
    rts <- rts[ord]
    ints <- ints[ord]

    rt_range <- range(rts, na.rm = TRUE)

    # Define cut points for 4 quartiles
    cuts <- seq(rt_range[1], rt_range[2], length.out = 5)[2:4]

    new_rts <- rts
    new_ints <- ints

    # Interpolate intensities at cut points
    for (ct in cuts) {
        interp_val <- approx(rts, ints, xout = ct)$y
        if (!is.na(interp_val)) {
            new_rts <- c(new_rts, ct)
            new_ints <- c(new_ints, interp_val)
        }
    }

    # Re-order after adding interpolated points
    ord_new <- order(new_rts)
    final_rts <- new_rts[ord_new]
    final_ints <- new_ints[ord_new]

    # Calculate trapezoidal area segments
    areas <- (final_ints[-1] + final_ints[-length(final_ints)]) / 2 * diff(final_rts)

    # Midpoints for binning
    midpoints <- (final_rts[-1] + final_rts[-length(final_rts)]) / 2

    # Bin areas into quartiles
    breaks <- seq(rt_range[1], rt_range[2], length.out = 5)
    bins <- cut(midpoints, breaks = breaks, include.lowest = TRUE, labels = FALSE)

    res <- numeric(4)
    for (i in 1:4) {
        res[i] <- sum(areas[which(bins == i)], na.rm = TRUE)
    }

    names(res) <- c("25%", "50%", "75%", "100%")
    attr(res, "areaUnderTicRtQuantiles") <- "MS:4000156"
    res
}

#' @rdname areaUnderTicRtQuantiles
setMethod("areaUnderTicRtQuantiles", "Chromatograms", function(object, msLevel = NULL, ...) {
    .areaUnderTicRtQuantiles_chromatograms(object, msLevel = msLevel, ...)
})

#' @title Area under TIC
#'
#' @rdname areaUnderTic
#'
#' @aliases areaUnderTIC,Spectra-method
#'
#' @description
#' The function `areaUnderTic` calculates the total area under the total ion
#' current chromatogram across all chromatograms.
#'
#' @details
#' Returns the sum of all intensity values across all chromatograms. Differences
#' between samples may indicate differences in dynamic range or sample content.
#'
#' id: MS:4000155
#' name: area under TIC in MS1
#' def: "The area under the total ion current chromatogram (MS:1000235) of all
#' MS1 spectra." [PSI:MS]
#' comment: The metric informs about the dynamic range of the acquisition.
#' Differences between samples of an experiment may indicate differences in the
#' dynamic range and/or in the sample content. The signal can be affected by
#' contamination, retention time shifts, or loss of hydrophilic/hydrophobic
#' peptides.
#' is_a: MS:4000003 ! single value
#' relationship: has_metric_category MS:4000009 ! ID free metric
#' relationship: has_metric_category MS:4000017 ! chromatogram metric
#'
#' @param na.rm `logical(1)` whether to remove `NA` values (default `TRUE`)
#'
#' @return `numeric(1)`
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L),
#'     mz = c(112.2, 123.3),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Returns total area under TIC
#' areaUnderTic(chr)
setMethod("areaUnderTic", "Chromatograms", function(object, na.rm = TRUE, ...) {
    all_ints <- unlist(intensity(object), use.names = FALSE)
    res <- sum(all_ints, na.rm = na.rm)
    attr(res, "areaUnderTic") <- "MS:4000155"
    res
})
