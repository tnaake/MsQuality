################################################################################
################## PART 1: FULL CHROMATOGRAMS METRICS ##########################
################################################################################

#' @title Calculate the duration of the chromatographic run per chromatogram
#'
#' @description
#' The function `chromatogramDuration` calculates The range of retention time
#' over which the data was acquired of the chromatographic run, defined as the
#' difference between the maximum and minimum retention time.
#'
#' @details
#' The function returns the duration of the chromatographic run.
#'
#' id: MS:4000055
#'
#' @param chromatograms `Chromatograms` object

#' @param ... further arguments passed to `min` and `max` (e.g. `na.rm`)
#'
#' @return `numeric(1)`
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @importFrom Chromatograms rtime
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
#' chromatogramDuration(chr)
chromatogramDuration <- function(chromatograms, ...) {
    rts <- unlist(rtime(chromatograms), use.names = FALSE)
    rts <- as.numeric(rts)

    if (length(rts) == 0 || all(is.na(rts))) {
        res <- NA_real_
    } else {
        res <- max(rts, na.rm = TRUE) - min(rts, na.rm = TRUE)
    }

    attr(res, "chromatogramDuration") <- "MS:4000055"
    res
}

#' @title Number of Chromatograms
#'
#' @description
#' The function `chromatogramCount` computes the number of chromatograms in the
#' `Chromatograms` object.
#'
#' @details
#' This metric corresponds to the PSI:MS term:
#' MS:4000056 (number of chromatograms).
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
    attr(res, "chromatogramCount") <- "MS:4000056"
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



#' @title Retention time acquisition range
#'
#' @description
#' The function `rtAcquisitionRangeChromatograms` returns the minimum and
#' maximum retention time across all chromatograms.
#'
#' @details
#' The function returns a named vector of length 2 with the minimum and maximum
#' retention times.
#'
#' id: MS:4000070
#' name: retention time acquisition range
#' def: "Upper and lower limit of retention time at which spectra are recorded." [PSI:MS]
#' is_a: MS:4000004 ! n-tuple
#' relationship: has_metric_category MS:4000009 ! ID free metric
#' relationship: has_metric_category MS:4000012 ! single run based metric
#' relationship: has_metric_category MS:4000016 ! retention time metric
#' relationship: has_units UO:0000010 ! second
#'
#' @param chromatograms `Chromatograms` object
#' @param ... further arguments passed to `range`
#'
#' @return `numeric(2)` named "min" and "max"
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
#' rtAcquisitionRangeChromatograms(chr)
rtAcquisitionRangeChromatograms <- function(chromatograms, ...) {
    rts <- unlist(rtime(chromatograms), use.names = FALSE)
    if (length(rts) == 0) {
        res <- c(min = NA_real_, max = NA_real_)
    } else {
        res <- range(rts, ...)
        names(res) <- c("min", "max")
    }
    attr(res, "rtAcquisitionRangeChromatograms") <- "MS:4000070"
    res
}

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
maxIntensity <- function(chromatograms, ...) {
    ints <- intensity(chromatograms)
    # Ensure we properly extract numeric values even in parallel context
    ints <- as.numeric(unlist(ints, use.names = FALSE))
    res <- max(ints, na.rm = TRUE, ...)
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
intensityMean <- function(chromatograms, ...) {
    ints <- intensity(chromatograms)
    # Ensure we properly extract numeric values even in parallel context
    ints <- as.numeric(unlist(ints, use.names = FALSE))
    res <- mean(ints, na.rm = TRUE, ...)
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
intensitySd <- function(chromatograms, ...) {
    ints <- intensity(chromatograms)
    # Ensure we properly extract numeric values even in parallel context
    ints <- as.numeric(unlist(ints, use.names = FALSE))
    res <- sd(ints, na.rm = TRUE, ...)
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
    qt <- summary(res, ...)
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
intensityRange <- function(chromatograms, ...) {
    res <- unlist(intensity(chromatograms), use.names = FALSE)
    min_val <- min(res, ...)
    max_val <- max(res, ...)
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
#' Note: MS:4000061 is "MS1 density quantiles", not peak count.
#'
#' @param chromatograms `Chromatograms` object
#' @param ... further arguments (currently ignored)
#'
#' @return `integer` vector of length equal to number of chromatograms
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
#' ## Returns number of data points per chromatogram: 5, 0, 5
#' peakCount(chr)
peakCount <- function(chromatograms, ...) {
    res <- lengths(chromatograms)
    attr(res, "peakCount") <- "custom_metric:peak_count"
    res
}

#' @title Retention Time IQR across chromatograms
#'
#' @description
#' The function `rtIqrChromatograms` calculates the interquartile range (IQR)
#' of the retention times across all chromatograms.
#'
#' @details
#' The function returns the IQR of the retention times across all chromatograms.
#'
#' No specific PSI:MS term exists. It should be created.
#'
#' @param chromatograms `Chromatograms` object
#' @param ... further arguments passed to `IQR`
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
#' ## Returns IQR of retention times across all chromatograms
#' rtIqrChromatograms(chr)
rtIqrChromatograms <- function(chromatograms, ...) {
    res <- IQR(unlist(rtime(chromatograms), use.names = FALSE), ...)
    attr(res, "rtIqr") <- "custom_metric:rt_iqr"
    res
}

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
baselineIntensity <- function(chromatograms, probs = 0.05, ...) {
    res <- quantile(
        unlist(intensity(chromatograms), use.names = FALSE),
        probs = probs,
        ...
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

#' @title Fraction of Intensity in Equal RT Intervals
#'
#' @description
#' The function `intensityQuantileRtFraction` calculates the fraction of the
#' total intensity that falls within 4 equal time intervals (retention time
#' quartiles) across the entire run.
#'
#' @details
#' The retention time range (min to max) is divided into 4 equal bins. The sum
#' of intensities for data points falling into each bin is calculated and
#' divided by the total intensity.
#'
#' @param chromatograms `Chromatograms` object
#' @param ... further arguments passed to `range`
#'
#' @return `numeric(4)` named vector with fractions for Q1, Q2, Q3, Q4
#'
#' @author Philippine Louail
#'
#' @export
#'
intensityQuantileRtFraction <- function(chromatograms, ...) {
    ints <- unlist(intensity(chromatograms), use.names = FALSE)
    rts <- unlist(rtime(chromatograms), use.names = FALSE)

    if (length(rts) == 0 || all(is.na(rts)) || sum(ints, na.rm = TRUE) == 0) {
        res <- setNames(rep(NA_real_, 4), c("Q1", "Q2", "Q3", "Q4"))
        attr(res, "intensityQuantileRtFraction") <- "custom_metric:intensity_rt_quantile_fraction"
        return(res)
    }

    rt_range <- range(rts, na.rm = TRUE, ...)

    if (rt_range[2] == rt_range[1]) {
        res <- setNames(c(1, 0, 0, 0), c("Q1", "Q2", "Q3", "Q4"))
        attr(res, "intensityQuantileRtFraction") <- "custom_metric:intensity_rt_quantile_fraction"
        return(res)
    }

    ## Create 4 bins: [Start, 25%], (25%, 50%], (50%, 75%], (75%, End]
    breaks <- seq(rt_range[1], rt_range[2], length.out = 5)

    ## include.lowest = TRUE ensures the first point is included in Q1
    bins <- cut(rts, breaks = breaks, include.lowest = TRUE, labels = FALSE)

    bin_sums <- numeric(4)
    for (i in 1:4) {
        bin_sums[i] <- sum(ints[which(bins == i)], na.rm = TRUE)
    }

    total_intensity <- sum(ints, na.rm = TRUE)
    res <- bin_sums / total_intensity

    names(res) <- c("Q1", "Q2", "Q3", "Q4")
    attr(res, "intensityQuantileRtFraction") <- "custom_metric:intensity_rt_quantile_fraction"
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
medianIntensityRtIqr <- function(chromatograms, ...) {
    all_rts <- unlist(rtime(chromatograms), use.names = FALSE)
    all_ints <- unlist(intensity(chromatograms), use.names = FALSE)

    if (length(all_rts) == 0 || all(is.na(all_rts))) {
        res <- NA_real_
        attr(res, "medianIntensityRtIqr") <- "custom_metric:median_intensity_rt_iqr"
        return(res)
    }

    q_rt <- quantile(all_rts, probs = c(0.25, 0.75), na.rm = TRUE, ...)
    mask <- all_rts >= q_rt[1] & all_rts <= q_rt[2] & !is.na(all_rts)

    if (!any(mask)) {
        res <- NA_real_
        attr(res, "medianIntensityRtIqr") <- "custom_metric:median_intensity_rt_iqr"
        return(res)
    }
    res <- median(all_ints[mask], na.rm = TRUE, ...)
    attr(res, "medianIntensityRtIqr") <- "custom_metric:median_intensity_rt_iqr"
    res
}

#' @title Extent of Intensity across chromatograms
#'
#' @description
#' The function `extentIntensity` calculates the range of intensity values
#' (max - min) across all chromatograms.
#'
#' @details
#' The function returns `max - min` intensity across all chromatograms.
#'
#' No specific PSI:MS term exists. It should be created.
#'
#' @param chromatograms `Chromatograms` object
#' @param ... further arguments passed to `min` and `max`
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
#' ## Returns max - min intensity across all chromatograms: 1120
#' extentIntensity(chr)
extentIntensity <- function(chromatograms, ...) {
    all_ints <- unlist(intensity(chromatograms), use.names = FALSE)
    if (length(all_ints) == 0 || all(is.na(all_ints))) {
        res <- NA_real_
    } else {
        res <- max(all_ints, ...) - min(all_ints, ...)
    }
    attr(res, "extentIntensity") <- "custom_metric:extent_intensity"
    res
}

#' @title Log Ratio of Intensity Quartiles across chromatograms
#'
#' @description
#' The function `intensityQuartileToQuartileLogRatio` calculates the log2 ratio
#' between specific intensity quartiles across all chromatograms.
#'
#' @details
#' It calculates `log2(Q_a / Q_b)` across all chromatograms.
#'
#' No specific PSI:MS term exists. It should be created.
#'
#' @param chromatograms `Chromatograms` object
#' @param probs `numeric(2)` probabilities for the quartiles (default c(0.25, 0.75))
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
#' ## Returns log2(Q3/Q1) ratio across all chromatograms
#' intensityQuartileToQuartileLogRatio(chr)
#' ## Use different quartiles
#' intensityQuartileToQuartileLogRatio(chr, probs = c(0.10, 0.90))
intensityQuartileToQuartileLogRatio <- function(
    chromatograms,
    probs = c(0.25, 0.75),
    ...
) {
    all_ints <- unlist(intensity(chromatograms), use.names = FALSE)
    if (length(all_ints) == 0 || all(is.na(all_ints))) {
        res <- NA_real_
        attr(res, "intensityQuartileToQuartileLogRatio") <- "custom_metric:intensity_quartile_to_quartile_log_ratio"
        return(res)
    }
    qs <- quantile(all_ints, probs = probs, na.rm = TRUE, ...)
    if (any(qs <= 0)) {
        res <- NA_real_
        attr(res, "intensityQuartileToQuartileLogRatio") <- "custom_metric:intensity_quartile_to_quartile_log_ratio"
        return(res)
    }
    res <- log2(qs[2] / qs[1])
    attr(res, "intensityQuartileToQuartileLogRatio") <- "custom_metric:intensity_quartile_to_quartile_log_ratio"
    res

}

#' @title Area Under Intensity-RT Quantiles
#'
#' @description
#' The function `areaUnderIntensityRtQuantiles` calculates the Area Under the
#' Curve (AUC) for 4 equal time intervals (retention time quartiles) across
#' the run.
#'
#' @details
#' The function splits the chromatogram into 4 equal time bins. To ensure accurate
#' area calculation at the boundaries, intensity values are linearly interpolated
#' exactly at the cut-off times.
#'
#' @param chromatograms `Chromatograms` object
#' @param ... further arguments
#'
#' @return `numeric(4)` named vector with Q1-Q4 area values
#'
#' @author Philippine Louail
#'
#' @importFrom stats approx
#' @export
#'
areaUnderIntensityRtQuantiles <- function(chromatograms, ...) {
    ints <- unlist(intensity(chromatograms), use.names = FALSE)
    rts <- unlist(rtime(chromatograms), use.names = FALSE)

    if (length(rts) < 2 || all(is.na(rts))) {
        res <- setNames(rep(NA_real_, 4), c("Q1", "Q2", "Q3", "Q4"))
        attr(res, "areaUnderIntensityRtQuantiles") <- "custom_metric:area_under_intensity_rt_quantiles"
        return(res)
    }

    ord <- order(rts)
    rts <- rts[ord]
    ints <- ints[ord]

    rt_range <- range(rts, na.rm = TRUE)

    cuts <- seq(rt_range[1], rt_range[2], length.out = 5)[2:4]

    new_rts <- rts
    new_ints <- ints

    for (ct in cuts) {
        interp_val <- approx(rts, ints, xout = ct)$y
        if (!is.na(interp_val)) {
            new_rts <- c(new_rts, ct)
            new_ints <- c(new_ints, interp_val)
        }
    }
        ord_new <- order(new_rts)
    final_rts <- new_rts[ord_new]
    final_ints <- new_ints[ord_new]


    areas <- (final_ints[-1] + final_ints[-length(final_ints)]) / 2 * diff(final_rts)

    midpoints <- (final_rts[-1] + final_rts[-length(final_rts)]) / 2

    breaks <- seq(rt_range[1], rt_range[2], length.out = 5)
    bins <- cut(midpoints, breaks = breaks, include.lowest = TRUE, labels = FALSE)

    res <- numeric(4)
    for (i in 1:4) {
        res[i] <- sum(areas[which(bins == i)], na.rm = TRUE)
    }

    names(res) <- c("Q1", "Q2", "Q3", "Q4")
    attr(res, "areaUnderIntensityRtQuantiles") <- "custom_metric:area_under_intensity_rt_quantiles"
    res
}

################################################################################
#################### PART 2: XIC/EIC METRICS FROM PSI-MS #######################
################################################################################

#' @title Distribution of Peak Widths (FWHM)
#'
#' @description
#' The function `xicFwhmQuantiles` calculates the Full Width at Half Maximum (FWHM)
#' for each chromatogram and returns the quantiles (min, 25%, median, 75%, max).
#'
#' @details
#' The FWHM is calculated by finding the maximum intensity peak, determining
#' 50% of that intensity, and linearly interpolating the time difference between
#' the left and right crossing points.
#'
#' This is a custom metric analogous to MS:4000120 (median retention time width).
#'
#' @param chromatograms `Chromatograms` object
#' @param ... further arguments
#'
#' @return `numeric(5)` Named vector of quantiles.
#'
#' @author Philippine Louail
#'
#' @importFrom stats approx quantile
#' @importFrom utils tail
#' @export
xicFwhmQuantiles <- function(chromatograms, ...) {

    fwhms <- vapply(
        seq_along(chromatograms),
        function(i) {
            rts <- rtime(chromatograms)[[i]]
            ints <- intensity(chromatograms)[[i]]

            ## Basic checks
            if (length(ints) < 3 || all(is.na(ints))) return(NA_real_)
            max_int <- max(ints, na.rm = TRUE)
            if (max_int == 0) return(NA_real_)

            max_idx <- which.max(ints)
            half_max <- max_int / 2

            ## Left side (rising): Find last point below half_max before peak
            left_candidates <- which(ints[1:max_idx] < half_max)
            if (length(left_candidates) == 0) return(NA_real_)
            left_idx <- tail(left_candidates, 1)

            ## Right side (falling): Find first point below half_max after peak
            right_candidates <- which(ints[max_idx:length(ints)] < half_max)
            if (length(right_candidates) == 0) return(NA_real_)
            right_idx <- max_idx + right_candidates[1] - 1

            ## Interpolation Logic
            ## Left: ints[left_idx] < half_max < ints[left_idx+1] (Increasing)
            rt_left <- approx(x = ints[c(left_idx, left_idx + 1)],
                              y = rts[c(left_idx, left_idx + 1)],
                              xout = half_max)$y

            ## Right: ints[right_idx-1] > half_max > ints[right_idx] (Decreasing)
            ## We MUST reverse the vectors so 'x' (intensity) is increasing for approx()
            rt_right <- approx(x = ints[c(right_idx, right_idx - 1)],
                               y = rts[c(right_idx, right_idx - 1)],
                               xout = half_max)$y

            return(rt_right - rt_left)
        },
        numeric(1)
    )

    res <- quantile(fwhms, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    attr(res, "xicFwhmQuantiles") <- "custom_metric:xic_fwhm_distribution"
    res
}

#' @title XIC50 fraction
#'
#' @description
#' The function `xic50Fraction` calculates the fraction of XIC peaks that account
#' for the top half of all XIC-FWHM values.
#'
#' @details
#' This metric characterizes the distribution of precursor peak widths, similar
#' to how N50 characterizes contig sizes in genome assembly.
#'
#' id: MS:4000050
#' name: XIC50 fraction
#' def: "The number of XIC that account for the top half of all XIC-FWHM divided
#' by the number of all XIC." [PSI:MS]
#' is_a: MS:4000003 ! single value
#' relationship: has_metric_category MS:4000009 ! ID free metric
#' relationship: has_metric_category MS:4000012 ! single run based metric
#' relationship: has_units UO:0000191 ! fraction
#'
#' @param chromatograms `Chromatograms` object
#' @param ... further arguments
#'
#' @return `numeric(1)` fraction value
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
#'     data.frame(rtime = c(4.0, 4.5, 5.0, 5.5, 6.0),
#'                intensity = c(50, 150, 300, 200, 75)),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Returns fraction of XICs accounting for top half of FWHM
#' xic50Fraction(chr)
xic50Fraction <- function(chromatograms, ...) {
    fwhm_vals <- vapply(
        seq_along(chromatograms),
        function(i) {
            rts <- rtime(chromatograms)[[i]]
            ints <- intensity(chromatograms)[[i]]

            if (length(rts) < 3 || all(is.na(ints))) {
                return(NA_real_)
            }

            max_idx <- which.max(ints)
            max_int <- ints[max_idx]
            half_max <- max_int / 2

            if (max_idx > 1 && max_idx < length(ints)) {
                left_idx <- max_idx - 1
                right_idx <- max_idx + 1

                while (left_idx > 1 && ints[left_idx] > half_max) {
                    left_idx <- left_idx - 1
                }
                while (right_idx < length(ints) && ints[right_idx] > half_max) {
                    right_idx <- right_idx + 1
                }

                if (left_idx < max_idx && right_idx > max_idx) {
                    rt_left <- approx(ints[c(left_idx, left_idx + 1)],
                                      rts[c(left_idx, left_idx + 1)],
                                      half_max)$y
                    rt_right <- approx(x = ints[c(right_idx, right_idx - 1)],
                        y = rts[c(right_idx, right_idx - 1)],
                        xout = half_max
                    )$y
                    return(rt_right - rt_left)
                }
            }
            NA_real_
        },
        numeric(1)
    )

    fwhm_vals <- fwhm_vals[!is.na(fwhm_vals)]
    if (length(fwhm_vals) < 2) {
        res <- NA_real_
    } else {
        ## Sort FWHM values in descending order
        sorted_fwhm <- sort(fwhm_vals, decreasing = TRUE)
        cumsum_fwhm <- cumsum(sorted_fwhm)
        total_fwhm <- sum(fwhm_vals)

        ## Find number of XICs accounting for top half
        n_top_half <- sum(cumsum_fwhm <= (total_fwhm / 2)) + 1
        res <- n_top_half / length(fwhm_vals)
    }
    attr(res, "xic50Fraction") <- "MS:4000050"
    res
}

#' @title XIC height quantile ratios
#'
#' @description
#' The function `xicHeightQuantileRatios` calculates the log ratios of successive
#' XIC height quantiles across all chromatograms.
#'
#' @details
#' Returns log ratios of Q2/Q1, Q3/Q2, and max/Q3 of peak heights.
#'
#' id: MS:4000182
#' name: XIC-Height quantile ratios
#' def: "The log ratio of successive XIC height quantiles. A value triplet
#' represents the original QuaMeter metrics, the log ratios of XIC-Height-Q2 to
#' XIC-Height-Q1, XIC-Height-Q3 to XIC-Height-Q2, XIC-Height max to XIC-Height-Q3." [PSI:MS]
#' is_a: MS:4000004 ! n-tuple
#' relationship: has_metric_category MS:4000009 ! ID free metric
#' relationship: has_metric_category MS:4000012 ! single run based metric
#' relationship: has_metric_category MS:4000018 ! XIC metric
#' relationship: has_units UO:0000191 ! fraction
#'
#' @param chromatograms `Chromatograms` object
#' @param probs `numeric` vector of probabilities (default c(0.25, 0.5, 0.75))
#' @param ... further arguments passed to `quantile`
#'
#' @return `numeric(3)` named vector with log ratios
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
#' ## Returns log ratios of Q2/Q1, Q3/Q2, max/Q3
#' xicHeightQuantileRatios(chr)
xicHeightQuantileRatios <- function(
    chromatograms,
    probs = c(0.25, 0.5, 0.75),
    ...
) {
    all_max_ints <- vapply(
        intensity(chromatograms),
        function(x) {
            if (length(x) == 0 || all(is.na(x))) {
                NA_real_
            } else {
                max(x, na.rm = TRUE)
            }
        },
        numeric(1)
    )

    all_max_ints <- all_max_ints[!is.na(all_max_ints)]

    if (length(all_max_ints) < 2) {
        res <- setNames(rep(NA_real_, 3), c("Q2/Q1", "Q3/Q2", "max/Q3"))
    } else {
        q_vals <- quantile(all_max_ints, probs = c(probs, 1), na.rm = TRUE, ...)

        ## Calculate log ratios
        res <- c(
            log2(q_vals[2] / q_vals[1]),
            log2(q_vals[3] / q_vals[2]),
            log2(q_vals[4] / q_vals[3])
        )
        names(res) <- c("Q2/Q1", "Q3/Q2", "max/Q3")
    }
    attr(res, "xicHeightQuantileRatios") <- "MS:4000182"
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
#' @aliases ticQuantileRtFraction, Spectra-method
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
setMethod("areaUnderTic", "Chromatograms", function(object, ...) {
    all_ints <- unlist(intensity(object), use.names = FALSE)
    res <- sum(all_ints, ...)
    attr(res, "areaUnderTic") <- "MS:4000155"
    res
})
