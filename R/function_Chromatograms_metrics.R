#' @title Maximum intensity per chromatogram
#'
#' @description
#' The function `maxIntensity` returns the maximum intensity value observed
#' within each chromatogram. \cr
#'
#' The metric is calculated as follows: \cr
#' (1) for each chromatogram in the `Chromatograms` object, the intensity
#' values are extracted, \cr
#' (2) the maximum value per chromatogram is obtained and returned.
#'
#' @param chromatograms `Chromatograms` object
#' @param na.rm `logical(1)` whether to remove `NA` values (default `TRUE`)
#' @param ... further arguments passed to `max`
#'
#' @return `numeric` of length equal to `length(chromatograms)`, one maximum
#'   intensity value per chromatogram
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
#' ## Returns max intensity per chromatogram: c(400, 1200)
#' maxIntensity(chr)
maxIntensity <- function(chromatograms, na.rm = TRUE, ...) {
    ints_list <- intensity(chromatograms)
    res <- vapply(ints_list, function(ints) {
        if (length(ints) == 0 || all(is.na(ints))) return(NA_real_)
        max(ints, na.rm = na.rm)
    }, numeric(1))
    attr(res, "maxIntensity") <- "custom_metric:max_intensity"
    res
}

#' @title Intensity statistics (Mean) per chromatogram
#'
#' @description
#' The function `intensityMean` calculates the mean of intensity values
#' within each chromatogram. \cr
#'
#' The metric is calculated as follows: \cr
#' (1) for each chromatogram in the `Chromatograms` object, the intensity
#' values are extracted, \cr
#' (2) the arithmetic mean per chromatogram is calculated and returned.
#'
#' @param chromatograms `Chromatograms` object
#' @param na.rm `logical(1)` whether to remove `NA` values (default `TRUE`)
#' @param ... further arguments passed to `mean`
#'
#' @return `numeric` of length equal to `length(chromatograms)`
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
#' ## Returns mean intensity per chromatogram
#' intensityMean(chr)
intensityMean <- function(chromatograms, na.rm = TRUE, ...) {
    ints_list <- intensity(chromatograms)
    res <- vapply(ints_list, function(ints) {
        if (length(ints) == 0) return(NA_real_)
        mean(ints, na.rm = na.rm)
    }, numeric(1))
    attr(res, "intensityMean") <- "custom_metric:intensity_mean"
    res
}

#' @title Intensity statistics (Standard Deviation) per chromatogram
#'
#' @description
#' The function `intensitySd` calculates the standard deviation of intensity
#' values within each chromatogram. \cr
#'
#' The metric is calculated as follows: \cr
#' (1) for each chromatogram in the `Chromatograms` object, the intensity
#' values are extracted, \cr
#' (2) the standard deviation per chromatogram is calculated and returned.
#'
#' @param chromatograms `Chromatograms` object
#' @param na.rm `logical(1)` whether to remove `NA` values (default `TRUE`)
#' @param ... further arguments passed to `sd`
#'
#' @return `numeric` of length equal to `length(chromatograms)`
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
#' ## Returns standard deviation of intensities per chromatogram
#' intensitySd(chr)
intensitySd <- function(chromatograms, na.rm = TRUE, ...) {
    ints_list <- intensity(chromatograms)
    res <- vapply(ints_list, function(ints) {
        if (length(ints) < 2) return(NA_real_)
        sd(ints, na.rm = na.rm)
    }, numeric(1))
    attr(res, "intensitySd") <- "custom_metric:intensity_sd"
    res
}

#' @title Intensity Quartiles per chromatogram
#'
#' @description
#' The function `intensityQuartiles` calculates the minimum, 1st quartile,
#' median, 3rd quartile, and maximum of intensity values within each
#' chromatogram. \cr
#'
#' The metric is calculated as follows: \cr
#' (1) for each chromatogram, the intensity values are extracted, \cr
#' (2) the summary statistics (Min, 1st Qu., Median, Mean, 3rd Qu., Max)
#' are calculated and returned as a row of the result matrix.
#'
#' @param chromatograms `Chromatograms` object
#' @param ... further arguments passed to `summary`
#'
#' @return `matrix` with `nrow` equal to `length(chromatograms)` and 6 columns
#'   (Min, 1st Qu., Median, Mean, 3rd Qu., Max)
#'
#' @author Philippine Louail
#'
#' @importFrom stats setNames
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
#' ## Returns a 2x6 matrix with summary statistics per chromatogram
#' intensityQuartiles(chr)
intensityQuartiles <- function(chromatograms, ...) {
    ints_list <- intensity(chromatograms)
    col_names <- c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max")
    res <- t(vapply(ints_list, function(ints) {
        if (length(ints) == 0 || all(is.na(ints)))
            return(setNames(rep(NA_real_, 6), col_names))
        q <- quantile(ints, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
        setNames(c(q[1], q[2], q[3], mean(ints, na.rm = TRUE), q[4], q[5]),
                 col_names)
    }, numeric(6)))
    colnames(res) <- col_names
    attr(res, "intensityQuartiles") <- "custom_metric:intensity_quartiles"
    res
}

#' @title Intensity Range per chromatogram
#'
#' @description
#' The function `intensityRange` calculates the range (min, max) of intensity
#' values within each chromatogram. \cr
#'
#' The metric is calculated as follows: \cr
#' (1) for each chromatogram, the intensity values are extracted, \cr
#' (2) the minimum and maximum values per chromatogram are obtained and
#' returned as a row of the result matrix.
#'
#' @param chromatograms `Chromatograms` object
#' @param na.rm `logical(1)` whether to remove `NA` values (default `TRUE`)
#' @param ... further arguments passed to `range`
#'
#' @return `matrix` with `nrow` equal to `length(chromatograms)` and
#'   2 columns (min, max)
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
#' ## Returns a 2x2 matrix with (min, max) per chromatogram
#' intensityRange(chr)
intensityRange <- function(chromatograms, na.rm = TRUE, ...) {
    ints_list <- intensity(chromatograms)
    res <- t(vapply(ints_list, function(ints) {
        if (length(ints) == 0 || all(is.na(ints)))
            return(c(min = NA_real_, max = NA_real_))
        c(min = min(ints, na.rm = na.rm), max = max(ints, na.rm = na.rm))
    }, numeric(2)))
    colnames(res) <- c("min", "max")
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
#'     removed before counting (default `FALSE`)
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
#' peakCount(chr)
peakCount <- function(chromatograms, na.rm = FALSE, ...) {
    if (length(chromatograms) == 0) return(integer(0))
    if (na.rm) {
        res <- vapply(intensity(chromatograms),
                      function(x) sum(!is.na(x)), integer(1L))
    } else {
        res <- lengths(intensity(chromatograms))
    }
    attr(res, "peakCount") <- "custom_metric:peak_count"
    res
}

#' @title Baseline Intensity across chromatograms
#'
#' @description
#' The function `baselineIntensity` estimates the baseline intensity as the
#' 5th percentile of intensity values across all chromatograms. \cr
#'
#' The metric is calculated as follows: \cr
#' (1) the intensity values from all chromatograms are extracted and
#' unlisted into a single numeric vector, \cr
#' (2) the quantile at the specified probability (default 5th percentile)
#' is calculated and returned.
#'
#' @details
#' In chromatograms with many zero-intensity points the chosen quantile may
#' fall within the zero-valued portion of the distribution, causing the
#' function to return \code{0}. This can propagate to downstream metrics that
#' rely on a non-zero baseline (e.g. \code{peakProminence} divides by the
#' baseline and would return \code{Inf}). Consider increasing \code{probs}
#' (e.g. \code{0.10} or higher) or pre-filtering zero-intensity data points
#' when working with zero-dominated chromatograms.
#'
#' This function returns a single summary value per chromatogram. For a
#' point-wise baseline estimate (one value per retention-time point),
#' see \code{\link[MsCoreUtils]{estimateBaseline}} which offers SNIP,
#' TopHat, ConvexHull, and median methods.
#'
#' @param chromatograms `Chromatograms` object
#' @param probs `numeric(1)` quantile probability (default 0.05)
#' @param na.rm `logical(1)` whether to remove `NA` values (default `TRUE`)
#' @param ... further arguments passed to `quantile`
#'
#' @return `numeric` of length equal to `length(chromatograms)`, one baseline
#' intensity per chromatogram
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
    ints_list <- intensity(chromatograms)
    res <- vapply(ints_list, function(ints) {
        if (length(ints) == 0) return(NA_real_)
        unname(quantile(ints, probs = probs, na.rm = na.rm))
    }, numeric(1))
    attr(res, "baselineIntensity") <- "custom_metric:baseline_intensity"
    res
}

#' @title Signal-to-Noise Ratio per chromatogram
#'
#' @description
#' The function `signalToNoiseRatio` estimates a classical signal-to-noise
#' ratio for each chromatogram using the MAD-based noise estimate from
#' \code{MsCoreUtils::noise()}. \cr
#'
#' The metric is calculated as follows (per chromatogram): \cr
#' (1) the retention time and intensity values are extracted, \cr
#' (2) a noise level is estimated at every data point using the Median
#' Absolute Deviation (MAD) via \code{MsCoreUtils::noise()}, \cr
#' (3) the signal is taken as the maximum intensity, \cr
#' (4) the noise is summarised as the median of the MAD-estimated noise
#' vector, \cr
#' (5) the ratio \code{signal / noise} is returned.
#'
#' @details
#' The MAD method (\code{MsCoreUtils::noise(x, y, method = "MAD")}) computes
#' a single global Median Absolute Deviation across all intensity values in
#' the chromatogram and replicates it to every data point. Because the
#' estimate is global, large peaks inflate the MAD, leading to a
#' **conservative** (under-estimated) S/N ratio. This is acceptable for
#' quality-control purposes where consistency across samples matters more
#' than an absolute noise floor.
#'
#' The MAD estimator is statistically unstable with very few data points
#' (roughly < 5); results from very short chromatograms should be
#' interpreted with caution.
#'
#' Returns \code{NA_real_} when the noise estimate is zero (e.g. constant
#' signal) or when the chromatogram is empty.
#'
#' @param chromatograms `Chromatograms` object
#' @param ... currently not used but included for consistency with other
#'     metric functions
#'
#' @return `numeric` of length equal to `length(chromatograms)`, one
#' signal-to-noise ratio per chromatogram
#'
#' @author Philippine Louail
#'
#' @importFrom MsCoreUtils noise
#'
#' @export
#'
signalToNoiseRatio <- function(chromatograms, ...) {
    rts_list <- rtime(chromatograms)
    ints_list <- intensity(chromatograms)
    res <- vapply(seq_along(ints_list), function(i) {
        rts <- rts_list[[i]]
        ints <- ints_list[[i]]
        if (length(ints) == 0 || all(is.na(ints))) return(NA_real_)
        noise_est <- noise(rts, ints, method = "MAD")
        noise_summary <- median(noise_est, na.rm = TRUE)
        if (is.na(noise_summary) || noise_summary == 0) return(NA_real_)
        max(ints, na.rm = TRUE) / noise_summary
    }, numeric(1))
    attr(res, "signalToNoiseRatio") <- "custom_metric:signal_to_noise_ratio"
    res
}

#' @name xicFwhm
#' @rdname xicFwhm
#' @title Full Width at Half Maximum (FWHM) per Chromatogram
#'
#' @description
#' The function `xicFwhm` calculates the Full Width at Half Maximum (FWHM)
#' for each chromatogram in the `Chromatograms` object. \cr
#'
#' The metric is calculated as follows (per chromatogram): \cr
#' (1) the retention time and intensity values are extracted, \cr
#' (2) if \code{peakBoundary} is provided, the data is subset to the peak
#' region, \cr
#' (3) the maximum intensity and its index are determined, \cr
#' (4) the half-maximum intensity (50\% of maximum) is calculated, \cr
#' (5) the left and right crossing points where intensity equals the
#' half-maximum are found by linear interpolation, \cr
#' (6) the FWHM is returned as the difference between the right and left
#' crossing retention times.
#'
#' @details
#' This metric is analogous to MS:4000051 (XIC-FWHM quantiles). \cr
#'
#' If \code{peakBoundary} is provided, the calculation is restricted to the peak
#' region defined by those boundaries. \code{peakBoundary} can be a `matrix`
#' (one row per chromatogram, as returned by `peakBoundary()`) or a
#' `numeric(2)` vector applied to all chromatograms.
#'
#' @param chromatograms `Chromatograms` object
#' @param peakBoundary optional peak boundary: either a `matrix` with columns
#'   `left_boundary` and `right_boundary` (one row per chromatogram), or a
#'   `numeric(2)` named vector applied to all chromatograms.
#' @param ... further arguments
#'
#' @return `numeric` of length equal to `length(chromatograms)`. Returns `NA`
#'   for chromatograms where FWHM cannot be calculated.
#'
#' @author Philippine Louail
#'
#' @importFrom stats approx
#' @importFrom utils tail
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
#' xicFwhm(chr)
#'
#' ## Use pre-computed peak boundaries for efficiency
#' pb <- peakBoundary(chr)
#' xicFwhm(chr, peakBoundary = pb)
NULL

#' Internal helper: compute FWHM for a single chromatogram
#' @noRd
.xicFwhmSingle <- function(rts, ints, peakBoundary = NULL) {
    if (!is.null(peakBoundary) && !any(is.na(peakBoundary))) {
        left_rt <- peakBoundary["left_boundary"]
        right_rt <- peakBoundary["right_boundary"]
        mask <- rts >= left_rt & rts <= right_rt
        rts <- rts[mask]
        ints <- ints[mask]
    }
    if (length(ints) < 3 || all(is.na(ints))) return(NA_real_)
    max_int <- max(ints, na.rm = TRUE)
    if (max_int == 0) return(NA_real_)
    max_idx <- which.max(ints)
    half_max <- max_int / 2
    left_candidates <- which(ints[1:max_idx] < half_max)
    if (length(left_candidates) == 0) return(NA_real_)
    left_idx <- tail(left_candidates, 1)
    right_candidates <- which(ints[max_idx:length(ints)] < half_max)
    if (length(right_candidates) == 0) return(NA_real_)
    right_idx <- max_idx + right_candidates[1] - 1
    left_ints <- ints[c(left_idx, left_idx + 1)]
    left_rts <- rts[c(left_idx, left_idx + 1)]
    right_ints <- ints[c(right_idx - 1, right_idx)]
    right_rts <- rts[c(right_idx - 1, right_idx)]
    if (any(is.na(left_ints)) || any(is.na(left_rts)) ||
        any(is.na(right_ints)) || any(is.na(right_rts))) return(NA_real_)
    rt_left <- approx(x = left_ints, y = left_rts, xout = half_max)$y
    rt_right <- approx(x = right_ints, y = right_rts, xout = half_max)$y
    rt_right - rt_left
}

#' @rdname xicFwhm
#' @export
xicFwhm <- function(chromatograms, peakBoundary = NULL, ...) {
    n <- length(chromatograms)
    rts_list <- rtime(chromatograms)
    ints_list <- intensity(chromatograms)
    res <- vapply(seq_len(n), function(i) {
        pb_i <- NULL
        if (!is.null(peakBoundary)) {
            pb_i <- if (is.matrix(peakBoundary)) peakBoundary[i, ] else peakBoundary
        }
        .xicFwhmSingle(rts_list[[i]], ints_list[[i]], pb_i)
    }, numeric(1))
    attr(res, "xicFwhm") <- "custom_metric:xic_fwhm"
    res
}

#' @name peakBoundary
#' @rdname peakBoundary
#' @title Peak Boundary per Chromatogram
#'
#' @description
#' The function `peakBoundary` finds the left and right retention time boundaries
#' of the main peak in each chromatogram (EIC). \cr
#'
#' The metric is calculated as follows (per chromatogram): \cr
#' (1) the retention time and intensity values are extracted, \cr
#' (2) the peak apex (maximum intensity) is identified, \cr
#' (3) the baseline intensity is estimated from the \code{baselineQuantile}
#' of all intensities, \cr
#' (4) local minima (valleys) on each side of the apex are found via
#' \code{MsCoreUtils::valleys()}, \cr
#' (5) boundaries are validated by checking if intensities are near
#' baseline level, \cr
#' (6) if valley-based boundaries are not at baseline, a fallback
#' threshold method is used, \cr
#' (7) the left and right boundary retention times are returned.
#'
#' @details
#' The baseline is estimated from a lower quantile of intensities (default 10th
#' percentile), and the threshold is calculated relative to peak height above
#' this baseline.
#'
#' @param chromatograms `Chromatograms` object
#' @param threshold `numeric(1)` fraction of peak height above baseline used
#'   as fallback threshold (default 0.1 = 10\%).
#' @param baselineThreshold `numeric(1)` maximum acceptable intensity at
#'   boundaries as a fraction of peak height above baseline. Default is 0.1.
#' @param baselineQuantile `numeric(1)` quantile used to estimate the baseline
#'   intensity (default 0.1 = 10th percentile of intensities).
#' @param ... further arguments (currently unused)
#'
#' @return `matrix` with `nrow` equal to `length(chromatograms)` and 2 columns
#'   (`left_boundary`, `right_boundary`). Returns `NA` values for chromatograms
#'   where boundaries cannot be determined.
#'
#' @author Philippine Louail
#'
#' @importFrom MsCoreUtils valleys
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L),
#'     mz = c(100.0, 200.0),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = c(1, 2, 3, 4, 5, 6, 7),
#'                intensity = c(0, 10, 50, 100, 50, 10, 0))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' ## Find peak boundaries for each chromatogram
#' peakBoundary(chr)
NULL

#' Internal helper: compute peak boundaries for a single chromatogram
#' @noRd
.peakBoundarySingle <- function(rts, ints, threshold = 0.1,
                                 baselineThreshold = 0.1,
                                 baselineQuantile = 0.1) {
    n <- length(ints)
    na_result <- c(left_boundary = NA_real_, right_boundary = NA_real_)
    if (n < 3 || all(is.na(ints))) return(na_result)
    max_int <- max(ints, na.rm = TRUE)
    if (max_int == 0) return(na_result)
    max_idx <- which.max(ints)
    baseline_int <- quantile(ints, probs = baselineQuantile, na.rm = TRUE)
    peak_height <- max_int - baseline_int
    baseline_thresh <- baseline_int + peak_height * baselineThreshold
    v <- valleys(ints, max_idx)
    left_idx <- if ("left" %in% colnames(v)) v[1L, "left"] else 1L
    right_idx <- if ("right" %in% colnames(v)) v[1L, "right"] else n
    left_ok <- !is.na(ints[left_idx]) && ints[left_idx] <= baseline_thresh &&
               !(left_idx > 1 && is.na(ints[left_idx - 1]))
    right_ok <- !is.na(ints[right_idx]) && ints[right_idx] <= baseline_thresh &&
                !(right_idx < n && is.na(ints[right_idx + 1]))
    if (!left_ok || !right_ok) {
        thresh_val <- baseline_int + peak_height * threshold
        left_cand <- which(ints[seq_len(max_idx)] <= thresh_val)
        right_cand <- which(ints[max_idx:n] <= thresh_val)
        left_idx <- if (length(left_cand)) max(left_cand) else 1L
        right_idx <- if (length(right_cand)) max_idx + min(right_cand) - 1L else n
    }
    c(left_boundary = rts[left_idx], right_boundary = rts[right_idx])
}

#' @rdname peakBoundary
#' @export
peakBoundary <- function(chromatograms,
                         threshold = 0.1,
                         baselineThreshold = 0.1,
                         baselineQuantile = 0.1,
                         ...) {
    n <- length(chromatograms)
    rts_list <- rtime(chromatograms)
    ints_list <- intensity(chromatograms)
    res <- t(vapply(seq_len(n), function(i) {
        .peakBoundarySingle(rts_list[[i]], ints_list[[i]],
                            threshold, baselineThreshold, baselineQuantile)
    }, numeric(2)))
    colnames(res) <- c("left_boundary", "right_boundary")
    attr(res, "peakBoundary") <- "custom_metric:peak_boundary"
    res
}

#' @title Peak Width per Chromatogram
#'
#' @description
#' The function `peakWidth` calculates the width of the main peak in each
#' chromatogram. \cr
#'
#' The metric is calculated as follows (per chromatogram): \cr
#' (1) if not provided, peak boundaries are determined using
#' \code{peakBoundary()}, \cr
#' (2) the peak width is calculated as the difference between the right
#' and left boundary retention times and returned.
#'
#' @details
#' This is different from FWHM which measures width at 50\% of max intensity.
#'
#' @param chromatograms `Chromatograms` object
#' @param peakBoundary optional peak boundary: either a `matrix` with columns
#'   `left_boundary` and `right_boundary` (one row per chromatogram, as returned
#'   by `peakBoundary()`), or a `numeric(2)` named vector applied to all
#'   chromatograms. If not provided, boundaries are calculated automatically.
#' @param ... further arguments passed to `peakBoundary()`
#'
#' @return `numeric` of length equal to `length(chromatograms)`. Returns `NA`
#'   for chromatograms where width cannot be calculated.
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L),
#'     mz = c(100.0, 200.0),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = c(1, 2, 3, 4, 5, 6, 7),
#'                intensity = c(0, 10, 50, 100, 50, 10, 0))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' peakWidth(chr)
#'
#' ## Use pre-computed peak boundaries for efficiency
#' pb <- peakBoundary(chr)
#' peakWidth(chr, peakBoundary = pb)
peakWidth <- function(chromatograms, peakBoundary = NULL, ...) {
    if (is.null(peakBoundary)) {
        peakBoundary <- peakBoundary(chromatograms, ...)
    }
    if (is.matrix(peakBoundary)) {
        res <- unname(peakBoundary[, "right_boundary"] -
                      peakBoundary[, "left_boundary"])
    } else {
        ## backward compat: numeric(2) vector
        res <- unname(peakBoundary["right_boundary"] -
                      peakBoundary["left_boundary"])
    }
    attr(res, "peakWidth") <- "custom_metric:peak_width"
    res
}

#' @name gaussianSimilarity
#' @rdname gaussianSimilarity
#' @title Gaussian Similarity (Peak Shape Quality) per Chromatogram
#'
#' @description
#' The function `gaussianSimilarity` assesses the shape quality of each
#' chromatographic peak by fitting a beta-distribution curve via
#' \code{MetaboCoreUtils::betaValues()} and returning two diagnostics. \cr
#'
#' The metric is calculated as follows (per chromatogram): \cr
#' (1) the retention time and intensity values are extracted, \cr
#' (2) if not provided, peak boundaries are determined using
#' \code{peakBoundary()}, \cr
#' (3) the data is subset to the peak region defined by the
#' boundaries, \cr
#' (4) \code{MetaboCoreUtils::betaValues()} is called on the peak region
#' with a set of skew parameters (\code{shape1 = 3, 3.5, 4, 4.5, 5},
#' \code{shape2 = 5}), \cr
#' (5) two values are returned: \code{gaussian_similarity} (maximum
#' correlation between the observed peak and the best-fitting beta
#' curve; values close to 1 indicate a symmetric, bell-shaped peak)
#' and \code{gaussian_residuals} (standard deviation of the residuals
#' after normalising and subtracting the best-fit curve; lower values
#' indicate a cleaner peak shape).
#'
#' @details
#' The underlying \code{betaValues()} function (Kumler et al. 2023) compares
#' the observed chromatographic peak to a family of beta-distribution curves
#' of varying skew.  The first returned value is the Pearson correlation
#' between the observed and best-fit curve, and the second is the standard
#' deviation of the normalised residuals.
#'
#' Requires at least 5 data points within the peak region.
#'
#' @references
#' Kumler W, Hazelton B J and Ingalls A E (2023) "Picky with peakpicking:
#' assessing chromatographic peak quality with simple metrics in metabolomics"
#' BMC Bioinformatics 24(1):404. doi: 10.1186/s12859-023-05533-4
#'
#' @param chromatograms `Chromatograms` object
#' @param peakBoundary optional peak boundary: either a `matrix` with columns
#'   `left_boundary` and `right_boundary` (one row per chromatogram), or a
#'   `numeric(2)` named vector applied to all chromatograms.
#' @param ... further arguments (currently unused)
#'
#' @return `matrix` with `nrow` equal to `length(chromatograms)` and 2 columns
#'   (`gaussian_similarity`, `gaussian_residuals`). Returns `NA` for
#'   chromatograms where values cannot be calculated.
#'
#' @author  William Kumler, Philippine Louail
#'
#' @importFrom MetaboCoreUtils betaValues
#'
#' @examples
#' library(Chromatograms)
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L),
#'     mz = c(100.0, 200.0),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' pdata <- list(
#'     data.frame(rtime = seq(1, 20, by = 0.5),
#'                intensity = c(10, 15, 25, 50, 100, 200, 350, 500, 650,
#'                              700, 650, 500, 350, 200, 100, 50, 25, 15,
#'                              10, 8, 6, 5, 4, 3, 2, 2, 1, 1, 1, 1, 1,
#'                              1, 1, 1, 1, 1, 1, 1, 1)),
#'     data.frame(rtime = c(1, 2, 3, 4, 5, 6, 7),
#'                intensity = c(0, 10, 50, 100, 50, 10, 0))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' gaussianSimilarity(chr)
NULL

#' Internal helper: compute gaussian similarity for a single chromatogram
#' @noRd
.gaussianSimilaritySingle <- function(rts, ints, peakBoundary = NULL) {
    na_result <- c(gaussian_similarity = NA_real_,
                   gaussian_residuals = NA_real_)
    if (is.null(peakBoundary)) {
        peakBoundary <- .peakBoundarySingle(rts, ints)
    }
    if (any(is.na(peakBoundary))) return(na_result)
    left_rt <- peakBoundary["left_boundary"]
    right_rt <- peakBoundary["right_boundary"]
    mask <- rts >= left_rt & rts <= right_rt
    peak_rts <- rts[mask]
    peak_ints <- ints[mask]
    if (length(peak_ints) < 5) return(na_result)
    beta_vals <- betaValues(intensity = peak_ints, rtime = peak_rts)
    c(gaussian_similarity = unname(beta_vals[1]),
      gaussian_residuals = unname(beta_vals[2]))
}

#' @rdname gaussianSimilarity
#' @export
gaussianSimilarity <- function(chromatograms, peakBoundary = NULL, ...) {
    n <- length(chromatograms)
    rts_list <- rtime(chromatograms)
    ints_list <- intensity(chromatograms)
    res <- t(vapply(seq_len(n), function(i) {
        pb_i <- NULL
        if (!is.null(peakBoundary)) {
            pb_i <- if (is.matrix(peakBoundary)) peakBoundary[i, ] else peakBoundary
        }
        .gaussianSimilaritySingle(rts_list[[i]], ints_list[[i]], pb_i)
    }, numeric(2)))
    colnames(res) <- c("gaussian_similarity", "gaussian_residuals")
    attr(res, "gaussianSimilarity") <- "custom_metric:gaussian_similarity"
    res
}

#' @name peakProminence
#' @rdname peakProminence
#' @title Peak Prominence (Peak-to-Baseline Ratio) per Chromatogram
#'
#' @description
#' The function `peakProminence` calculates the prominence of the main
#' chromatographic peak relative to its baseline in each chromatogram. \cr
#'
#' The metric is calculated as follows (per chromatogram): \cr
#' (1) the retention time and intensity values are extracted, \cr
#' (2) if \code{peakBoundary} is provided, the data is subset to the peak
#' region, \cr
#' (3) the maximum intensity is determined, \cr
#' (4) the baseline intensity is estimated from the
#' \code{baselineQuantile} of intensities, \cr
#' (5) the prominence is calculated as
#' \code{(max - baseline) / baseline} and returned.
#'
#' @details
#' Higher values indicate more prominent peaks that stand out clearly from
#' the baseline.
#'
#' @param chromatograms `Chromatograms` object
#' @param peakBoundary optional peak boundary: either a `matrix` with columns
#'   `left_boundary` and `right_boundary` (one row per chromatogram), or a
#'   `numeric(2)` named vector applied to all chromatograms.
#' @param baselineQuantile `numeric(1)` quantile used to estimate the baseline
#'   intensity (default 0.1 = 10th percentile of intensities).
#' @param ... further arguments (currently unused)
#'
#' @return `numeric` of length equal to `length(chromatograms)`. Returns `NA`
#'   for chromatograms where prominence cannot be calculated.
#'
#' @author Philippine Louail
#'
#' @examples
#' library(Chromatograms)
#'
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L),
#'     mz = c(100.0, 200.0),
#'     dataOrigin = c("mem1", "mem1")
#' )
#' ## Good peak with high prominence and noisy plateau with low prominence
#' pdata <- list(
#'     data.frame(rtime = 1:20,
#'                intensity = c(100, 100, 100, 200, 500, 1000, 2000, 5000,
#'                              10000, 15000, 10000, 5000, 2000, 1000, 500,
#'                              200, 100, 100, 100, 100)),
#'     data.frame(rtime = 1:20,
#'                intensity = c(3000, 3500, 4000, 5000, 8000, 10000, 12000,
#'                              11000, 10000, 11000, 12000, 10000, 9000,
#'                              8000, 7000, 6000, 5000, 4000, 3500, 3000))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' peakProminence(chr)
NULL

#' Internal helper: compute peak prominence for a single chromatogram
#' @noRd
.peakProminenceSingle <- function(rts, ints, peakBoundary = NULL,
                                   baselineQuantile = 0.1) {
    if (!is.null(peakBoundary) && !any(is.na(peakBoundary))) {
        left_rt <- peakBoundary["left_boundary"]
        right_rt <- peakBoundary["right_boundary"]
        mask <- rts >= left_rt & rts <= right_rt
        ints <- ints[mask]
    }
    if (length(ints) < 3 || all(is.na(ints))) return(NA_real_)
    max_int <- max(ints, na.rm = TRUE)
    baseline_int <- quantile(ints, probs = baselineQuantile, na.rm = TRUE)
    if (is.na(baseline_int) || baseline_int <= 0) return(NA_real_)
    (max_int - baseline_int) / baseline_int
}

#' @rdname peakProminence
#' @export
peakProminence <- function(chromatograms, peakBoundary = NULL,
                           baselineQuantile = 0.1, ...) {
    n <- length(chromatograms)
    rts_list <- rtime(chromatograms)
    ints_list <- intensity(chromatograms)
    res <- vapply(seq_len(n), function(i) {
        pb_i <- NULL
        if (!is.null(peakBoundary)) {
            pb_i <- if (is.matrix(peakBoundary)) peakBoundary[i, ] else peakBoundary
        }
        .peakProminenceSingle(rts_list[[i]], ints_list[[i]], pb_i,
                               baselineQuantile)
    }, numeric(1))
    attr(res, "peakProminence") <- "custom_metric:peak_prominence"
    res
}
