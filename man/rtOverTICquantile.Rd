% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_Spectra_metrics.R
\name{rtOverTicQuantile}
\alias{rtOverTicQuantile}
\title{RT over TIC quantile (QC:4000054)}
\usage{
rtOverTicQuantile(
  spectra,
  probs = seq(0, 1, 0.25),
  msLevel = 1L,
  relative = TRUE,
  ...
)
}
\arguments{
\item{spectra}{`Spectra` object}

\item{probs}{`numeric` defining the quantiles. See [stats::quantile()] for
details. Defaults to `probs = seq(0, 1, 0.25)`.}

\item{msLevel}{`integer`}

\item{relative}{`logical`, if set to `TRUE` the relative retention time 
will be returned instead of the abolute retention time}

\item{...}{additional arguments passed to [stats::quantile()].}
}
\value{
`numeric` of length equal to length `probs` with the relative
   duration (duration divided by the total run time) after which the TIC
   exceeds the respective quantile of the TIC.
}
\description{
"The interval when the respective quantile of the TIC accumulates divided by 
retention time duration. The number of quantiles observed is given by the 
size of the tuple." [PSI:QC]
id: QC:4000054

The metric is calculated as follows:
(1) the `Spectra` object is ordered according to the retention time,

(2) the cumulative sum of the ion count is calculated (TIC),

(3) the quantiles are calculated according to the `probs` argument, e.g.
when `probs` is set to `c(0, 0.25, 0.5, 0.75, 1)` the 0%, 25%, 50%, 75% and
100% quantile is calculated,

(4) the retention time/relative retention time (retention time divided by 
the total run time taking into account the minimum retention time) is 
calculated,

(5) the (relative) duration of the LC run after which the cumulative
}
\details{
is_a: QC:4000004 ! n-tuple
is_a: QC:4000010 ! ID free
is_a: QC:4000021 ! retention time metric
is_a: QC:4000022 ! chromatogram metric
}
\examples{
library(S4Vectors)
library(Spectra)

spd <- DataFrame(
    msLevel = c(2L, 2L, 2L),
    polarity = c(1L, 1L, 1L),
    id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
    name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
## Assign m/z and intensity values
spd$mz <- list(
    c(109.2, 124.2, 124.5, 170.16, 170.52),
    c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
    c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
        111.0551, 123.0429, 138.0662, 195.0876))
spd$intensity <- list(
    c(3.407, 47.494, 3.094, 100.0, 13.240),
    c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
    c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
spd$rtime <- c(9.44, 9.44, 15.84)
sps <- Spectra(spd)
rtOverTicQuantile(spectra = sps, msLevel = 2L)
}
\author{
Thomas Naake, \email{thomasnaake@googlemail.com}, Johannes Rainer
}
