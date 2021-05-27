## MSnExp

library(msdata)
library(xcms)
swath_file <- system.file("TripleTOF-SWATH",
                          "PestMix1_SWATH.mzML",
                          package = "msdata")
swath_data <- readMSData(swath_file, mode = "onDisk")
cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10,
                     peakwidth = c(3, 30))
swath_data <- findChromPeaks(swath_data, param = cwp)
cwp <- CentWaveParam(snthresh = 3, noise = 10, ppm = 10,
                     peakwidth = c(3, 30))
swath_data <- findChromPeaksIsolationWindow(swath_data, param = cwp)
swath_spectra <- reconstructChromPeakSpectra(swath_data, minCor = 0.9,
                                             return.type = "Spectra")


dda_file <- system.file("TripleTOF-SWATH", "PestMix1_DDA.mzML",
                        package = "msdata")
dda_data <- readMSData(dda_file, mode = "onDisk")

dda_data %>% filterMsLevel(1L) %>% precursorMz()
dda_data %>% filterMsLevel(2L) %>% precursorIntensity()
estimatePrecursorIntensity(dda_data) ## for SciEx data

bpis <- chromatogram(dda_data)
plot(bpis, col = 1)
tic(dda_data)

## Spectra

library(Spectra)

spd <- DataFrame(
    msLevel = c(2L, 2L, 2L, 2L, 2L, 2L),
    polarity = c(1L, 1L, 1L, 1L, 1L, 1L),
    id = c("HMDB0000001", "HMDB0000001", "HMDB0001847", "HMDB0000517)", "HMDB0000177", "HMDB0000162"),
    name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine", "Arginine", "Histidine", "Proline"))

## Assign m/z and intensity values.
spd$mz <- list(
    c(109.2, 124.2, 124.5, 170.16, 170.52),
    c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
    c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
      111.0551, 123.0429, 138.0662, 195.0876),
    c(60.00, 70.00, 116),
    c(109.12, 156.07),
    c(69.91, 116.07))
spd$intensity <- list(
    c(3.407, 47.494, 3.094, 100.0, 13.240),
    c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
    c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994),
    c(200, 1000, 290),
    c(300, 1100),
    c(202, 471))
sps <- Spectra(spd)
sps$rtime <- c(1, 3, 7, 9, 10, 13)
sps
spectra <- sps