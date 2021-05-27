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
