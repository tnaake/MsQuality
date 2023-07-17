# MsQuality

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
![R-CMD-check-bioc](https://github.com/tnaake/MsQuality/workflows/R-CMD-check-bioc/badge.svg)
[![codecov.io](http://codecov.io/github/tnaake/MsQuality/coverage.svg?branch=master)](http://codecov.io/github/tnaake/MsQuality?branch=main)
[![license](http://img.shields.io/badge/license-GPL%20%28%3E=%203%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)

Calculation of QC metrics from mass spectrometry data

## Description
Data quality assessment is an integral part of preparatory data analysis 
to ensure sound biological information retrieval. 

We present here the `MsQuality` package, which provides functionality to calculate
quality metrics for mass spectrometry-derived, spectral data at the per-sample
level. `MsQuality` relies on the [`mzQC`](https://github.com/HUPO-PSI/mzQC) 
framework of quality metrics defined by the Human Proteome 
Organization-Proteomics Standards Intitiative (HUPO-PSI). These metrics 
quantify the quality of spectral raw files using a controlled vocabulary. 

The package is especially addressed towards users that acquire 
mass spectrometry data on a large scale (e.g. data sets from clinical settings 
consisting of several thousands of samples): while it is easier to control 
for high-quality data acquisition in small-scale experiments, typically run
in one or few batches, clinical data sets are often acquired over longer 
time frames and are prone to higher technical variation that is often
unnoticed. `MsQuality` tries to address this problem by calculating metrics that
can be stored along the spectral data sets (raw files or feature-extracted 
data sets). `MsQuality`, thus, facilitates the tracking of shifts in data quality
and quantifies the quality using multiple metrics. It should be thus easier
to identify samples that are of low quality (high-number of missing values,
termination of chromatographic runs, low instrument sensitivity, etc.).

The `MsQuality` package allows to calculate low-level quality metrics that require
minimum information on mass spectrometry data: retention time, m/z values, 
and associated intensities.
The list included in the `mzQC` framework is excessive, also including 
metrics that rely on more high-level information, that might not be readily
accessible from .raw or .mzML files, e.g. pump pressure mean, or rely
on alignment results, e.g. retention time mean shift, signal-to-noise ratio,
precursor errors (ppm). 

The `MsQuality` package is built upon the `Spectra` and the `MsExperiment` package.
Metrics will be calculated based on the information stored in a 
`Spectra` object, thus, the spectral data of each sample should be stored
in one `Spectra` object. The `MsExperiment` serves as a container to 
store the mass spectral data of multiple samples. `MsQuality` enables the user
to calculate quality metrics both on `Spectra` and `MsExperiment` objects. 


## Contact 

You are welcome to 

 * write a mail to <thomasnaake (at) googlemail (dot) com> 
 * submit suggestions and issues: <https://github.com/tnaake/MsQuality/issues>
 * send a pull request: <https://github.com/tnaake/MsQuality/issues> 

## Install

`MsQuality` is available via Bioconductor. To install the package, users can
either install from the 
[devel branch](https://bioconductor.org/packages/devel/bioc/html/MsQuality.html) 
or from the current 
[RELEASE branch](https://bioconductor.org/packages/release/bioc/html/MsQuality.html).

To install `MsQuality`, you have first to install the 
[`BiocManager`](https://www.bioconductor.org/install/) and
[remotes](http://cran.r-project.org/web/packages/cran/index.html) 
package: 

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
```

Install the `MsQuality` package then via

```r
## to install from Bioconductor
BiocManager::install("MsQuality")

## to install the development version from GitHub
BiocManager::install("tnaake/MsQuality")
```


