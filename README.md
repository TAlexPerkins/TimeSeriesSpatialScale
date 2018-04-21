# TimeSeriesSpatialScale
========================

This repository contains code related to the following preprint:

TA Perkins<sup>&#10033;</sup>, I Rodriguez-Barraquer<sup>&#10033;</sup>, C Manore<sup>&#10033;</sup>, AS Siraj, G España, CM Barker, MA Johansson, RC Reiner<sup>&#10033;</sup>. 2018. **Heterogeneous local dynamics revealed by classification analysis of spatially disaggregated time series data.** *bioRxiv* doi:[10.1101/276006](https://www.biorxiv.org/content/early/2018/03/05/276006)

All code contained within this repository is released under the [CRAPL v0.1 License](http://matt.might.net/articles/crapl/).

========================

### Overview

Contents of this repository are organized according to `code`, `data`, and `output`. Within each of those folders, there are three subfolders divided according to the three components of the analyses performed in this study, which are referred to consistently throughout the Abstract, Introduction, Methods, and Results sections of the preprint.

All code was written in the R language (version 3.3.1) and was executed on a Mac with OS X. At the time the research was done, all R packages used in this analysis were available on CRAN and were straightforward to install and load, with one exception. For those interested in running code used to generate cartograms, see [https://www.r-bloggers.com/cartograms-with-r/](https://www.r-bloggers.com/cartograms-with-r/) for guidance on how to set up this capability on your machine. We provide no further advice or troubleshooting on that topic.

In general, we anticipate that the structure of this repository will be relatively easy to navigate for someone who has first read the preprint. However, there are a few comments we make to clear up possible confusion about a few things.

* In `code/2_classification/`, the script ending in muni should be run before the one ending in dept, as there is a dependency in the latter on an output file from the former.
* In `code/3_elucidation/`, the script `simulate_ensemble_colombia.R` should be run before `analyze_curves_simulated_muni_reps.R`, as there is a dependency in the latter on an output file from the former.
* All data and output files required are available here, with one exception. If one were to run all the code, an output file called `simulations.RData` would be produced in `output/3_elucidation/`. We excluded this file from the repository due to its large size (459.2 MB). If you would like this file, you will need to either run the code or contact the authors.
