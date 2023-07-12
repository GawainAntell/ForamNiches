[![DOI](https://zenodo.org/badge/198804848.svg)](https://zenodo.org/badge/latestdoi/198804848)

# Temperature niche stability across glacial cycles
This repository contains the code to analyze fossil occurrence and paleoclimate data across the last 700 ka, as described in the article:

Gawain T. Antell, Isabel S. Fenton, Paul J. Valdes, and Erin E. Saupe (2021). "Thermal niches of planktonic foraminifera are static throughout glacial–interglacial climate change." _Proceedings of the National Academy of Sciences_, 118 (DOI:10.1073/pnas.2017105118).

Files in run order:
* BRIDGE GCM data prep
* foram occ data prep
* foram niche construction
* analysis time series
* analysis niches vs climate
* analysis type 2 error
* analysis phylo eco
* exploratory plot all spp KDEs
* movies sampling through time

These rely on several functions:
* raster brick import fcn
* paleocoords fcn
* species kde buildr

**Data files** are archived at Zenodo along with release v3.0 of this repo: [https://doi.org/10.5281/zenodo.4658884](https://doi.org/10.5281/zenodo.4658884). For file-path compatibility when running scripts, unzip the /Data folder within the same top-level location as the R files. The fossil occurrence data are a subset of the dataset published by Fenton and others in _Scientific Data_, "Triton, a new species-level database of Cenozoic planktonic foraminiferal occurrences."

The first two scripts format the paleoclimate and occurrence data, respectively. The formatted output files are already saved with the data, so it is not necessary to source these scripts before running any of the analysis files. The file named "foram niche construction" does kernel density estimation on species' occupied temperature data. The relevant niche values are saved with the other data files, so this script can also be skipped. Each analysis file is independent and can be run in isolation. The "time series" analysis file conducts both the "trait evolution models" and "time series correlation" analyses. The "niches vs climate" analysis file conducts both the "niche lability vs climate change magnitude" and "niche lability among extreme climate intervals" analyses. The "type 2 error" and "phylo eco" files each do their singular, namesake analysis. The "exploratory plot all spp KDEs" file produces Figure 1 and SI Appendix Figure S1. The "movies sampling through time" file produces a series of png images that were combined in ImageMagick to make SI Appendix Movie S1. 

Several of the files require functions from the [kerneval](https://github.com/GawainAntell/kerneval) package, which can be installed with the command: devtools::install_github('GawainAntell/kerneval')
