# Repo overview: foraminiferal climatic niches project
Files in run order:
* BRIDGE GCM data prep
* foram occ data prep
* foram niche construction
* analysis time series
* analysis niches vs climate
* analysis type 2 error
* analysis phylo eco
* exploratory plot all spp KDEs

These rely on several functions:
* raster brick import fcn
* paleocoords fcn
* species kde buildr

This repository contains the R scripts to analyze fossil occurrence and paleoclimate data for planktonic foraminifera across the last 700 ka. The data (initial files as well as intermediate products from computationally-intensive analyses) are available in a separate storage space for large files; all data will be archived at Zenodo upon manuscript acceptance and are available to reviewers in the meantime. For file-path compatibility, save the data to a folder named /Data within the same top-level folder as the R files.

The first two scripts format the paleoclimate and occurrence data, respectively. The formatted output files are already saved with the data, so it is not necessary to source these scripts before running any of the analysis files. The file named "foram niche construction" does kernel density estimation on species' occupied temperature data. The relevant niche values are saved with the other data files, so this script can also be skipped. Each analysis file is independent and can be run in isolation. 
The "time series" analysis file conducts both the "trait evolution models" and "time series correlation" analyses. The "niches vs climate" analysis file conducts both the "niche lability vs climate change magnitude" and "niche lability among extreme climates" analyses. The "exploratory plot all spp KDEs" file produces Figure 1 and Supplemental Figure 1. Several of the files require functions from the kerneval package, which can be installed with the command: devtools::install_github('GwenAntell/kerneval')
