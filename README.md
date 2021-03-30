# Quantification of thermal niche change across glacial cycles
This repository contains the R scripts to necessary to run the analyses and produce the final results of the article:

Gwen S. Antell, Isabel S. Fenton, Paul J. Valdes, and Erin E. Saupe (2021). "Thermal niches of planktonic foraminifera are static throughout glacial–interglacial climate change." _Proceedings of the National Academy of Sciences_.

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

The code is designed to analyze fossil occurrence and paleoclimate data across the last 700 ka. For file-path compatibility, save the data to a folder named /Data within the same top-level folder as the R files. The fossil occurrence data are a subset of the dataset published by Fenton et al. in _Scientific Data_, "Triton: a new database of Cenozoic planktonic foraminifera."

The first two scripts format the paleoclimate and occurrence data, respectively. The formatted output files are already saved with the data, so it is not necessary to source these scripts before running any of the analysis files. The file named "foram niche construction" does kernel density estimation on species' occupied temperature data. The relevant niche values are saved with the other data files, so this script can also be skipped. Each analysis file is independent and can be run in isolation. 
The "time series" analysis file conducts both the "trait evolution models" and "time series correlation" analyses. The "niches vs climate" analysis file conducts both the "niche lability vs climate change magnitude" and "niche lability among extreme climate intervals" analyses. The "exploratory plot all spp KDEs" file produces Figure 1 and Supplemental Figure 1. Several of the files require functions from the [kerneval](https://github.com/GwenAntell/kerneval) package, which can be installed with the command: devtools::install_github('GwenAntell/kerneval')
