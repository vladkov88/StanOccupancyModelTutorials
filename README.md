# Occupancy model tutorials in Stan

This repository contains a tutorial for occupancy modeling using Stan and R.
The tutorials is written as R Markdown files.
The purpose of this tutorial is to help someone who is familiar with R and Stan to build and run occupancy models in with Stan in R. 
This tutorial is necessary because occupancy models in Stan can be difficult and require some programming tricks due to discrete latent variables. 

This document is a work in progress and will have additional materials added to it.
Additionally, I am open to having people contribute examples or improve my examples.

## File structure

- `README.md` is this file
- `OccupancyStanIntroduciton.Rmd` is the source file for the tutorial 
- `OccupancyStanIntroduciton.html` and `OccupancyStanIntroduciton.pdf` are complied versions of this tutorial
- `./Chapters/` are the markdown source files for each chapters and are called by `OccupancyStanIntroduction.Rmd` 
- `./ChaptersCode/` are the R and Stan source files for each chapters and are called by `OccupancyStanIntroduction.Rmd` 

## Contact

Richard A. Erickson (rerickson@usgs.gov)

## Disclaimer

This software is in the public domain because it contains materials that originally came from the U.S. Geological Survey, an agency of the United States Department of Interior. For more information, see the [official USGS copyright policy](https://www2.usgs.gov/visual-id/credit_usgs.html#copyright/).


This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.

This software is provided "AS IS".
