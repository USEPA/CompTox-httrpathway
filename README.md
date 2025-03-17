httrpathway R package
===========================

System Requirements
-------------------

### R

All code in this package has been tested on R 3.6.0 and 4.4.1, but any version of R 3.x or 4.x should be compatible. Other package versions may have impact the use of *httrpathway* due to changes in functionality, function parameter updates (e.g. GSVA), or changes in API.


Installing httrpathway
-------------------

httrpathway is currently implemented as an in-development R package and can be installed using the *devtools* R package (http://devtools.r-lib.org/):

**Local installation:**

1. Clone or download the package from github (TBD):
```bash
git clone (TBD)
```
2. Install the package using *devtools* functions in R/Rstudio:
```r
#load devtools library
library(devtools)

#use devtools functions to load package and documentation
devtools::install_local(path = "/path/to/httrpathway/", force = TRUE)

library(httrpathway)
```

**Using devtools::install_github():** *Note: This is still under development* 

1. Run devtools function to install directly from Github
```r
#load devtools library
library(devtools)

#install from Github
devtools::install_github(repo = "TBD")
```

Using httrpathway
-------------------

Please see the *httrpathway* vignette (`/httrpathway/vignettes/httrpathway-vignette.Rmd`) for instructions for how to execute the core functions of *httrpathway*. This includes an overview of the structure of the log2 fold-change tables *httrpathway* expects as input, as well as how to run the core concentration-response analysis of high-throughput transcriptomics (HTTr) data for either signature/pathways or genes/probes. 

### Version History

**v0.2 (2/6/2025)** -- in progress

+ Various code improvements
+ Improved roxygen2 documentation
+ devtools package load functions work as intended
+ Refined documentation including manual and vignette
+ Deprecated many non-core functions

**v0.1 (12/13/2021)**

+ Initial commit from public Github repo (12/13/2021)


Contributors
------------

+ **[Derik E. Haggard](mailto:haggard.derik@epa.gov)**
+ **[Joseph Bundy](mailto:bundy.joseph@epa.gov)**
+ **[Logan J. Everett](mailto:everett.logan@epa.gov)**
+ **Lionel Girardin**


