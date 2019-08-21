# plotRegulome

## Description
`plotRegulome` is an R package integrates different types of genomic and epigenomic data into a characteristic type of plot.

## How to install
Clone the repository with `git clone`. 
Move to the folder containing the repository, open `R` and use `devtools::install()` to install this package.

### Needed packages
- CRAN
    - [scales](https://CRAN.R-project.org/package=scales)
    - [ggrepel](https://cran.r-project.org/package=ggrepel)
    - [dplyr](https://cran.r-project.org/package=dplyr)
    - [cowplot](https://cran.r-project.org/package=cowplot)
    - [ggplot2](https://cran.r-project.org/package=ggplot2)

- Bioconductor
    - [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
    - [regioneR](https://bioconductor.org/packages/release/bioc/html/regioneR.html)

## Pending to implement  
2. Change proportions of plot depending on number of genes (width constant, height variable).

## Current Bugs
- ~~Package does not load necessary packages (???) GGplot, for example.~~
- ~~Generating spurious Rplots object when running Rscript ([more info](https://github.com/STAT545-UBC/Discussion/issues/59)).~~
- ~~Error with ideogram where it's completely grey.~~

## Possible improvements
- Add possibility of plotting several TFs data.
