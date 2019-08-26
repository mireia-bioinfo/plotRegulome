# plotRegulome

![](http://isletregulome.com/isletregulome/favicon.png)

## Description
`plotRegulome` is an R package that integrates different types of genomic and epigenomic data into the Islet Regulome Browser characteristic plot.

For a detailed description of all plot elements you can visit [isletregulome.com](isletregulome.com) and go to Info > Description of the plot.

## How to install
Clone the repository with `git clone`. 
Move to the folder containing the repository, open `R` and use `devtools` to install this package.

```
devtools::install(build_vignettes=TRUE)
```

### Required packages
- CRAN  
    - [scales](https://CRAN.R-project.org/package=scales)
    - [ggrepel](https://cran.r-project.org/package=ggrepel)
    - [dplyr](https://cran.r-project.org/package=dplyr)
    - [cowplot](https://cran.r-project.org/package=cowplot)
    - [ggplot2](https://cran.r-project.org/package=ggplot2)
    - [Hmisc](https://cran.r-project.org/package=Hmisc)

- Bioconductor  
    - [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
    - [regioneR](https://bioconductor.org/packages/release/bioc/html/regioneR.html)

## Quick start

```
library(plotRegulome)

## Download IRB database
dir <- "~/data/"
path <- downloadIRB(output_dir=dir) # Download IRB data, returns path for database

## Define parameters for plotting
coordinates <- "chr5:50678921-51318155"
snps_dataset <- "diagram"
contacts_dataset <- "ISL1"
maps_dataset <- "openChromatinClasses"
cluster_dataset <- "enhancerClusters"
tfs_dataset <- "adult"

## Call ploRegulome function
plotRegulome(coordinates=coordinates,
             snps_dataset=snps_dataset,
             contacts_dataset=contacts_dataset,
             maps_dataset=maps_dataset,
             cluster_dataset=cluster_dataset,
             tfs_dataset=tfs_dataset,
             path=path)
```

More information on datasets and plotting with `plotRegulome` can be found in the package vignette:

``` 
vignette("using_plotRegulome", package="plotRegulome")
``` 
