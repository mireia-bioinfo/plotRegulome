# plotRegulome 0.2.2

* Fixed bug in which gene exons were not filled.

# plotRegulome 0.2.1

* Added legend for RNA-seq data (islet-specific genes).

# plotRegulome 0.2.0

* `create_xxxRegulome()` functions now accept characters and dataframes as `coordintes` argument, not only `GRanges` objects.
* Added examples for `create_xxxRegulome()`.
* Renamed `cluster_dataset` to `clusters_dataset` for consistency reasons.
* Renamed generic S3 method from `plot()` to `plotR()`.
* Added legend for gene annotation. It is only included when plotting using the S3 method `plotR`, but not when unsing `plotRegulome` function.


# plotRegulome 0.1.0

* Changed whole database structure and file loading.
* Added `downloadIRB()` function for downloading required IRB datasets.
* Added `randomIRB` parameter in `plotRegulome()` to randomly generate Islet Regulome plots.
* Fixed great amount of warnings and messages when producing the plot.
* Fixed error `get_legend()` when no chromatin maps in view.
* Fixed bug with legend when `  tfs_dataset=""`.
* Added vignette `using_plotRegulome`.
* Added a `NEWS.md` file to track changes to the package.
