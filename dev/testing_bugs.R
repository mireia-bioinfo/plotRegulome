library(plotRegulome)

## Download IRB database
# dir <- "~/data/"
# path <- downloadIRB(output_dir=dir) # Download IRB data, returns path for database

## Define parameters for plotting
coordinates <- "chr13:28444157-28550368"
snps_dataset <- "diagram"
contacts_dataset <- "PDX1"
maps_dataset <- "openChromatinClasses"
clusters_dataset <- "enhancerClusters"
tfs_dataset <- "adult"
path <- "../PIRB_database/"

## Call ploRegulome function
plotRegulome(coordinates=coordinates,
             snps_dataset=snps_dataset,
             contacts_dataset=contacts_dataset,
             maps_dataset=maps_dataset,
             clusters_dataset=clusters_dataset,
             tfs_dataset=tfs_dataset,
             path=path)


## Genes
genes_col=c("gene"="dark grey",
            "spec"="darkorchid3",
            "lnc"="black")

genesObject <- create_genesRegulome(coordinates=coordinates,
                                    genes_col = genes_col,
                                    showLongestTranscript=TRUE,
                                    genome="hg19",
                                    path=path)
ggplot2::ggplot() + plotR(genesObject)
