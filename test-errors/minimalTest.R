## Minimal example for testing
devtools::load_all()

coordinates <- #GRanges("chr5:50678921-51318155") #GRanges("chr13:28444157-28550368")
snps_dataset <- ""#"70KforT2D"
contacts_dataset <- 570519 #183446
maps_dataset <- "chromatinClassesReduced"
cluster_dataset <- "enhancerClusters"
tfs_dataset <- "adult"
path <- "~/data/IRB/"

files.contacts <- paste0(path,
                         "hg19/",
                         "new_Virtual4C/bw/",
                         "bg", c(0,3,5), "_score/",
                         substr(contacts_dataset, 1,1),
                         "/",
                         "PI_Merged_Digest_Human_HindIII_BaitID", contacts_dataset, "_bg", c(0,3,5), "_score.bw")

# files.contacts <- ""

devtools::load_all()
## General plot --------------------------
plotRegulome(coordinates=coordinates,
             snps_dataset=snps_dataset,
             contacts_dataset=contacts_dataset,#iles.contacts,
             maps_dataset=maps_dataset,
             cluster_dataset=cluster_dataset,
             tfs_dataset=tfs_dataset,
             path=path)

## Create different objects --------------
contactsObject <- create_contactsRegulome(coordinates=coordinates,
                                          contacts_dataset=570519,
                                          genome="hg19",
                                          path=path)

snpsObject <- create_snpsRegulome(coordinates=coordinates,
                                  snps_dataset=snps_dataset,
                                  genome="hg19",
                                  snp_scaling=contactsObject$moreArgs$maxContact,
                                  path=path)

mapsObject <- create_mapsRegulome(coordinates=coordinates,
                                  maps_dataset=maps_dataset,
                                  path=path)

clustersObject <- create_clustersRegulome(coordinates=coordinates,
                                          cluster_dataset=cluster_dataset,
                                          path=path)
tfsObject <- create_tfsRegulome(coordinates=coordinates,
                                tfs_dataset=tfs_dataset,
                                path=path)

genesObject <- create_genesRegulome(coordinates=coordinates)

ggplot2::ggplot() + plot(contactsObject)
ggplot2::ggplot() + plot(snpsObject)



####################---------------------
snps_dataset=""
snps_col="dark red"
# Contacts -------
contacts_dataset=""
contacts_col=c("0"="grey", "3"="blue", "5"="dark orange")
# Maps -------
maps_dataset=""
# Clusters -------
cluster_dataset=""
cluster_col="dark green"
# TFs -------
tfs_dataset=""
tfs_col="dark blue"
# Genes -------
genes_col=c("gene"="dark grey",
            "spec"="darkorchid3",
            "lnc"="black")
showLongestTranscript=TRUE
# General -------
genome="hg19"
path="~/data/IRB/"


coordinates <- "chr5:131817301-131826490"#GRanges("chr5:50678921-51318155") #GRanges("chr13:28444157-28550368")
snps_dataset <- "70KforT2D"
contacts_dataset <- 570519 #183446
maps_dataset <- "progenitors"
cluster_dataset <- ""
tfs_dataset <- ""
path <- "~/data/IRB/"


plotRegulome(coordinates=coordinates,
             snps_dataset=snps_dataset,
             contacts_dataset=contacts_dataset,#iles.contacts,
             maps_dataset=maps_dataset,
             cluster_dataset=cluster_dataset,
             tfs_dataset=tfs_dataset,
             path=path)
