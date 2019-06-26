## Minimal example for testing
devtools::load_all()

coordinates <- GRanges("chr13:28444157-28550368")
snps.type <- "70KforT2D"
contacts.type <- 183446
maps.type <- "chromatinClassesReduced"
cluster.type <- "enhancerClusters"
tfs.type <- "adult"
path <- "~/tools/IRB/isletregulome_shiny/static_data/RData/"

files.contacts <- paste0(path,
                         "hg19/",
                         "new_Virtual4C/bw/",
                         "bg", c(0,3,5), "_score/",
                         substr(contacts.type, 1,1),
                         "/",
                         "PI_Merged_Digest_Human_HindIII_BaitID", contacts.type, "_bg", c(0,3,5), "_score.bw")

# files.contacts <- ""

devtools::load_all()
## General plot --------------------------
plotRegulome(coordinates=coordinates,
             snps.type=snps.type,
             contacts.type=files.contacts,
             maps.type=maps.type,
             cluster.type=cluster.type,
             tfs.type=tfs.type,
             path=path)

## Create different objects --------------
contactsObject <- create_contactsRegulome(coordinates=coordinates,
                                          contacts.type=files.contacts,
                                          genome="hg19",
                                          path=path)

snpsObject <- create_snpsRegulome(coordinates=coordinates,
                                          snps.type=snps.type,
                                          genome="hg19",
                                  scaling=contactsObject$moreArgs$maxContact,
                                          path=path)

mapsObject <- create_mapsRegulome(coordinates=coordinates,
                                  maps.type=maps.type,
                                  path=path)

clustersObject <- create_clustersRegulome(coordinates=coordinates,
                                          cluster.type=cluster.type,
                                          path=path)
tfsObject <- create_tfsRegulome(coordinates=coordinates,
                                tfs.type=tfs.type,
                                path=path)

ggplot2::ggplot() + plot(tfsObject)
