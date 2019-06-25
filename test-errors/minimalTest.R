## Minimal example for testing
devtools::load_all(".")

coordinates <- GRanges("chr13:28444157-28550368")
snps.type <- ""
contacts.type <- ""#183446
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

files.contacts <- ""

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
