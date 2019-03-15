
### test error
path = "static_data/RData/"

## Load viartual 4c
load(paste0(path, "shared/baitID_and_name_virtual4C.rda"))
ids <-  unique(ids[order(ids$baitName), 1:2])
list.art4C <- ids$baitID
names(list.art4C) <- ids$baitName
rm(ids)

## Input parameters
input <- list(genome="hg19",
              gene="DEXI",
              ranges=5e4,
              maps.type="openChromatinClasses",
              clusters.type="enhancerClusters",
              tfs.type="adult",
              snps.type="",
              contacts.type=list.art4C[grep("DEXI", names(list.art4C))])

## Load gene name
load(paste0(path,
            input$genome, "/",input$genome, "_genesCoding_ensemblv75.rda"))
genes <- genes[genes$gene_name==input$gene,]
coord <- regioneR::extendRegions(genes, 
                                 extend.start=as.numeric(input$ranges),
                                 extend.end=as.numeric(input$ranges))

## Select contact files
files.contacts <- paste0(path,
                         input$genome, "/",
                         "new_Virtual4C/bw/",
                         "bg", c(0,3,5), "_score/",
                         substr(input$contacts.type, 1,1),
                         "/",
                         "PI_Merged_Digest_Human_HindIII_BaitID", input$contacts.type, "_bg", c(0,3,5), "_score.bw")


suppressWarnings(plotRegulome(coord,
             snps.type=gsub("-", "", input$snps.type),
             contacts.type=files.contacts,
             maps.type=gsub("-", "", input$maps.type),
             cluster.type=gsub("-", "", input$clusters.type),
             tfs.type=gsub("-", "", input$tfs.type),
             showLongestTranscript=TRUE,
             genome=input$genome,
             path=path))


####
contactsObject <- create_contactsRegulome(coordinates=coord,
                        contacts.type=files.contacts,
                        genome=input$genome,
                        path=path)
