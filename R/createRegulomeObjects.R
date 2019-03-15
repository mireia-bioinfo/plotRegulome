#' Create snpsRegulome Object
#'
#' This function loads the GWAS information for SNPs from the specified dataset, selects those
#' present in \code{coordinates} and adds all the elements needed for plotting.
#' @param coordinates GRanges object with coordinates you want to plot.
#' @param snps.type GWAS SNPs dataset name to use for the analysis. The value can be: diagram, magic, 70KforT2D or ""
#' (default, no SNPs will be plotted).
#' @param col Color for the GWAS SNPs. Default: "dark red"
#' @param genome Character string indicating the genome for the coordinates. Default: hg19.
#' @param top Number of most significant SNPs for showing ID. Default: 3.
#' @param maxContacts Maximum score for contacts in contactsRegulome object. If no contacts present, should be NULL.
#' @export
#' @import GenomicRanges
create_snpsRegulome <- function(coordinates,
                                snps.type,
                                col="dark red",
                                genome,
                                top=3,
                                maxContacts=NULL,
                                path="~/data/IsletRegulome/Rdata/") {
  ## Load SNPs data
  if (snps.type!="") {
    load(paste0(path,
                    genome, "/", genome, "_snps_", snps.type, "_",
                    as.character(seqnames(coordinates)), ".rda"))
  } else {
    snps <- GRanges()
  }

  ## Create snpsRegulome object
  snps.sel <- subsetByOverlaps(snps, coordinates)
  if (length(snps.sel)!=0) max <- max(-log10(snps.sel$PVAL)) else max <- NULL

  snps <- list("name"=snps.type,
               "col"=col,
               "genome"=genome,
               "coordinates"=coordinates,
               "value"=snps.sel,
               "moreArgs"=list("top"=top,
                               "maxContact"=maxContacts,
                               "maxLogPVAL"=max))
  snps <- structure(snps, class="snpsRegulome")

  return(snps)
}

#' Create contactsRegulome Object
#'
#' This function uses virtual 4C data (3 bigWig files, scores 0, 3 and 5) to generate an object containing
#' the data present in \code{coordinates} and all the elements needed for its plotting method.
#' @param coordinates GRanges object with coordinates you want to plot.
#' @param contacts.type Character vector including the coordinates for the bw files.
#' @param col Colors for plotting virtual 4C data. Defaults to c("grey", "blue", "dark orange").
#' @param genome Character string indicating the genome for the coordinates. Default: hg19.
#' @param smooth Arbitrary number for smoothing contact signal. Default: 5.
#' @export
#' @import GenomicRanges
create_contactsRegulome <- function(coordinates,
                                    contacts.type,
                                    col=c("0"="grey", "3"="blue", "5"="dark orange"),
                                    genome,
                                    smooth=5,
                                    path="~/data/IsletRegulome/Rdata/") {
  ## Load coverage data
  if (contacts.type != "") {
    cov <- lapply(contacts.type,
                  rtracklayer::import,
                  which=coordinates)
    cov <- GRangesList(cov)
    names(cov) <- c("0", "3", "5")
  } else {
    cov <- GenomicRanges::GRanges()
    col <- NULL
  }

  ## Smooth coverage
  idx <- sapply(cov, function(x) length(x)>0)
  if (length(idx)==0) idx=0
  cov.smooth <- lapply(cov[idx],
                       smoothCoverage,
                       smooth=smooth,
                       coordinates=coordinates)
  cov.smooth <- GRangesList(cov.smooth)

  max <- max(unlist(lapply(cov.smooth, function(x) max(mcols(x)[,1]))))
  if(is.infinite(max) | max==0) max <- NULL

  ## Control for ranges with score 0
  if (sum(cov.smooth$`0`$meanScore)==0) cov.smooth <- GenomicRanges::GRanges()

  ## Load info of coverage
  load(paste0(path,
              "shared/baitID_and_name_virtual4C.rda"))

  ## Obtain ID from path
  id.type <- unlist(strsplit(contacts.type[1], "/"))
  id.type <- gsub("PI_Merged_Digest_Human_HindIII_BaitID", "", id.type[length(id.type)])
  id.type <- as.numeric(gsub("_bg0_score.bw", "", id.type))

  ## Create snpsRegulome object
  contacts <- list("name"=ids$baitName[ids$baitID==id.type],
                   "col"=col,
                   "genome"=genome,
                   "coordinates"=coordinates,
                   "value"=cov.smooth,
                   "moreArgs"=list("smooth"=smooth,
                                   "viewpoint"=ids$position[ids$baitID==id.type],
                                   "maxContact"=max))
  if (length(contacts$name)==0) contacts$name <- ""
  contacts <- structure(contacts, class="contactsRegulome")

  return(contacts)
}

#' Create mapsRegulome Object
#'
#' This function loads chromatin maps, obtains the classes present in \code{coordinates} and adds
#' the needed parameters for plotting.
#' @param coordinates GRanges object with coordinates you want to plot.
#' @param maps.type Name of the chromatin maps to plot. The value can be: chromatinClasses, chromatinClassesReduced, chromatinStates, openChromatinClasses, progenitors or "" (default, no map).
#' @param genome Character string indicating the genome for the coordinates. Default: hg19.
#' @export
#' @import GenomicRanges
create_mapsRegulome <- function(coordinates,
                                maps.type,
                                genome,
                                path="~/data/IsletRegulome/Rdata/") {
  ## Load Maps data
  if(maps.type!="") {
    load(paste0(path,
                    genome, "/", genome, "_map_", maps.type, ".rda"))
  } else {
    map.name=""
    map_col=NULL
    map=GenomicRanges::GRanges()
  }

  ## Create mapsRegulome object
  maps <- list("name"=map.name,
               "col"=map_col,
               "genome"=genome,
               "coordinates"=coordinates,
               "value"=subsetByOverlaps(map, coordinates),
               "moreArgs"=list())
  maps <- structure(maps, class="mapsRegulome")

  return(maps)
}

#' Create clustersRegulome Object
#'
#' This function loads the selected cluster classification, selects the clusters present in \code{coordinates}
#' and adds the needed parameters for plotting.
#' @param coordinates GRanges object with coordinates you want to plot.
#' @param cluster.type Name of the cluster to plot. The value can be: enhancerClusters, enhancerHubs, stretchEnhancers, superEnhancers or "" (no clusters).
#' @param col Color for plotting clusters. Default: "dark green".
#' @param genome Character string indicating the genome for the coordinates. Default: hg19.
#' @export
#' @import GenomicRanges
create_clustersRegulome <- function(coordinates,
                                    cluster.type,
                                    col="dark green",
                                    genome,
                                    path="~/data/IsletRegulome/Rdata/") {
  ## Load Cluster data
  if(cluster.type!="") {
    load(paste0(path,
                    genome, "/", genome, "_cluster_", cluster.type, ".rda"))
  } else {
    clusters.n=""
    clusters <- GRanges()
  }

  ## Select clusters in coordinates
  clusters.sel <- subsetByOverlaps(clusters, coordinates)

  ## Adjust start and end to match start and end of coordinates
  start(ranges(clusters.sel))[start(ranges(clusters.sel)) < start(ranges(coordinates))] <-
    start(ranges(coordinates))
  end(ranges(clusters.sel))[end(ranges(clusters.sel)) > end(ranges(coordinates))] <-
    end(ranges(coordinates))

  ## Create clusterRegulome object
  clusters <- list("name"=clusters.n,
               "col"=col,
               "genome"=genome,
               "coordinates"=coordinates,
               "value"=clusters.sel,
               "moreArgs"=list("lwd"=c("CDS"=2,
                                       "intron"=1,
                                       "exon"=3)))
  clusters <- structure(clusters, class="clustersRegulome")

  return(clusters)
}

#' Create tfsRegulome Object
#'
#' This function load the TFBS data, selects those binding sites present in \code{coordinates} and
#' adds the needed parameters for plotting.
#' @param coordinates GRanges object with coordinates you want to plot.
#' @param tfs.type Name of the TF dataset to plot. The value can be: adult, progenitors, structure or "" (don't plot TFs).
#' @param col Color for the tfs circles border.
#' @param position Position of the center of TFs circle, should be between c(0,1).
#' @param genome Character string indicating the genome for the coordinates. Default: hg19.
#' @export
create_tfsRegulome <- function(coordinates,
                               tfs.type,
                               col="dark blue",
                               position.y=0.5, # from 0 to 1
                               genome,
                               path="~/data/IsletRegulome/Rdata/") {
  ## Load TFs data
  if (tfs.type!="") {
    load(paste0(path,
                    genome, "/", genome, "_tfs_", tfs.type, ".rda"))
  } else {
    tfs.name <- ""
    tfs <- GRanges()
  }

  ## Create tfsRegulome object
  tfs <- list("name"=tfs.name,
               "col"=col,
               "genome"=genome,
               "coordinates"=coordinates,
               "value"=subsetByOverlaps(tfs, coordinates),
               "moreArgs"=list("position.y"=position.y))
  tfs <- structure(tfs, class="tfsRegulome")

  return(tfs)
}

#' Create genesRegulome Object
#'
#' This funciton loads all genes in chr, selects those transcripts present in \code{coordinates}
#' and adds the needed parameters for plotting.
#' @param coordinates GRanges object with coordinates you want to plot.
#' @param col Color for the genes.
#' @param showLongestTranscript When plotting gene data, set to TRUE (default) if you want to reduce the number of transcripts by only plotting the longest transcript per gene. If set to FALSE, will plot all the transcripts.
#' @param genome Character string indicating the genome for the coordinates. Default: hg19.
#' @export
#' @import GenomicRanges
create_genesRegulome <- function(coordinates,
                                 col=c("gene"="dark grey",
                                       "spec"="darkorchid3",
                                       "lnc"="black"),
                                 showLongestTranscript=TRUE,
                                 genome,
                                 path="~/tools/isletregulomebrowser_shiny/static_data/RData/",
                                 specPath=paste0(path, "shared/specific_genes.rda")) {
  ## Load gene and non-coding data
  load(paste0(path,
              genome, "/", genome, "_gene_annotation_allNotReduced_ensemblv75_",
              as.character(seqnames(coordinates)), ".rda"))
  # load("../isletregulomebrowser_shiny/data/hg19/hg19_gene_annotation_codingNotReduced_ensemblv75_chr20.rda")
  load(paste0(path,
                  genome, "/", genome, "_lncRNA.rda"))
  lncRNA$tx_id <- lncRNA$id
  mcols(lncRNA) <- mcols(lncRNA)[,c(2,3,4,1)]

  ## Select genes of interest
  genes <- subsetByOverlaps(genes, coordinates, ignore.strand=TRUE)
  if (length(genes)>0) genes$group <- "gene"
  lnc <- subsetByOverlaps(lncRNA, coordinates)
  if (length(lnc)>0) lnc$group <- "lnc"
  genes <- c(genes, lnc)

  load(specPath)
  genes$group[genes$gene_name %in% spec_genes] <- "spec"

  ## If showLongestTranscript is TRUE, select longest transcripts for each gene
  genes.split <- split(genes, genes$gene_name) ## split according to gene name
  longest <- unique(unlist(lapply(genes.split,
                                  function(x) x[which(width(x)==max(width(x))),]$tx_id)))
  genes <- genes[genes$tx_id %in% longest,]

  ## Adjust start and end to match start and end of coordinates
  start(ranges(genes))[start(ranges(genes)) < start(ranges(coordinates))] <-
    start(ranges(coordinates))
  end(ranges(genes))[end(ranges(genes)) > end(ranges(coordinates))] <-
    end(ranges(coordinates))

  ## Create genesRegulome object
  genes <- list("name"="geneAnnotation",
                "col"=col,
                "genome"=genome,
                "coordinates"=coordinates,
                "value"=genes,
                "moreArgs"=list())
  genes <- structure(genes, class="genesRegulome")

  return(genes)
}
