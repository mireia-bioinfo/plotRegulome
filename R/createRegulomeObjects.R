#' Create snpsRegulome Object
#'
#' This function loads the GWAS information for SNPs from the specified dataset, selects those
#' present in \code{coordinates} and adds all the elements needed for plotting.
#' @param top Number indicating how many of the most significant SNPs in the region will be shown
#' as labels.
#' @param snps_scaling Number indicating the value to which scale the maximum -log10 PVAL for SNPs.
#' Defaults to NULL, meaning that no scaling will be performed.
#' @inheritParams plotRegulome
#' @return An object of class \code{snpsRegulome}, which is basically a list with the following
#'         attributes:
#'         \itemize{
#'             \item{\strong{name}: Name of the selected SNP dataset.}
#'             \item{\strong{col}: Color for plotting the SNPs.}
#'             \item{\strong{genome}: Genome build.}
#'             \item{\strong{coordinates}: Selected coordinates to plot as a \code{GRanges} object.}
#'             \item{\strong{value}: \code{GRanges} object containing all SNP values in \code{coordinates}.}
#'             \item{\strong{moreArgs}: Additional arguments, which are:
#'                   \itemize{
#'                       \item{\strong{top}: contains the
#'                             number of SNPs with lower P-values whose name will be plot}
#'                       \item{\strong{scaleTo}: parameter that includes a scaling parameter in case
#'                             virtual4C is also plot along SNPs (it is the maximum virtual4C value and
#'                             if no virtual4C is plot, is the maximum SNP -log10 p-value).}
#'                       }
#'                  }
#'             }
#' @export
#' @import GenomicRanges
create_snpsRegulome <- function(coordinates,
                                snps_dataset="",
                                snps_col="dark red",
                                genome="hg19",
                                top=3,
                                snps_scaling=NULL,
                                path="~/data/IRB/") {
  ## Convert coordinates to GRanges
  if (class(coordinates)=="GRanges") {
    coordinates <- coordinates
  } else if (class(coordinates)=="data.frame") {
    coordinates <- GRanges(seqnames=coordinates[1,1],
                           ranges=IRanges(start=coordinates[1,2],
                                          end=coordinates[1,3]))
  } else if (class(coordinates)=="character") {
    coordinates <- GRanges(coordinates)
  }

  ## Load SNPs data
  if (snps_dataset!="") {
    file <- paste0(path,
                   genome, "/snps/", genome, "_snps_", snps_dataset, "_",
                   as.character(seqnames(coordinates)), ".rda")
    if (file.exists(file)) {
      load(file)
    } else {
      stop("Can't find file in ", file,
           ". Please check if path is correctly defined and if your dataset is listed among: ",
           paste0(c(IRB$snps, '""'), collapse=", "),
           ".\nCheck ?plotRegulome for more informtion on available datasets.")
    }
  } else {
    snps <- GRanges() # create empty GRanges if no SNPs to plot
  }

  ## Create snpsRegulome object
  snps_sel <- subsetByOverlaps(snps, coordinates)

  ## Return paramters needed for scaling
  if (length(snps_sel)>0) max <- max(-log10(snps_sel$PVAL)) else max <- NULL
  if (is.null(snps_scaling)) scale_to <- max else scale_to <- snps_scaling # scale to max pval if no scaling

  snps <- list("name"=snps_dataset,
               "col"=snps_col,
               "genome"=genome,
               "coordinates"=coordinates,
               "value"=snps_sel,
               "moreArgs"=list("top"=top,
                               "scaleTo"=scale_to))
  snps <- structure(snps, class="snpsRegulome")

  return(snps)
}

#' Create contactsRegulome Object
#'
#' This function uses virtual 4C data (3 bigWig files, scores 0, 3 and 5) to generate an object containing
#' the data present in \code{coordinates} and all the elements needed for its plotting method.
#' @inheritParams plotRegulome
#' @return An object of class \code{contactsRegulome}, which is basically a list with the following
#'         attributes:
#'         \itemize{
#'             \item{\strong{name}: Name of bait for the selected virtual4C data.}
#'             \item{\strong{col}: Color for plotting virtual4C data, grouping by chicago scores (0=0-3; 3=3-5; 5=>5).}
#'             \item{\strong{genome}: Genome build.}
#'             \item{\strong{coordinates}: Selected coordinates to plot as a \code{GRanges} object.}
#'             \item{\strong{value}: \code{GRanges} object containing all CHiCAGO scores in \code{coordinates}.}
#'             \item{\strong{moreArgs}: Additional arguments, which are:
#'                   \itemize{
#'                       \item{\strong{viewpoint}: Position for the bait, will drow a triangle in this position.}
#'                       \item{\strong{maxContact}: Maximum CHiCAGO score in coordinates that will be used for
#'                             scaling SNP data in case it is also plot.}
#'                       }
#'                  }
#'             }
#' @export
#' @import GenomicRanges
create_contactsRegulome <- function(coordinates,
                                    contacts_dataset,
                                    contacts_col=c("0"="grey", "3"="blue", "5"="dark orange"),
                                    genome="hg19",
                                    path="~/data/IRB/") {
  ## Convert coordinates to GRanges
  if (class(coordinates)=="GRanges") {
    coordinates <- coordinates
  } else if (class(coordinates)=="data.frame") {
    coordinates <- GRanges(seqnames=coordinates[1,1],
                           ranges=IRanges(start=coordinates[1,2],
                                          end=coordinates[1,3]))
  } else if (class(coordinates)=="character") {
    coordinates <- GRanges(coordinates)
  }

  ## Load coverage data
  if (contacts_dataset != "") {
    load(paste0(path, genome, "/virtual4c/baitID_keyTable.rda"))

    if (contacts_dataset %in% ids$baitID) {
      sel <- ids[ids$baitID==contacts_dataset,]
    } else if (any(grepl(contacts_dataset, ids$baitName))) {
      sel <- ids[grep(contacts_dataset, ids$baitName),]

      sp <- strsplit(sel$baitName, ";")

      pos <- lapply(sp, grep, pattern=contacts_dataset)
      logi <- as.logical(sapply(pos, length))

      sel <- sel[which(logi),]
    }

    load(paste0(path, sel$file))

    virtual4c <- subsetByOverlaps(virtual4c, coordinates)
  } else {
    virtual4c <- GenomicRanges::GRanges()
    contacts_col <- NULL
    sel <- data.frame(baitName="",
                      position=NA)
  }

  ## Get maximum value
  if (length(virtual4c)>1) max <- max(virtual4c$score) else max <- NULL

  ## Create snpsRegulome object
  contacts <- list("name"=sel$baitName,
                   "col"=contacts_col,
                   "genome"=genome,
                   "coordinates"=coordinates,
                   "value"=virtual4c,
                   "moreArgs"=list("viewpoint"=sel$position,
                                   "maxContact"=max))
  contacts <- structure(contacts, class="contactsRegulome")

  return(contacts)
}

#' Create mapsRegulome Object
#'
#' This function loads chromatin maps, obtains the classes present in \code{coordinates} and adds
#' the needed parameters for plotting.
#' @inheritParams plotRegulome
#' @return An object of class \code{mapsRegulome}, which is basically a list with the following
#'         attributes:
#'         \itemize{
#'             \item{\strong{name}: Name of selected chromatin map dataset.}
#'             \item{\strong{col}: Named vector with the color of each type of chromatin class.}
#'             \item{\strong{genome}: Genome build.}
#'             \item{\strong{coordinates}: Selected coordinates to plot as a \code{GRanges} object.}
#'             \item{\strong{value}: \code{GRanges} object containing coordinates and type of chromatin class.}
#'             \item{\strong{moreArgs}: Additional arguments, which are empty.}
#'             }
#' @export
#' @import GenomicRanges
create_mapsRegulome <- function(coordinates,
                                maps_dataset,
                                genome="hg19",
                                path="~/data/IRB/") {
  ## Convert coordinates to GRanges
  if (class(coordinates)=="GRanges") {
    coordinates <- coordinates
  } else if (class(coordinates)=="data.frame") {
    coordinates <- GRanges(seqnames=coordinates[1,1],
                           ranges=IRanges(start=coordinates[1,2],
                                          end=coordinates[1,3]))
  } else if (class(coordinates)=="character") {
    coordinates <- GRanges(coordinates)
  }

  ## Load Maps data
  if(maps_dataset!="") {
    file <- paste0(path, genome, "/maps/",
                genome, "_map_", maps_dataset, ".rda")

    if (file.exists(file)) {
      load(file)
    } else {
      stop("Can't find file in ", file,
           ". Please check if path is correctly defined and if your dataset is listed among: ",
           paste0(c(IRB$maps, '""'), collapse=", "),
           ".\nCheck ?plotRegulome for more informtion on available datasets.")
    }
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
#' @inheritParams plotRegulome
#' @return An object of class \code{clustersRegulome}, which is basically a list with the following
#'         attributes:
#'         \itemize{
#'             \item{\strong{name}: Name of selected enhancer cluster dataset.}
#'             \item{\strong{col}: Character vector with the color of the enhancer clusters.}
#'             \item{\strong{genome}: Genome build.}
#'             \item{\strong{coordinates}: Selected coordinates to plot as a \code{GRanges} object.}
#'             \item{\strong{value}: \code{GRanges} object containing coordinates for the enhancer clusters and
#'                 the class of the line that will be drawn.}
#'             \item{\strong{moreArgs}: Additional arguments, which include:
#'                 \itemize{
#'                     \item{\strong{lwd}: Information on the linewidth for each class of enhancer cluster
#'                         feature.}
#'                     }
#'                  }
#'             }
#' @export
#' @import GenomicRanges
create_clustersRegulome <- function(coordinates,
                                    clusters_dataset,
                                    cluster_col="dark green",
                                    genome="hg19",
                                    path="~/data/IRB/") {
  ## Convert coordinates to GRanges
  if (class(coordinates)=="GRanges") {
    coordinates <- coordinates
  } else if (class(coordinates)=="data.frame") {
    coordinates <- GRanges(seqnames=coordinates[1,1],
                           ranges=IRanges(start=coordinates[1,2],
                                          end=coordinates[1,3]))
  } else if (class(coordinates)=="character") {
    coordinates <- GRanges(coordinates)
  }

  ## Load Cluster data
  if(clusters_dataset!="") {
    file <- paste0(path, genome, "/clusters/",
                   genome, "_cluster_", clusters_dataset, ".rda")

    if (file.exists(file)) {
      load(file)
    } else {
      stop("Can't find file in ", file,
           ". Please check if path is correctly defined and if your dataset is listed among: ",
           paste0(c(IRB$clusters, '""'), collapse=", "),
           ".\nCheck ?plotRegulome for more informtion on available datasets.")
    }
  } else {
    clusters.n=""
    clusters <- GRanges()
  }

  ## Select clusters in coordinates
  clust_sel <- subsetByOverlaps(clusters, coordinates)

  ## Adjust start and end to match start and end of coordinates
  start(clust_sel)[start(clust_sel) < start(coordinates)] <- start(coordinates)
  end(clust_sel)[end(clust_sel) > end(coordinates)] <- end(coordinates)

  ## Create clusterRegulome object
  clusters <- list("name"=clusters.n,
               "col"=cluster_col,
               "genome"=genome,
               "coordinates"=coordinates,
               "value"=clust_sel,
               "moreArgs"=list("lwd"=c("CDS"=2, # Info for enhancer hubs
                                       "intron"=1,
                                       "exon"=3)))
  clusters <- structure(clusters, class="clustersRegulome")

  return(clusters)
}

#' Create tfsRegulome Object
#'
#' This function load the TFBS data, selects those binding sites present in \code{coordinates} and
#' adds the needed parameters for plotting.
#' @param position_y Number from 0 (bottom) to 1 (top) indicating the position at which the circles should
#' be plotted. Default = 0.5.
#' @inheritParams plotRegulome
#' @return An object of class \code{tfsRegulome}, which is basically a list with the following
#'         attributes:
#'         \itemize{
#'             \item{\strong{name}: Name of selected transcription factor dataset.}
#'             \item{\strong{col}: Character vector with the color of the enhancer clusters.}
#'             \item{\strong{genome}: Genome build.}
#'             \item{\strong{coordinates}: Selected coordinates to plot as a \code{GRanges} object.}
#'             \item{\strong{value}: \code{GRanges} object containing coordinates for the transcription factor binding
#'                 sites and the name of the TF that is binding.}
#'             \item{\strong{moreArgs}: Additional arguments, which include:
#'                 \itemize{
#'                     \item{\strong{positionY}: Information on the linewidth for each class of enhancer cluster
#'                         feature.}
#'                     }
#'                  }
#'             }
#' @export
create_tfsRegulome <- function(coordinates,
                               tfs_dataset,
                               tfs_col="dark blue",
                               position_y=0.5, # from 0 to 1
                               genome="hg19",
                               path="~/data/IRB/") {
  ## Convert coordinates to GRanges
  if (class(coordinates)=="GRanges") {
    coordinates <- coordinates
  } else if (class(coordinates)=="data.frame") {
    coordinates <- GRanges(seqnames=coordinates[1,1],
                           ranges=IRanges(start=coordinates[1,2],
                                          end=coordinates[1,3]))
  } else if (class(coordinates)=="character") {
    coordinates <- GRanges(coordinates)
  }

  ## Load TFs data
  if (tfs_dataset!="") {
    file <- paste0(path, genome, "/tfs/",
                genome, "_tfs_", tfs_dataset, ".rda")

    if (file.exists(file)) {
      load(file)
    } else {
      stop("Can't find file in ", file,
           ". Please check if path is correctly defined and if your dataset is listed among: ",
           paste0(c(IRB$tfs, '""'), collapse=", "),
           ".\nCheck ?plotRegulome for more informtion on available datasets.")
    }
  } else {
    tfs.name <- ""
    tfs <- GRanges()
  }

  ## Create tfsRegulome object
  tfs <- list("name"=tfs.name,
               "col"=tfs_col,
               "genome"=genome,
               "coordinates"=coordinates,
               "value"=subsetByOverlaps(tfs, coordinates),
               "moreArgs"=list("positionY"=position_y))
  tfs <- structure(tfs, class="tfsRegulome")

  return(tfs)
}

#' Create genesRegulome Object
#'
#' This funciton loads all genes in chr, selects those transcripts present in \code{coordinates}
#' and adds the needed parameters for plotting.
#' @inheritParams plotRegulome
#' @return An object of class \code{mapsRegulome}, which is basically a list with the following
#'         attributes:
#'         \itemize{
#'             \item{\strong{name}: Gene Annotation.}
#'             \item{\strong{col}: Named vector with the color of each group of gene annotation, currently
#'                 "gene", "spec", and "lnc".}
#'             \item{\strong{genome}: Genome build.}
#'             \item{\strong{coordinates}: Selected coordinates to plot as a \code{GRanges} object.}
#'             \item{\strong{value}: \code{GRanges} object containing coordinates for the gene annotations, together
#'                 with the type of gene annotation (EXON or GENE) and the group..}
#'             \item{\strong{moreArgs}: Additional arguments, which are empty.}
#'             }
#' @export
#' @import GenomicRanges
create_genesRegulome <- function(coordinates,
                                 genes_col=c("gene"="dark grey",
                                             "spec"="darkorchid3",
                                             "lnc"="black"),
                                 showLongestTranscript=TRUE,
                                 genome="hg19",
                                 path="~/data/IRB/") {
  ## Convert coordinates to GRanges
  if (class(coordinates)=="GRanges") {
    coordinates <- coordinates
  } else if (class(coordinates)=="data.frame") {
    coordinates <- GRanges(seqnames=coordinates[1,1],
                           ranges=IRanges(start=coordinates[1,2],
                                          end=coordinates[1,3]))
  } else if (class(coordinates)=="character") {
    coordinates <- GRanges(coordinates)
  }

  ## Load gene and non-coding data
  load(paste0(path, genome, "/genes/",
              genome, "_gene_annotation_ensemblv75_",
              as.character(seqnames(coordinates)), ".rda"))
  genes <- subsetByOverlaps(genes, coordinates, ignore.strand=TRUE)

  load(paste0(path, genome, "/genes/",
              genome, "_lncRNA.rda"))
  lnc <- subsetByOverlaps(lncRNA, coordinates)
  lnc$tx_id <- lnc$id
  mcols(lnc) <- mcols(lnc)[,c(2,3,4,1)]

  ## Select genes of interest
  if (length(genes)>0) genes$group <- "gene"

  if (length(lnc)>0) {
    lnc$gene_biotype <- "lncRNA"
    lnc$longest <- T
    lnc$group <- "lnc"
    }
  genes <- c(genes, lnc)

  load(paste0(path, "shared/specific_genes.rda"))
  genes$group[genes$gene_name %in% spec_genes] <- "spec"

  ## If showLongestTranscript is TRUE, select longest transcripts for each gene
  if (showLongestTranscript) genes <- genes[genes$longest,]

  ## Adjust start and end to match start and end of coordinates
  start(genes)[start(genes) < start(coordinates)] <- start(coordinates)
  end(genes)[end(genes) > end(coordinates)] <- end(coordinates)

  ## Create genesRegulome object
  genes <- list("name"="Gene Annotation",
                "col"=genes_col,
                "genome"=genome,
                "coordinates"=coordinates,
                "value"=genes,
                "moreArgs"=list())
  genes <- structure(genes, class="genesRegulome")

  return(genes)
}
