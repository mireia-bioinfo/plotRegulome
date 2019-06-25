#' Generic plotting methof for Regulome data
#'
#' This is a generic method that depending on the class of the passed object will plot according
#' to the regulome Plot.
#' @param regulomeObject An object created by \code{create_xxRegulome}, depending on the data you want to show.
#' @export
#' @import GenomicRanges

plot <- function(regulomeObject) {
  UseMethod("plot")
}

plot.snpsRegulome <- function(snpsObject) {
  if (length(snpsObject$value)!=0) {
    ## Add ids for top snps
    top <- snpsObject$moreArgs$top
    snpsObject$value <- snpsObject$value[order(-log10(snpsObject$value$PVAL),
                                               decreasing=TRUE),]
    snpsObject$value$id.show <- ""

    if(length(snpsObject$value)>top) {
      snpsObject$value$id.show[1:top] <- snpsObject$value$SNP[1:top]
    } else {
      snpsObject$value$id.show <- snpsObject$value$SNP
    }

    ## Convert SNPs to data.frame
    df <- data.frame(snpsObject$value)

    ## Rescale if coverage are present
    if (!is.null(snpsObject$moreArgs$maxContact)) {
      df$rescPVAL <- scales::rescale(-log10(df$PVAL),
                                     c(0, snpsObject$moreArgs$maxContact))
    } else {
      df$rescPVAL <- -log10(df$PVAL)
    }

    ## List with arguments to be used after ggplot()
    snpsPlot <- list(# SNP points
      ggplot2::geom_point(data=df,
                 ggplot2::aes(start, rescPVAL,
                     fill=-log10(PVAL),
                     color=-log10(PVAL),
                     alpha=-log10(PVAL)),
                 size=1, pch=21),
      ## SNP labels
      ggrepel::geom_text_repel(data=df,
                               ggplot2::aes(start, rescPVAL,
                                   label=id.show),
                               size=3),
      ## SNP scale
      ggplot2::scale_fill_continuous(high = snpsObject$col,
                            low = "white",
                            name="-log10\nPvalue",
                            limits=c(0, NA)),
      ggplot2::scale_color_continuous(high = snpsObject$col,
                             low = "white",
                             name="-log10\nPvalue",
                             limits=c(0, NA)),
      ## SNP legends
      ggplot2::guides(fill=ggplot2::guide_colorbar(title.position="bottom",
                                 barwidth=0.5,
                                 barheight=2),
             alpha=FALSE),
      ## SNP y label
      ggplot2::ylab(""),
      ## Limits for plot
      ggplot2::xlim(start(ranges(snpsObject$coordinates)),
           end(ranges(snpsObject$coordinates)))
    )
  } else {
    ## Return empty plot when GRanges is empty
    snpsPlot <- emptyPlot(snpsObject$coordinates, 0)
  }

  return(snpsPlot)
  # ggplot() + snpsPlot
}

plot.contactsRegulome <- function(contactsObject) {
  ## Check if GRanges contains elements
  if (length(contactsObject$value)!=0) {
    ## Add colors as variables
    len <- sapply(contactsObject$value, length)
    if (length(len)<3) {
      len <- len[c("0", "3", "5")]
      names(len) <- c("0", "3", "5")
      len[is.na(len)] <- 0
    }


    contacts.smooth <- lapply(contactsObject$value,
                              as.data.frame)
    contacts.smooth <- do.call(rbind, contacts.smooth)
    contacts.smooth$color <- unlist(mapply(rep, contactsObject$col, each=len))

    contactPlot <- list(## Plot rects (histogram-like) coverage
      ggplot2::geom_rect(data=contacts.smooth,
                ggplot2::aes(xmin=start, xmax=end,
                    ymin=0, ymax=meanScore),
                fill=contacts.smooth$color),
      ## Add triange for viewpoint (bottom)
      ggplot2::annotate("point", x=contactsObject$moreArgs$viewpoint,
               y=-0, pch=25, color="black", fill="black",
               size=3),
      ## Subtitle of axis y name
      # annotate("text", x=-Inf, y=0,
      #          label=paste0("Viewpoint: ",
      #                       contactsObject$name),
      #          angle=90, vjust=-3.5, size=6),
      ## Y axis label
      ggplot2::ylab("Virtual 4C"),
      ## Limits for plot
      ggplot2::xlim(start(ranges(contactsObject$coordinates)),
           end(ranges(contactsObject$coordinates))))
  } else {
    ## If GRanges is empty, return an empty plot
    contactPlot <- emptyPlot(contactsObject$coordinates, 0)
  }
  return(contactPlot)
  # ggplot() + contactPlot
}

plot.mapsRegulome <- function(mapsObject) {
  if (length(mapsObject$value)>0) {
    ## Convert RE type to character
    mapsObject$value$type <- as.character(mapsObject$value$type)
    maps.df <- data.frame(mapsObject$value)

    ## Add colors to data
    colors <- data.frame(type=names(mapsObject$col),
                         color=mapsObject$col,
                         stringsAsFactors=FALSE)
    maps.df <- dplyr::left_join(maps.df, colors)

    ## Plot maps
    mapsPlot <- list(## Plot chromatin maps
      ggplot2::geom_rect(data=maps.df,
                mapping=ggplot2::aes(xmin=start, xmax=end,
                            ymin=0.85, ymax=1),
                fill=maps.df$color),
      ## Limits for plot
      ggplot2::xlim(start(ranges(mapsObject$coordinates)),
           end(ranges(mapsObject$coordinates))))
  } else {
    mapsPlot <- emptyPlot(mapsObject$coordinates, 1)
  }

  return(mapsPlot)
  # ggplot() + mapsPlot
}

plot.clustersRegulome <- function(clustersObject) {
  if (length(clustersObject$value)>0) {
    ## Convert to data.frame
    clusters.df <- data.frame(clustersObject$value)

    ## Plot clusters
    clustersPlot <- list(ggplot2::geom_segment(data=clusters.df, ## Plot chromatin clusters
                                      ggplot2::aes(x=start, xend=end,
                                          y=1.05, yend=1.05, lwd=class),
                                      color=clustersObject$col),
                         ggplot2::scale_size_manual(values=clustersObject$moreArgs$lwd))
  } else {
    clustersPlot <- emptyPlot(clustersObject$coordinates, 1)
  }

  return(clustersPlot)
  # ggplot() + clustersPlot
}

plot.tfsRegulome <- function(tfsObject) {

  if(length(tfsObject$value)!=0) {
    ## Add number of TFs binding in same sites
    tfsObject$value$coloc <- countOverlaps(tfsObject$value, tfsObject$value)

    ## Select names of TFs present in data.frame
    tfNames <- unique(tfsObject$value$TF)

    ## Create data.frame with position of circles representing TFs
    points <- tile(tfsObject$coordinates, n=length(tfNames)+1)[[1]]
    points <- end(ranges(points))[-(length(tfNames)+1)]
    points.df <- data.frame("position"=points,
                            "height"=rep(tfsObject$moreArgs$position.y, length(points)),
                            "TF"=tfNames)

    ## Add value for circle sizes
    sizes <- seq(27, 29, by=0.5)
    points.df <- do.call(rbind,
                         replicate(length(sizes), points.df, simplify=FALSE))
    points.df$sizes <- rep(rev(sizes), each=length(tfNames))

    ## Add midpoint for binding sites
    tfs.df <- data.frame(tfsObject$value)
    tfs.df <- dplyr::left_join(tfs.df, points.df[,-4])
    tfs.df$midpoint <- tfs.df$start + tfs.df$width/2

    ## Plot TFs
    tfsPlot <- list(## Plot unions TF to chromatin
      ggplot2::geom_segment(data=tfs.df,
                                 ggplot2::aes(x=position, xend=midpoint,
                                     y=height, yend=0.84, color=coloc),
                                 lwd=0.3),
      ## Add scale for union lines
      ggplot2::scale_colour_gradient(low="light grey", high="black",
                            guide=FALSE),
      ## Plot lines for TF binding
      ggplot2::geom_segment(data=tfs.df,
                   ggplot2::aes(x=start, xend=end, y=0.84, yend=0.84),
                   color=tfsObject$col, lwd=2),
      ## Plot circles with TF
      ggplot2::geom_point(data=points.df,
                 ggplot2::aes(x=position, y=height),
                 pch=21, size=points.df$sizes, fill="grey",
                 color=tfsObject$col),
      ## Plot TF name in circles
      ggplot2::geom_text(data=unique(points.df[,-4]), ## Plot TF labels
                ggplot2::aes(x=position, y=height, label=TF),
                size=4),
      ## Label for Y axis name
      ggplot2::ylab(paste0(tfsObject$name,"\n")),
      ## Limits for plot
      ggplot2::xlim(start(ranges(tfsObject$coordinates)),
           end(ranges(tfsObject$coordinates))))
  } else {
    tfsPlot <- emptyPlot(tfsObject$coordinates, 1)
  }

  return(tfsPlot)
  # ggplot() + tfsPlot
}

plot.genesRegulome <- function(genesObject) {
  if (length(genesObject$value)!=0) {
    ## Add stepping to avoid overlapping genes
    genes.step <- addStepping(genesObject$value,
                              genesObject$coordinates)

    genes.df <- data.frame(genes.step)
    color <- data.frame(group=names(genesObject$col),
                        color=genesObject$col)
    genes.df <- dplyr::left_join(genes.df, color)

    distance <- width(genesObject$coordinates)*0.01
    genes.df$modEnd <- genes.df$end + distance

    genesPlot <- list(ggplot2::geom_segment(data=genes.df,
                                   ggplot2::aes(x=start, y=stepping,
                                       xend=end, yend=stepping)),
                      ggplot2::geom_rect(data=genes.df[genes.df$type=="EXON",],
                                ggplot2::aes(xmin=start, xmax=end,
                                    ymin=(stepping-0.3), ymax=(stepping+0.3)),
                                fill=genes.df$color[genes.df$type=="EXON"],
                                color="black"),
                      ggplot2::geom_text(data=genes.df[genes.df$type=="GENE",],
                                ggplot2::aes(x=modEnd, y=stepping,
                                    label=gene_name),
                                hjust=0, fontface=3,
                                size=3)
    )
  } else {
    genesPlot <- emptyPlot(genesObject$coordinates, 1)
  }
  return(genesPlot)
  # ggplot() + genesPlot
}

#' Smooth Coverage/Contacts
#'
#' Given a GRanges dataset that contains scores in the first mcol, it will smooth the data.
#' @param coverage.gr GRanges object with the first mcol containing the score of the contact.
#' @param smooth Arbitary number with smoothing value. Default: 5
#' @param coordinates GRanges object with coordinates you want to plot.
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
smoothCoverage <- function(coverage.gr,
                           smooth=5,
                           coordinates) {
  ## Create bins for smoothing data
  pixels <- 713
  factor <- pixels/smooth
  bin <- (max(end(ranges(coordinates))) - min(start(ranges(coordinates))))/factor

  ## GRangesList with binned regions
  regions.bin <- GenomicRanges::tile(coordinates, width=bin)[[1]]
  names(regions.bin) <- NULL

  hits <- findOverlaps(regions.bin, coverage.gr)
  hits.df <- data.frame(hits)

  score_list <- split(score(coverage.gr[subjectHits(hits)]), queryHits(hits))
  score_mean <- sapply(score_list, mean, na.rm=T)

  regions.score <- regions.bin[unique(queryHits(hits)),]
  regions.score$meanScore <- score_mean
  return(regions.score)
}

#' Add stepping for plotting genes
#'
#' Given a GRanges dataset representing genes, will add an arbitrary value for them to be plotted in
#' the Y axis without overlapping each other.
#' @param genesDat GRanges object containing gene information.
#' @param coordinates GRanges object with coordinates you want to plot.
#' @import GenomicRanges
addStepping <- function(genesDat, coordinates) {
  ## Create extension for avoiding overlap with gene names
  ext <- sapply(genesDat$gene_name, nchar) * width(coordinates)/50
  genesDat.ext <- regioneR::extendRegions(genesDat, extend.end=ext)

  ## Add stepping to data
  dict_stepping <- data.frame("tx_id"=unique(genesDat.ext$tx_id),
                              "stepping"=disjointBins(genesDat.ext[genesDat.ext$type=="GENE"],
                                                      ignore.strand=TRUE))
  mcols(genesDat) <- dplyr::left_join(data.frame(mcols(genesDat)),
                                      dict_stepping)
  return(genesDat)
}

#' Generate emtpy plot
#'
#' In those cases where there's no output, it generates an empty plot.
#' @param coordinates GRanges object with coordinates you want to plot.
#' @param y Value for the empty point in Y axis. Usually set to either 0 or 1.
#' @import GenomicRanges
emptyPlot <- function(coordinates, y) {
  df.empty <- data.frame(x=start(ranges(coordinates)),
                         y=y)
  emptyPlot <- list(ggplot2::geom_point(data=df.empty,
                              ggplot2::aes(x=x, y=y),
                              color=NA),
                    ggplot2::xlab(""), ggplot2::ylab(""))
  return(emptyPlot)
}
