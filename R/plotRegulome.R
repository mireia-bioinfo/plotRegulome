#' Plot Regulome Data
#'
#' @param coordinates Either a GRanges object, data.frame or character string (chr11:17227951-17589050) indicating
#' the coordinates for the region the user wants to plot.
#' @param snps.type GWAS SNPs dataset name to use for the analysis. The value can be: diagram, magic, 70KforT2D or ""
#' (default, no SNPs will be plotted).
#' @param snps.col Color for the GWAS SNPs. Default: "dark red"
#' @param contacts.type ID for the bait you want to plot. See correspondance in \url{http://gattaca.imppc.org/genome_browser/lplab/IsletRegulome/Rdata/hg19/baitID_and_name_virtual4C.rda}
#' @param maps.type Name of the chromatin maps to plot. The value can be: chromatinClasses, chromatinClassesReduced, chromatinStates, openChromatinClasses, progenitors or "" (default, no map).
#' @param cluster.type Name of the cluster to plot. The value can be: enhancerClusters, enhancerHubs, stretchEnhancers, superEnhancers or "" (no clusters).
#' @param cluster.col Color for the clusters. Default: "dark green"
#' @param tfs.type Name of the TF dataset to plot. The value can be: adult, progenitors, structure or "" (don't plot TFs).
#' @param tfs.col Color for the tfs circles border.
#' @param showLongestTranscript When plotting gene data, set to TRUE (default) if you want to reduce the number of transcripts by only plotting the longest transcript per gene. If set to FALSE, will plot all the transcripts.
#' @param genome Character string indicating the genome for the coordinates. Default: hg19.
#' @param path Path containing the genomes folder (for example "hg19"). Default: "http://gattaca.imppc.org/genome_browser/lplab/IsletRegulome/Rdata/"
#' @return A ggplot2 object that can be plotted or saved with ggsave.
#' @export
#' @import GenomicRanges

plotRegulome <- function(coordinates,
                         snps.type="",
                         snps.col="dark red",
                         contacts.type="",
                         maps.type="",
                         cluster.type="",
                         cluster.col="dark green",
                         tfs.type="",
                         tfs.col="dark blue",
                         showLongestTranscript=TRUE,
                         genome="hg19",
                         path="http://gattaca.imppc.org/genome_browser/lplab/IsletRegulome/Rdata/") {

  if (class(coordinates)=="GRanges") {
    coordinates <- coordinates
  } else if (class(coordinates)=="data.frame") {
    coordinates <- GRanges(seqnames=coordinates[1,1],
                           ranges=IRanges(start=coordinates[1,2],
                                          end=coordinates[1,3]))
  } else if (class(coordinates)=="character") {
    coordinates <- GRanges(coordinates)
  }

  ###### Define Objects
  contactsObject <- create_contactsRegulome(coordinates=coordinates,
                                            contacts.type=contacts.type,
                                            genome=genome,
                                            path=path)

  snpsObject <- create_snpsRegulome(coordinates=coordinates,
                                    snps.type=snps.type,
                                    maxContacts=contactsObject$moreArgs$maxContact,
                                    genome=genome,
                                    path=path)

  mapsObject <- create_mapsRegulome(coordinates=coordinates,
                                    maps.type = maps.type,
                                    genome=genome,
                                    path=path)

  clustersObject <- create_clustersRegulome(coordinates=coordinates,
                                            cluster.type=cluster.type,
                                            genome=genome,
                                            path=path)

  tfsObject <- create_tfsRegulome(coordinates=coordinates,
                                  tfs.type=tfs.type,
                                  genome=genome,
                                  col="dark blue",
                                  position.y=0.5,
                                  path=path)

  genesObject <- create_genesRegulome(coordinates=coordinates,
                                      showLongestTranscript=TRUE,
                                      genome="hg19",
                                      path=path)

  #### Full plot
  scale <- snpsObject$moreArgs$maxLogPVAL/contactsObject$moreArgs$maxContact
  if(length(scale)==0) scale <- 1 else if(is.infinite(scale)) scale <- 1

  if(length(snpsObject$value)!=0) snpsName <- snpsObject$name else snpsName <- ""

  p1 <- ggplot2::ggplot() + plot(snpsObject) + plot(contactsObject) + themes$panel1 +
    ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~.*scale,
                                           name=snpsName)) +
    ggplot2::xlim(start(ranges(coordinates)), end(ranges(coordinates)))

  p2 <- ggplot2::ggplot() + plot(mapsObject) +
    ggplot2::geom_rect(ggplot2::aes(xmin=start(ranges(coordinates)),
                  xmax=end(ranges(coordinates)),
                  ymin=0.85, ymax=1), color="black", fill=NA) +
    plot(clustersObject) +
    plot(tfsObject) +
    themes$panel2 +
    ggplot2::ylim(0.25,1.1) +
    ggplot2::xlim(start(ranges(coordinates)), end(ranges(coordinates)))

  p3 <- ggplot2::ggplot() + plot(genesObject) + themes$panel3 +
    ggplot2::xlab("Genomic Position (bp)") +
    ggplot2::xlim(start(ranges(coordinates)), end(ranges(coordinates)))

  main <- cowplot::plot_grid(p1, p2, p3,
                    nrow=3, align="v", rel_heights = c(0.4,0.3,0.3))

  legend <- generateLegendGG(contactsObject, mapsObject,
                           clustersObject, tfsObject)

  RegulomePlot <- cowplot::plot_grid(
                                plotIdeogramRegulome(coordinates,
                                 genome),
                            main,
                            ncol=1,
                            rel_heights = c(0.10,0.9))
  RegulomePlot <- cowplot::plot_grid(RegulomePlot,
                                     legend,
                                     ncol=2,
                                     rel_widths = c(0.75,0.25))
  return(RegulomePlot)
}
