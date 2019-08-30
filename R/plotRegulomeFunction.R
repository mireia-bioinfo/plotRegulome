#' Plot Regulome Data
#'
#' By providing the coordinates and IRB datasets of interest, it will produce the characteristic
#'  Islet Regulome Broser plot.
#' @param coordinates Either a \code{GRanges} object, data.frame or character string
#' (i.e. \code{chr11:17227951-17589050}) indicating the coordinates for the region the user wants
#' to plot.
#' @param snps_dataset GWAS SNPs dataset name to use for the analysis. The value should be one
#' of the defined below (column \emph{Dataset}). Defaults to \code{""}, which will produce an
#' empty plot.
#' \tabular{lll}{
#' \strong{Name}                     \tab \strong{Dataset}        \tab \strong{Reference}  \cr
#' 70KforT2D                         \tab 70KforT2D               \tab S. Bonàs-Guarch et al. (2018)   \cr
#' Diagram                           \tab diagram                 \tab A. P. Morris et al. (2012)         \cr
#' Magic                             \tab magic                   \tab R. A. Scott et al. (2012)          \cr
#' }
#' @param snps_col Color for the GWAS SNPs. Default: "dark red"
#' @param contacts_dataset BaitID or bait gene name for the virtual 4C data (Miguel-Escalada
#' et al. (2019)) to plot.
#' @param contacts_col Named vector ("0", "3" and "5"), with the colors for each group of CHiCAGO
#' contacts.
#' @param maps_dataset Name of the chromatin maps to plot. The value should be one
#' of the defined below (column \emph{Dataset}). Defaults to \code{""}, which will produce an
#' empty plot.
#' \tabular{lll}{
#' \strong{Name}                     \tab \strong{Dataset}        \tab \strong{Reference}  \cr
#' Adult Islets - Chromatin Classes  \tab chromatinClassesReduced \tab Miguel-Escalada et al. (2019)\cr
#' Adult Islets -  Chromatin Classes \tab openChromatinClasses    \tab Pasquali et al. (2014)       \cr
#' Pancreatic Progenitors            \tab progenitors             \tab Cebola et al. (2015)         \cr
#' Adult Islets - Chromatin States   \tab chromatinStates         \tab Parker et al. (2013)         \cr
#' }
#' @param clusters_dataset Name of the enhncer cluster data to plot. The value should be one
#' of the defined below (column \emph{Dataset}). Defaults to \code{""}, which will produce an
#' empty plot.
#' \tabular{lll}{
#' \strong{Name}                     \tab \strong{Dataset}        \tab \strong{Reference}  \cr
#' Enhancer Hubs                     \tab enhancerHubs            \tab Miguel-Escalada et al. (2019)\cr
#' Super Enhancers                   \tab superEnhancers          \tab Miguel-Escalada et al. (2019)\cr
#' Enhancer Clusters                 \tab enhancerClusters        \tab Pasquali et al. (2014)       \cr
#' Stretch Enhancers                 \tab stretchEnhancers        \tab Parker et al. (2013)         \cr
#' COREs                             \tab cores                   \tab K. J. Gaulton et al. (2010)        \cr
#' }
#' @param cluster_col Color for the clusters. Default: "dark green"
#' @param tfs_dataset Name of the TF dataset to plot. The value should be one
#' of the defined below (column \emph{Dataset}). Defaults to \code{""}, which will produce an
#' empty plot.
#' \tabular{lll}{
#' \strong{Name}                     \tab \strong{Dataset}        \tab \strong{Reference}  \cr
#' Adult Islets – Structural         \tab structure               \tab Miguel-Escalada et al. (2019)\cr
#' Pancreatic Progenitors            \tab progenitors             \tab Cebola et al. (2015)         \cr
#' Adult Islets – Tissue-specific    \tab adult                   \tab Pasquali et al. (2014)
#' }
#' @param tfs_col Color for the TFs circle border.
#' @param genes_col Named character vector ("gene", "lnc" and "spec") with the colors for each type of feture plotted in the gene
#' annotation track.
#' @param showLongestTranscript When plotting gene data, set to TRUE (default) if you want to reduce the
#' number of transcripts by only plotting the longest transcript per gene. If set to FALSE, will plot all
#' available transcripts.
#' @param randomIRB Generate random combinations of IRB datasets.
#' @param genome Character string indicating the genome for the coordinates. Default: hg19.
#' @param path Path containing the genomes folder (for example "hg19").
#' Default: "~/data/IRB/"
#' @return A ggplot2 object that can be plotted or saved with ggsave. Recomended device size in
#' inches is: \code{width=11, height=6}
#' @export
#' @import GenomicRanges
#' @examples \dontrun{
#'     coordinates <- "chr5:131817301-131826490"
#'     snps_dataset <- "70KforT2D"
#'     contacts_dataset <- 570519 #183446
#'     maps_dataset <- "progenitors"
#'     clusters_dataset <- ""
#'     tfs_dataset <- ""
#'     path <- "~/data/IRB/"
#'
#'     plotRegulome(coordinates=coordinates,
#'                  snps_dataset=snps_dataset,
#'                  contacts_dataset=contacts_dataset,
#'                  maps_dataset=maps_dataset,
#'                  clusters_dataset=clusters_dataset,
#'                  tfs_dataset=tfs_dataset,
#'                  path=path)
#' }

plotRegulome <- function(coordinates,
                         # SNPs -------
                         snps_dataset="",
                         snps_col="dark red",
                         # Contacts -------
                         contacts_dataset="",
                         contacts_col=c("0"="grey", "3"="blue", "5"="dark orange"),
                         # Maps -------
                         maps_dataset="",
                         # Clusters -------
                         clusters_dataset="",
                         cluster_col="dark green",
                         # TFs -------
                         tfs_dataset="",
                         tfs_col="dark blue",
                         # Genes -------
                         genes_col=c("gene"="dark grey",
                                     "spec"="darkorchid3",
                                     "lnc"="black"),
                         showLongestTranscript=TRUE,
                         # General -------
                         randomIRB=FALSE,
                         genome="hg19",
                         path="~/data/IRB/") {

  if (randomIRB) {
    rnd <- generateRandomIRB()
    coordinates <- rnd$coordinates
    snps_dataset <- rnd$snps
    maps_dataset <- rnd$maps
    clusters_dataset <- rnd$clusters
    tfs_dataset <- rnd$tfs
  }

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
                                            contacts_dataset=contacts_dataset,
                                            contacts_col=contacts_col,
                                            genome=genome,
                                            path=path)

  snpsObject <- create_snpsRegulome(coordinates=coordinates,
                                    snps_dataset=snps_dataset,
                                    snps_col=snps_col,
                                    snps_scaling=contactsObject$moreArgs$maxContact,
                                    genome=genome,
                                    path=path)

  mapsObject <- create_mapsRegulome(coordinates=coordinates,
                                    maps_dataset = maps_dataset,
                                    genome=genome,
                                    path=path)

  clustersObject <- create_clustersRegulome(coordinates=coordinates,
                                            clusters_dataset=clusters_dataset,
                                            cluster_col=cluster_col,
                                            genome=genome,
                                            path=path)

  tfsObject <- create_tfsRegulome(coordinates=coordinates,
                                  tfs_dataset=tfs_dataset,
                                  genome=genome,
                                  tfs_col=tfs_col,
                                  position_y=0.5,
                                  path=path)

  genesObject <- create_genesRegulome(coordinates=coordinates,
                                      genes_col = genes_col,
                                      showLongestTranscript=TRUE,
                                      genome="hg19",
                                      path=path)

  #### Full plot ---------------------------
  ## SNPs & contacts -----------------------
  plotSNPS <- plotR(snpsObject)
  plotContacts <- plotR(contactsObject)

  p1 <-
      ggplot2::ggplot() +
        plotSNPS[-length(plotSNPS)] +
        plotContacts[-length(plotContacts)] +
        scaleXCoordinates(chr=as.character(seqnames(snpsObject$coordinates)),
                                   limits=c(start(coordinates),
                                            end(coordinates)))


  if (length(snpsObject$value)>0 &
      length(contactsObject$value)>0) {
    scale <- max(-log10(snpsObject$value$PVAL), na.rm=T)/contactsObject$moreArgs$maxContact
    p1 <- p1 + ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~.*scale,
                                                                        name=expression("SNPs "*-log[10]*" P-value")))
  } else if (length(snpsObject$value)>0 &
             length(contactsObject$value)==0) {
    p1 <- p1 + ggplot2::scale_y_continuous(name=expression("SNPs "*-log[10]*" P-value"))
  }

  ## Maps & TFs -----------------------
  plotMaps <- plotR(mapsObject)
  plotClusters <- plotR(clustersObject)
  plotTFS <- plotR(tfsObject)

  p2 <-
    ggplot2::ggplot() +
      plotMaps[-length(plotMaps)] +
      ggplot2::geom_rect(ggplot2::aes(xmin=start(ranges(coordinates)),
                    xmax=end(ranges(coordinates)),
                    ymin=0.85, ymax=1), color="black", fill=NA) +
      plotClusters[-length(plotClusters)] +
      plotTFS[-length(plotTFS)] +
      ggplot2::ylim(0.25,1.1) +
      ggplot2::ylab(tfsObject$name)

  ## Genes -----------------------
  p3 <-
    ggplot2::ggplot() +
    plotR(genesObject)

  ## Compose main plot --------------
  main <- cowplot::plot_grid(p1 + themeXblank(legend.position="none"),
                             p2 + themeXYblank(title=T, legend.position="none"),
                             p3 + ggplot2::theme(legend.position="none"),
                             nrow=3, align="v", rel_heights = c(0.4,0.3,0.3))

  if (snpsObject$name=="" & contactsObject$name=="") {
    main$layers <- main$layers[-1]
  }

  legend <- generateLegendGG(contactsObject,
                             mapsObject,
                             clustersObject,
                             tfsObject,
                             snpsObject)

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

ggsave <- function(...) {
  ggplot2::ggsave(...)
}
