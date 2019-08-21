#' Plot Regulome Data
#'
#' @param coordinates Either a GRanges object, data.frame or character string (chr11:17227951-17589050) indicating
#' the coordinates for the region the user wants to plot.
#' @param snps_dataset GWAS SNPs dataset name to use for the analysis. The value can be: diagram, magic,
#' 70KforT2D or "".
#' (default, no SNPs will be plotted).
#' @param snps_col Color for the GWAS SNPs. Default: "dark red"
#' @param contacts_dataset BaitID or bait gene name for the virtual 4C data you want to plot.
#' @param maps_dataset Name of the chromatin maps to plot. The value can be: chromatinClasses,
#' chromatinClassesReduced, chromatinStates, openChromatinClasses, progenitors or "" (default, no map).
#' @param cluster_dataset Name of the cluster to plot. The value can be: enhancerClusters, enhancerHubs,
#' stretchEnhancers, superEnhancers or "" (no clusters).
#' @param cluster_col Color for the clusters. Default: "dark green"
#' @param tfs_dataset Name of the TF dataset to plot. The value can be: adult, progenitors, structure or ""
#' (don't plot TFs).
#' @param tfs_col Color for the tfs circles border.
#' @param genes_col Named character vector with the colors for each type of feture plotted in the gene
#' annotation track.
#' @param showLongestTranscript When plotting gene data, set to TRUE (default) if you want to reduce the
#' number of transcripts by only plotting the longest transcript per gene. If set to FALSE, will plot all
#' the transcripts.
#' @param genome Character string indicating the genome for the coordinates. Default: hg19.
#' @param path Path containing the genomes folder (for example "hg19").
#' Default: "http://gattaca.imppc.org/genome_browser/lplab/IsletRegulome/Rdata/"
#' @return A ggplot2 object that can be plotted or saved with ggsave.
#' @export
#' @import GenomicRanges

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
                         cluster_dataset="",
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
                         genome="hg19",
                         path="~/data/IRB/") {

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
                                            cluster_dataset=cluster_dataset,
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
                                      showLongestTranscript=TRUE,
                                      genome="hg19",
                                      path=path)

  #### Full plot ---------------------------
  ## SNPs & contacts -----------------------
  p1 <-
      ggplot2::ggplot() +
        plot(snpsObject) +
        plot(contactsObject)


  if (length(snpsObject$value)>0 &
      length(contactsObject$value)>0) {
    scale <- max(-log10(snpsObject$value$PVAL), na.rm=T)/contactsObject$moreArgs$maxContact
    p1 <- p1 + ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~.*scale,
                                                                        name=expression("SNPs "*-log[10]*" P-value")))
  }

  ## Maps & TFs -----------------------
  p2 <-
    ggplot2::ggplot() +
      plot(mapsObject) +
      ggplot2::geom_rect(ggplot2::aes(xmin=start(ranges(coordinates)),
                    xmax=end(ranges(coordinates)),
                    ymin=0.85, ymax=1), color="black", fill=NA) +
      plot(clustersObject) +
      plot(tfsObject) +
      ggplot2::ylim(0.25,1.1)

  ## Genes -----------------------
  p3 <-
    ggplot2::ggplot() +
    plot(genesObject)

  ## Compose main plot --------------
  main <- cowplot::plot_grid(p1 + themeXblank(legend.position="none"),
                             p2 + themeXYblank(title=T, legend.position="none"),
                             p3,
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
