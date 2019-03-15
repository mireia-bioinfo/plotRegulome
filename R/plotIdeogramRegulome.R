#' Plot Ideogram for Regulome Plot
#'
#' This function creates an ideogram representation of the chromosome in \code{coordinates}.
#' @param coordinates GRanges object with coordinates you want to plot.
#' @param genome Character string indicating the genome for the coordinates. Default: hg19.
#' @import GenomicRanges

plotIdeogramRegulome <- function(coordinates,
                                 genome) {
  # load("data/ideoCyto.rda")
  hg19 <- ideoCyto$hg19

  dat <- hg19[as.character(seqnames(hg19)) %in% as.character(seqnames(coordinates))]

  dat <- data.frame(dat)

  max <- max(dat$end)
  space <- 0.3
  ymin=0
  ymax.ideo=10
  ymax=20

  ## scale
  sc <- width(coordinates)/10
  units <- c("bp"=1,
             "Kb"=1e3,
             "Mb"=1e6)
  cont=0
  while(sc>1e3) {
    cont <- cont + 1
    sc <- sc/units[cont]
    names(sc) <- names(units)[cont]
  }

  lab.scale <- paste0(round(sc, 0), names(sc))

  ideo <-
  ggplot2::ggplot() +
    ggplot2::geom_rect(data=dat,
              ggplot2::aes(xmin=start, xmax=end,
                  ymin=0, ymax=ymax.ideo,
                  fill=gieStain)) +
    ggplot2::geom_rect(ggplot2::aes(xmin=1,
                  xmax=max(dat$end),
                  ymin=0, ymax=ymax.ideo),
              fill=NA, color="black") +
    ggplot2::scale_fill_manual(values=cytobandColor) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position="none") +
    ggplot2::annotate("text",
             max+max*0.1,
             ymax.ideo/2, label=genome,
             size=4) +
    ggplot2::annotate("text",
             max/2,
             ((ymax-ymax.ideo)/1.9+ymax.ideo),
             label=as.character(seqnames(coordinates)),
             size=4) +
    ggplot2::geom_rect(data=data.frame(coordinates),
              ggplot2::aes(xmin=start, xmax=end, ymin=0-1, ymax=ymax.ideo+1),
              color="dark green", fill="dark green", alpha=0.8) +
    ggplot2::geom_segment(data=data.frame(coordinates),
                 ggplot2::aes(x=-0.2*max, xend=-0.06*max,
                     y=ymax.ideo/2, yend=ymax.ideo/2)) +
    ggplot2::geom_segment(data=data.frame(coordinates),
                 ggplot2::aes(x=-0.2*max, xend=-0.2*max,
                     y=ymax.ideo/2+2, yend=ymax.ideo/2-2)) +
    ggplot2::geom_segment(data=data.frame(coordinates),
                 ggplot2::aes(x=-0.06*max, xend=-0.06*max,
                     y=ymax.ideo/2+2, yend=ymax.ideo/2-2)) +
    ggplot2::annotate("text",
             -0.2*max+((0.2*max - 0.06*max)/2),
             ymax.ideo/2+8, label=lab.scale,
             size=4) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position="none",
          strip.text.y = ggplot2::element_blank()) +
    ggplot2::xlim(-space*max, max+max*space) +
    ggplot2::ylim(-10, ymax)

  return(ideo)
}
