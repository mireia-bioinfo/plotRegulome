#' Blank X axis ggplot2 theme
#'
#' When added to a ggplot, will set all elements in the x axis to blank.
#' @param ... Additional parameters to pass to \link[ggplot2]{theme} function.
#' @export
themeXblank <- function(...) {
  theme <- cowplot::theme_cowplot() + ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                          axis.title.x=ggplot2::element_blank(),
                          axis.ticks.x=ggplot2::element_blank(),
                          axis.line.x=ggplot2::element_blank(),
                          ...)
  return(theme)
}

#' Blank Y axis ggplot2 theme
#'
#' When added to a ggplot, will set all elements in the y axis to blank.
#' @param title Logical indicating if the title of the axis should be included. Default = FALSE.
#' @param ... Additional parameters to pass to \link[ggplot2]{theme} function.
#' @export
themeYblank <- function(title=FALSE, ...) {
  if(title) {
    theme <- cowplot::theme_cowplot() + ggplot2::theme(axis.text.y=ggplot2::element_blank(),
                                                       axis.ticks.y=ggplot2::element_blank(),
                                                       axis.line.y=ggplot2::element_blank(),
                                                       ...)
  } else {
    theme <- cowplot::theme_cowplot() + ggplot2::theme(axis.text.y=ggplot2::element_blank(),
                                                       axis.title.y=ggplot2::element_blank(),
                                                       axis.ticks.y=ggplot2::element_blank(),
                                                       axis.line.y=ggplot2::element_blank(),
                                                       ...)
  }

  return(theme)
}

#' Blank X and Y axis ggplot2 theme
#'
#' When added to a ggplot, will set all elements in the x and y axis to blank.
#' @param title Logical indicating if the title of the axis should be included. Default = FALSE.
#' @param ... Additional parameters to pass to \link[ggplot2]{theme} function.
#' @export
themeXYblank <- function(title=FALSE, ...) {
  if(title) {
    theme <- cowplot::theme_cowplot() + ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                                                       axis.title.x=ggplot2::element_blank(),
                                                       axis.ticks.x=ggplot2::element_blank(),
                                                       axis.line.x=ggplot2::element_blank(),
                                                       axis.text.y=ggplot2::element_blank(),
                                                       axis.ticks.y=ggplot2::element_blank(),
                                                       axis.line.y=ggplot2::element_blank(),
                                                       ...)
  } else {
    theme <- cowplot::theme_cowplot() + ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                                                       axis.title.x=ggplot2::element_blank(),
                                                       axis.ticks.x=ggplot2::element_blank(),
                                                       axis.line.x=ggplot2::element_blank(),
                                                       axis.text.y=ggplot2::element_blank(),
                                                       axis.title.y=ggplot2::element_blank(),
                                                       axis.ticks.y=ggplot2::element_blank(),
                                                       axis.line.y=ggplot2::element_blank(),
                                                       ...)
  }
}

#' X axis for genomic coordinates
#'
#' When added to a ggplot, will set the X axis to represent genomic coordintes data, by rounding
#' X axis labels to Mb and showing the chromosome number as an X axis title.
#' @param chr Logical indicating if the title of the axis should be included. Default = FALSE.
#' @param ... Additional parameters to pass to \link[ggplot2]{scale_x_continuous} function.
#' @export
scaleXCoordinates <- function(chr,
                              ...) {
  scale <- ggplot2::scale_x_continuous(name=paste0("Coordinates ", chr, " (Mb)"),
                              labels=function(x) round(x/1e6, 2),
                              ...)
  return(scale)
}
