generateLegendGG <- function(contactsObject,
                             mapsObject,
                             clustersObject,
                             tfsObject) {
  ## Maps Object
  df2 <- data.frame("x"=rep(1, length(mapsObject$col)),
                    "y"=rep(1, length(mapsObject$col)),
                    "class"=names(mapsObject$col))

  ## Contacts Object
  if (contactsObject$name!="") {
    df1 <- data.frame("x"=rep(1, length(contactsObject$col)+1),
                      "y"=rep(1, length(contactsObject$col)+1),
                      "class"=c(paste("CHiCAGO score", c("<3", "3-5", ">5")),
                                paste0("Viewpoint (", contactsObject$name, ")")))
    contactsObject$col <- c(contactsObject$col, "black")
    names(contactsObject$col) <- df1$class
  } else {
    df1 <- data.frame()
  }


  ## Other Object
  if (tfsObject$name!="") tfs.name="TF binding" else tfs.name=""

  df3 <- data.frame("x"=1:2,
                    "y"=1:2,
                    "class"=c(clustersObject$name,
                              tfs.name),
                    "color"=c(clustersObject$col,
                              tfsObject$col))
  df3 <- df3[df3$class!="",]

  ## Create legends
  mapLegend <- list(ggplot2::geom_bar(data=df2,
                                      ggplot2::aes(x, y, fill=class),
                                      stat="identity"),
                    ggplot2::scale_fill_manual(values=mapsObject$col,
                                                 name=gsub("/", "\n", mapsObject$name)))

  contactLegend <- list(ggplot2::geom_point(data=df1,
                                            ggplot2::aes(x,y, color=class, shape=class),
                                            size=4, fill="black"),
                          ggplot2::scale_color_manual(values=contactsObject$col,
                                                      name="Virtual 4C"),
                          ggplot2::scale_shape_manual(values=c(rep(15, 3), 25),
                                                      name="Virtual 4C"))

  otherLegend <- list(ggplot2::geom_line(data=df3,
                                         ggplot2::aes(x, y, size=class),
                                         color=df3$color),
                      ggplot2::scale_size_manual(values=c(2,1),
                                                   name="Other"))

  legends <- list(contactLegend,
                  mapLegend,
                  otherLegend)

  dfs <- list(df1, df2, df3)
  idx <- sapply(dfs, nrow)
  idx <- idx>0

  legends <- legends[idx]

  if (length(legends)>0) {
    # tmp <- ggplot2::ggplot() +
    #   legends +
    #   ggplot2::guides(fill=ggplot2::guide_legend(ncol=1,order=2),
    #                   color=ggplot2::guide_legend(ncol=1,order=1),
    #                   size=ggplot2::guide_legend(ncol=1,order=3)) +
    #   ggplot2::theme(legend.margin=ggplot2::margin(10,10,10,10,unit="pt"),
    #                  legend.text=ggplot2::element_text(size=12),
    #                  legend.title=ggplot2::element_text(size=14, face="bold"),
    #                  legend.key	= ggplot2::element_rect(fill=NA),
    #                  legend.key.size=ggplot2::unit(12,"pt"),
    #                  legend.box="vertical")
    #
    # legend <- cowplot::get_legend(tmp)
    # plot_grid(legend)
    theme_legend <- ggplot2::theme(legend.margin=ggplot2::margin(10,10,10,10,unit="pt"),
                                   legend.text=ggplot2::element_text(size=12),
                                   legend.title=ggplot2::element_text(size=14, face="bold"),
                                   legend.key	= ggplot2::element_rect(fill=NA),
                                   legend.key.size=ggplot2::unit(12,"pt"),
                                   legend.box="vertical")

    leg.cont <- cowplot::get_legend(ggplot2::ggplot() + contactLegend + theme_legend)
    leg.map <- cowplot::get_legend(ggplot2::ggplot() + mapLegend + theme_legend)
    leg.other <- cowplot::get_legend(ggplot2::ggplot() + otherLegend + theme_legend)

    legend <- cowplot::plot_grid(leg.cont,
                                 leg.map,
                                 leg.other,
                                 ncol=1)
  } else {
    t <- emptyPlot(mapsObject$coordinates, 1)
    legend <- ggplot2::ggplot() + t + ggplot2::theme_void()
  }

  return(legend)
}
