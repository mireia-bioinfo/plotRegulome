generateLegendGG <- function(contactsObject,
                             mapsObject,
                             clustersObject,
                             tfsObject,
                             snpsObject) {

  ## Automatic legends ------------------
  leg.list <- list(snps=snpsObject,
                   maps=mapsObject)
  leg.plot <- list()
  for (i in 1:length(leg.list)) {
    if (leg.list[[i]]$name!="") {
      leg <- ggplot2::ggplot() + plot(leg.list[[i]])
      leg.plot[[names(leg.list[i])]] <- leg
    }
  }

  ## Manual legends ----------------------
  ## Contacts Object
  if (contactsObject$name!="") {
    df1 <- data.frame("x"=rep(1, length(contactsObject$col)+1),
                      "y"=rep(1, length(contactsObject$col)+1),
                      "class"=c("0", "3", "5", "Viewpoint"),
                      stringsAsFactors = F)
    contactsObject$col <- c(contactsObject$col, "black")
    names(contactsObject$col) <- df1$class

  } else {
    df1 <- data.frame()
  }

  ## Other Object
  if (tfsObject$name!="") tfs.name="TF binding" else tfs.name=""

  df3 <- data.frame(x=1:2,
                    y=1:2,
                    "class"=c(clustersObject$name,
                              tfs.name),
                    "color"=c(clustersObject$col,
                              tfsObject$col))
  df3 <- df3[df3$class!="",]

  df3.plot <- rbind(df3, df3)

  if (nrow(df3.plot)>0) {
    df3.plot$x <- seq(1, nrow(df3.plot))
    df3.plot$y <- seq(1, nrow(df3.plot))
  }


  ## Manually create legends
  contactLegend <- list(ggplot2::geom_point(data=df1,
                                            ggplot2::aes(x,y, color=class, shape=class),
                                            size=4, fill="black"),
                          ggplot2::scale_color_manual(values=contactsObject$col,
                                                      name="Virtual 4C",
                                                      labels=c("0"="CHiCAGO score <3",
                                                               "3"="CHiCAGO score 3-5",
                                                               "5"="CHiCAGO score >5",
                                                               "Viewpoint"=bquote('Viewpoint ('~italic(.(contactsObject$name))~")"))),
                          ggplot2::scale_shape_manual(values=c(rep(15, 3), 25),
                                                      name="Virtual 4C",
                                                      labels=c("0"="CHiCAGO score <3",
                                                               "3"="CHiCAGO score 3-5",
                                                               "5"="CHiCAGO score >5",
                                                               "Viewpoint"=bquote('Viewpoint ('~italic(.(contactsObject$name))~")"))))

  otherLegend <- list(ggplot2::geom_line(data=df3.plot,
                                         ggplot2::aes(x,
                                                      y,
                                                      size=class,
                                                      color=class)),
                      ggplot2::scale_size_manual(values=c(2,1),
                                                   name="Other"),
                      ggplot2::scale_color_manual(values=unique(as.character(df3$color)),
                                                  name="Other"))

  legends <- list(ggplot2::ggplot() + contactLegend,
                  ggplot2::ggplot() + otherLegend)
  dfs <- list(df1, df3)
  idx <- sapply(dfs, nrow)
  idx <- idx>0
  legends <- legends[idx]

  ## Add to other legends
  if (length(legends)==2) {
    leg.plot <- c(legends[1],
                  leg.plot,
                  legends[2])
  } else {
    leg.plot <- c(leg.plot, legends)
  }

  if (length(leg.plot)>0) {
    theme_legend <- ggplot2::theme(legend.margin=ggplot2::margin(10,10,10,10,unit="pt"),
                                   legend.text=ggplot2::element_text(size=12),
                                   legend.title=ggplot2::element_text(size=14, face="bold"),
                                   legend.key	= ggplot2::element_rect(fill=NA),
                                   legend.key.size=ggplot2::unit(12,"pt"),
                                   legend.box="vertical")

    ## Check if a plot has no legend
    grobs <- lapply(leg.plot, function(x) cowplot::plot_to_gtable(x)$grobs)
    legendIndex <- as.logical(sapply(lapply(grobs,
                                 function(x) which(sapply(x, function(x) x$name) == "guide-box")),
                          length))

    leg.all <- lapply(leg.plot[legendIndex],
                      function(x) cowplot::get_legend(x + theme_legend))


    legend <- cowplot::plot_grid(plotlist=leg.all,
                                 ncol=1,
                                 align="v")
  } else {
    t <- emptyPlot(mapsObject$coordinates, 1)
    legend <- ggplot2::ggplot() + t + ggplot2::theme_void()
  }

  return(legend)
}
