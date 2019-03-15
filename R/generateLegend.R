generateLegend <- function(contactsObject,
                           mapsObject,
                           clustersObject,
                           tfsObject) {
  lab <- c("Virtual 4C",
           mapsObject$name,
           "Other")
  pad <- 0.01
    ## Always plot Maps and Other (clusters+tfs), only plot contacts when present
  if (contactsObject$name=="" &
      mapsObject$name=="") {
    widths <- c("maps"=0.7,
                "other"=0.3)
    lab <- lab[-1]
    red <- c(pad/2, pad/2)

    ## Calculate text positions
  } else {
    widths <- c("contacts"=0.22,
                "maps"=0.58,
                "other"=0.2)
    red <- c(pad/2, pad, pad/2)
  }

  ## Calculate X positions for boxes

  widths <- widths - red

  positions.x.all <- c(0)
  for (i in 1:length(widths)) {
    if (i < length(widths)) {
      ele1 <- pad + sum(widths[1:i])-pad/2
      ele2 <- pad + sum(widths[1:i])+pad/2
      positions.x.all <- c(positions.x.all,
                           ele1, ele2)
    } else {
      ele <- pad + sum(widths) + pad
      positions.x.all <- c(positions.x.all,
                           ele)
    }
  }

  # Create indexes for odd positions and real widths after scaling
  idx <- 1:length(positions.x.all)
  idx <- idx[idx %% 2==1]
  widths <- positions.x.all[idx+1] - positions.x.all[idx]

  ## Maximum y position (text)
  ymax <- 0.9

  ## Virtual 4C legend
  if (contactsObject$name!="") {
    virt4c <- returnPositions(nele=length(contactsObject$col),
                              ncol=1, colors=contactsObject$col,
                              start=positions.x.all[1],
                              end=positions.x.all[2],
                              pad=pad, ymax=ymax)
    virt4c$label <- paste0("CHiCAGO score ", c("<3", "3-5", ">5"))
    virt4c$type <- 22
    virt4c$size <- 4
    positions.noVirt <- positions.x.all[-c(1:2)]
  } else {
    positions.noVirt <- positions.x.all
  }

  ## Main legend
  main <- returnPositions(nele=length(mapsObject$col),
                            ncol=3, colors=mapsObject$col,
                            start=positions.noVirt[1],
                            end=positions.noVirt[2],
                            pad=pad, ymax=ymax)
  main$label <- names(mapsObject$col)
  main$type <- 22
  main$size <- 4

  ## Other legend
  if (contactsObject$name=="") {
    nele=2
    cols <- c(tfsObject$col,
              clustersObject$col)
    otlab <- c("TF binding", clustersObject$name)
    ottype <- c(20, 20)
    ottsize <- c(0.2,0.2)
  } else {
    nele=3
    cols <- c(tfsObject$col,
              clustersObject$col,
              "black")
    otlab <- c("TF binding", clustersObject$name, "Viewpoint")
    ottype <- c(20, 20, 17)
    ottsize <- c(NA,NA,4)
  }

  other <- returnPositions(nele=nele,
                           ncol=1, colors=cols,
                           start=positions.noVirt[3],
                           end=positions.noVirt[4],
                           pad=pad, ymax=ymax)
  other$label <-otlab
  other$type <- ottype
  other$size <- ottsize

  elementsLegend <- rbind(main, other)

  if (contactsObject$name!="") elementsLegend <- rbind(elementsLegend, virt4c)

  legend <- ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(xmin=positions.x.all[idx],
                  xmax=positions.x.all[idx+1],
                  ymin=0, ymax=ymax),
              color="black", fill=NA) +
    ggplot2::geom_tile(ggplot2::aes(x=positions.x.all[idx]+widths/2,
                  y=ymax,
                  width=widths/1.5, height=0.2),
              color="white", fill="white") +
    ggplot2::geom_text(ggplot2::aes(x=positions.x.all[idx]+widths/2,
                  y=ymax,
                  label=lab),
              fontface=2, size=4) +
    ggplot2::geom_point(data=elementsLegend,
                        ggplot2::aes(posX+pad, posY),
               pch=elementsLegend$type, fill=elementsLegend$col,
               color="black", size=elementsLegend$size) +
    ggplot2::geom_text(data=elementsLegend,
                       ggplot2::aes(posX+pad*2, posY, label=label),
            hjust=0, size=3.5)+
    ggplot2::geom_segment(data=elementsLegend[elementsLegend$type==20,],
                          ggplot2::aes(x=posX-pad/2, xend=posX+pad+pad/2,
                                       y=posY, yend=posY),
                          color=elementsLegend$col[elementsLegend$type==20]) +
    ggplot2::ylim(0,1) +
    ggplot2::theme_void()

  return(legend)
}

returnPositions <- function(nele, ncol,
                            colors,
                            start, end,
                            pad,
                            ymax=0.9) {
  positions.x <- seq(start+pad, end-pad, length=ncol+1)
  positions.x <- positions.x[-length(positions.x)]

  if (ceiling(nele/ncol)>2) {
    positions.y <- rev(seq(0+pad*20, ymax-pad*20, length=ceiling(nele/ncol)))
  } else {
    positions.y <- rev(seq(0+pad*20, ymax-pad*20, length=ceiling(nele/ncol)+1))
    positions.y <- positions.y[-length(positions.y)]
  }

  loc <- data.frame(posX=rep(positions.x, length(positions.y))[1:nele],
                    posY=rep(positions.y, each=length(positions.x))[1:nele],
                    col=colors)
  return(loc)
}

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
                                "Viewpoint"),
                      "pch"=c(rep(15, 3), 25))
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
                                            ggplot2::aes(x,y, color=class),
                                            pch=df1$pch, size=4, fill="black"),
                          ggplot2::scale_color_manual(values=contactsObject$col,
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
    tmp <- ggplot2::ggplot() +
      legends +
      ggplot2::guides(fill=ggplot2::guide_legend(ncol=1,order=2),
                      color=ggplot2::guide_legend(ncol=1,order=1),
                      size=ggplot2::guide_legend(ncol=1,order=3)) +
      ggplot2::theme(legend.margin=ggplot2::margin(10,10,10,10,unit="pt"),
                     legend.text=ggplot2::element_text(size=12),
                     legend.title=ggplot2::element_text(size=14, face="bold"),
                     legend.key	= ggplot2::element_rect(fill=NA),
                     legend.key.size=ggplot2::unit(12,"pt"),
                     legend.box="vertical")

    legend <- cowplot::get_legend(tmp)
    # plot_grid(legend)
  } else {
    t <- emptyPlot(mapsObject$coordinates, 1)
    legend <- ggplot2::ggplot() + t + ggplot2::theme_void()
  }

  return(legend)
}
