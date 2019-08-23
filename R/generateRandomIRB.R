
generateRandomIRB <- function() {
  randomIRB <- list()

  width <- floor(runif(1, 100, 5e5)) # randomly select region width

  # Randomize coordinates to show
  coordinates <- GRanges(paste0("chr1:1-", width))
  coordinates <- regioneR::randomizeRegions(coordinates, mask=T)
  randomIRB[["coordinates"]] <- coordinates

  # Randomize datasets
  for (i in 1:length(IRB)) {
    sel <- IRB[[i]]
    n <- round(runif(1, 1, length(sel)))
    randomIRB[[names(IRB)[i]]] <- sel[n]
  }

  message(rep("-",
              nchar(paste0("-- Random set generated at ",
                           as.character(coordinates), " --"))),
          "\n-- Random set generated at ",
          as.character(coordinates),
          " --\n", rep("-",
                       nchar(paste0("-- Random set generated at ",
                                         as.character(coordinates), " --"))),
          "\nClusters: ", randomIRB$clusters,
          "\nMaps: ", randomIRB$maps,
          "\nSNPs: ", randomIRB$snps,
          "\nTFs: ", randomIRB$tfs, "\n",
          rep("-",
              nchar(paste0("-- Random set generated at ",
                           as.character(coordinates), " --"))))

  return(randomIRB)
}
