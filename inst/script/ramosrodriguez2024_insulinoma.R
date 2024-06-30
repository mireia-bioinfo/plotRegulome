tab_irds <- "~/Downloads/Table S6 - IRD.xlsx"
ird <- readxl::read_excel(tab_irds, skip=1)
download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
              destfile = "hg38ToHg19.over.chain.gz")



library(GenomicRanges)
library(rtracklayer)
library(dplyr)
ird <- GRanges(ird)
ird_19 <- rtracklayer::liftOver(
  ird,
  chain = rtracklayer::import.chain("hg38ToHg19.over.chain")
)
ird_19 <- unlist(ird_19)


load("PIRB/IRB_database/hg19/maps/hg19_map_chromatinClasses.rda")

clusters <- ird_19
mcols(clusters) <- NULL
clusters$class <- "CDS"

clusters.n <- "Insulinoma Regulatory Domains"
clusters.l <- TRUE
save(clusters, clusters.n, clusters.l, file="PIRB/IRB_database/hg19/clusters/hg19_cluster_irds.rda")


## Regulatory Elements
tab_res <- "~/Downloads/Table S3 - Differential H3K27ac.xlsx"
res <- readxl::read_excel(tab_res, skip=1)
res <- GRanges(res)
res_19 <- rtracklayer::liftOver(
  res,
  chain = rtracklayer::import.chain("hg38ToHg19.over.chain")
)
res_19 <- unlist(res_19) %>% 
  plyranges::filter(Type == "Gained") 

type <- data.frame(res_19) %>% 
  mutate(type = case_when(
    IRD == "Yes" & H3K27me3.in.controls == "Yes" ~ "DeIRD",
    IRD == "Yes" & H3K27me3.in.controls == "No" ~ "IRD",
    IRD == "No" & H3K27me3.in.controls == "No" ~ "Other gained",
    IRD == "No" & H3K27me3.in.controls == "Yes" ~ "Other gained"
  )) %>% 
  pull(type)

mcols(res_19) <- NULL
res_19$type <- factor(type, levels=c("DeIRD", "IRD", "Other gained"))
map <- res_19
map.name <- "Regulatory Elements/Insulinoma"
map_col <- c("#7F2704", "#05684B", "#51C398")
names(map_col) <- levels(map$type)
save(map, map.name, map_col, file="PIRB/IRB_database/hg19/maps/hg19_map_insREs.rda")
