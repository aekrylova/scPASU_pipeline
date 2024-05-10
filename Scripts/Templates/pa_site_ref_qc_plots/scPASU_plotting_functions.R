# This code is adapted from Ting lab's bulk APA pipeline to repurpose it for scRNA APA analysis
# Author: Surbhi Sona

# original script from bulkAPA pipeline: '08_pas_stats_figure.r'
# hg38 polya db: atlas.clusters.2.0.GRCh38.96.bed  (downloaded from :https://polyasite.unibas.ch/atlas#2)
# hg19 db: polyAsite_v1_hg19_clusters.bed


library(data.table)
library(goldmine)
#library(handy)
library(IRanges)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm10)

#### Hexamer plotting function ####

countHexes <- function(gr,us=50)
{
  hexes <- c("AATAAA","ATTAAA","AGTAAA","TATAAA","CATAAA","GATAAA","AATATA","AATACA","AATAGA","AATGAA","ACTAAA","AACAAA","TTTAAA")
  
  message("Generating flanks")
  
  # Also including the PR itself in addition to the flank here
  fl <- resize(gr,us+width(gr),fix="end")
  #fl <- flank(re,width=us)
  
  
  message("Pulling sequence data from hg38")
  seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,names=fl)
  
  message("Searching hexamer occurances")
  pd <- PDict(x=hexes,max.mismatch=0)
  vc <- vcountPDict(pd, seq)
  vc <- t(vc)
  
  
  message("Selecting matches")
  ap <- apply(vc,1,function(x) which(x>=1)[1])
  hp <- hexes[ap]
  gr$hex <- hp
  gr[rowSums(vc)==0]$hex <- "NONE"
  #table(gr$hex)
  
  gr$hex <- factor(gr$hex,levels=c(hexes,"NONE"))
  
  return(gr)
}

# Find the most up/downstream hexamer 

findHex <- function(gr,us=50,direction='upstream',most.downstream=TRUE,
                    hexes=c("AATAAA","ATTAAA","AGTAAA","TATAAA",
                            "CATAAA","GATAAA","AATATA","AATACA",
                            "AATAGA","AATGAA","ACTAAA","AACAAA","TTTAAA")) {
  
  cat("Generating flanks",us,'bp',direction,'\n')
  
  # Also including the PR itself in addition to the flank here
  if (direction == 'upstream') {
    fl <- resize(gr,us+width(gr),fix="end")
  } else if (direction == 'downstream') {
    fl <- resize(gr,us+width(gr),fix="start")  
  } else {stop('Direction should be either upstream or downstream \n')}
  
  cat("Pulling sequence data from hg38 \n")
  seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,names=fl)
  
  cat("Locate most",if (most.downstream) {paste0('downstream')} else {paste0('upstream')},"hexamer if any \n")
  pd <- PDict(x=hexes,max.mismatch=0)
  mc <- lapply(seq,function(x) {
    mcoord <- matchPDict(pd,x) %>% as.data.frame
    if (nrow(mcoord) == 0) {
      return('NONE')
    } else {
      if (most.downstream) {
        most_ds_hex_idx <- mcoord$group[which(mcoord$start == max(mcoord$start))]
      } else {
        most_ds_hex_idx <- mcoord$group[which(mcoord$start == min(mcoord$start))]
      }
      most_ds_hex <- hexes[most_ds_hex_idx]
      return(most_ds_hex)
    }
  })
  mc <- do.call('c',mc)
  
  gr$hex <- mc
  
  gr$hex <- factor(gr$hex,levels=c(hexes,"NONE"))
  
  return(gr)
}
