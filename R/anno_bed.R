#' Assign a position to genomic features by ChIPseeker
#'
#' filter SNP or INDELs in a vcf or maf file
#'
#' @param pos  three columns of the positions, chromosome,start,end
#' @param reference Choose hg19,hg38,mm10 ,defaults: hg38
#'
#' @return assigned results
#' @examples
#' anno_bed
#' anno_bed(pos,reference='hg38' )

anno_bed <- function(pos,reference='hg38' ){
  require(ChIPseeker)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(GenomicRanges)
  peak <- GRanges(seqnames=Rle(pos[,1]),
                 ranges=IRanges(pos[,2], pos[,3]), strand=rep(c("*"), nrow(pos)))
  peak

  if(reference=='hg38'){
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    txdb=TxDb.Hsapiens.UCSC.hg38.knownGene
    peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")
  }
  if(reference=='hg19'){
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
    peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")
  }
  if(reference=='mm10'){
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    txdb=TxDb.Mmusculus.UCSC.mm10.knownGene
    peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                             TxDb=txdb, annoDb="org.Mm.eg.db")
  }

  return(as.data.frame(peakAnno))
}
