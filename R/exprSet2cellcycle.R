#' Assign cell cycle status based on expression matrix
#'
#' The rownames of expression matrix, should be ensembl IDs or gene symbols
#' One should set the species as human or mouse.
#'
#' @param exprSet The expression matrix(which shoud be normalized,like log2(cpm+1) or log2(tpm+1))
#' @param geneType Choose ensembl or symbol,defaults: ensembl
#' @param species Choose human or mouse,defaults: human
#'
#' @return assigned A list of (phases,scores,normalized.scores) return from cyclone(scran)
#' @examples
#' exprSet2cellcycle
#' exprSet2cellcycle(exprSet,geneType='ensembl',species='human')

exprSet2cellcycle <- function(exprSet,geneType='ensembl',species='human' ){
  library(scran)
  sce <- SingleCellExperiment(list(counts=exprSet))
  if(species=='human'){
    library(org.Hs.eg.db)
    mm.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
    if(geneType=='ensembl'){
      assigned <- cyclone(sce, pairs=mm.pairs )
    }
    if(geneType=='symbol'){
      ensembl <- mapIds(org.Hs.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")
      assigned <- cyclone(sce, pairs=mm.pairs, gene.names=ensembl)
    }
  }
  if(species=='mouse'){
    library(org.Mm.eg.db)
    mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
    if(geneType=='ensembl'){
      assigned <- cyclone(sce, pairs=mm.pairs )
    }
    if(geneType=='symbol'){
      ensembl <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")
      assigned <- cyclone(sce, pairs=mm.pairs, gene.names=ensembl)
    }
  }
  # head(cycles$scores)
  # table(cycles$phases)
  # dat=cbind(cycles$score,cycles$phases)
  # colnames(dat)
  # attach(dat)
  # library(scatterplot3d)
  # scatterplot3d(G1, S, G2M, angle=20,color = rainbow(3)[as.numeric(as.factor(cycles$phases))],
  #               grid=TRUE, box=FALSE)
  # detach(dat)
  return(assigned)
}
