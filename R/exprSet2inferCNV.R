#' Create two input files for inferCNV based on expression matrix
#'
#' expression matrix should be numeric values matrix in which  each row is a gene, and each column is a cell or samples.
#' The rownames for this expression matrix  should be ensembl IDs or gene symbols
#' One should set the species as human or mouse.
#' Check the document for input format at : https://github.com/broadinstitute/inferCNV/wiki
#' The genomic positions for airway data is : human_geneInfo_genecode_v25
#'
#' @param exprSet The expression matrix(which shoud be normalized,like log2(cpm+1) or log2(tpm+1))
#' @param geneType Choose ensembl or symbol,defaults: ensembl
#' @param species Choose human or mouse,defaults: human
#' @param prefix The prefix for files,defaults:example
#' @param dir Choose where to put the files,defaults: ./
#'
#' @return Two files \code{example_inferCNV_pos.txt} and \code{example_inferCNV_exprSet.txt}
#' @examples
#' exprSet2inferCNV
#' exprSet2inferCNV(exprSet,geneType='ensembl',species='human',prefix='airwary',dir=dir)

exprSet2inferCNV <- function(exprSet,geneType='ensembl',species='human',prefix='example',dir='./'){

  if(species=='human'){
    head(human_geneInfo_genecode_v25)
    if(geneType=='ensembl'){
      pos=human_geneInfo_genecode_v25
      exprSet=exprSet[rownames(exprSet) %in% pos$ensembl,]
      dim(exprSet)
      pos=pos[match(rownames(exprSet),pos$ensembl),c(4,1:3)]
    }
    if(geneType=='symbol'){
      pos=human_geneInfo_genecode_v25
      exprSet=exprSet[rownames(exprSet) %in% pos$symbol,]
      dim(exprSet)
      pos=pos[match(rownames(exprSet),pos$symbol),c(6,1:3)]
    }
    new_chr=gsub('chr','',pos$chr)
    table(new_chr)
    new_chr[new_chr=='X']=23
    new_chr[new_chr=='Y']=24
    new_chr=as.numeric(new_chr)
    pos$chr=new_chr
    pos=pos[order(pos$chr,pos$start),]
  }
  if(species=='mouse'){
    head(mouse_geneInfo_genecode_vM12)
    if(geneType=='ensembl'){
      pos=mouse_geneInfo_genecode_vM12
      exprSet=exprSet[rownames(exprSet) %in% pos$ensembl,]
      dim(exprSet)
      pos=pos[match(rownames(exprSet),pos$ensembl),c(4,1:3)]
    }
    if(geneType=='symbol'){
      pos=mouse_geneInfo_genecode_vM12
      exprSet=exprSet[rownames(exprSet) %in% pos$symbol,]
      dim(exprSet)
      pos=pos[match(rownames(exprSet),pos$symbol),c(5,1:3)]
    }
    new_chr=gsub('chr','',pos$chr)
    table(new_chr)
    new_chr[new_chr=='X']=20
    new_chr[new_chr=='Y']=21
    new_chr=as.numeric(new_chr)
    pos$chr=new_chr
    pos=pos[order(pos$chr,pos$start),]
  }

  write.table(pos,file.path(dir,paste0(prefix,'_inferCNV_pos.txt')),row.names = F,col.names = F,sep = '\t',quote = F)
  write.table(exprSet,file.path(dir,paste0(prefix,'_inferCNV_exprSet.txt')),quote = F,sep = '\t')

  # Rscript ~/biosoft/scRNA_cnv/inferCNV/scripts/inferCNV.R    \
  # --output_dir  test airwary_inferCNV_exprSet.txt airwary_inferCNV_pos.txt

}

