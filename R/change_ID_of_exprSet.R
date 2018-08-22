#' Change the colnames of expression matrix between ensembl and symbol
#'
#' expression matrix should be numeric values matrix in which  each row is a gene, and each column is a cell or samples.
#' The rownames for this expression matrix  should be ensembl IDs or gene symbols
#' One should set the species as human or mouse.
#' The expression value will be sum if two or more genes to one gene
#'
#' @param exprSet The expression matrix
#' @param geneType Choose ensembl or symbol,defaults: ensembl
#' @param species Choose human or mouse,defaults: human
#'
#' @return exprSet The expression matrix which had been changed.
#' @examples
#' change_ID_of_expreSet
#' change_ID_of_expreSet(exprSet,geneType='ensembl',species='human' )

change_ID_of_expreSet <- function(exprSet,geneType='ensembl',species='human' ){
  print(dim(exprSet))
  pos=data.frame()
  if(species=='human'){
    pos=human_geneInfo_genecode_v25[,c('ensembl','symbol')]
  }else{
    pos=mouse_geneInfo_genecode_vM12[,c('ensembl','symbol')]
  }
  head(pos)
  if(geneType=='ensembl'){
    exprSet=exprSet[rownames(exprSet) %in% pos$ensembl,]
    print(dim(exprSet))
    pos=pos[match(rownames(exprSet),pos$ensembl), ]
    #tmp=split(as.data.frame(exprSet),pos$symbol);x=tmp$ZNF385C
    tmp=lapply(split(as.data.frame(exprSet),pos$symbol), function(x){
      if(class(x)== "data.frame"){
        return(colSums(x))
      }else(return(x))
    })
    newExprSet <- do.call(rbind,tmp)
    # head(newExprSet);dim(newExprSet)
    colnames(newExprSet)=colnames(exprSet)
    print(dim(newExprSet))
  }
  if(geneType=='symbol'){
    exprSet=exprSet[rownames(exprSet) %in% pos$symbol,]
    print(dim(exprSet))
    pos=pos[match(rownames(exprSet),pos$symbol), ]
    #tmp=split(as.data.frame(exprSet),pos$symbol);x=tmp$ZNF385C
    tmp=lapply(split(as.data.frame(exprSet),pos$ensembl), function(x){
      if(class(x)== "data.frame"){
        return(colSums(x))
      }else(return(x))
    })
    newExprSet <- do.call(rbind,tmp)
    # head(newExprSet);dim(newExprSet)
    colnames(newExprSet)=colnames(exprSet)
    print(dim(newExprSet))
  }
  return(newExprSet)
}

