#' Calculate CNV based on normalized expression matrix
#'
#' expression matrix should be numeric values matrix in which  each row is a gene, and each column is a cell or sample.
#' The rownames for this expression matrix  should be ensembl IDs or gene symbols
#' One should set the species as human or mouse.
#' The genomic positions is human_geneInfo_genecode_v25 or mouse_geneInfo_genecode_vM12
#'
#' @param exprSet The expression matrix(which shoud be normalized,like log2(cpm+1) or log2(tpm+1))
#' @param geneType Choose ensembl or symbol,defaults: ensembl
#' @param species Choose human or mouse,defaults: human
#'
#' @return all_cnv CNV value matrix for each gene in each cell
#' @examples
#' exprSet2CNV
#' exprSet2CNV(exprSet,geneType='ensembl',species='human')

exprSet2CNV <- function(exprSet,geneType='ensembl',species='human' ){

  pos = get_genomic_positions(rownames(exprSet),geneType ,species )
  exprSet=exprSet[rownames(exprSet) %in% pos[,1],];dim(exprSet)
  pos=pos[pos[,1] %in% rownames(exprSet),];dim(pos)
  exprSet=exprSet[pos[,1],]


  res=cbind(pos,exprSet)
  table(res$chr)
  all_cnv <- lapply(split(res,res$chr), function(x){
    # x=split(res,res$chr)[[1]]
    anno=x[,1:4]
    ## the expression matrix for each chromosome
    dat=x[,5:ncol(x)]
    # At first, expression matrix is log2(cpm+1), we need to scale it by gene.
    # Then, we defined relative expression by centering the expression levels, Er[i,j]=E[i,j]-average(E[i,1...n]).
    dat=apply(dat, 1, function(x) x-mean(x))
    dat=t(dat)
    ## Then, To avoid considerable impact of any particular gene on the moving average
    ## we limited the relative expression values to [-3,3] by replacing all values above 3 by 3,
    ## and replacing values below -3 by -3.
    dat[dat>3]=3
    dat[dat < -3 ] = -3
    if(nrow(dat)>100){
      cnv <- lapply(51:(nrow(dat)-50), function(i){
        this_cnv <- unlist( lapply(1:ncol(dat), function(j){
          sum(dat[(i-50):(i+50),j])/101
        }))
        return(this_cnv)
      })
      cnv=do.call(rbind,cnv)
      cnv=cbind(anno[51:(nrow(x)-50),],cnv)
      # cnv[1:4,1:8]
    }else{
      return(NULL)
    }
  })
  all_cnv=do.call(rbind,all_cnv)
  head(all_cnv[1:4,1:8])
  table(all_cnv$chr)
  colnames(all_cnv)[5:ncol(all_cnv)]=colnames(exprSet)
  return(all_cnv)

}

