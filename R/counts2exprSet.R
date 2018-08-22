#' transform and normlize the raw counts matrix
#'
#' row counts matrix will be filtered and log2(cpm+1) tranformed.
#'
#' @param counts The raw counts matrix from featureCounts or other tools
#'
#' @return exprSet, the normalized expression matrix, log2(cpm+1)
#' @examples
#' counts2exprSet
#' counts2exprSet(counts)

counts2exprSet <- function(counts){
  library(edgeR)
  library(DESeq2)
  exprSet=counts
  geneLists=rownames(exprSet)
  keepGene=rowSums(cpm(exprSet)>0) >=2
  table(keepGene);dim(exprSet)
  dim(exprSet[keepGene,])
  exprSet=exprSet[keepGene,]
  rownames(exprSet)=geneLists[keepGene]

 # boxplot(exprSet,las=2)
  # CPM normalized counts.
  exprSet=log2(cpm(exprSet)+1)
  # boxplot(exprSet,las=2)
  exprSet[1:4,1:4]
  return(exprSet)
}
