#' Draw heatmap for CNV matrix
#'
#' CNV matrix calculated base on expression matrix by exprSet2CNV
#'
#' @param cnv The CNV results from exprSet2CNV
#' @param meta Choose ensembl or symbol,defaults: ensembl
#' @param prefix The prefix of the filename of the PDF heatmap
#' @param noise_filter A value must be atleast this much more or less than the reference to be plotted [Default 0.2]
#' @param upper The maximum value for heatmap  [Default 2]
#' @param species Choose human or mouse,defaults: human
#'
#' @return heatmap
#' @examples
#' ht_cnv
#' exprSet2CNV(exprSet,geneType='ensembl',species='human')

ht_cnv <- function(cnv,meta,prefix='test',noise_filter=0.2,upper=2){
  library(pheatmap)
  all_cnv=cnv
  D=t(scale(all_cnv[,5:ncol(all_cnv)] ))
  apply(D, 1, summary)
  apply(D[,1:10], 2, summary)
  D[D> upper]=upper
  D[D< -upper] = -upper
  D[abs(D) < noise_filter]=0

  dim(D)
  colnames(D)=paste0('genes_',1:ncol(D))
  rownames(D)=colnames(exprSet)

  require(RColorBrewer)
  cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

  library(stringr)
  annotation_row = data.frame(
    patients=str_split(rownames(D),'_',simplify = T)[,1]
  )

  rownames(annotation_row) =  rownames(D)

  annotation_col = data.frame(
    chr= factor(all_cnv$chr,levels =  unique(all_cnv$chr))
  )
  rownames(annotation_col) = colnames(D)
  pheatmap(D,cluster_rows = T,col=rev(cols),
           annotation_col=annotation_col,
           annotation_row = annotation_row,
           cluster_cols = F,show_rownames=F,show_colnames=F,filename=paste0(prefix,'_cnv.pdf'))
}
