#' Calculate DM(distance to median) values based on cpm(counts per million) expression matrix
#'
#' expression matrix should be numeric values matrix in which  each row is a gene, and each column is a cell or sample.
#' The rownames for this expression matrix  should be ensembl IDs or gene symbols
#' One should set the species as human or mouse.
#' The gene length information calculated base onf GENCODE database.
#'
#' @param exprSet cpm(counts per million) expression matrix
#' @param geneType Choose ensembl or symbol,defaults: ensembl
#' @param species Choose human or mouse,defaults: human
#'
#' @return log10cv2_adj DM(distance to median) values for each gene
#' @examples
#' cal_DM
#' cal_DM(exprSet,geneType='ensembl',species='human')
#'
cal_DM <- function(exprSet,geneType='ensembl',species='human' ){
  # In paper : normalised read counts (reads per million)
  if(species=='human'){
    head(human_geneLength_genecode_v25)
    len=human_geneLength_genecode_v25
    if(geneType=='ensembl'){
      exprSet=exprSet[rownames(exprSet) %in% len$ensembl,]
      length_per_gene=len[match(rownames(exprSet) ,  len$ensembl ),2]
    }
    if(geneType=='symbol'){
      ## TODO
      len=merge(human_geneLength_genecode_v25,human_geneInfo_genecode_v25,by='ensembl')
      len=len[,c('symbol','length')]
      exprSet=exprSet[rownames(exprSet) %in% len$symbol,]
      length_per_gene=len[match(rownames(exprSet) ,  len$symbol ),2]
    }
  }
  if(species=='mouse'){
    head(mouse_geneLength_genecode_vM12)
    len=mouse_geneLength_genecode_vM12
    if(geneType=='ensembl'){
      exprSet=exprSet[rownames(exprSet) %in% len$ensembl,]
      length_per_gene=len[match(rownames(exprSet) ,  len$ensembl ),2]
    }
    if(geneType=='symbol'){
      ## TODO
      ## TODO
      len=merge(mouse_geneLength_genecode_vM12,mouse_geneInfo_genecode_vM12,by='ensembl')
      len=len[,c('symbol','length')]
      exprSet=exprSet[rownames(exprSet) %in% len$symbol,]
      length_per_gene=len[match(rownames(exprSet) ,  len$symbol ),2]

    }
  }

  ## step1: check the correlations among different characteristics of a expression matrix
  # exprSet should be a  log2(cpm+1) expression matrix.

  mean_per_gene <- apply(exprSet, 1, mean, na.rm = TRUE)
  sd_per_gene <- apply(exprSet, 1, sd, na.rm = TRUE)
  mad_perl_gene <-   apply(exprSet, 1, mad, na.rm = TRUE)
  cv_per_gene <- sd_per_gene/mean_per_gene
  cha <- data.frame(mean = log10(mean_per_gene),
                            sd = sd_per_gene,
                            mad=mad_perl_gene,
                            cv = cv_per_gene,
                    len=log10(length_per_gene))
  rownames(cha) <- rownames(exprSet)
  # pairs(cha)
  # It's clear that these characteristics are related with each.
  # plot(cha[,c(1,4)])
  # Squared coefficient of variation (CV2) vs. average normalized read count of genes
  # As gene expression levels increase, genes are more likely to show lower levels of variation.

  # step 2 :Compute rolling medians of CV2 across all samples.

  # https://jdblischak.github.io/singleCellSeq/analysis/cv-adjusted-wo-19098-r2.html
  library(zoo)
  # Order of genes by mean expression levels
  order_gene <- order( mean_per_gene )
  cv=cv_per_gene
  # Rolling medians of log10 squared CV by mean expression levels
  roll_medians_mean <- rollapply(log10(cv^2)[order_gene], width = 50, by = 25,
                                 FUN = median, fill = list("extend", "extend", "NA") )
  ## then change the NA values in the roll_medians_mean
  table(is.na(roll_medians_mean))
  ii_na <- which( is.na(roll_medians_mean) )
  roll_medians_mean[ii_na] <- median( log10(cv^2)[order_gene][ii_na] )
  names(roll_medians_mean) <- rownames(exprSet)[order_gene]

  # re-order rolling medians according to the expression matrix
  roll_medians_mean <- roll_medians_mean[  match(rownames(exprSet), names(roll_medians_mean) ) ]
  stopifnot( all.equal(names(roll_medians_mean), rownames(exprSet) ) )

  # adjusted coefficient of variation on log10 scale
  log10cv2_adj <-  log10( cv^2) - roll_medians_mean

  if(F){
    plot(log10cv2_adj,log10(mean_per_gene))
    #install.packages("basicTrendline")
    library(basicTrendline)
    trendline(log10cv2_adj,log10(mean_per_gene),model="line2P")
  }

  # step 2 :Compute rolling medians of gene length(log10) across all samples.

  order_gene <- order( log10(length_per_gene) )
  cv=log10cv2_adj
  roll_medians_length <- rollapply(cv[order_gene], width = 50, by = 25,
                                   FUN = median, fill = list("extend", "extend", "NA") )
  ## then change the NA values in the roll_medians_length
  table(is.na(roll_medians_length))
  ii_na <- which( is.na(roll_medians_length) )
  roll_medians_length[ii_na] <- median( cv[order_gene][ii_na] )
  names(roll_medians_length) <- rownames(exprSet)[order_gene]
  roll_medians_length <- roll_medians_length[  match(rownames(exprSet), names(roll_medians_length) ) ]
  stopifnot( all.equal(names(roll_medians_length), rownames(exprSet) ) )
  log10cv2_adj <-  cv -  roll_medians_length

  if(F){
    pheatmap::pheatmap(log2(exprSet[names(head(sort(log10cv2_adj),50)),]+1))
    pheatmap::pheatmap(log2(exprSet[names(tail(sort(log10cv2_adj),50)),]+1))
  }
  return(log10cv2_adj)
}
