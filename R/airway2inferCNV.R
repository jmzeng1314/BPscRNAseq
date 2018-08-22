#' Create two input files for inferCNV for airway data
#'
#' Airway is the expression matrix  from public bulk RNA-seq data.
#' Check the document for input format at : https://github.com/broadinstitute/inferCNV/wiki
#' The genomic positions for airway data is : human_geneInfo_genecode_v25
#'
#' @param run Choose TRUE or FALSE,defaults: FALSE
#' @param dir Choose where to put the files,defaults: ./
#'
#' @return Two files \code{airwary_inferCNV_pos.txt} and \code{airwary_inferCNV_exprSet.txt}
#' @examples
#' airway2inferCNV
#' airway2inferCNV(TRUE)
#' airway2inferCNV(TRUE, '~/biosoft/scRNA_cnv/project/airway')

airway2inferCNV <- function(run=FALSE,dir='./'){

  if(run){
    library(airway)
    library(edgeR)
    library(DESeq2)
    data(airway)
    airway
    counts=assay(airway)
    counts[1:4,1:4];dim(counts)
    exprSet=counts2exprSet(counts)
    exprSet[1:4,1:4];dim(exprSet)
    exprSet2inferCNV(exprSet,geneType='ensembl',species='human',prefix='airwary',dir=dir)

    # Rscript ~/biosoft/scRNA_cnv/inferCNV/scripts/inferCNV.R    \
    # --output_dir  test airwary_inferCNV_exprSet.txt airwary_inferCNV_pos.txt

  }

}

