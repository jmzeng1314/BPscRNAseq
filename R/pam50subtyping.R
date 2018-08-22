#' Pam50 subtyping based on normalized expression matrix
#'
#' expression matrix should be numeric values matrix in which  each row is a gene, and each column is a cell or sample.
#' The rownames for this expression matrix  should be ensembl IDs or gene symbols
#' One should set the species as human or mouse.
#'
#' @param exprSet The expression matrix(which shoud be normalized,like log2(cpm+1) or log2(tpm+1))
#' @param geneType Choose ensembl or symbol,defaults: ensembl
#' @param species Choose human or mouse,defaults: human
#'
#' @return all_cnv molecular.subtyping results
#' @examples
#' pam50subtyping
#' pam50subtyping(exprSet,geneType='ensembl',species='human')

pam50subtyping <- function(exprSet,geneType='ensembl',species='human' ){
  suppressPackageStartupMessages(library(genefu))
  data(pam50)
  pam50genes=pam50$centroids.map[c(1,3)]
  pam50genes[pam50genes$probe=='CDCA1',1]='NUF2'
  pam50genes[pam50genes$probe=='KNTC2',1]='NDC80'
  pam50genes[pam50genes$probe=='ORC6L',1]='ORC6'
  rownames(pam50genes)=pam50genes$probe



  # CDCA1 -->  NUF2	NUF2, NDC80 Kinetochore Complex Component
  # KNTC2 --> NDC80
  # ORC6L --> ORC6 Origin Recognition Complex Subunit 6
  if(species=='human'){
    head(human_geneInfo_genecode_v25)

    if(geneType=='ensembl'){
      pos=human_geneInfo_genecode_v25
      pam50genes$ensembl=pos[match(rownames(pam50genes),pos$symbol),'ensembl']
      exprSet=exprSet[rownames(exprSet) %in% pam50genes$ensembl,]
      ddata=t(exprSet)
      dannot=pam50genes[match(colnames(ddata),pam50genes$ensembl),]
      dannot$probe=dannot$ensembl
      rownames(dannot)=dannot$probe

    }
    if(geneType=='symbol'){
      exprSet=exprSet[rownames(exprSet) %in% pam50genes$probe,]
      ddata=t(exprSet)
      dannot=pam50genes[match(colnames(ddata),pam50genes$probe),]
    }

  }
  if(species=='mouse'){
    pam50genes$probe=human2mouse_symbols[match(pam50genes$probe , human2mouse_symbols[,1]),2]
    rownames(pam50genes)=pam50genes$probe

    head(mouse_geneInfo_genecode_vM12)
    if(geneType=='ensembl'){
      pos=mouse_geneInfo_genecode_vM12
      pam50genes$ensembl=pos[match(rownames(pam50genes),pos$symbol),'ensembl']
      exprSet=exprSet[rownames(exprSet) %in% pam50genes$ensembl,]
      ddata=t(exprSet)
      dannot=pam50genes[match(colnames(ddata),pam50genes$ensembl),]
      dannot$probe=dannot$ensembl
      rownames(dannot)=dannot$probe
    }
    if(geneType=='symbol'){
      exprSet=exprSet[rownames(exprSet) %in% pam50genes$probe,]
      ddata=t(exprSet)
      dannot=pam50genes[match(colnames(ddata),pam50genes$probe),]
    }
  }

  message(paste0(nrow(dannot), ' of 50 genes are used to subtype'))
  PAM50Preds<-molecular.subtyping(sbt.model = "pam50",data=ddata,
                                  annot=dannot,do.mapping=TRUE)
  table(PAM50Preds$subtype)
  return(PAM50Preds)

}

