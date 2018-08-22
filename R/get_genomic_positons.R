#' Get the genomic positions based on a list of gene
#'
#' A list of gene, should be ensembl IDs or gene symbols
#' One should set the species as human or mouse.
#' The genomic positions is human_geneInfo_genecode_v25 or mouse_geneInfo_genecode_vM12
#'
#' @param geneList Should be the rownames of expression matrix
#' @param geneType Choose ensembl or symbol,defaults: ensembl
#' @param species Choose human or mouse,defaults: human
#'
#' @return pos the genomic positions, columns should be : gene/chr/start/end
#' @examples
#' get_genomic_positions
#' get_genomic_positions(rownames(exprSet),geneType='ensembl',species='human')

get_genomic_positions <- function(geneList,geneType='ensembl',species='human'){

  if(species=='human'){
    head(human_geneInfo_genecode_v25)
    if(geneType=='ensembl'){
      pos=human_geneInfo_genecode_v25
      pos=pos[match(geneList,pos$ensembl),c(4,1:3)]
    }
    if(geneType=='symbol'){
      pos=human_geneInfo_genecode_v25
      pos=pos[match(geneList,pos$symbol),c(6,1:3)]
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
      pos=pos[match(geneList,pos$ensembl),c(4,1:3)]
    }
    if(geneType=='symbol'){
      pos=mouse_geneInfo_genecode_vM12
      pos=pos[match(geneList,pos$symbol),c(5,1:3)]
    }
    new_chr=gsub('chr','',pos$chr)
    table(new_chr)
    new_chr[new_chr=='X']=20
    new_chr[new_chr=='Y']=21
    new_chr=as.numeric(new_chr)
    pos$chr=new_chr
    pos=pos[order(pos$chr,pos$start),]
  }
  return(pos)
}
