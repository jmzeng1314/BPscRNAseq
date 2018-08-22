options(stringsAsFactors = F)

if(F){
  ls('package:BPscRNAseq')
}

## for human gencode.v25.annotation.gtf
if(F){

  a=read.table('data/human.gene.positions')[,c(2:4,1,6,7)]
  colnames(a)=c('chr','start','end','ensembl','type','symbol')
  length(unique(a$symbol))
  length(unique(a$ensembl))
  head(a)
  human_geneInfo_genecode_v25=a
  devtools::use_data(human_geneInfo_genecode_v25, overwrite = T)

  a=read.table('data/human_ENSG_length')
  colnames(a)=c( 'ensembl','length' )
  head(a)
  human_geneLength_genecode_v25=a
  devtools::use_data(human_geneLength_genecode_v25, overwrite = T)
}


## for mouse gencode.vM12.annotation.gtf.gz
if(F){
  options(stringsAsFactors = F)
  a=read.table('data/mouse.gene.positions')[,c(2:4,1,7,6)]
  colnames(a)=c('chr','start','end','ensembl','symbol','type')
  length(unique(a$symbol))
  length(unique(a$ensembl))
  head(a)
  mouse_geneInfo_genecode_vM12=a
  devtools::use_data(mouse_geneInfo_genecode_vM12, overwrite = T)

  a=read.table('data/mouse_ENSG_length')
  colnames(a)=c( 'ensembl','length' )
  head(a)
  mouse_geneLength_genecode_vM12=a
  devtools::use_data(mouse_geneLength_genecode_vM12, overwrite = T)

}

## The orthologous genes between human and mouse

if(F){
  a=read.csv('data/human2mouse.csv',header = F )
  table(a[,1] == toupper(a[,2]))
  colnames(a)=c('human','mouse')
  rmGenes=apply(a,1,function(x) sum(x=='') >0)
  a=a[!rmGenes,]
  table(a[,1] == toupper(a[,2]))
  human2mouse_symbols=a
  devtools::use_data(human2mouse_symbols, overwrite = T)
}

if(F){
  library(BPscRNAseq)
  library(airway)
  library(edgeR)
  library(DESeq2)
  data(airway)
  airway
  counts=assay(airway)
  counts[1:4,1:4];dim(counts)
  exprSet=counts2exprSet(counts)
  exprSet[1:4,1:4];dim(exprSet)
  airway_exprSet =exprSet
  devtools::use_data(airway_exprSet, overwrite = T)
}









