# I hope It will be the best practice for scRNA seq data(downstream analysis)

### Upstream workflow 

I prefer `STAR+FeatureCounts` to generate the  raw counts expression matrix .

One can also download it from published paper, such as : [Cell Rep.](https://www.ncbi.nlm.nih.gov/pubmed/29091775#) 2017 [Single-Cell RNA-Seq Analysis of Infiltrating Neoplastic Cells at the Migrating Front of Human Glioblastoma.](https://www.ncbi.nlm.nih.gov/pubmed/29091775)  

- The raw sequence data can be found at GEO:   [GSE84465](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465)   
- The expression matrix can be dowload from:  [3,589 cells in a cohort of four patients  ](http://gbmseq.org/) 

### Create input files for inferCNV 

Once we generate the raw counts expression matrix based on scRNA-seq data, such as `'GBM_raw_gene_counts.csv'` , we can use the code below to create  input files for inferCNV

```R
options(stringsAsFactors = F)
dir='/Users/jmzeng/biosoft/scRNA_cnv/project/gbm_2017/'
## download from http://gbmseq.org/#downloadData 
counts=read.table(file.path(dir,'GBM_raw_gene_counts.csv'))
#counts=counts[,1:100]
counts[1:4,1:4];dim(counts)
library(BPscRNAseq)
exprSet=counts2exprSet(counts)
exprSet[1:4,1:4];dim(exprSet)
exprSet2inferCNV(exprSet,geneType='symbol',species='human',prefix='gbm_2017',dir=dir)
```

Then we can run inferCNV as https://github.com/broadinstitute/inferCNV/wiki 

```Shell
Rscript ~/biosoft/scRNA_cnv/inferCNV/scripts/inferCNV.R  --output_dir   test  gbm_inferCNV_exprSet.txt gbm_inferCNV_pos.txt
```

There's also a more convenient way to get the example input files for inferCNV, by just run `airway2inferCNV(TRUE)`   

### Calculate CNV values and draw heatmap 

If you don't want use the inferCNV based on the scripts from Broad institute, you can also calculate CNV in R and draw heatmap as below:

```R
library(airway)
library(edgeR)
library(DESeq2)
data(airway)
airway
counts=assay(airway)

counts[1:4,1:4];dim(counts)
exprSet=counts2exprSet(counts)
exprSet[1:4,1:4];dim(exprSet)
cnv=exprSet2CNV(exprSet,geneType='ensembl',species='human') 
cnv[1:4,1:8];dim(cnv)
ht_cnv(cnv)
```



### Compare the results between inferCNV and my function



```R
rm(list=ls())
dir='/Users/jmzeng/biosoft/scRNA_cnv/project/gbm_2014'
load(file.path(dir,'GBM_for_CNV_input.Rdata'))
exprSet[1:4,1:4];dim(exprSet)
cnv=exprSet2CNV(exprSet,geneType='symbol',species='human') 
cnv[1:4,1:8];dim(cnv)
ht_cnv(cnv,prefix = 'gbm_2014_jimmy')
exprSet2inferCNV(exprSet,geneType='symbol',species='human',prefix='gbm_2014',dir=dir)
# Rscript ~/biosoft/scRNA_cnv/inferCNV/scripts/inferCNV.R  --output_dir   test  
# gbm_2014_inferCNV_exprSet.txt gbm_2014_inferCNV_pos.txt
```



### Calculate the DM (distance to median)  values 

```R
library(BPscRNAseq)
library(airway)
library(edgeR)
library(DESeq2)
data(airway)
airway
counts=assay(airway)
counts[1:4,1:4];dim(counts)
geneLists=rownames(counts)
# removed lowly expressed genes whose mean normalised read counts (reads per million) are less than10, 
# since we cannot distinguish biological noise from technical noise for these genes.
keepGene=rowMeans(cpm(counts) ) >=10
table(keepGene);dim(counts)
dim(counts[keepGene,])
exprSet=counts[keepGene,]
rownames(exprSet)=geneLists[keepGene]
exprSet=cpm(exprSet)
exprSet[1:4,1:4];dim(exprSet)
DM=cal_DM(exprSet,geneType='ensembl',species='human')
pheatmap::pheatmap(log2(exprSet[names(head(sort(DM),50)),]+1))
pheatmap::pheatmap(log2(exprSet[names(tail(sort(DM),50)),]+1))

```

