% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/airway2inferCNV.R
\name{airway2inferCNV}
\alias{airway2inferCNV}
\title{Create two input files for inferCNV for airway data}
\usage{
airway2inferCNV(run = FALSE, dir = "./")
}
\arguments{
\item{run}{Choose TRUE or FALSE,defaults: FALSE}

\item{dir}{Choose where to put the files,defaults: ./}
}
\value{
Two files \code{airwary_inferCNV_pos.txt} and \code{airwary_inferCNV_exprSet.txt}
}
\description{
Airway is the expression matrix  from public bulk RNA-seq data.
Check the document for input format at : https://github.com/broadinstitute/inferCNV/wiki
The genomic positions for airway data is : human_geneInfo_genecode_v25
}
\examples{
airway2inferCNV
airway2inferCNV(TRUE)
airway2inferCNV(TRUE, '~/biosoft/scRNA_cnv/project/airway')
}
