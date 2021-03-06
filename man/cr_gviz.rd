\name{cr_gviz}
\alias{cr_gviz}
\title{
Customized Gviz plot for a gene model
}
\description{
Customized Gviz plot for a gene model
}
\usage{
cr_gviz(cr, gi, expr, txdb, gene_start = NULL, gene_end = NULL,
    species = "hg19", gf_list = NULL, hm_list = NULL, symbol = NULL)
}
\arguments{

  \item{cr}{correlated regions generated by \code{\link{filter_correlated_regions}}}
  \item{gi}{gene id}
  \item{expr}{the expression matrix which is same as in \code{\link{correlated_regions}}}
  \item{txdb}{a \code{GenomicFeatures::GRanges} object.}
  \item{gene_start}{start position of gene}
  \item{gene_end}{end position of the gene}
  \item{species}{species}
  \item{gf_list}{a list of \code{\link[GenomicRanges]{GRanges}} objects which contains additional annotations}
  \item{hm_list}{a list of \code{\link[GenomicRanges]{GRanges}}}
  \item{symbol}{symbol of the gene}

}
\details{
Several information on the correlated regions in an extended gene model are visualized by Gviz package:

\itemize{
  \item gene models. Multiple transcripts will also be plotted.
  \item correlation for every CpG window
  \item heatmap for methylation
  \item significant correlated regions
  \item CpG density
  \item annotation to other genomic features if provided
  \item annotation to other ChIP-Seq peak data if provided
}
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
