\name{cr_hilbert}
\alias{cr_hilbert}
\title{
Hilbert curve visualization of correlated regions
}
\description{
Hilbert curve visualization of correlated regions
}
\usage{
cr_hilbert(cr, template, txdb, chromosome = paste0("chr", 1:22), merge_chr = TRUE)
}
\arguments{

  \item{cr}{filtered correlated regions from \code{\link{filter_correlated_regions}}}
  \item{template}{template path of cr objects, chromosome name should be marked by \code{@{chr}}}
  \item{txdb}{a \code{\link[GenomicFeatures]{TxDb}} object}
  \item{chromosome}{a vector of chromosome names}
  \item{merge_chr}{whether put chromsomes in one single plot}

}
\details{
It can both visualize unfiltered correlated regions (generated by \code{\link{correlated_regions}}) and significant correlated
regions (generated by \code{\link{filter_correlated_regions}}). To visualize unfiltered cr, \code{template} should be specified and to
visualize filtered correlated regions, \code{cr} should be specified.
}
\value{
No value is returned
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
