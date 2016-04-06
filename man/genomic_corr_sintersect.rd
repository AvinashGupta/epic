\name{genomic_corr_sintersect}
\alias{genomic_corr_sintersect}
\title{
Intersection between two sets of genomic regions
}
\description{
Intersection between two sets of genomic regions
}
\usage{
genomic_corr_sintersect(query, reference, background = NULL)
}
\arguments{

  \item{query}{genomic region 1, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{reference}{genomic region 2, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{background}{subset of sites that should be only looked into, a \code{\link[GenomicRanges]{GRanges}} object}

}
\details{
It calculates the total length of overlapped regions in \code{query}.

If the interest is e.g. the number of CpG sites both in \code{query} and in \code{reference}
\code{background} can be set with a GRanges object which contains positions of CpG sites.

Be careful with the \code{strand} in your GRanges object!!
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
