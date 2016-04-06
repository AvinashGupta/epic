\name{genomic_corr_pintersect}
\alias{genomic_corr_pintersect}
\title{
Intersection between two sets of genomic regions
}
\description{
Intersection between two sets of genomic regions
}
\usage{
genomic_corr_pintersect(query, reference, ...)
}
\arguments{

  \item{query}{genomic region 1, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{reference}{genomic region 2, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{...}{pass to \code{\link{percentOverlaps}}}

}
\details{
For each region in \code{query}, it calculates the percent that is covered by \code{reference}.

The returned value is percent which is how much \code{query} is covered by \code{reference} (by default).

Be careful with the \code{strand} in your GRanges object!!
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
