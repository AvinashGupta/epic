\name{percentOverlaps}
\alias{percentOverlaps}
\title{
For every interval in \code{query}, it calculates the percent that is covered by \code{subject}
}
\description{
For every interval in \code{query}, it calculates the percent that is covered by \code{subject}
}
\usage{
percentOverlaps(query, subject, ...)
}
\arguments{

  \item{query}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{subject}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{...}{pass to \code{\link[GenomicRanges]{findOverlaps}}}

}
\value{
a numeric vector which is same as the length of \code{query}.

Be careful with \code{strand} in your \code{\link[GenomicRanges]{GRanges}} object!!
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
gr1 = GRanges(seqname = "chr1", ranges = IRanges(start = c(4, 10), end = c(6, 16)))
gr2 = GRanges(seqname = "chr1", ranges = IRanges(start = c(7, 13), end = c(8, 20)))
percentOverlaps(gr1, gr2)
percentOverlaps(gr2, gr1)
}
