\name{common_regions}
\alias{common_regions}
\title{
Find common genomic regions across several samples
}
\description{
Find common genomic regions across several samples
}
\usage{
common_regions(gr_list, min_width = 0, min_coverage = 1, gap = bp(0))
}
\arguments{

  \item{gr_list}{a list of \code{\link[GenomicRanges]{GRanges}}, better have names}
  \item{min_width}{minimum width of the common regions}
  \item{min_coverage}{minimum cross-sample coverage for the common regions}
  \item{gap}{gap to merge common regions, pass to \code{\link{reduce2}}}

}
\details{
A common region is defined a region which at least i samples overlap with.
The output can be sent to \code{\link{subgroup_specific_genomic_regions}} to find subgroup
specific regions.
}
\value{
a GRanges object contains coordinates of common regions. The columns in meta data
are percent of the common region which is covered by regions in every sample.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
