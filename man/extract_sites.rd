\name{extract_sites}
\alias{extract_sites}
\title{
Extract subset of sites which are in a set of intervals
}
\description{
Extract subset of sites which are in a set of intervals
}
\usage{
extract_sites(start, end, site, return_index = FALSE, min_sites = 0)
}
\arguments{

  \item{start}{start position, a vector}
  \item{end}{end position, a vector. Note there should be no overlap between all \code{[start, end]} (You may use \code{\link[IRanges]{reduce}} to merge the overlapping intervals.)}
  \item{site}{positions of all sites, should be sorted.}
  \item{return_index}{whether return the index in the position vector or just the position itself?}
  \item{min_sites}{minimal number of sites in aninterval}

}
\details{
Providing a huge vector of genomic positions, we want to extract subset of sites which
locate in a specific region (e.g. extract CpG sites in DMRs). Normally, we will use:

  \preformatted{
	site = sort(sample(10000000, 1000000))
	start = 123456
	end = 654321
	subsite = site[site >= start & site <= end]  }

Unfortunately, in above code, the whole vector \code{site} will be scaned for four times
(\code{>=}, \code{<=}, \code{&} and \code{[}).
If you want to look for sites in more than one regions (e.g. 1000 regions), in every
loop, the whole \code{site} vector will be re-scanned again and again which is very time-consuming.

Here we have \code{\link{extract_sites}} function which uses binary search to do subsetting.
Of course, \code{site} should be sorted non-decreasing beforehand.

  \preformatted{
	subsite = extract_sites(start, end, site, index = FALSE)  }

Not only for single interval, you can also extract sites in a list of genomic regins,
by setting \code{start} and \code{end} as a vector.

  \preformatted{
	start = c(123456, 234567, 345678)
	end = c(133456, 244567, 355678)
	subsite = extract_sites(start, end, site)  }

You can choose to return index only or positions.

  \preformatted{
	subsite = extract_sites(start, end, site, return_index = FALSE)
	head(subsite)
	subsite_index = extract_sites(start, end, site, return_index = TRUE)
	head(subsite_index)
	head(site[subsite_index])  }
}
\value{
a vector of positions or index.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
