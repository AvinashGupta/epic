\name{logical_segment}
\alias{logical_segment}
\title{
Segmentation by continuous logical values
}
\description{
Segmentation by continuous logical values
}
\usage{
logical_segment(l)
}
\arguments{

  \item{l}{a logical vector}

}
\details{
the logical vector will be segmented according to their values.
It returns intervals for continuous \code{\link{TRUE}} values
}
\value{
a data frame in which the first column is the index of start sites
the second column is the index of end sites.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
l = c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE)
logical_segment(l)
}
