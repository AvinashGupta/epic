\name{make_transition_matrix_from_chromHMM}
\alias{make_transition_matrix_from_chromHMM}
\title{
Generate transition matrix from chromHMM results
}
\description{
Generate transition matrix from chromHMM results
}
\usage{
make_transition_matrix_from_chromHMM(gr_list_1, gr_list_2, window = NULL,
    min_1 = round(length(gr_list_1)/2), min_2 = round(length(gr_list_2)/3))
}
\arguments{

  \item{gr_list_1}{a list of \code{\link[GenomicRanges]{GRanges}} objects which contain chromatin states in group 1. The first column in meta columns should be the states. please note the start position when converting bed format to \code{GRanges} format (0-based and 1-based).}
  \item{gr_list_2}{a list of \code{\link[GenomicRanges]{GRanges}} objects which contains chromatin states in group 2.}
  \item{window}{window size which was used to do chromHMM states prediction. If it is not specified, the greatest common divisor of all region width is used.}
  \item{min_1}{minimal recurrency in \code{gr_list_1}}
  \item{min_2}{minimal recurrency in \code{gr_list_2}}

}
\value{
A transition matrix in which values represent width of regions that transite from one state to the other.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
