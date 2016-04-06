\name{cor_cols}
\alias{cor_cols}
\title{
pair-wise correlation of rows in a matrix with huge number of columns
}
\description{
pair-wise correlation of rows in a matrix with huge number of columns
}
\usage{
cor_cols(x, abs_cutoff = 0.5, size = 1000, mc = 1, ...)
}
\arguments{

  \item{x}{a matrix, correlation is calculated by columns}
  \item{abs_cutoff}{cutoff of absolute correlation}
  \item{size}{size of blocks}
  \item{mc}{multiple cores}
  \item{...}{pass to \code{\link[stats]{cor}}}

}
\details{
The code is adapted from \url{https://rmazing.wordpress.com/2013/02/22/bigcor-large-correlation-matrices-in-r/}
}
\value{
a vector that represents how many other columns correlate to current column under the correlation cutoff
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
