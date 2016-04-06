\name{available_gencode_field}
\alias{available_gencode_field}
\title{
Returns all supported fields in gtf data
}
\description{
Returns all supported fields in gtf data
}
\usage{
available_gencode_field(file, level = "gene")
}
\arguments{

  \item{file}{the input gtf file}
  \item{level}{level of the annotation (e.g. gene, transcript, exon, ...)}

}
\details{
These fields are stored in the 9th column in the gtf file.

This function only works under Linux-like OS.
}
\value{
A vector of available fields
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
