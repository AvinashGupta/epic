\name{methylation_subtype_classfication}
\alias{methylation_subtype_classfication}
\title{
Classify subtypes by methylation data
}
\description{
Classify subtypes by methylation data
}
\usage{
methylation_subtype_classfication(gr, sample_id, chromosome, p_cutoff, corr_cutoff, k, cgi, n_class, ha = NULL)
}
\arguments{

  \item{gr}{genomic regions}
  \item{sample_id}{sample id}
  \item{chromosome}{chromosomes}
  \item{p_cutoff}{cutoff for p-values}
  \item{corr_cutoff}{cutoff for absolute correlation}
  \item{k}{number}
  \item{cgi}{CpG islands}
  \item{n_class}{number of classes}
  \item{ha}{additional annotation}

}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
