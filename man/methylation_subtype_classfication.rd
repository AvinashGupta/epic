\name{methylation_subtype_classfication}
\alias{methylation_subtype_classfication}
\title{
Classify subtypes by methylation data
}
\description{
Classify subtypes by methylation data
}
\usage{
methylation_subtype_classfication(gr, n_class, pct_cutoff, corr_cutoff,
    k, ha = NULL, cgi = NULL, shore = NULL)
}
\arguments{

  \item{gr}{genomic regions which contains mean methylation, should be generated by \code{\link{get_mean_methylation_in_genomic_features}}}
  \item{n_class}{number of classes expected}
  \item{pct_cutoff}{percent of most variable rows}
  \item{corr_cutoff}{cutoff for absolute correlation}
  \item{k}{number of correlated windows}
  \item{ha}{additional annotation}
  \item{cgi}{cpg island}
  \item{shore}{cgi shores}

}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
