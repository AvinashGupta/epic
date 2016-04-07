\name{run}
\alias{run}
\title{
run pipeline
}
\description{
run pipeline
}
\usage{
run()
}
\details{
\preformatted{
 Usage: Rscript -e "epic::run()" cmd [options]

 Available cmd:

     chromatin_states_transitions
     correlated_enriched
     correlated_regions
     correlated_regions_downstream
     correlated_regions_filter
     correlated_regions_gviz
     correlated_regions_reduce
     differential_methylation_in_cgi_and_shore
     differential_methylation_in_genomic_features
     general_methylation_distribution
     methylation_subtype_classification_in_cgi_and_shore  }

For each cmd, use \code{\link[Rscript -e "epic]{run()" cmd --help}} to get help
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
