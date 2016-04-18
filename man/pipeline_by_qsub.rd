\name{pipeline_by_qsub}
\alias{pipeline_by_qsub}
\title{
run pipeline through qsub system
}
\description{
run pipeline through qsub system
}
\usage{
pipeline_by_qsub(config_file, prefix = "", email = NULL, enforce = FALSE, Rscript_binary = "Rscript")
}
\arguments{

  \item{config_file}{path of configuration script}
  \item{prefix}{prefix of the job name}
  \item{email}{email}
  \item{enforce}{enforce run all the steps}
  \item{Rscript_binary}{path of Rscript binary}

}
\details{
Automatically run
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
