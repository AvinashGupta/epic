\name{set_counter}
\alias{set_counter}
\title{
Set a counter for a loop
}
\description{
Set a counter for a loop
}
\usage{
set_counter(n, fmt = "\%s")
}
\arguments{

  \item{n}{maximum number of looping}
  \item{fmt}{format of message}

}
\value{
a function which should be called in the loop
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
counter = set_counter(1000)
for(i in 1:1000) {counter()}
counter = set_counter(1000, fmt = "processing \%s")
for(i in 1:1000) {counter()}
}
