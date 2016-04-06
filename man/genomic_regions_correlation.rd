\name{genomic_regions_correlation}
\alias{genomic_regions_correlation}
\title{
Correlation between two sets of genomic regions
}
\description{
Correlation between two sets of genomic regions
}
\usage{
genomic_regions_correlation(gr_list_1, gr_list_2, background = NULL,
    chromosome = paste0("chr", 1:22), species = "hg19",
    nperm = 0, mc.cores = 1, stat_fun = genomic_corr_jaccard, ...,
    bedtools_binary = Sys.which("bedtools"), tmpdir = tempdir())
}
\arguments{

  \item{gr_list_1}{a list of \code{\link[GenomicRanges]{GRanges}} objects, should be a named list.}
  \item{gr_list_2}{a list of \code{\link[GenomicRanges]{GRanges}} objects, should be a named list.}
  \item{background}{a \code{\link[GenomicRanges]{GRanges}} object, the genomic background to be restricted in}
  \item{chromosome}{chromosomes}
  \item{species}{species, used for random shuffling genomic regions}
  \item{nperm}{number of permutations. If it is set to 0, no permutation will be performed.}
  \item{mc.cores}{number of cores for parallel calculation}
  \item{stat_fun}{method to calculate correlations. There are some pre-defined functions: \code{\link{genomic_corr_reldist}}, \code{\link{genomic_corr_absdist}}, \code{\link{genomic_corr_jaccard}}, \code{\link{genomic_corr_nintersect}}, \code{\link{genomic_corr_pintersect}}, \code{\link{genomic_corr_sintersect}}. The self-defined function should accept at least two arguments which are two GRanges object. The third argument is \code{...} which is passed from the main function. The function should only return a numeric value.}
  \item{...}{pass to \code{stat_fun}}
  \item{bedtools_binary}{file for __bedtools__}
  \item{tmpdir}{tempoary dir}

}
\details{
The correlation between two sets of genomic regions basicly means how much the first type of genomic regions
are overlapped or close to the second type of genomic regions.

The significance of the correlation is calculated by random shuffling the regions. 
In random shuffling, regions in \code{gr_list_1} will be shuffled. So if you want to shuffle \code{gr_list_2},
just switch the first two arguments.

Pleast note random shuffling is done by _bedtools_, so _bedtools_ should be installed and exists in \code{PATH}
and should support \code{-i -g -incl} options.
}
\value{
A list containing:

\describe{
  \item{foldChange}{stat/E(stat), stat divided by expected value which is generated from random shuffling}
  \item{p.value}{p-value for over correlated. So, 1 - p.value is the significance for being no correlation}
  \item{stat}{statistic value}
  \item{stat_random_mean}{mean value of stat in random shuffling}
  \item{stat_random_sd}{standard deviation in random shuffling}
}
}
\seealso{
\code{\link{genomic_corr_reldist}}, \code{\link{genomic_corr_jaccard}}, \code{\link{genomic_corr_absdist}}, \code{\link{genomic_corr_nintersect}}, \code{\link{genomic_corr_pintersect}}, \code{\link{genomic_corr_sintersect}}
}
\examples{
# There is no example
NULL

}
