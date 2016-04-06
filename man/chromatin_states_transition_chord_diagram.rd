\name{chromatin_states_transition_chord_diagram}
\alias{chromatin_states_transition_chord_diagram}
\title{
Chord diagram for chromatin states transistion
}
\description{
Chord diagram for chromatin states transistion
}
\usage{
chromatin_states_transition_chord_diagram(mat, max_mat = mat, remove_unchanged_transition = FALSE,
    state_col = NULL, ...)
}
\arguments{

  \item{mat}{the transition matrix}
  \item{max_mat}{if there are several matrix to be compared, set it to the matrix with maximum sum}
  \item{remove_unchanged_transition}{whether to remove regions that states are not changed}
  \item{state_col}{color for grids. It should be a vector of which names correspond to states}
  \item{...}{pass to \code{\link[circlize]{chordDiagram}}}

}
\details{
Rows of \code{mat} always locates at the bottom.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
