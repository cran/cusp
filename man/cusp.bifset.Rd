\name{cusp.bifset}
\alias{cusp.bifset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ compute normal/symmetry factor borders of bifurcation set of cusp catastrophe }
\description{
  Given bifurcation/splitting factor values this function computes the border values of the normal/symmetry factor for the bifurcation set of the cusp catastrophe.
}
\usage{
cusp.bifset(beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{ values of the bifurcation/splitting factor at which the border values of the normal/symmetry factor is computed }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  %~Describe the value returned
  Matrix with columns named \code{beta}, \code{alpha.l}, \code{alpha.u}. The latter two columns give respectively the lower and upper border values of the normal/symmetry factor. Negative values of \code{beta} give \code{NaN} values for the normal factor.
}
\references{ See \code{\link{cusp-package}} }
\author{ Raoul Grasman }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{cusp-package}} }
\examples{
    cusp.bifset(-3:3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
\keyword{ models }% __ONLY ONE__ keyword per line
