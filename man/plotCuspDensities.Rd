\name{plotCuspDensities}
\alias{plotCuspDensities}
\title{ Plot Cusp State Variable Densities Conditioned on Control Parameter Values}
\description{
  Plot density of state variables conditioned on their location on the cusp control surface.
}
\usage{
plotCuspDensities(object, main = "Conditional density", ...)
}
\arguments{
  \item{object}{ cusp fit object returned by \code{cusp} }
  \item{main}{ title of plot }
  \item{\dots}{ named arguments that are passed to \code{plot} and \code{draw.cusp.bifset} }
}
\details{
  This function is mainly intended for internal use by \code{plot.cusp}.
}
\value{
No return value. Called for its side effect.
}
\author{ Raoul Grasman }
%\note{ Any further notes to the user?
%
% Organize in sections like Warning with \section{Warning }{....} 
%}
\seealso{ \code{\link{plot.cusp}} }
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
%\keyword{ models }% __ONLY ONE__ keyword per line
