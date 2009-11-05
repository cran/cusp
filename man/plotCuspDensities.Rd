\name{plotCuspDensities}
\alias{plotCuspDensities}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot Cusp State Variable Densities Conditioned on Control Parameter Values}
\description{
  Plot density of state variables conditioned on their location on the cusp control surface.
}
\usage{
plotCuspDensities(object, main = "Conditional density", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ cusp fit object returned by \code{cusp} }
  \item{main}{ title of plot }
  \item{\dots}{ named arguments that are passed to \code{plot} and \code{draw.cusp.bifset} }
}
\details{
  This function is mainly intended for internal use by \code{plot.cusp}.
}
%\value{
%  ~Describe the value returned
%  If it is a LIST, use
%  \item{comp1 }{Description of 'comp1'}
%  \item{comp2 }{Description of 'comp2'}
%  ...
%}
%\references{ ~put references to the literature/web site here ~ }
\author{ Raoul Grasman }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{plot.cusp}} }
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
%\keyword{ models }% __ONLY ONE__ keyword per line
