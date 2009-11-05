\name{plotCuspResidfitted}
\alias{plotCuspResidfitted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Residuals against Fitted Plot for Cusp Model Fit }
\description{
  Plot Residuals against Fitted Values for a Cusp Model Fit.
}
\usage{
plotCuspResidfitted(object, caption = "Residual vs Fitted", 
    xlab = paste("Fitted (", colnames(fitted(object))[1], " convention)", sep = ""), 
    ylab = "Residual", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ cusp fit object returned by \code{cusp} }
  \item{caption}{ plot caption }
  \item{xlab}{ label for x-axis }
  \item{ylab}{ label for y-axis }
  \item{\dots}{ named arguments that are passed to \code{\link{plot}} }
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
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
