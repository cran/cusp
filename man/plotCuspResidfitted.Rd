\name{plotCuspResidfitted}
\alias{plotCuspResidfitted}
\title{ Residuals against Fitted Plot for Cusp Model Fit }
\description{
  Plot Residuals against Fitted Values for a Cusp Model Fit.
}
\usage{
plotCuspResidfitted(object, caption = "Residual vs Fitted", 
    xlab = paste("Fitted (", colnames(fitted(object))[1], " convention)", sep = ""), 
    ylab = "Residual", ...)
}
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
\value{
No return value. Called for its side effect.
}
%\references{ 
%}
\author{ Raoul Grasman }
%\note{ further instruction
% Use sections like Warning with \section{Warning }{....} to organize
%}
\seealso{ \code{\link{plot.cusp}} }
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }% __ONLY ONE__ keyword per line
