\name{plot.cusp}
\alias{plot.cusp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Graphical Diagnostic Display of Cusp Catastrophe Data Fit }
\description{
  This function generates diagnostic graphical displays of fits of a cusp catastrophe model to data obtained with \code{\link{cusp}}
}
\usage{
\method{plot}{cusp}(x, what = c("all", "bifurcation", "residual", "densities"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Object returned by \code{\link{cusp}} }
  \item{what}{ 1-character string giving the type of plot desired. The following values are possible: \code{"all"} for a panel plot with all diagnostic plots, \code{"bifurcation"} for a plot of the bifurcation surface with estimated control parameter locations superimposed, "residual" for a plot of the residuals against fitted values, "densities" for a plot of density estimates conditioned on the estimated location on the bifurcation surface. }
  \item{\dots}{ named arguments that are passed to lower level plotting function }
}
\details{
  These diagnostic plots help to identify problems with the fitted model. In optimal
  cases the fitted locations in the parameter plane are dispersed over regions of
  qualitatively different behavior. Within each region the fitted dependent
  values have a density of the appropriate shape (e.g., bimodal in the bifurcation
  set).
}
\value{
No return value. Called for its side effect.
}
\references{ See \code{\link{cusp-package}} }
\author{ Raoul Grasman }
%\note{ Any notes to user?
%
% You can organize into other sections like Warning with \section{Warning }{....} 
%}
\seealso{ \code{\link{plotCuspBifurcation}}, \code{\link{plotCuspResidfitted}}, \code{\link{plotCuspDensities}}  }
\examples{
set.seed(20)
x1 = runif(150)
x2 = runif(150)
z = Vectorize(rcusp)(1, 4*x1-2, 4*x2-1)
data <- data.frame(x1, x2, z)
fit <- cusp(y ~ z, alpha ~ x1+x2, beta ~ x1+x2, data)
\dontrun{
plot(fit)

# just densities
layout(matrix(1:4,2))
plot(fit, what="densities")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ models }% __ONLY ONE__ keyword per line
