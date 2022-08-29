\name{plotCuspBifurcation}
\alias{plotCuspBifurcation}
\title{ Display Fitted Data on Control Plane of Cusp Catastrophe. }
\description{
  Displays fitted data points on the control plane of cusp catastrophe. The function takes a fit object obtained with \code{cusp} and generates a plot. Different diagnostic plots may be chosen, or all can be combined in a single plot (the default).
}
\usage{
plotCuspBifurcation(object, xlim = a + c(-0.3, 0.3), ylim = b + c(-0.1,
                 0.1), xlab = expression(alpha), ylab =
                 expression(beta), hue = 0.5 + 0.25 * tanh(object$y),
                 col = hsv(h = hue, s = 1, alpha = 0.4), cex.xlab =
                 1.55, cex.ylab = cex.xlab, axes = TRUE, box = TRUE,
                 add = FALSE, bifurcation.set.fill = gray(0.8),
                 cex.scale = 15, cex = (cex.scale/log(NROW(ab))) *
                 dens/max(dens), pch = 20)}

\arguments{
  \item{object}{ object returned by \code{cusp} }
  \item{xlim}{ the x limits (x1, x2) of the plot. }
  \item{ylim}{ the y limits of the plot. }
  \item{xlab}{ a label for the x axis.  }
  \item{ylab}{ a label for the x axis.  }
  \item{hue}{ hue of points (see \code{hsv}) }
  \item{col}{ color used in plots }
  \item{cex.xlab, cex.ylab}{ see \code{\link{par}} }
  \item{axes}{ logical. Should the axes be displayed? }
  \item{box}{ logical. Should a box be drawn around the plot? }
  \item{add}{ logical. Add to current plot? }
  \item{bifurcation.set.fill}{ 1-character string. Color used to fill the bifurcation set (see \code{\link{colors}}).}
  \item{cex.scale, cex, pch}{ see \code{\link{par}} }
}
\details{
  The default hue of each dot is a function of the height of the cusp surface to which it is closest. This is especially usefull in the bifurcation set. Purple dots are higher than green dots.

  The size of the dots depends on the density of dots at its location. The higher the density the larger the dot.
}
\value{
No return value. Called for its side effect.
}
\references{ See \code{\link{cusp-package}} }
\author{ Raoul Grasman }
%\note{ Notes can be added here
%
% Organize in sections like Warning with \section{Warning }{....} 
%}
\seealso{ \code{\link{plot.cusp}}, \code{\link{cusp3d}} }
\examples{
set.seed(20)
# example with regressors
x1 = runif(150)
x2 = runif(150)
z = Vectorize(rcusp)(1, 4*x1-2, 4*x2-1)
data <- data.frame(x1, x2, z)
fit <- cusp(y ~ z, alpha ~ x1+x2, beta ~ x1+x2, data)
\dontrun{
plot(fit, what='bifurcation', box=TRUE, axes=FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
%\keyword{ models }% __ONLY ONE__ keyword per line
