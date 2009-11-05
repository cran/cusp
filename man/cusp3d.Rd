\name{cusp3d}
\alias{cusp3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate 3D plot of Cusp Catatrophe Model Fit }
\description{
  This function generates a 3D display of the fit (object) of a cusp model.
}
\usage{
cusp3d(y, alpha = if (!missing(y) && is.list(y)) y$lin[, "alpha"], 
    beta = if (!missing(y) && is.list(y)) y$lin[, "beta"], w = 0.03, 
    theta = 170, phi = 35, B = 4, Y = 3, Yfloor = -15, 
    np = 180, n.surface = 30, surface.plot = TRUE, 
    surf.alpha = 0.75, surf.gamma = 1.5, surf.chroma = 35, surf.hue = 240, 
    surf.ltheta = 0, surf.lphi = 45, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ object returned by \code{\link{cusp}} or a vector of observed state values }
  \item{alpha}{ vector of normal/symmetry factor values corresponding to the state values in \code{y} }
  \item{beta}{ vector of bifurcation/splitting factor values corresponding to the state values in \code{y} }
  \item{w}{ number that specifies the size of the data points plotted on the cusp surface }
  \item{theta,
  phi}{ angles defining the viewing direction. \code{theta} gives the azimuthal direction and \code{phi} the colatitude. }
  \item{B}{ range of the splitting factor axis }
  \item{Y}{ range of the state variable axis }
  \item{Yfloor}{ location on state variable axis where the control surface is plotted }
  \item{np}{ factor that determines the fineness of the drawing }
  \item{n.surface}{ factor that determines the fineness of the rendered surface }
  \item{surface.plot}{ plot the surface? }
  \item{surf.alpha}{ transparency level of rendered surface }
  \item{surf.gamma}{ factor that determines the shading of surface facets (\code{surf.gamma}<1 deminishes shading, \code{surf.gamma}>1 exagerates shading) }
  \item{surf.chroma,
  surf.hue}{ chroma and hue of surface color (see \code{\link{hcl}}) }
  \item{surf.ltheta,
  surf.lphi}{ the surface is shaded as though it was being illuminated from the direction specified by azimuth \code{surf.ltheta} and colatitude \code{surf.lphi} }
  \item{\dots}{ named parameters that are pased to \code{\link{persp}} }
}
\details{
  This function is experimental.
}
\value{
    \code{cusp3d} returns the viewing transformation matrix, say \code{VT}, a 4 x 4 matrix suitable for projecting 3D coordinates (x,y,z) into the 2D plane using homogeneous 4D coordinates (x,y,z,t). It can be used to superimpose additional graphical elements on the 3D plot, by lines() or points(), using the simple function trans3d().
}
\references{ See \code{\link{cusp-package}} }
\author{ Raoul Grasman }
\note{ Currently still somewhat buggy.
}
\seealso{ \code{\link{persp}}, \code{\link{plot.cusp}}, \code{\link{cusp3d.surface}} }
\examples{
x1 = runif(150)
x2 = runif(150)
z = Vectorize(rcusp)(1, 4*x1-2, 4*x2-1)
data <- data.frame(x1, x2, z)
fit <- cusp(y ~ z, alpha ~ x1+x2, beta ~ x1+x2, data)
\dontrun{
cusp3d(fit)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ models }% __ONLY ONE__ keyword per line
