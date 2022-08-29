\name{cusp3d.surface}
\alias{cusp3d.surface}
\title{ Generate 3D plot of the Cusp surface }
\description{
  This function generates a 3D display of the cusp equilibrium surface.
}
\usage{
cusp3d.surface(alpha = c(-5, 5), beta = c(-3, 3), y = 41, 
xlim = range(alpha), ylim = range(beta), zlim = c(-5, 4), 
xlab = expression(alpha), ylab = expression(beta), zlab = "equilibrium states", 
main = NULL, sub = NULL, phi = 20, theta = 160, 
r = sqrt(3), d = 1, scale = TRUE, expand = 1, hue = 240, 
chroma = 35, surf.alpha = 0.75, gamma = 1.5, bcol = NA, 
lcol = "gray", ltheta = 90, lphi = 70, box = TRUE, 
axes = FALSE, nticks = 5, ticktype = "simple", floor.lines = TRUE, ...)
}

\arguments{
  \item{alpha}{ numeric 2-vector specifying the normal/symmetry factor axis range }
  \item{beta}{ numeric 2-vector specifying the bifurcation/splitting factor axis range }
  \item{y}{ numeric specifying the iso contours used to render the surface (see detais below) }
  \item{xlim,ylim,zlim}{ numeric 2-vectors (see \code{\link{persp}}) }
  \item{xlab,ylab,zlab,main,sub}{ strings (see \code{\link{persp}}) }
  \item{phi,theta}{ numeric, determine viewing direction (see \code{\link{persp}}) }
  \item{r}{ numeric, distance to center of the plotting box (see \code{\link{persp}}) }
  \item{d}{ numeric, strength of perspective transformation (see \code{\link{persp}}) }
  \item{scale,expand}{ logical, see \code{\link{persp}} }
  \item{hue,chroma,surf.alpha}{ hue, chroma and alpha (transparency) of the surface segments (see \code{\link{hcl}})}
  \item{gamma}{ gamma for shading of surface (see \code{\link{cusp3d}}) }
  \item{bcol}{ color, \code{NA}, or string \code{"surface"}. Color of the border of each surface element; \code{NA} gives transparent borders; \code{"surface"} tries to hide the border as much as possible by giving it the same color as the surface segment. }
  \item{lcol}{ color of the lines on the floor of the plotting cube }
  \item{ltheta,lphi}{ numeric, direction of illumination of the surface (similar to \code{\link{persp}})}
  \item{box,axes,nticks,ticktype}{ (see \code{\link{persp}}) }
  \item{floor.lines}{ logical, if \code{TRUE} (default) iso-contours are projected on the floor of the plotting cube (revealing the bifurcation set) }
  \item{\dots}{ extra arguments that are passed to \code{\link{lines}} and \code{\link{polygon}}}
}
\details{
  If \code{y} has length 1, it is interpreted as the number of contours. Otherwise it is interpreted as a vector of contour levels from which the surface must be determined. If \code{y} is a number, the exact range of \code{y} is determined by the ranges of \code{alpha} and \code{beta} through the cusp equilibrium equation below.
  
	The surface is constructed from the iso-contours of the cusp equilibrium surface that makes up the solutions to \deqn{\alpha + \beta*y - y^3 = 0}
	as a (multi-)function of the asymmetry variable \eqn{\alpha} and bifurcation variable \eqn{\beta}. For each possible solution \eqn{y} the iso-contours are given by the equation \deqn{\alpha = (\beta*y - y^3)/y,} which are linear in \eqn{\beta}. For each value of \eqn{y} the values of \eqn{alpha} are determined for the end points of the \eqn{beta} range specified by \code{beta}. The two 3D coordinates (\eqn{\alpha}, \eqn{\beta}, \eqn{y}) are projected onto the 2D canvas using the \code{\link{persp}} transformation matrix and used for drawing the lines and polygons.
}
\value{
    \code{cusp3d.surface} returns the viewing transformation matrix, say \code{VT}, a 4 x 4 matrix suitable for projecting 3D coordinates (x,y,z) into the 2D plane using homogeneous 4D coordinates (x,y,z,t). It can be used to superimpose additional graphical elements on the 3D plot, by lines() or points(), using the simple function trans3d().
}
\references{ See \code{\link{cusp-package}}, \code{\link{cusp3d}} }
\author{ Raoul Grasman }
\note{   This function is an alternative to \code{\link{cusp3d}} which uses a different method of rendering and also plots fitted points on the surface.
}
\seealso{ \code{\link{persp}}, \code{\link{plot.cusp}} }
\examples{
\dontrun{
p = cusp3d.surface(chroma=40,lcol=1,surf.alpha=.95,phi=30,theta=150,
	bcol="surface",axes=TRUE,main="Cusp Equilibrium Surface")
lines(trans3d(c(5,5), c(3,3), c(-5,4), p), lty=3) # replot some of the box outlines 
lines(trans3d(c(-5,5), c(3,3), c(4,4), p), lty=3)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ models }% __ONLY ONE__ keyword per line
