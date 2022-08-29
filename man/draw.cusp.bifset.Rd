\name{draw.cusp.bifset}
\alias{draw.cusp.bifset}
\title{ Add Cusp Bifurcation Set Diagram to Existing Plot }
\description{
  Add a miniature bifurcation set for the cusp catastrophe to an existing plot. 
}
\usage{
draw.cusp.bifset(rx = par("usr")[1:2], ry = par("usr")[3:4], xpos = min(rx) + 
    0.01 * diff(rx)[1], ypos = max(ry) - 0.01 * diff(ry)[1], 
    xscale =  0.1 * diff(rx), yscale =  0.1 * diff(ry) / xscale, 
    aspect = 1, mark = 1, col = hsv(0.7, s = 0.8, alpha = 0.5), 
    border = NA, density = NA, bifurcation.set.fill = gray(0.8), 
    background = hsv(0.1, s = 0.1, alpha = 0.5), ..., X)  
}
\arguments{
  \item{rx}{ x-axis range of the plot window }
  \item{ry}{ y-axis range of the plot window }
  \item{xpos}{ x-axis position of drawing }
  \item{ypos}{ y-axis position of drawing }
  \item{xscale}{ scaling applied to drawing along x-axis }
  \item{yscale}{ scaling applied to drawing along y-axis }
  \item{aspect}{ aspect ratio  }
  \item{mark}{ 0, 1, 2, 3, or 4; indicates which part of the cusp surface should be marked }
  \item{col}{ color used for marking a part of the cusp surface }
  \item{border}{ color used for the marked part of the cusp surface. See \code{polygon} for details. }
  \item{density}{ the density of shading lines of the marked part of the cusp surface, in lines per inch. The default value of NULL means that no shading lines are drawn. See \code{polygon} for details. }
  \item{bifurcation.set.fill}{ color for marking the bifurcation set }
  \item{background}{ background color of the cusp surface }
  \item{\dots}{ arguments passed to \code{rect} and \code{polygon} }
  \item{X}{ \code{data.frame}, depricated }
}
\details{
  This function is mainly intended for internal use by \code{cusp.plot}.
}
\value{
No return value. Called for its side effect.
}
\references{ \url{http://www.scholarpedia.org/article/Cusp_bifurcation} }
\author{ Raoul Grasman }
%\note{ More notes?
% Make other sections within this like Warning with \section{Warning }{....}
%}
\seealso{ \code{\link{plot.cusp}}, \code{\link{polygon}} }
\examples{
\dontrun{
plot(1:10)
draw.cusp.bifset(mark=0) # no marking
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }
\keyword{ models }% __ONLY ONE__ keyword per line
