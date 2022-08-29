\name{cusp.extrema}
\alias{cusp.extrema}
\title{ Locate Extrema of Cusp Catastrophe Potential Function }
\description{
  This function computes the locations of the extrema of the cusp catastrophe potential function.
}
\usage{
cusp.extrema(alpha, beta)
}
\arguments{
  \item{alpha}{ (single) value of normal/symmetry factor }
  \item{beta}{ (single) value of bifurcation/splitting factor }
}
\details{
  The locations are determined by computing the solutions to the equation \deqn{\alpha+\beta\,X - X^3 = 0.}{\alpha+\beta * X - X^3 = 0.}
}
\value{ Ordered vector with locations of extremes.
}
\references{ \url{http://www.scholarpedia.org/article/Cusp_bifurcation} }
\author{ Raoul Grasman }
\note{ Use \code{\link{Vectorize}} to allow for array input. 
}
\seealso{ \code{\link{cusp.bifset}} }
\examples{
    # simple use
    cusp.extrema(2,3)
    
    # using vectorize to allow for array input; 
    # returns a matrix with locations in each column
    Vectorize(cusp.extrema)(-3:3, 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
\keyword{ utilities }% __ONLY ONE__ keyword per line
