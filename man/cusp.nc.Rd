\name{cusp.nc}
\alias{cusp.nc}
\alias{cusp.nc.C}
\alias{cusp.nc.c}
\alias{cusp.nc.vec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate the Normalizing Constant of Cobb's Cusp Density}
\description{
  A family of functions that return the normalization constant for the cusp density given the values of the bifurcation and asymmetry parameters (default), or returns the moment of a specified order (\code{cusp.nc}).
}
\usage{
cusp.nc(alpha, beta, mom.order = 0, ...)
cusp.nc.c(alpha, beta, ..., keep.order = TRUE) 
cusp.nc.C(alpha, beta, subdivisions = 100, rel.tol = .Machine$double.eps^0.25, 
    abs.tol = rel.tol, stop.on.error = TRUE, aux = NULL, keep.order = TRUE) 
cusp.nc.vec(alpha, beta, ..., keep.order = FALSE) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{ the asymmetry parameter in Cobb's cusp density (see \code{\link{cusp}}) }
  \item{beta}{ the bifurcation parameter in Cobb's cusp density (see \code{\link{cusp}}) }
  \item{mom.order}{ the moment order to be computed (see details below) }
  \item{subdivisions,
  rel.tol, 
  abs.tol,
  stop.on.error,
  aux}{ arguments used by the internal integration routine of R (see \code{\link{integrate}})}
  \item{keep.order}{logical, that indicates wether the order of the output should be the same as the order of the input}
  \item{\dots}{ extra arguments in \code{cusp.nc.c} that are passed to \code{cusp.nc.C} }
}
\details{
	The function \code{cusp.nc} returns \eqn{\Psi} if \code{mom.order = 0} and \eqn{\Psi} times the moment of order \code{mom.order} otherwise.
	
	The function \code{cusp.nc} is internally used if the C-routine symbol \code{"cuspnc"} is not loaded.
	The functions \code{cusp.nc.c} and \code{cusp.nc.C} call this C routine, which is considerably faster than
	\code{cusp.nc}. 
	
	These functions are not intended to be called directly by the user.
}
\value{
    \code{cusp.nc, cusp.nc.c, cusp.nc.vec} return a numeric vector of the same length as \code{alpha} and \code{beta} with normalizing constants, or the indicated moments times the normalization constant (\code{cusp.nc} only).
    
    \code{cusp.nc.C} returns a list with vectors with the results obtained from \code{\link{integrate}}. 
    \code{cusp.nc.c} first sorts the input in such a way that the numerical integrals can be evaluated more quickly than in arbitrary order
}
\author{ Raoul Grasman }
%\note{ further notes 
%
% Make other sections like Warning with \section{Warning }{....} 
%}
\seealso{ \code{\link{pcusp}}, \code{\link{dcusp}} }
%\examples{
%    # dcusp(x,2,1) = exp(2*x + x^2 - x^4/4) / cusp.nc(2,1), hence
%    # evaluating this at x=0, cusp.nc(2,1) is equal to
%    1/dcusp(0,2,1)
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
\keyword{ models }% __ONLY ONE__ keyword per line
