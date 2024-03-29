\name{dcusp}
\alias{dcusp}
\alias{rcusp}
\alias{pcusp}
\alias{qcusp}
\title{ Cobb's Cusp Distribution }
\description{
  Functions for the cusp distribution.
}
\usage{
dcusp(y, alpha, beta)
pcusp(y, alpha, beta, subdivisions = 100, rel.tol = .Machine$double.eps^0.25, 
            abs.tol = rel.tol, stop.on.error = TRUE, aux = NULL, keep.order = TRUE)
qcusp(p, alpha, beta)
rcusp(n, alpha, beta) 
}

\arguments{
  \item{y}{ vector of quantiles }
  \item{p}{ vector of probabilities }
  \item{n}{ number of observations. }
  \item{alpha}{ normal/asymmetry factor value of cusp density }
  \item{beta}{ bifurcation/splitting factor value of cusp density }
  \item{subdivisions}{See \code{\link{cusp-package}}.}
  \item{rel.tol}{See \code{\link{cusp-package}}.}
  \item{abs.tol}{See \code{\link{cusp-package}}.}
  \item{stop.on.error}{See \code{\link{cusp-package}}.}
  \item{aux}{See \code{\link{cusp-package}}.}
  \item{keep.order}{ logical. If true the order of the output values is the same as those of the input values \code{y}}
}
\details{
  The cusp distribution is defined by \deqn{ f(y) = \Psi \exp(\alpha y + \beta y^2/2 - y^4/4), } where \eqn{\Psi} is the normalizing constant.
  
  \code{rcusp} uses rejection sampling to generate samples.
  
  \code{qcusp} implements binary search and is rather slow.
}
\value{ \code{dcusp} gives the density function, \code{pcusp} gives the distribution function, \code{qcusp} gives the quantile function, and \code{rcusp} generates observations.
}
\references{ See \code{\link{cusp-package}}, \code{\link{integrate}} }
\author{ Raoul Grasman }
%\note{ any other notes go here
%
% Make other sections like Warning with \section{Warning }{....}
%}
\seealso{ \code{\link{cusp-package}} } 
\examples{
    # evaluate density and distribution
    dcusp(0,2,3)
    pcusp(0,2,3)
    pcusp(qcusp(0.125,2,3),2,3) # = 0.125
    
    # generate cusp variates
    rcusp(100, 2, 3) 
    
    # generate cusp variates for random normal and splitting factor values
    alpha = runif(20, -3, 3)
    beta =  runif(20, -3, 3)
    Vectorize(rcusp)(1, alpha, beta)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
\keyword{ distribution }% __ONLY ONE__ keyword per line
\keyword{ univar }
