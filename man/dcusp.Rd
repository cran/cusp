\name{dcusp}
\alias{dcusp}
\alias{rcusp}
\alias{pcusp}
\alias{qcusp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Cobb's Cusp Distribution }
\description{
  Functions for the cusp distribution.
}
\usage{
dcusp(y, alpha, beta)
pcusp(y, alpha, beta)
qcusp(p, alpha, beta)
rcusp(n, alpha, beta) %, scale = 4.2, ..., ignore.warnings = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ vector of quantiles }
  \item{p}{ vector of probabilities }
  \item{n}{ number of observations. }
  \item{alpha}{ normal/asymmetry factor value of cusp density }
  \item{beta}{ bifurcation/splitting factor value of cusp density }
  %\item{scale}{ depricated } %\code{scale} parameter of \code{metrop}, see Details below }
  %\item{ignore.warnings}{ depricated } %ignore warnings from \code{metrop}? }
}
\details{
  The cusp distribution is defined by \deqn{ f(y) = \Psi \exp(\alpha y + \beta y^2/2 - y^4/4), } where \eqn{\Psi} is the normalizing constant.
  
  \code{rcusp} uses rejection sampling to generate samples.
  
  \code{qcusp} implements binary search and is rather slow.
}
\value{ \code{dcusp} gives the density function, \code{pcusp} gives the distribution function, \code{qcusp} gives the quantile function, and \code{rcusp} generates observations.
}
\references{ See \code{\link{cusp-package}} }
\author{ Raoul Grasman }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{cusp-package}} } %, \code{\link{metrop}} }
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
