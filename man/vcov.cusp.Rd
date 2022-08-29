\name{vcov.cusp}
\alias{vcov.cusp}
\alias{confint.cusp}
\title{ Calculate Variance-Covariance Matrix for a Fitted Cusp Model Object }
\description{
  Returns an estimate of the variance-covariance matrix of the main parameters of a fitted cusp model object.
}
\usage{
\method{vcov}{cusp}(object, ...)
\method{confint}{cusp}(object, parm, level = 0.95, ...)
}
\arguments{
  \item{object}{ a fitted cusp model object. }
  \item{parm}{ a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered. }
  \item{level}{ the confidence level required.}  
  \item{\dots}{ additional arguments for method functions.  }
}
\details{
  The variance-covariance matrix is estimated by the inverse of the Hessian matrix of the log-likelhood at the maximum likelihood estimate (\code{vcov}).

  Normal theory confidence intervals are computed for all parameters in the cusp model object using \code{vcov} to obtain the standard errors (\code{confint}).
}
\value{
    The variance-covariance matrix (\code{vcov}).
    
    A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2 in % (by default 2.5% and 97.5%) (\code{confint}).
}
\references{ Seber, Wild (2005) Nonlinear regression. New York: Wiley}
\author{ Raoul Grasman }
%\note{ Some other notes
%
% Use sections like Warning with \section{Warning }{....} to organize notes.
%}
\seealso{ \code{\link{vcov}}, \code{\link{cusp}} }
\examples{
set.seed(123)
x1 = runif(150)
x2 = runif(150)
z = Vectorize(rcusp)(1, 4*x1-2, 4*x2-1)
data <- data.frame(x1, x2, z)
fit <- cusp(y ~ z, alpha ~ x1+x2, beta ~ x1+x2, data)

vcov(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
\keyword{ models }% __ONLY ONE__ keyword per line
