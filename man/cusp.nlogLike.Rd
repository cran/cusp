\name{cusp.nlogLike}
\alias{cusp.nlogLike}
\alias{cusp.nlogLike.c}
\alias{cusp.logLike}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Negative log-likelihood for Cobb's cusp density}
\description{
  (Negative) log-likelihood for Cobb's cusp probability density function used by \code{cusp}. This function is not to be called by the user. See \code{help(cusp)}.
}
\usage{
cusp.nlogLike(p, y, X.alpha, X.beta = X.alpha, ..., verbose = FALSE)
cusp.nlogLike.c(p, y, X.alpha, X.beta = X.alpha, ..., verbose = FALSE)
cusp.logLike(p, x, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ parameter vector }
  \item{x}{ vector of observed values for the state variable in the cusp (\code{cusp.nlogLike} only) }
  \item{y}{ design matrix predicting state values at which the likelihood is evaluated }
  \item{X.alpha}{ design matrix predicting alpha in the model }
  \item{X.beta}{ design matrix predicting beta in the model }
  \item{\dots}{ unused extra arguments }
  \item{verbose}{ logical, if \code{TRUE} the value of the parameters are printed to the console }
}
\details{
	\code{cusp.nlogLike} is the R version of the corresponding C function wrapped by \code{cusp.nlogLike.c}
	
	These functions are not intended to be called directly by the user.
}
\value{
	The value of the negative log-likelihood function (\code{cusp.nlogLike}, \code{cusp.nlogLike.c}), 
	the value of the log-likelihood function (\code{cusp.logLike}).
}
\references{ See \link{cusp-package} }
\author{ Raoul Grasman }
\note{
	The functions are not to be called by the user directly.
}
\seealso{ \code{\link{cusp}}, \link{cusp-package} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\examples{
	y = rcusp(100,1,0)
    cusp:::cusp.logLike(c(a=1,b=0,w0=0,w1=1),y)
}
\keyword{ utilities }
\keyword{ models }% __ONLY ONE__ keyword per line
