\name{oliva}
\alias{oliva}
\docType{data}
\title{ Synthetic cusp data set}
\description{
  Synthetic \sQuote{multivariate} data from the cusp catastrophe as generated from the equations specified by Oliva et al. (1987).
}
\usage{data(oliva)}
\format{
  A data frame with 50 observations on the following 12 variables.
  \describe{
    \item{\code{x1}}{splitting factor predictor}
    \item{\code{x2}}{splitting factor predictor}
    \item{\code{x3}}{splitting factor predictor}
    \item{\code{y1}}{the bifurcation  factor predictor}
    \item{\code{y2}}{the bifurcation  factor predictor}
    \item{\code{y3}}{the bifurcation  factor predictor}
    \item{\code{y4}}{the bifurcation  factor predictor}
    \item{\code{z1}}{the state factor predictor}
    \item{\code{z2}}{the state factor predictor}
    \item{\code{alpha}}{the true \eqn{alpha}'s}
    \item{\code{beta}}{the true \eqn{beta}'s}
    \item{\code{y}}{the true state variable values}
  }
}
\details{
	The data in Oliva et al. (1987) are obtained from the equations
	\deqn{\alpha_i = X_{i1} - .969\,X_{i2} - .201\,X_{i3}, }
	\deqn{\beta_i  = .44\,Y_{i1} + 0.08\,Y_{i2} + .67\,Y_{i3} + .19\,Y_{i4}, }
	\deqn{y_i      = -0.52\,Z_{i1} - 1.60\,Z_{i2}.}
	Here the \eqn{X_{ij}}'s are uniformly distributed on (-2,2), and the \eqn{Y_{ij}}'s and \eqn{Z_{i1}} are 
	uniform on (-3,3).
	The states \eqn{y_i} were then generated from the cusp density, using \code{\link{rcusp}}, with their respective 
	\eqn{\alpha_i}'s and \eqn{\beta_i}'s as normal and splitting factors, and then \eqn{Z_2} was computed as 
	\deqn{Z_{i2} = (y_i + 0.52 Z_{i1} )/( 1.60).}
}
\source{
  Oliva T, Desarbo W, Day D, Jedidi K (1987). GEMCAT: A general multivariate 
methodology for estimating catastrophe models. Behavioral Science, 32(2), 121137.
}
\references{
  Oliva T, Desarbo W, Day D, Jedidi K (1987). GEMCAT: A general multivariate 
methodology for estimating catastrophe models. Behavioral Science, 32(2), 121137.
}
\examples{
data(oliva)
set.seed(121)
fit <- cusp(y ~ z1 + z2 - 1, 
	alpha ~ x1 + x2 + x3 - 1, ~ y1 + y2 + y3 + y4 - 1, 
	data = oliva, start = rnorm(9))
summary(fit)
\dontrun{
cusp3d(fit, B=5.25, n.surf=50, theta=150) 
# B modifies the range of beta (is set here to 5.25 to make 
# sure all points lie on the surface)
}
}
\keyword{datasets}
