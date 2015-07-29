\name{predict.cusp}
\alias{predict.cusp}
\title{
Predict method for Cusp Model Fits
}
\description{
Predicted values based on a cusp model object.
}
\usage{
\method{predict}{cusp}(object, newdata, se.fit = FALSE, interval =
  c("none", "confidence", "prediction"), level = 0.95, type = c("response", "terms"),
  terms = NULL, na.action = na.pass, pred.var = res.var/weights, weights = 1,
  method = c("delay", "maxwell", "expected"), keep.linear.predictors = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object of class "\code{cusp}"
}
  \item{newdata}{
An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
}
  \item{se.fit}{
See \code{\link{predict.lm}}. Not yet used.
}
  \item{interval}{
See \code{\link{predict.lm}}. Not yet used.
}
  \item{level}{
See \code{\link{predict.lm}}. Not yet used.
}
  \item{type}{
See \code{\link{predict.lm}}. Not yet used.
}
  \item{terms}{
See \code{\link{predict.lm}}. Not yet used.
}
  \item{na.action}{
See \code{\link{predict.lm}}. Not yet used.
}
  \item{pred.var}{
See \code{\link{predict.lm}}. Not yet used.
}
  \item{weights}{
See \code{\link{predict.lm}}. Not yet used.
}
  \item{method}{
Type of prediction convention to use. Can be abbreviated. (\code{expected} should currently not be trusted).
}
  \item{keep.linear.predictors}{
Logical. Should the linear predictors (alpha, beta, and y) be returned?
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
\details{
\code{predict.cusp} produces predicted values, obtained by evaluating the regression functions from the
cusp \code{object} in the frame \code{newdata} using \code{predict.lm}. This results in linear
predictors for the cusp control variables \code{alpha}, and \code{beta}, and, if \code{method = "delay"},
for the behavioral cusp variable \code{y}. These are then used to compute predicted values: If
\code{method = "delay"} these are the points \eqn{y*} on the cusp surface defined by
\deqn{V'(y*) = \alpha + \beta y* - y*^3 = 0} that are closest to \code{y}. If \code{method = "maxwell"} they are
the points on the cusp surface corresponding to the minimum of the associated potential function
\eqn{V(y*) = \alpha y* + 0.5 y*^2 - 0.25 y*^4}.
}
\value{
A vector of predictions. If \code{keep.linear.predictors} the return value has a \code{"data"} attribute
which links to \code{newdata} augmented with the linear predictors \code{alpha}, \code{beta}, and, if
\code{method = "delay"}, \code{y}. If \code{method = "expected"}, the expected value from the equilibrium
distribution of the stochastic process \deqn{dY_t = V'(Y_t;\alpha, \beta)dt + dW_t,} where \eqn{W_t} is
a Wiener proces (aka Brownian motion) is returned. (This distribution is implemented in
\code{\link{dcusp}}.)
}
\references{
See \code{\link{cusp-package}}.
}
\author{
Raoul Grasman
}
\note{
Currently \code{method = "expected"} should not be trusted.
}


\seealso{
\code{\link{cusp-package}}, \code{\link{predict.lm}}.
}
\examples{
set.seed(123)
# example with regressors
x1 = runif(150)
x2 = runif(150)
z = Vectorize(rcusp)(1, 4*x1-2, 4*x2-1)
data <- data.frame(x1, x2, z)
fit <- cusp(y ~ z, alpha ~ x1+x2, beta ~ x1+x2, data)

newdata = data.frame(x1 = runif(10), x2 = runif(10), z = 0)
predict(fit, newdata)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
