\name{cusp}
\alias{cusp}
\alias{print.cusp}
\title{ Fit a Cusp Catatrophe Model to Data }
\description{
  This function fits a cusp catatrophe model to data using the maximum likelihood method of Cobb. Both the state variable may be modelled by a linear combination of variables and design factors, as well as the normal/asymmetry factor \code{alpha} and bifurction/splitting factor \code{beta}.
}
\usage{
cusp(formula, alpha, beta, data, weights, offset, ..., control =
    glm.control(), method = "cusp.fit", optim.method = "L-BFGS-B", model = TRUE,
    contrasts = NULL)
}
\arguments{
  \item{formula}{ \code{\link{formula}} that models the canonical state variable }
  \item{alpha}{ \code{\link{formula}} that models the canonical normal/asymmetry factor }
  \item{beta}{ \code{\link{formula}} that models the canonical bifurcation/splitting factor }
  \item{data}{ \code{\link{data.frame}} that contains all the variables named in the formulas }
  \item{weights}{ vector of weights by which each data point is weighted (experimental) }
  \item{offset}{ vector of offsets for the data (experimental) }
  \item{\dots}{ named arguments that are passed to \code{\link{optim}} }
  \item{control}{ \code{\link{glm.control}} object, currently unused }
  \item{method}{ string, currently unused }
  \item{optim.method}{ string passed to \code{optim} to choose the optimization algorithm }
  \item{model}{ should the model matrix be returned? }
  \item{contrasts}{ matrix of \code{\link{contrasts}}, experimental }
}
\details{
  \code{cusp} fits a cusp catastrophe model to data. Cobb's definition for the canonical form of the stochastic cusp catastrophe is the stochastic differential equation \deqn{dY_t = (\alpha + \beta Y_t - Y_t^3)dt + dW_t.}{dY(t) = (\alpha + \beta Y(t_ - Y(t)^3)dt + dW(t).} The stationary distribution of the \sQuote{behavioral}, or \sQuote{state} variable \eqn{Y}, given the control parameters \eqn{\alpha} (\sQuote{asymmetry} or \sQuote{normal} factor) and \eqn{\beta} (\sQuote{bifurcation} or \sQuote{splitting} factor) is \deqn{ f(y) = \Psi \exp(\alpha y + \beta y^2/2 - y^4/4), } where \eqn{\Psi} is a normalizing constant.

  The behavioral variable and the asymmetry and bifurcation factors are usually not directly related to the dependent and independent variables in the data set. These are therefore used to predict the state variable and control parameters:
     \deqn{y_i = w_0 + w_1 Y_{i1} + \cdots + w_p Y_{ip}}{y[i] = w[0] + w[1] * Y[i,1] + \cdots + w[p] * Y[i,p],}
     \deqn{\alpha_i = a_0 + a_1 X_{i1} + \cdots + a_p X_{ip}}{\alpha[i] = a[0] + a[1] * X[i,1] + \cdots + a[p] * X[i,p],}
     \deqn{\beta_i = b_0 + b_1 X_{i1} + \cdots + b_q X_{iq}}{\beta[i] = b[0] + b[1] * X[i,1] + \cdots + b[p] * X[i,p],}
     in which the \eqn{a_j}{a[j]}'s, \eqn{b_j}{b[j]}'s, and \eqn{w_j}{w[j]}'s are estimated by means of maximum likelihood.
     Here, the \eqn{Y_{ij}}{Y[i,j]}'s and \eqn{X_{ij}}{X[i,j]}'s are variables constructed from variables in the data set. Variables predicting the \eqn{\alpha}'s and \eqn{\beta}'s need not be the same.

  The state variable and control parameters can be modelled by specifying a model \code{\link{formula}}: \deqn{\code{y ~ model},} \deqn{\code{alpha ~ model},} \deqn{\code{beta ~ model},} in which \code{model} can be any valid \code{\link{formula}} specified in terms of variables that are present in the \code{data.frame}.
}
\value{
  % Describe the structure of value returned
  List with components
   \item{ coefficients }{Estimated coefficients }
 %  \item{ residuals }{ vector of residuals (observed minus predicted assuming Delay convention) }
 %  \item{ fitted.values }{  }
 %  \item{ effects }{description }
 %  \item{ R }{?? of Hessian of objective function evaluated at minimum }
   \item{ rank }{rank of Hessian matrix }
   \item{ qr }{\code{\link{qr}} decomposition of the Hessian matrix }
 %  \item{ family }{description }
   \item{ linear.predictors }{two column matrix containing the \eqn{\alpha_i}{\alpha[i]}'s and \eqn{\beta_i}{\beta[i]}'s for each case }
   \item{ deviance }{sum of squared errors using Delay convention }
   \item{ aic }{ AIC }
   \item{ null.deviance }{variance of canonical state variable }
   \item{ iter }{ number of optimization iterations }
   \item{ weights }{ weights provided through weights argument }
 %  \item{ prior.weights }{description }
   \item{ df.residual }{ residual degrees of freedom }
   \item{ df.null }{ degrees of freedom of constant model for state variable }
   \item{ y }{ predicted values of state variable }
   \item{ converged }{ convergence status given by \code{optim} }
 %  \item{ boundary }{description }
   \item{ par }{ parameter estimates for \code{qr} standardized data }
   \item{ Hessian }{Hessian matrix of negative log likelihood function at minimum }
   \item{ hessian.untransformed }{ Hessian matrix of negative log likelihood for \code{qr} standardized data }
 %  \item{ qr.transform }{description }
   \item{ code }{\code{\link{optim}} convergence indicator }
   \item{ model }{list with model design matrices }
 %  \item{ x }{description }
   \item{ call }{function call that created the object }
   \item{ formula }{list with the formulas }
 %  \item{ terms }{description }
   \item{ OK }{ logical. \code{TRUE} if Hessian matrix is positive definite at the minimum obtained }
   \item{ data }{original data.frame }
 %  \item{ offset }{description }
 %  \item{ control }{description }
 %  \item{ method }{description }
 %  \item{ contrasts }{description }
 %  \item{ xlevels }{description }"
}
\references{ See \link{cusp-package}  }
\author{ Raoul Grasman }
%\note{ notes to user
%
% Make other sections like Warning with \section{Warning }{....}
%}
\seealso{ \code{\link{cusp-package}}.

  \code{\link{summary.cusp}} for summaries and model assessment.

  The generic functions \code{\link{coef}}, \code{\link{effects}}, \code{\link{residuals}}, \code{\link{fitted}}, \code{\link{vcov}}.

  \code{\link{predict}} for estimated values of the control parameters \eqn{\alpha[i]} and \eqn{\beta[i]},
}
\examples{
set.seed(123)
# example with regressors
x1 = runif(150)
x2 = runif(150)
z = Vectorize(rcusp)(1, 4*x1-2, 4*x2-1)
data <- data.frame(x1, x2, z)
fit <- cusp(y ~ z, alpha ~ x1+x2, beta ~ x1+x2, data)
print(fit)
summary(fit)
\dontrun{
plot(fit)
cusp3d(fit)
}

# useful use of OK
\dontrun{
while(!fit$OK)
    fit <- cusp(y ~ z, alpha ~ x1+x2, beta ~ x1+x2, data,
            start=rnorm(fit$par)) # use different starting values
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
