\name{cusp.logist}
\alias{cusp.logist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit a Logistic Surface Model to Data }
\description{
  This function fits a logistic curve model to data using maximum likelihood under the assumption of normal errors (i.e., nonlinear least squares). Both the response variable may be modelled by a linear combination of variables and design factors, as well as the normal/asymmetry factor \code{alpha} and bifurction/splitting factor \code{beta}.
}
\usage{
cusp.logist(formula, alpha, beta, data, ..., model = TRUE, x =
                 FALSE, y = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula, alpha, beta}{ \code{\link{formula}}s for the response variable and the regression variables (see below) }
  \item{data}{ \code{\link{data.frame}} containing \eqn{n} observations of all the variables named in the formulas }
  \item{\dots}{ named arguments that are passed to \code{\link{nlm}} }
  \item{model, x, y}{ logicals. If \code{TRUE} the corresponding components of the fit (the model frame, the model matrix, and the response are returned. }
}
\details{
    A nonlinear regression is carried out of the model
    \deqn{y_i = \frac{1}{1+\exp(-\alpha_i/\beta_i^2)} + \epsilon_i}{%
             y[i] = 1/(1+exp(-\alpha[i]/\beta[i]^2)) + \epsilon[i]}
     for \eqn{i = 1, 2, \ldots, n},
     where
     \deqn{y_i = w_0 + w_1 Y_{i1} + \cdots + w_p Y_{ip}}{y[i] = w[0] + w[1] * Y[i,1] + \cdots + w[p] * Y[i,p],}
     \deqn{\alpha_i = a_0 + a_1 X_{i1} + \cdots + a_p X_{ip}}{\alpha[i] = a[0] + a[1] * X[i,1] + ... + a[p] * X[i,p],}
     \deqn{\beta_i = b_0 + b_1 X_{i1} + \cdots + b_q X_{iq}}{\beta[i] = b[0] + b[1] * X[i,1] + ... + b[p] * X[i,p],}
     in which the \eqn{a_j}{a[j]}'s, and \eqn{b_j}{b[j]}'s, are estimated. The \eqn{Y_{ij}}{Y[i,j]}'s are variables in the data set
     and specified by \code{formula}; the \eqn{X_{ij}}{X[i,j]}'s are variables in the data set and are specified in \code{alpha}
     and \code{beta}. Variables in \code{alpha} and \code{beta} need not be the same. The \eqn{w_j}{w[j]}'s are estimated implicitely
     using concentrated likelihood methods, and are not returned explicitely.
}
\value{
  List with components
    \item{minimum}{Objective function value at minimum}
    \item{estimate}{Coordinates of objective function minimum}
    \item{gradient}{Gradient of objective function at minimum.}
    \item{code}{Convergence \code{code} returned by \code{optim}}
    \item{iterations}{Number of iterations used by \code{optim}}
    \item{coefficients}{A named vector of estimates of \eqn{a_j, b_j}{a[j], b[j]}'s}
    \item{linear.predictors}{Estimates of \eqn{\alpha_i}{\alpha[i]}'s and \eqn{\beta_i}{\beta[i]}'s.}
    \item{fitted.values}{Predicted values of \eqn{y_i}{y[i]}'s as determined from the \code{linear.predictors}}
    \item{residuals}{Residuals}
    \item{rank}{Numerical rank of matrix of predictors for \eqn{\alpha_i}{\alpha[i]}'s plus rank of matrix of predictors for \eqn{\beta_i}{\beta[i]}'s plus rank of matrix of predictors for the \eqn{y_i}{y[i]}'s.}
    \item{deviance}{Residual sum of squares.}
    \item{logLik}{Log of the likelihood at the minimum.}
    \item{aic}{Akaike's information criterion}
    \item{rsq}{R Squared (proportion of explained variance)}
    \item{df.residual}{Degrees of freedom for the residual}
    \item{df.null}{Degrees of freedom for the Null residual}
    \item{rss}{Residual sum of squares}
    \item{hessian}{Hessian matrix of objective function at the minimum if \code{hessian=TRUE}.}
    \item{Hessian}{Hessian matrix of log-likelihood function at the minimum (currently unavailable)}
    \item{qr}{QR decomposition of the \code{hessian} matrix}
    \item{converged}{Boolean indicating if optimization convergence is proper (based on exit code \code{optim}, gradient, and, if \code{hessian=TRUE} eigen values of the hessian).}
    \item{weights}{\code{weights} (currently unused)}
    \item{call}{the matched call}
    \item{y}{If requested (the default), the matrix of response variables used.}
    \item{x}{If requested, the model matrix  used.}
    \item{null.deviance}{The sum of squared deviations from the mean of the estimated \eqn{y_i}{y[i]}'s.}
}
\references{
Hartelman PAI (1997). \emph{Stochastic Catastrophe Theory.} Amsterdam: University of Amsterdam, PhDthesis.
}
\author{ Raoul Grasman }
%\note{
%}
\seealso{ \code{\link{summary.cusp}} }
%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ nonlinear }% __ONLY ONE__ keyword per line
\keyword{ multivariate }
