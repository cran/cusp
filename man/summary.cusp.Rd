\name{summary.cusp}
\alias{summary.cusp}
\alias{print.summary.cusp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Summarizing Cusp Catastrophe Model Fits }
\description{
  \code{summary} method for class \dQuote{cusp}
}
\usage{
\method{summary}{cusp}(object, correlation = FALSE, symbolic.cor = FALSE, logist = FALSE, ...)

\method{print}{summary.cusp}(x, digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ Object returned by \code{cusp} }
  \item{x}{ \sQuote{\code{summary.cusp}} object }
  \item{correlation}{ logical; if \code{TRUE} the correlation matrix is returned }
  \item{symbolic.cor}{ logical; currently unused }
  \item{logist}{ logical. If \code{TRUE} a logistic model is fitted for cusp model assesment (see \code{\link{cusp.logist}} for details).}
  \item{digits}{ numeric; the number of significant digits to use when printing. }
  \item{signif.stars}{ logical. If \code{TRUE}, significance stars are printed for each coefficient. }
  \item{\dots}{ further arguments passed to or from other methods. }
}
\details{
\code{print.summary.cusp} tries to be smart about formatting the coefficients, standard errors, etc. and additionally gives significance stars if \code{signif.stars} is \code{TRUE}.

Correlations are printed to two decimal places (or symbolically): to see the actual correlations print \code{summary(object)$correlation} directly.}
\value{
The function \code{summary.cusp} computes and returns a list of summary statistics of the fitted linear model given in object, using the components (list elements) \dQuote{\code{call}} and \dQuote{\code{terms}} from its argument, plus 
    \item{ call }{the matched call}
    \item{ terms }{ the \code{\link{terms}} object used.}
%    \item{ family }{}
    \item{ deviance }{ sum of squared residuals of cusp model fit }
    \item{ aic }{ Akaike Information Criterion for cusp model fit }
    \item{ contrasts }{ contrasts used}
    \item{ df.residual }{ degrees of freedom for the residuals of the cusp model fit }
    \item{ null.deviance }{ variance of canonical state variable}
    \item{ df.null }{degrees of freedom of constant model for state variable}
    \item{ iter }{number of optimization iterations}
    \item{ deviance.resid }{ residuals computed by \code{\link{residuals.glm}} using \code{type="deviance"}}
    \item{ coefficients }{ a \eqn{p \times 4} matrix with columns for the estimated coefficient, its standard error, t-statistic and corresponding (two-sided) p-value. Aliased coefficients are omitted.}
    \item{ aliased }{ named logical vector showing if the original coefficients are aliased.}
    \item{ dispersion }{ always 1 }
    \item{ df }{ 3-vector containing the rank of the model matrix, residual degrees of freedom, and model degrees of freedom. }
    \item{ resid.name }{ string specifying the convention used in determining the residuals (i.e., "Delay" or "Maxwell").}
    \item{ cov.unscaled }{ the unscaled (dispersion = 1) estimated covariance matrix of the estimated coefficients.}
%    \item{ cov.scaled }{}
    \item{ r2lin.r.squared }{\eqn{R^2}, the \sQuote{fraction of variance explained} by the linear regression model \deqn{w_0+w_1 Y_{i1} + \cdots + w_p Y_{ip} = \beta_0 + \beta_1 X_{i1} + \cdots + \beta_q X_{iq} + \epsilon_i,}{w[0]+w[1]*Y[i,1] + ... + w[p]*Y[i,p] = \beta[0]+\beta[1]*X[i,1] + ... + \beta[q]*X[i,q] + \epsilon[i],} where \eqn{Y} containes all explanatory variables for the behavioral states in the cusp model, and \eqn{X} containes all explanatory variables for the control parameters of the cusp model. This is computed from the largest canonical correlation.}
    \item{ r2lin.dev }{ residual sums of squares of the linear model }
    \item{ r2lin.df }{ degrees of freedom for the linear model }
    \item{ r2lin.logLik }{ value of the log-likelihood for the linear model assuming normal errors }
    \item{ r2lin.npar }{ number of parameters in the linear model }
    \item{ r2lin.aic }{ AIC for the linear model }
    \item{ r2lin.aicc }{ corrected AIC for the linear model }
    \item{ r2lin.bic }{ BIC for the linear model}
    \item{ r2log.r.squared }{\eqn{R^2}, the \sQuote{fraction of variance explained} by the logistic model. See \code{\link{cusp.logist}} for details. }
    \item{ r2log.dev }{ if \code{logist = TRUE} residual sums of square for the logistic  model}
    \item{ r2log.df }{ ditto, degrees of freedom for the logistic model}
    \item{ r2log.logLik }{ ditto, value of log-likelihood function for the logistic model assuming normal errors. }
    \item{ r2log.npar }{ ditto, number of parameters for the logistic model }
    \item{ r2log.aic }{ ditto, AIC for logistic model }
    \item{ r2log.aicc }{ ditto, corrected AIC for logistic model }
    \item{ r2log.bic }{ ditto, BIC for logistic model }
    \item{ r2cusp.r.squared }{ pseudo-\eqn{R^2}, the \sQuote{fraction of variance explained by the cusp model}, \deqn{R^2 = 1 - \frac{Var(residuals_i)}{Var(y_i)}.}{R^2 = 1 - Var(residuals[i])/Var(y[i]).} This value can be negative.}
    \item{ r2cusp.dev }{residual sums of squares for cusp model}
    \item{ r2cusp.df }{residual degrees of freedom for cusp model}
    \item{ r2cusp.logLik }{value of the log-likelihood function for the cusp model}
    \item{ r2cusp.npar }{number of parameters in the cusp model}
    \item{ r2cusp.aic }{AIC for cusp model fit}
    \item{ r2cusp.aicc }{corrected AIC for cusp model fit}
    \item{ r2cusp.bic }{BIC for cusp model fit.}
}
\references{ 
    Cobb L, Zacks S (1985). \emph{Applications of Catastrophe Theory for Statistical Modeling in the Biosciences.} Journal of the American Statistical Association, 80(392), 793--802.
    
    Hartelman PAI (1997). \emph{Stochastic Catastrophe Theory.} Amsterdam: University of Amsterdam, PhDthesis. 
    
    Cobb L (1998). \emph{An Introduction to Cusp Surface Analysis.} 
        \cr\url{http://www.aetheling.com/models/cusp/Intro.htm}. 
 }
\author{ Raoul Grasman }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{cusp}}, \code{\link{cusp.logist}} }
\examples{
x1 = runif(150)
x2 = runif(150)
z = Vectorize(rcusp)(1, 4*x1-2, 4*x2-1)
data <- data.frame(x1, x2, z)
fit <- cusp(y ~ z, alpha ~ x1+x2, beta ~ x1+x2, data)
print(fit)
summary(fit, logist=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
\keyword{ models }% __ONLY ONE__ keyword per line
