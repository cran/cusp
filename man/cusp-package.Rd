\name{cusp-package}
\alias{cusp-package}
%\alias{cusp}
\docType{package}
\title{
Cusp Catastrophe Modeling
}
\description{
Fits cusp catastrophe to data using Cobb's maximum likelihood method with a different algorithm. The package contains  utility functions for plotting, and for comparing the model to linear regression and logistic curve models. The package allows for multivariate response subspace modelling in the sense of the GEMCAT software of Oliva et al.
}
\details{
\tabular{ll}{
Package: \tab cusp\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2008-02-14\cr
License: \tab GNU GPL v2 (or higher)\cr
}
This package helps fitting Cusp catastrophy models to data, as advanced in Cobb et al. (1985). The main functions are 
\tabular{ll}{
    \command{cusp} \tab Fit Cobb's Cusp catastrophe model; see example below.\cr
    \command{summary.cusp} \tab Summary statistics of cusp model fit.\cr
    \command{confint.cusp} \tab Confidence intervals for parameter estimates\cr
    \command{plot.cusp} \tab Diagnostic plots for cusp model fit\cr
    \command{cusp3d} \tab 3D graphical display of cusp model fit (experimental).\cr
    \command{dcusp} \tab Density of Cobb's cusp distribtution\cr
    \command{pcusp} \tab Cumulative probability function of Cobb's cusp distribution\cr
    \command{qcusp} \tab Quantile function of Cobb's cusp distribution\cr
    \command{rcusp} \tab Sample from Cobb's cusp distribution.\cr
    \command{cusp.logist} \tab Fit logistic model for bifurctation testing (experimental)\cr
}
}
\author{
Raoul Grasman <rgrasman@uva.nl>
}
\references{
L. Cobb and S. Zacks (1985) \emph{Applications of Catastrophe Theory for Statistical Modeling in the Biosciences (article)}, Journal of the American Statistical Association, 392:793--802.

P. Hartelman (1996). \emph{Stochastic Catastrophy Theory.} Unpublished PhD-thesis.

H. L. J. van der Maas, R. Kolstein, J van der Pligt (2003). \emph{Sudden Transitions in Attitudes}, Sociological Methods \& Research, 32:125-152.

Oliva, DeSarbo, Day, \& Jedidi. (1987) \emph{GEMCAT : A General Multivariate Methodology for Estimating Catastrophe Models}, Behavioral Science, 32:121-137.

R. P. P. P. Grasman, H. L. J. van der Maas, \& E-J. Wagenmakers (2009). \emph{Fitting the Cusp
  Catastrophe in R: A cusp Package Primer.} Journal of Statistical Software 32(8), 1-28. URL
 \url{ http://www.jstatsoft.org/v32/i08/}.
}
\keyword{ package }
%\seealso{
%}
\examples{
# fitting cusp to cusp data
x <- rcusp(100, alpha=0, beta=1) 
fit <- cusp(y ~ x, alpha ~ 1, beta ~ 1)
print(fit)

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

# use of OK
npar <- length(fit$par)
\dontrun{
while(!fit$OK) # refit if necessary until convergence is OK
    fit <- cusp(y ~ z, alpha ~ x1+x2, beta ~ x1+x2, data, start=rnorm(npar))
}

# example 1 from paper
data(attitudes)
data(attitudeStartingValues)
fit.attitudes <- cusp(y ~ Attitude, alpha ~ Orient + Involv, beta ~ Involv, 
 data = attitudes, start=attitudeStartingValues)

summary(fit.attitudes)
plot(fit.attitudes)
cusp3d(fit.attitudes, B = 0.75, Y = 1.35, theta = 170, phi = 30, Yfloor = -9)
}
