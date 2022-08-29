\name{zeeman}
\alias{zeeman}
\alias{zeeman1}
\alias{zeeman2}
\alias{zeeman3}
\docType{data}
\title{ Measurements from Zeeman's Catastrophe Machine}
\description{ Data sets with measurements from different physical instances of Zeeman's Catastrophe Machine}
\usage{
data(zeeman1)
data(zeeman2)
data(zeeman3)
}
\format{
  A data frame with 150/198/282 observations on the following 3 variables.
  \describe{
    \item{\code{x}}{a control plane variable that are manipulable by the experimentalist}
    \item{\code{y}}{a control plane variable that are manipulable by the experimentalist}
    \item{\code{z}}{the state variable of the machine: the shortest distance to the longitudinal axis of the machine}
  }
}
\details{
  The behavior Zeeman's catastrophe machine is archetypal for the Cusp catastrophe. 
  This device consists of a wheel is tethered by an elastic chord to a fixed point. 
  Another elastic, also attached to the wheel is moved about in the \sQuote{control plane} 
  area opposite to the fixed point. The shortest distance between the strap point on the 
  wheel and the axis defined by the fixed point and the control plane is recorded as a 
  function of the position in the control plane. (In the original machine the angle 
  between this axis and the line through the wheel center and the strap point is used.) 
  See \url{http://www.math.stonybrook.edu/~tony/whatsnew/column/catastrophe-0600/cusp4.html} 
  for a vivid demonstration. These data sets were obtained from 3 different physical 
  instances of this machine, made by different people.

  Measurements were made by systematically sampling different points in the control plane.
  
  See vignette for example analysis with all three data sets.
  
  For pictures of the machines, see
  \describe{
  \item{Zeeman catastrophy machine 1}{\url{https://purl.oclc.org/net/rgrasman/cusp/zeeman1}}
  \item{Zeeman catastrophy machine 2}{\url{https://purl.oclc.org/net/rgrasman/cusp/zeeman2}}
  \item{Zeeman catastrophy machine 3}{\url{https://purl.oclc.org/net/rgrasman/cusp/zeeman3}}
  }
}
\source{
  \describe{
  \item{\code{zeeman1}}{ is due to Noemi Schuurman}
  \item{\code{zeeman2}}{ is due to Karin Visser }
  \item{\code{zeeman3}}{ is due to Mats Nagel & Joris ?}
  }
  See \url{https://sites.google.com/site/zeemanmachine/data-repository}
}
\references{
  Zeeman (1976).
}
\examples{
data(zeeman1)
data(zeeman2)
data(zeeman3)
\dontrun{
fit <- cusp(y~z, alpha~x+y, beta~x+y, data=zeeman1)
plot(fit)
cusp3d(fit, surf.hue = 40, theta=215, phi=37.5, B=5.25)
}
}
\keyword{datasets}
