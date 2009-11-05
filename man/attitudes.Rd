\name{attitudes}
\alias{attitudes}
\alias{attitudeStartingValues}
\docType{data}
\title{ Multistability in political attitudes  }
\description{ Data set reflecting bistability in political attitudes }
\usage{data(attitudes)
data(attitudeStartingValues)}
\format{
  A data frame with 1387 observations on the following 3 variables.
  \describe{
    \item{\code{Orient}}{a numeric vector}
    \item{\code{Involv}}{a numeric vector}
    \item{\code{Attitude}}{a numeric vector}
  }

  The format of attitudeStartingValues is:
 num [1:7]  0.153 -0.453 -0.097 -0.124 -0.227 ...
}
\details{
The data set was taken from (van der Maas, Kolstein, & van der Pligt, 2003). It concerns attitudinal response 
transitions with respect to the statement \dQuote{The government must force companies to let their 
workers benefit from the profit as much as the shareholders do}. Responses of some 1387 Dutch respondents are included who 
indicated their level of agreement with this statement on a 5 point scale (1 = total ly agree, 
5 = total ly disagree). As a normal factor political orientation (measures on a 10 point scale 
from 1 = left wing to 10 = right wing) was used. As a bifurcation factor the total score on a 
12 item political involvement scale was used. The theoretical social psychological details are 
discussed in (van der Maas et al. 2003). 

  The starting values provided here for a cusp analysis of the \code{attitude} data set give proper convergence in one run. They were found after many trial starting values that yielded improper convergence.
}
\source{
    van der Maas HLJ, Kolstein R, van der Pligt J (2003). Sudden Transitions in Attitudes. 
    Sociological Methods & Research, 23(2), 125152.
}
\references{
    van der Maas HLJ, Kolstein R, van der Pligt J (2003). Sudden Transitions in Attitudes. 
    Sociological Methods & Research, 23(2), 125152.
}
\examples{
data(attitudes)
data(attitudeStartingValues)
\dontrun{
fit <- cusp(y ~ Attitude, 
			alpha ~ Orient + Involv, 
			beta ~ Involv, 
			data = attitudes, start=attitudeStartingValues) 
}
## maybe str(attitudeStartingValues) ; plot(attitudeStartingValues) ...
}
\keyword{datasets}
