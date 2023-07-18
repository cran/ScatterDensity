\name{ScatterDensity-package}
\alias{ScatterDensity-package}
\alias{ScatterDensity}
\docType{package}
\title{
\packageTitle{ScatterDensity}
}
\description{
\packageDescription{ScatterDensity}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{ScatterDensity}
\packageIndices{ScatterDensity}
}
\author{
\packageAuthor{ScatterDensity}

Maintainer: \packageMaintainer{ScatterDensity}
}

\keyword{ package }

\references{
[Brinkmann et al., 2023]  Brinkmann, L., Stier, Q., & Thrun, M. C.: Computing Sensitive Color Transitions for the Identification of Two-Dimensional Structures, Proc. Data Science, Statistics & Visualisation (DSSV) and the European Conference on Data Analysis (ECDA), p.109, Antwerp, Belgium, July 5-7, 2023. 

}

\examples{
\donttest{
# Create two bimodial distributions
x1=rnorm(n = 7500,mean = 0,sd = 1)
y1=rnorm(n = 7500,mean = 0,sd = 1)
x2=rnorm(n = 7500,mean = 2.5,sd = 1)
y2=rnorm(n = 7500,mean = 2.5,sd = 1)
x=c(x1,x2)
y=c(y1,y2)

DensityScatter.DDCAL(x, y, Marginals = TRUE)
}
}
