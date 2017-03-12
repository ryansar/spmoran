\name{resf}
\alias{resf}
\title{Random effects eigenvector spatial filtering}
\usage{
resf( y, x = NULL, meig, method = "reml" )
}
\description{
This function estimates the random effects eigenvector spatial filtering (RE-ESF) model.
}
\arguments{
\item{y}{Vector of explained variables (\eqn{N} x 1)}
\item{x}{Matrix of explanatory variables (\eqn{N} x K). Default is NULL}
\item{meig}{Moran's eigenvectors and eigenvalues. Output from \code{\link{meigen}} or \code{\link{meigen_f}}}
\item{method}{Estimation method. Restricted maximum likelihood method ("reml") and maximum likelihood method ("ml") are available. Default is "reml"}
}
\value{
\item{b}{Matrix with columns for the estimated coefficients on \eqn{x}, their standard errors, t-values, and p-values (\eqn{K} x 4)}
\item{r}{Vector of estimated random coefficients on Moran's eigenvectors (\eqn{L} x 1)}
\item{sf}{Vector of estimated spatial dependent component (\eqn{N} x 1)}
\item{s}{Vector of estimated shrinkage parameters (2 x 1). The first and the second elements denote the standard error and the spatial scale of the estimated spatial dependent component, respectively (see Murakami and Griffith, 2015)}
\item{e}{Vector whose elements are residual standard error (resid_SE), adjusted conditional R2 (adjR2(cond)), restricted log-likelihood (rlogLik), Akaike information criterion (AIC), and Bayesian information criterion (BIC). When method = "ml", restricted log-likelihood (rlogLik) is replaced with log-likelihood (logLik)}
\item{pred}{Vector of predicted values (\eqn{N} x 1)}
\item{resid}{Vector of residuals (\eqn{N} x 1)}
}
\references{
Murakami, D. and Griffith, D.A. (2015) Random effects specifications in eigenvector spatial filtering: a simulation study. Journal of Geographical Systems, 17 (4), 311-331.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{meigen}}, \code{\link{meigen_f}}
}
\examples{
require(spdep)
data(boston)
y	<- boston.c[, "CMEDV" ]
x	<- boston.c[,c("CRIM","ZN","INDUS", "CHAS", "NOX","RM", "AGE",
                       "DIS" ,"RAD", "TAX", "PTRATIO", "B", "LSTAT")]
coords  <- boston.c[,c("LAT","LON")]
meig 	<- meigen(coords=coords)
res	<- resf(y=y,x=x,meig=meig, method = "reml")
res$b
res$s
res$e

#########Fast approximation
meig 	<- meigen_f(coords=coords)
res	<- resf(y=y,x=x,meig=meig, method = "reml")
}

