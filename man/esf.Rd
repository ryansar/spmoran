\name{esf}
\alias{esf}
\title{Eigenvector spatial filtering}
\usage{
esf( y, x = NULL, vif = NULL, meig, fn = "r2" )
}
\description{
This function estimates the linear eigenvector spatial filtering (ESF) model. Following many of ESF studies, Moran's eigenvectors are selected by a forward stepwise method.
}
\arguments{
\item{y}{Vector of explained variables (\eqn{N} x 1)}
\item{x}{Matrix of explanatory variables (\eqn{N} x \eqn{K}). Default is NULL}
\item{vif}{Maximum permissive value for the variance inflation factor (VIF) (scalar). For example, if \eqn{vif} = 10, eigenvectors are selected so that the maximum VIF value among explanatory variables and eigenvectors does not exceed 10. Default is NULL}
\item{meig}{Moran's eigenvectors and eigenvalues. Output from \code{\link{meigen}} or \code{\link{meigen_f}}}
\item{fn}{Objective function for the stepwise eigenvector selection. The adjusted R2 ("r2"), AIC ("aic"), or BIC ("bic") are available. It is possible to introduce all eigenvectors in meig by setting fn = "all". This is allowable for large samples. Default is "r2"}
}
\value{
\item{b}{Matrix with columns for the estimated coefficients on \eqn{x}, their standard errors, t-values, and p-values (\eqn{K} x 4)}
\item{r}{Matrix with columns for the estimated coefficients on Moran's eigenvectors, their standard errors, t-values, and p-values (\eqn{L} x 4)}
\item{sf}{Vector of estimated spatial dependent component (\eqn{E\gamma}) (\eqn{N} x 1)}
\item{e}{Vector whose elements are residual standard error (resid_SE), adjusted R2 (adjR2), log-likelihood (logLik), AIC, and BIC}
\item{pred}{Vector of predicted values (\eqn{N} x 1)}
\item{resid}{Vector of residuals (\eqn{N} x 1)}
}
\references{
Tiefelsdorf, M., and Griffith, D. A. (2007). Semiparametric filtering of spatial autocorrelation: the eigenvector approach. Environment and Planning A, 39 (5), 1193-1221.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{resf}}
}
\examples{
require(spdep)
data(boston)
y	<- boston.c[, "CMEDV" ]
x	<- boston.c[,c("CRIM","ZN","INDUS", "CHAS", "NOX","RM", "AGE")]
coords  <- boston.c[,c("LAT","LON")]

#########Distance-based ESF
meig 	<- meigen(coords=coords)
esfD	<- esf(y=y,x=x,meig=meig, vif=5)
esfD$vif
esfD$b
esfD$e

#########Fast approximation
meig 	<- meigen_f(coords=coords)
esfD	<- esf(y=y,x=x,meig=meig, vif=10, fn="all")
}
