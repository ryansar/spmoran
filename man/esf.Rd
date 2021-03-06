\name{esf}
\alias{esf}
\title{Spatial regression with eigenvector spatial filtering}
\usage{
esf( y, x = NULL, vif = NULL, meig, fn = "r2" )
}
\description{
This function estimates the linear eigenvector spatial filtering (ESF) model. The eigenvectors are selected by a forward stepwise method.
}
\arguments{
\item{y}{Vector of explained variables (N x 1)}
\item{x}{Matrix of explanatory variables (N x K). Default is NULL}
\item{vif}{Maximum acceptable value for the variance inflation factor (VIF) (scalar). For example, if vif = 10, eigenvectors are selected so that the maximum VIF value among explanatory variables and eigenvectors is equal to or less than 10. Default is NULL}
\item{meig}{Moran's eigenvectors and eigenvalues. Output from \code{\link{meigen}} or \code{\link{meigen_f}}}
\item{fn}{Objective function for the stepwise eigenvector selection. The adjusted R2 ("r2"), AIC ("aic"), or BIC ("bic") are available. Alternatively, all the eigenvectors in meig are use if fn = "all". This is acceptable for large samples (see Murakami and Griffith, 2018). Default is "r2"}
}
\value{
\item{b}{Matrix with columns for the estimated coefficients on x, their standard errors, t-values, and p-values (K x 4)}
\item{r}{Matrix with columns for the estimated coefficients on Moran's eigenvectors, their standard errors, t-values, and p-values (L x 4)}
\item{vif}{Vector of variance inflation factors of the explanatory variables (N x 1)}
\item{e}{Vector whose elements are residual standard error (resid_SE), adjusted R2 (adjR2), log-likelihood (logLik), AIC, and BIC}
\item{sf}{Vector of estimated spatial dependent component (E\eqn{\gamma}) (N x 1)}
\item{pred}{Vector of predicted values (N x 1)}
\item{resid}{Vector of residuals (N x 1)}
\item{other}{List of other outcomes, which are internally used}
}
\references{
Tiefelsdorf, M., and Griffith, D. A. (2007). Semiparametric filtering of spatial autocorrelation: the eigenvector approach. Environment and Planning A, 39 (5), 1193-1221.

Murakami, D. and Griffith, D.A. (2018) Eigenvector spatial filtering for large data sets: fixed and random effects approaches. Geographical Analysis, doi: 10.1111/gean.12156.
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
meig_f<- meigen_f(coords=coords)
esfD	<- esf(y=y,x=x,meig=meig_f, vif=10, fn="all")

############################Not run
#########Topoligy-based ESF (it is commonly used in regional science)
#
#cknn	<- knearneigh(coordinates(coords), k=4) #4-nearest neighbors
#cmat	<- nb2mat(knn2nb(cknn), style="B")
#meig <- meigen(cmat=cmat, threshold=0.25)
#esfT	<- esf(y=y,x=x,meig=meig)

}
