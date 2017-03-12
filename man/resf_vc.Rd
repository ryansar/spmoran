\name{resf_vc}
\alias{resf_vc}
\title{RE-ESF with spatially varying coefficients}
\usage{
resf_vc( y, x = NULL, xconst = NULL, se = TRUE, meig, method = "reml" )
}
\description{
This function estimates the random effects eigenvector spatial filtering (RE-ESF) model with varying coefficients.
}
\arguments{
\item{y}{Vector of explained variables (\eqn{N} x 1)}
\item{x}{Matrix of explanatory variables whose coefficients are allowed to vary across geographical space (\eqn{N} x \eqn{k}). Default is NULL}
\item{xconst}{Matrix of explanatory variables whose coefficients are assumed constant across space (\eqn{N} x \eqn{K_const}). Default is NULL}
\item{meig}{Moran's eigenvectors and eigenvalues. Output from \code{\link{meigen}} or \code{\link{meigen_f}}}
\item{se}{If it is TRUE, standard errors, t-values, and p-value sfor spatially varying coefficients are calculated. The computation can be slow for large samples. Default is TRUE}
\item{method}{Estimation method. Restricted maximum likelihood method ("reml") and maximum likelihood method ("ml") are available. Default is "reml"}
}
\value{
\item{b}{Matrix with columns for the estimated coefficients on \eqn{xconst}, their standard errors, t-values, and p-values (\eqn{K_const} x 4)}
\item{b_vc}{Matrix of estimated spatially varying coefficients (SVCs) on \eqn{x} (\eqn{N} x \eqn{K})}
\item{bse_vc}{Matrix of estimated standard errors for the SVCs (\eqn{N} x \eqn{k})}
\item{t_vc}{Matrix of estimated t-values for the SVCs (\eqn{N} x \eqn{K})}
\item{p_vc}{Matrix of estimated p-values for the SVCs (\eqn{N} x \eqn{K})}
\item{s}{Matrix of estimated shrinkage parameters (2 x \eqn{K}). The (1, \eqn{k})-th element denotes the standard error of the \eqn{k}-th SVCs, while the (2, \eqn{k})-th element denotes the spatial scale of the same SVCs (see Murakami et al., 2016)}
\item{e}{Vector whose elements are residual standard error (resid_SE), adjusted conditional R2 (adjR2(cond)), restricted log-likelihood (rlogLik), Akaike information criterion (AIC), and Bayesian information criterion (BIC). When method = "ml", restricted log-likelihood (rlogLik) is replaced with log-likelihood (logLik)}
\item{pred}{Vector of predicted values (\eqn{N} x 1)}
\item{resid}{Vector of residuals (\eqn{N} x 1)}
}
\note{
Spatially varying coefficients estimates can be unstable if too many explanatory variables are included in \eqn{x} (see Wheeler and Tiefelsdorf, 2005). In case with many explanatory variables, it is preferable to include at most four or five explanatory variables for \eqn{x}, and the other explanatory variables for \eqn{xconst}.
}
\references{
Murakami, D., Yoshida, T., Seya, H., Griffith, D.A., and Yamagata, Y. (2016) A Moran coefficient-based mixed effects approach to investigate spatially varying relationships. Spatial Statistics, 19, 68-89.

Wheeler, D. and Tiefelsdorf, M. (2005) Multicollinearity and correlation among local regression coefficients in geographically weighted regression. Journal of Geographical Systems, 7 (2), 161-187.
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
y	<- boston.c[, "CMEDV"]
x	<- boston.c[, "NOX"  ]#Coefficients on "NOX" are allowed to vary spatially
xconst  <- boston.c[,c("ZN","INDUS", "CHAS", "RM", "DIS", "AGE", "CRIM",
                       "RAD", "TAX", "PTRATIO", "B", "LSTAT")]
coords  <- boston.c[,c("LAT","LON")]
meig 	<- meigen(coords=coords)
res	<- resf_vc(y=y,x=x,xconst=xconst,meig=meig, method ="reml")
res$b
res$s
res$e
res$b_vc[1:10,]
res$bse_vc[1:10,]
res$t_vc[1:10,]
res$p_vc[1:10,]

#########Fast approximation
meig 	<- meigen_f(coords=coords)
res	<- resf_vc(y=y,x=x,xconst=xconst, meig=meig, method ="reml", se = FALSE)
}
