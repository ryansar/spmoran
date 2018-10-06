\name{resf_qr}
\alias{resf_qr}
\title{Spatial filter unconditional quantile regression}
\usage{
resf_qr( y, x = NULL, meig, tau = NULL, boot = TRUE, iter = 200 )
}
\description{
This function estimates the spatial filter unconditional quantile regression (SF-UQR) model.
}
\arguments{
\item{y}{Vector of explained variables (N x 1)}
\item{x}{Matrix of explanatory variables (N x K). Default is NULL}
\item{meig}{Moran's eigenvectors and eigenvalues. Output from \code{\link{meigen}} or \code{\link{meigen_f}}}
\item{tau}{The quantile(s) to be modeled. It must be a number (or a vector of numbers) strictly between 0 and 1. By default, \eqn{tau} = c(0.1, 0.2, ..., 0.9)}
\item{boot}{If it is TRUE, confidence intervals for regression coefficients are estimated by a semiparametric bootstrapping. Default is TRUE}
\item{iter}{The number of bootstrap replications. Default is 200}
}
\value{
\item{b}{Matrix of estimated regression coefficients (K x Q), where Q is the number of quantiles (i.e., the length of tau)}
\item{r}{Matrix of estimated random coefficients on Moran's eigenvectors (L x Q)}
\item{s}{Vector of estimated shrinkage parameters (2 x Q). The first element denotes the standard error of the spatially dependent component (shrink_sf_SE), and the second represents the spatial scale of the component (shrink_sf_alpha) (see Murakami and Griffith, 2015)}
\item{e}{Vector whose elements are residual standard error (resid_SE) and adjusted quasi conditional R2 (quasi_adjR2(cond))}
\item{B}{Q matrices (K x 4) summarizing bootstrapped estimates for the regression coefficients Columns of these matrices consist of the estimated coefficients, the lower and upper bounds for the 95 percent confidencial intervals, and p-values. It is returned if boot = TRUE}
\item{S}{Q matrices (2 x 3) summarizing bootstrapped estimates for the shrinkage parameters. Columns of these matrices consist of the estimated parameters, the lower and upper bounds for the 95 percent confidencial intervals. It is returned if boot = TRUE}
\item{B0}{List of Q matrices (K x iter) summarizing bootstrapped coefficients. The q-th matrix consists of the coefficients on the q-th quantile. It is returned if boot = TRUE}
\item{S0}{List of Q matrices (2 x iter) summarizing bootstrapped shrinkage parameters. The q-th matrix consists of the parameters on the q-th quantile. It is returned if boot = TRUE}
}
\seealso{
\code{\link{plot_qr}}
}
\references{
Murakami, D. and Seya, H. (2017) Spatially filtered unconditional quantile regression. ArXiv.
}
\author{
Daisuke Murakami
}
\examples{
require(spdep)
data(boston)
y	<- boston.c[, "CMEDV" ]
x	<- boston.c[,c("CRIM","ZN","INDUS", "CHAS", "NOX","RM", "AGE",
                       "DIS" ,"RAD", "TAX", "PTRATIO", "B", "LSTAT")]
coords  <- boston.c[,c("LAT","LON")]
meig 	<- meigen(coords=coords)
res	<- resf_qr(y=y,x=x,meig=meig, boot=FALSE)
res$b
res$s
res$e
plot_qr(res,1) #The first explanatory vafriable (intercept)
plot_qr(res,2) #The second explanatory variable
plot_qr(res,1,"s") #shrink_sf_SE
plot_qr(res,2,"s") #shrink_sf_alpha

###Not run
#res	<- resf_qr(y=y,x=x,meig=meig, boot=TRUE)
#res$B
#res$S
#res$e
#plot_qr(res,1) #The first explanatory vafriable (intercept)
#plot_qr(res,2) #The second explanatory variable
#plot_qr(res,1,"s") #shrink_sf_SE
#plot_qr(res,2,"s") #shrink_sf_alpha
}

