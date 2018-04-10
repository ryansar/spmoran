\name{predict0_vc}
\alias{predict0_vc}
\title{Prediction of explained variables and spatially varying coefficients}
\usage{
predict0_vc( mod, meig0, x0 = NULL, xconst0 = NULL )
}
\description{
This function predicts explained variables and spatially varying coefficients using a spatially varying coefficient model, which is based on the random effects ESF. The Nystrom extension is used to minimize the expected prediction error
}
\arguments{
\item{mod}{spatially varying coefficient model estimates. Output from \code{\link{resf_vc}}}
\item{meig0}{Moran's eigenvectors at predicted sites. Output from \code{\link{meigen0}}}
\item{x0}{Matrix of explanatory variables at predicted sites whose coefficients are allowed to vary across geographical space (\eqn{N_0} x \eqn{K}). Default is NULL}
\item{xconst0}{Matrix of explanatory variables at predicted sites whose coefficients are assumed constant across space (\eqn{N_0} x \eqn{K_const}). Default is NULL}
}
\value{
\item{pred}{Matrix with the first column for the predicted values. The second and the third columns are the trend component (i.e., component explained by x0 and xconst0) and the spatial component in the predicted values (\eqn{N} x 3)}
\item{b_vc}{Matrix of estimated spatially varying coefficients (SVCs) on x0 (\eqn{N_0} x \eqn{K})}
\item{bse_vc}{Matrix of estimated standard errors for the SVCs (\eqn{N_0} x \eqn{k})}
\item{t_vc}{Matrix of estimated t-values for the SVCs (\eqn{N_0} x \eqn{K})}
\item{p_vc}{Matrix of estimated p-values for the SVCs (\eqn{N_0} x \eqn{K})}

}

\references{
Drineas, P. and Mahoney, M.W. (2005) On the Nystrom method for approximating a gram matrix for improved kernel-based learning. Journal of Machine Learning Research, 6 (2005), 2153-2175.

Murakami, D., Yoshida, T., Seya, H., Griffith, D.A., and Yamagata, Y. (2017) A Moran coefficient-based mixed effects approach to investigate spatially varying relationships. Spatial Statistics, 19, 68-89.
}
\seealso{
\code{\link{meigen0}}, \code{\link{predict0}}
}

\examples{
require(spdep)
data(boston)
samp    <- sample( dim( boston.c )[ 1 ], 300)

d       <- boston.c[ samp, ]    ## Data at observed sites
y	      <- d[, "CMEDV"]
x       <- d[,c("CRIM", "ZN", "INDUS", "RM" ,"LSTAT")]
xconst  <- d[,c("NOX", "CHAS", "AGE", "DIS", "RAD", "TAX", "PTRATIO", "B" )]
coords  <- d[,c("LAT","LON")]

d0      <- boston.c[-samp, ]    ## Data at unobserved sites
x0      <- d0[,c("CRIM", "ZN", "INDUS", "RM" ,"LSTAT")]
xconst0 <- d0[,c("NOX", "CHAS", "AGE", "DIS", "RAD", "TAX", "PTRATIO", "B" )]
coords0 <- d0[,c("LAT","LON")]

############ Model estimation
meig 	  <- meigen( coords = coords )
mod	    <- resf_vc(y=y, x=x, xconst=xconst, meig=meig)

############ Spatial prediction of y and spatially varying coefficients
meig0 	<- meigen0( meig = meig, coords0 = coords0 )
pred0   <- predict0_vc( mod = mod, x0 = x0, xconst0=xconst0, meig0 = meig0 )
pred0$pred[1:10,]
pred0$b_vc[1:10,]
pred0$bse_vc[1:10,]
pred0$t_vc[1:10,]
pred0$p_vc[1:10,]

############ or spatial prediction of spatially varying coefficients
pred00  <- predict0_vc( mod = mod, meig0 = meig0 )
pred00$b_vc[1:10,]
pred00$bse_vc[1:10,]
pred00$t_vc[1:10,]
pred00$p_vc[1:10,]
}
