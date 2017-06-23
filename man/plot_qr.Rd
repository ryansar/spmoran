\name{plot_qr}
\alias{plot_qr}
\title{Plot quantile regression coefficients estimated by SF-UQR}
\usage{
plot_qr( mod, pnum = 1, par = "b", cex.main = 28, cex.lab = 26, cex.axis = 24, lwd = 1.5 )
}
\description{
This function plots regression coefficients estimated by the spatial filter unconditional quantile regression (SF-UQR) approach.
}
\arguments{
\item{mod}{Object produced by the \code{\link{resf_qr}} function}
\item{pnum}{A number specifying which parameter is plotted. If \eqn{par} = "b", coefficients on \eqn{pnum}-th explanatory variable are plotted (intercepts are plotted if \eqn{pnum} = 1). On the other hand, standard error for the reidual spatial process is plotted if \eqn{par} = "s" and \eqn{pnum} = 1, whereas the scale (degree) parameter for the process is plotted when \eqn{par} = "s" and \eqn{pnum} = 2}
\item{par}{If it is "b", regression coefficeints are plotted. If it is "s", shrinkage (variance) parameters for the residual spatial process are plotted. Default is "b"}
\item{cex.main}{Graphical parameter specifying the size of the main title}
\item{cex.lab}{Graphical parameter specifying the size of the x and y axis labels}
\item{cex.axis}{Graphical parameter specifying the size of the tick label numbers}
\item{lwd}{Graphical parameters specifying the width of the line drawing the coefficient estimates}
}
\note{
See \code{\link{par}} for the graphical parameters
}
\seealso{
\code{\link{resf_qr}}
}
