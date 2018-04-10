\name{meigen0}
\alias{meigen0}
\title{Nystrom extension of Moran's eigenvectors}
\usage{
meigen0( meig, coords0 )
}
\description{
This function estimates Moran's eigenvectors at predicted sites using the Nystrom extension.
}
\arguments{
\item{coords0}{Matrix of spatial point coordinates of predicted sites (\eqn{N_0} x 2)}
\item{meig}{Moran's eigenvectors and eigenvalues. Output from \code{\link{meigen}} or \code{\link{meigen_f}}}
}
\value{
\item{sf}{Matrix of the first \eqn{L} eigenvectors at predicted sites (\eqn{N_0} x \eqn{L})}
\item{ev}{Vector of the first \eqn{L} eigenvalues (\eqn{L} x 1)}
\item{ev_full}{Vector of all eigenvalues (\eqn{N} x 1)}
}
\references{
Drineas, P. and Mahoney, M.W. (2005) On the Nystrom method for approximating a gram matrix for improved kernel-based learning. Journal of Machine Learning Research, 6 (2005), 2153-2175.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{meigen}}, \code{\link{meigen_f}}
}
