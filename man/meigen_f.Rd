\name{meigen_f}
\alias{meigen_f}
\title{Fast approximation of Moran's eigenvectors}
\usage{
meigen_f( coords, model = "exp", enum = 200 )
}
\description{
This function performs a fast approximation of Moran's eigenvectors and their corresponding eigenvalues.
}
\arguments{
\item{coords}{Matrix of spatial point coordinates (\eqn{N} x 2)}
\item{model}{Type of kernel to model spatial dependence. The currently available options are "exp" for the exponential kernel, "gau" for the Gaussian kernel, and "sph" for the spherical kernel}
\item{enum}{Number of eigenvectors and eigenvalues to be extracted (scalar). Default is 200}
}
\details{
This function extracts approximated Moran's eigenvectors from MCM. M = I - 11'/\eqn{N} is a centering operator, and C is a spatial connectivity matrix whose (\eqn{i}, \eqn{j})-th element is given by exp( -\eqn{d_i,j} / \eqn{h} ), where \eqn{d_i,j} is the Euclidean distance between sample sites \eqn{i} and \eqn{j}, and \eqn{h} is a range parameter given by the maximum length of the minimum spanning tree connecting sample sites (see Dray et al., 2006).

Following a simulation result that 200 eigenvectors are sufficient for accurate approximation of ESF models (Murakami and Griffith, 2017), this function approximates the first 200 eigenvectors by default (i.e., \eqn{enum} = 200). If \eqn{enum} is given by a smaller value like 100, the computation time will be shorter, but with greater approximation error. Following \code{\link{meigen}}, eigenvectors corresponding to negative eigenvalues are omitted among the 200 eigenvectors.
}
\value{
\item{sf}{Matrix of the first \eqn{L} approximated eigenvectors (\eqn{N} x \eqn{L})}
\item{ev}{Vector of the first \eqn{L} approximated eigenvalues (\eqn{L} x 1)}
\item{ev_full}{Vector of all approximated eigenvalues (\eqn{enum} x 1)}
\item{other}{List of other outcomes, which are internally used}
}
\references{
Dray, S., Legendre, P., and Peres-Neto, P.R. (2006) Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbour matrices (PCNM). Ecological Modelling, 196 (3), 483-493.

Murakami, D. and Griffith, D.A. (2018) Eigenvector spatial filtering for large data sets: fixed and random effects approaches. Geographical Analysis, doi: 10.1111/gean.12156.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{meigen}}
}
