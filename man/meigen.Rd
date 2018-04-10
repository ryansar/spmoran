\name{meigen}
\alias{meigen}
\title{Extraction of Moran's eigenvectors}
\usage{
meigen( coords, model = "exp", threshold = 0, enum = NULL, cmat = NULL )
}
\description{
This function calculates Moran's eigenvectors and their corresponding eigenvalues.
}
\arguments{
\item{coords}{Matrix of spatial point coordinates (\eqn{N} x 2)}
\item{model}{Type of kernel to model spatial dependence. The currently available options are "exp" for the exponential kernel, "gau" for the Gaussian kernel, and "sph" for the spherical kernel}
\item{threshold}{Threshold for eigenvalues (scalar). Suppose that \eqn{\lambda}_1 is the first eigenvalue, eigenvectors whose eigenvalues that are equal or greater than [\eqn{threshold} x \eqn{\lambda}_1] are extracted. The value must lie between 0 and 1. Default is zero (see Details)}
\item{cmat}{Optional. A spatial connectivity matrix (\eqn{N} x \eqn{N}). It must be provided when the user wants to apply a spatial connectivity matrix other than the default exponential decay matrix (e.g. k-nearest neighbor, adjacency; see Detail)}
\item{enum}{Optional. The muxmum acceptable mumber of eigenvectors to be extracted (scalar)}
}
\details{
If \eqn{cmat} is not provided (default), this function extracts Moran's eigenvectors from MCM, where M = I - 11'/\eqn{N} is a centering operator. C is a \eqn{N} x \eqn{N} connectivity matrix whose (\eqn{i}, \eqn{j})-th element equals exp(\eqn{-d_i,j / h}), where \eqn{d_i,j} is the Euclidean distance between sample sites \eqn{i} and \eqn{j}, and \eqn{h} is a range parameter given by the maximum length of the minimum spanning tree connecting sample sites (see Dray et al., 2006). If \eqn{cmat} is provided, this function performs the same calculation after C is replaced with \eqn{cmat}.

If \eqn{threshold} is not provided (default), all eigenvectors corresponding to positive eigenvalues are extracted. It implies to consider all elements describing positive spatial dependence. If \eqn{threshold} is provided, eigenvectors whose corresponding eigenvalues are equal to or greater than [\eqn{threshold} x \eqn{\lambda}_1] are extracted. \eqn{threshold} = 0.00 or 0.25 are standard assumptions (see Griffith, 2003; Murakami and Griffith, 2015).
}
\value{
\item{sf}{Matrix of the first \eqn{L} eigenvectors (\eqn{N} x \eqn{L})}
\item{ev}{Vector of the first \eqn{L} eigenvalues (\eqn{L} x 1)}
\item{ev_full}{Vector of all eigenvalues (\eqn{N} x 1)}
\item{other}{List of other outcomes, which are internally used}
}
\references{
Dray, S., Legendre, P., and Peres-Neto, P.R. (2006) Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbour matrices (PCNM). Ecological Modelling, 196 (3), 483-493.

Griffith, D.A. (2003) Spatial autocorrelation and spatial filtering: gaining understanding through theory and scientific visualization. Springer Science & Business Media.

Murakami, D. and Griffith, D.A. (2015) Random effects specifications in eigenvector spatial filtering: a simulation study. Journal of Geographical Systems, 17 (4), 311-331.

}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{meigen_f}} for fast eigen-decomposition
}
