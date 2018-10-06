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
\item{coords}{Matrix of spatial point coordinates (N x 2)}
\item{model}{Type of kernel to model spatial dependence. The currently available options are "exp" for the exponential kernel, "gau" for the Gaussian kernel, and "sph" for the spherical kernel}
\item{threshold}{Threshold for the eigenvalues (scalar). Suppose that \eqn{\lambda}_1 is the maximum eigenvalue, eigenvectors whose corresponding eigenvalues that are equal or greater than [threshold x \eqn{\lambda}_1] are extracted. threshold must be a value between 0 and 1. Default is zero (see Details)}
\item{cmat}{Optional. A user-specified spatial connectivity matrix (N x N). It must be provided when the user wants to use a spatial connectivity matrix other than the default matrices}
\item{enum}{Optional. The muxmum acceptable mumber of eigenvectors to be extracted (scalar)}
}
\details{
If cmat is not provided and model = "exp" (default), this function extracts Moran's eigenvectors from MCM, where M = I - 11'/N is a centering operator. C is a N x N connectivity matrix whose (i, j)-th element equals exp(-d(i,j)/h), where d(i,j) is the Euclidean distance between the sample sites i and j, and h is given by the maximum length of the minimum spanning tree connecting sample sites (see Dray et al., 2006). If cmat is provided, this function performs the same calculation after C is replaced with cmat.

If threshold is not provided (default), all eigenvectors corresponding to positive eigenvalues are extracted. It implies to consider all the elements describing positive spatial dependence. If threshold is provided, eigenvectors whose corresponding eigenvalues are equal to or greater than [threshold x \eqn{\lambda}_1] are extracted. threshold = 0.00 or 0.25 are standard assumptions (see Griffith, 2003; Murakami and Griffith, 2015).
}
\value{
\item{sf}{Matrix of the first L eigenvectors (N x L)}
\item{ev}{Vector of the first L eigenvalues (L x 1)}
\item{ev_full}{Vector of all eigenvalues (N x 1)}
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
