\name{simage2}
\alias{simage2}
\title{Shortcut to sparse matrix image}
\description{\code{image} for sparsity.
}
\usage{
simage2(X, \dots)
}
\arguments{
 \item{X}{object of class \code{"matrix"}}
 \item{\dots}{Further arguments to be passed to the image function.}
}
\value{
  vector of image
 }
\author{Giovanny Covarrubias}
\examples{

simage2(matrix(0:8,3,3))

}
\seealso{\code{\link{image}}, \code{\link{lmeb}}}
\keyword{models}

