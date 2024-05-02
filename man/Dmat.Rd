\name{Dmat}
\alias{Dmat}
\title{vector of the diagonal for the D matrix from the decomposition A= TDT'}
\usage{Dmat(ped)}
\description{numeric vector that should be the diagonal elements of the diagonal matrix D}
\details{
  Determine the diagonal factor in the decomposition of the relationship
  matrix from a pedigree equal to TDT'. Where T is unit lower triangular 
  and D is a diagonal matrix.
  This function returns a numeric vector with the entries of  D
}
\value{a numeric vector}
\arguments{\item{ped}{an object that inherits from class \linkS4class{pedigree}}}
\examples{
ped <- pedigree(sire = c(NA,NA,1, 1,4,5),
                dam  = c(NA,NA,2,NA,3,2), label= 1:6)
Dmat(ped)
}
\keyword{misc}
