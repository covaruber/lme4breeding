\name{getAInv}
\alias{getAInv}
\title{Inverse of the relationship matrix}
\usage{getAInv(ped)}
\description{Inverse of the Relationship matrix from a pedigree}
\details{
  Determine the inverse of the relationship matrix from a pedigree
   \code{ped}.
}
\references{
   2010. A.I. Vazquez, D.M. Bates, G.J.M. Rosa, D. Gianola and K.A. Weigel.
   Technical Note: An R package for fitting generalized linear mixed models
   in animal breeding. Journal of Animal Science, 88:497-504.
}
\value{sparse matrix, inverse of the relationship matrix}
\arguments{
  \item{ped}{a pedigree that includes the individuals who occur in \code{labs}}
}
\examples{

## Example from chapter 2 of Mrode (2005)

ped <- pedigree(sire = c(NA,NA,1, 1,4,5),
                dam  = c(NA,NA,2,NA,3,2), label= 1:6)
getAInv(ped)
}
\keyword{array}
\keyword{algebra}
