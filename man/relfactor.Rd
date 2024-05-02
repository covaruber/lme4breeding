\name{relfactor}
\alias{relfactor}
\title{Relationship factor from a pedigree...}
\usage{relfactor(ped, labs)}
\description{Relationship factor from a pedigree}
\details{
  Determine the right Cholesky factor of the relationship matrix for the
  pedigree \code{ped}, possibly restricted to the specific labels that
  occur in \code{labs}.
}
\references{
   2010. A.I. Vazquez, D.M. Bates, G.J.M. Rosa, D. Gianola and K.A. Weigel.
   Technical Note: An R package for fitting generalized linear mixed models
   in animal breeding. Journal of Animal Science, 88:497-504.
}
\value{an upper triangular, sparse (right) Cholesky factor of the relationship matrix}
\arguments{
  \item{ped}{a pedigree that includes the individuals who occur in \code{labs}}
  \item{labs}{a character vector or a factor giving the labels to
    which to restrict the relationship matrix. If \code{labs} is a
    factor then the levels of the factor are used as the labels.
    Default is the complete set of labels in the pedigree.}
}
\examples{

## Example from chapter 2 of Mrode (2005)

ped <- pedigree(sire = c(NA,NA,1, 1,4,5),
                dam  = c(NA,NA,2,NA,3,2), label= 1:6)
(fac <- relfactor(ped))
crossprod(fac)    # the relationship matrix
getA(ped)         # the relationship matrix 
}
\keyword{array}
\keyword{algebra}
