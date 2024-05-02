\name{getA}
\alias{getA}
\title{Additive Relationship Matrix}
\usage{getA(ped)}
\description{Additive relationship matrix from a pedigree}
\details{
  Returns the additive relationship matrix for the pedigree \code{ped} .
}
%\references{}
\value{Sparse matrix}
\arguments{
  \item{ped}{a pedigree that includes the individuals who occur in \code{labs}}
}
\examples{

## Example from chapter 2 of Mrode (2005)

ped <- pedigree(sire = c(NA,NA,1, 1,4,5),
                dam  = c(NA,NA,2,NA,3,2), label= 1:6)
(getA(ped))
}
\keyword{array}
\keyword{algebra}
