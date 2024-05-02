\name{inbreeding}
\alias{inbreeding}
\title{Inbreeding coefficients from a pedigree...}
\usage{inbreeding(ped)}
\description{Inbreeding coefficients from a pedigree}
\details{
  Determine the inbreeding coefficients for all the individuals of a pedigree. 
  This function a numeric vector. 
}
\arguments{\item{ped}{an object that inherits from class \linkS4class{pedigree}}}
\value{a numeric vector}
\references{
   2010. A.I. Vazquez, D.M. Bates, G.J.M. Rosa, D. Gianola and K.A. Weigel.
   Technical Note: An R package for fitting generalized linear mixed models
   in animal breeding. Journal of Animal Science, 88:497-504.
}
\examples{
ped <- pedigree(sire = c(NA,NA,1, 1,4,5),
                dam  = c(NA,NA,2,NA,3,2), label= 1:6)
inbreeding(ped)
}
\source{
  Sargolzaei, M. and H. Iwaisaki, 2005. Comparison of four direct algorithms
  for computing the inbreeding coefficients. J. Anim. Sci, 76: 401-406. 
}
\keyword{misc}
