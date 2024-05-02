\name{mastitis}
\alias{mastitis}
\docType{data}
\title{Mastitis cases in dairy cattle}
\description{
  Records of the number of cases of clinical mastitis during the first
  lactation of 1,675 cows, primarily Holsteins.
  Cows belonged to 41 herds and were daughters of 38 sires. 
  There were 1,491 healthy cows, 134 had only one case of mastitis, 
  36 had 2 cases, and 14 had between 4 and  cases. 
  Overall, mastitis incidence was 0.11. 
  Calving years for these records were from 2000 through 2005. 
  The sire, herd and days in milk are also recorded for each
  cow.
}
\format{
  A data frame with 1675 observations on the following 8 variables.
  \describe{
    \item{\code{id}}{Identifier of the animal.}
    \item{\code{sire}}{Identifier of the animal's sire.}
    \item{\code{birth}}{year of birth of the animal (as a factor).}
    \item{\code{herd}}{herd id number (as a factor).}
    \item{\code{calvingYear}}{year of calving for this lactation.}
    \item{\code{DIM}}{total number of days in milk for the lactation.}
    \item{\code{mastitis}}{a factor indicating if the cow had any
      incidents of clinical mastitis during the lactation.}
    \item{\code{NCM}}{An ordered factor giving the number of clinical
      mastitis cases for the cow during this lactation.}
  }
}
\details{The pedigree of the sires is given in the companion
  \code{\link{pedSires}} data set.}
\source{
    Vazquez, A.I. 2007. Analysis of number of episodes of clinical mastitis
    in Norwegian Red and Holstein cows with Poisson and categorical data mixed models. 
    Master of Science Thesis. University of Wisconsin - Madison. 162 pp.

   2010. A.I. Vazquez, D.M. Bates, G.J.M. Rosa, D. Gianola and K.A. Weigel.
   Technical Note: An R package for fitting generalized linear mixed models
   in animal breeding. Journal of Animal Science, 88:497-504.
}
\references{
   2010. A.I. Vazquez, D.M. Bates, G.J.M. Rosa, D. Gianola and K.A. Weigel.
   Technical Note: An R package for fitting generalized linear mixed models
   in animal breeding. Journal of Animal Science, 88:497-504.
}
\seealso{\code{\link{pedSires}}, \code{\link{pedigree}}}
\examples{
str(mastitis)
summary(mastitis, maxsum = 10)
}
\keyword{datasets}
