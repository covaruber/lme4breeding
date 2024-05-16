\name{leg}
\alias{leg}
\title{Legendre polynomial matrix}
\description{
  Legendre polynomials of order 'n' are created given a vector 'x' and normalized to lay between values u and v.
}
\usage{
  leg(x,n=1,u=-1,v=1, intercept=TRUE, intercept1=FALSE)
}

\arguments{
  \item{x}{numeric vector to be used for the polynomial.}
  \item{n}{order of the Legendre polynomials.}
  \item{u}{lower bound for the polynomial.}
  \item{v}{upper bound for the polynomial.}
  \item{intercept}{a TRUE/FALSE value indicating if the intercept should be included.}
  \item{intercept1}{a TRUE/FALSE value indicating if the intercept should have value 1 (is multiplied by sqrt(2)).}
}
\value{
  \describe{
    
    \item{$S3}{ an Legendre polynomial matrix of order n.}
  }
}
\references{
 
Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
x <- sort(rep(1:3,100))

\donttest{

library(orthopolynom)
 leg(x, n=1)
 leg(x, n=2)
}
## see dataset data(DT_legendre) for a random regression modeling example

}
\seealso{The core function of the package \code{\link{lmebreed}}}