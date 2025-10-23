\name{smm}
\alias{smm}
\title{sparse model matrix}
\description{
  \code{smm} creates a sparse model matrix for the levels of the random effect  to be used with the \code{\link{lmeb}} solver.
}
\usage{
  smm(x)
}

\arguments{
  \item{x}{vector of observations for the random effect.}
}
\value{
  \describe{
    \item{$res}{ a model matrix for a given factor.}
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
x <- as.factor(c(1:5,1:5,1:5));x
smm(x)

}
\seealso{The \code{\link{lmeb}} solver.}