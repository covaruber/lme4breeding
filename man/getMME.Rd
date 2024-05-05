\name{getMME}
\alias{getMME}

\title{
Build the mixed model equations from a lmebreed object
}
\description{
Uses the lmebreed object and builds the coefficient matrix (C) and returns its inverse and the solutions to the mixed model equations.
}
\usage{
getMME(object, vc)
}

\arguments{
  \item{object}{
an object of class lmebreed.
}

  \item{vc}{
optional variance components to force in the mixed model equations. This this to be the outpur of the \code{\link{VarCorr}} function.
}

}
\details{

Uses the lmebreed object and builds the coefficient matrix (C) and returns its inverse and the solutions to the mixed model equations. It is internally used by the ranef function when the user wants standard errors for the BLUPs.

}
\value{

\describe{
\item{$Ci}{inverse of the coefficient matrix.}
\item{$bu}{solutions to the mixed model equations}
}

}

\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}

\examples{

data(DT_example)
DT <- DT_example
A <- A_example

ansMain <- lmebreed(Yield ~ Env + (1|Name),
                    relmat = list(Name = Matrix::chol(A) ),
                    data=DT)
mme <- getMME(ansMain)

}
\seealso{\code{\link{lmebreed}}-- the core function of the package }
