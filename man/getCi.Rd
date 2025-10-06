\name{getCi}
\alias{getCi}
\title{Get the inverse coefficient matrix}
\description{
Retrieve the inverse of the coefficient matrix from the mixed model equations.
}
\usage{
getCi(object)
}
\arguments{
 \item{object}{model object of class \code{\link{lmebreed}} }
}
\value{
  Matrix
 }
\author{Giovanny Covarrubias}
\examples{

data(DT_example)
DT <- DT_example
A <- A_example
head(DT)

## Compound simmetry (CS) model
ans1 <- lmebreed(Yield~Env + (1|Name) + (1|Env:Name),
                 data=DT)
vc <- VarCorr(ans1); print(vc,comp=c("Variance"))

PEV <- getCi(ans1)
# image(PEV)

}
\seealso{\code{\link{image}}, \code{\link{lmebreed}}}
\keyword{models}

