\name{condVarRotation}
\alias{condVarRotation}
\title{Get the inverse coefficient matrix}
\description{
Retrieve the conditional variance matrix from the random effects rotated by the 
Cholesky factors and Eigen factors if applicable.
}
\usage{
condVarRotation(object)
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

condVarMat <- condVarRotation(ans1)
# image(condVarMat)

}
\seealso{\code{\link{image}}, \code{\link{lmebreed}}}
\keyword{models}

