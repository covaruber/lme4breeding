\name{condVarRotated}
\alias{condVarRotated}
\title{Get the rotated conditional variance matrix}
\description{
Retrieve the conditional variance matrix from the random effects rotated by the 
Cholesky factors and Eigen factors if applicable.
}
\usage{
condVarRotated(object)
}
\arguments{
 \item{object}{model object of class \code{\link{lmeb}} }
}
\value{
  Matrix
 }
\author{Giovanny Covarrubias}
\examples{

data(DT_example, package="enhancer")
DT <- DT_example
A <- A_example
head(DT)

## Compound simmetry (CS) model
ans1 <- lmeb(Yield~Env + (1|Name) + (1|Env:Name),
                 data=DT)
vc <- VarCorr(ans1); print(vc,comp=c("Variance"))

condVarMat <- condVarRotated(ans1)
# image(condVarMat)

}
\seealso{\code{\link{image}}, \code{\link{lmeb}}}
\keyword{models}

