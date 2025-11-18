\name{mkMmeIndex}
\alias{mkMmeIndex}
\title{Table to track effects}
\description{
For an lmeb object it creates a table to track the indices of each effect for
a given slope and intercept.
}
\usage{
mkMmeIndex(object)
}
\arguments{
 \item{object}{model object of class \code{\link{lmeb}} }
}
\value{
  data frame with 4 columns: 
  
  1) index: indicates the position in the effects and PEV matrix (\code{\link{condVarRotated}})
  
  2) level: the name of the level for a given effect
  
  3) variable: the intercept effect in the model
  
  4) group: the slope effect in the model
  
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

head(mkMmeIndex(ans1))


}
\seealso{\code{\link{condVarRotated}}, \code{\link{lmeb}}}
\keyword{models}

