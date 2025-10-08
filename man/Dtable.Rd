\name{Dtable}
\alias{Dtable}
\title{Get the hypertable for predictions}
\description{
For an lmebreed object it creates a template for the hypertable indication.
}
\usage{
Dtable(object)
}
\arguments{
 \item{object}{model object of class \code{\link{lmebreed}} }
}
\value{
  
     A data frame with columns; variable, group, type, include, average.
     
     A pure "include" term means that the model matrices for that fixed or random effect is filled with 1's for the positions where column names and row names match.

     An "include and average" term means that the model matrices for that fixed or random effect is filled with 1/#1's in that row.

     An "average" term alone means that all rows for such fixed or random effect will be filled with 1/#levels in the effect.

     If a term is not considered "include" or "average" is then totally ignored in the BLUP and SE calculation.

     The default behavior when the user doesn't provide the hyperTable is to include and average any fixed effect that is not part of classify. Include any term making a perfect match with the classify argument and include and average any imperfect match with classify argument (e.g., interactions).
  
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
Dtable(ans1)


}
\seealso{\code{\link{image}}, \code{\link{lmebreed}}}
\keyword{models}

