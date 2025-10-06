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
  data frame
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

