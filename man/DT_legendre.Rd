\name{DT_legendre}
\alias{DT_legendre}
\docType{data}
\title{
Simulated data for random regression
}
\description{
A data frame with 4 columns; SUBJECT, X, Xf and Y to show how to use the Legendre polynomials in the lmebreed function using a numeric variable X and a response variable Y. 

}
\usage{data("DT_legendre")}
\format{
  The format is:
 chr "DT_legendre"
}
\source{
This data was simulated for fruit breeding applications.
}
\references{

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}
\examples{

data(DT_legendre)
DT <- DT_legendre
head(DT)

\donttest{

library(orthopolynom)
Z <- with(DT, smm(leg(X,1)) )
## diagonal random regression
form <- paste0( "Y ~ Xf + (0+", paste(colnames(Z), collapse = "+"), "|| SUBJECT)")
## unstructured random regression
form <- paste0( "Y ~ Xf + (0+", paste(colnames(Z), collapse = "+"), "| SUBJECT)")
for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
mRR2b<-lmebreed(as.formula(form),
           , data=DT)
vc <- VarCorr(mRR2b); print(vc,comp=c("Variance"))
sigma(mRR2b)^2 # error variance

}


}
\keyword{datasets}
\seealso{The core function of the package \code{\link{lmebreed}} }
