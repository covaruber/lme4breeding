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

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
 
}
\examples{
# ####=========================================####
# #### For CRAN time limitations most lines in the 
# #### examples are silenced with one '#' mark, 
# #### remove them and run the examples using
# #### command + shift + C |OR| control + shift + C
# ####=========================================####
# # you need to install the orthopolynom library to do random regression models
# library(orthopolynom)
data(DT_legendre)
# DT <- DT_legendre
# 
# Z <- with(DT, dsc(leg(X,1)) )$Z
# # diagonal random regression
# form <- paste0( "Y ~ Xf + (0+", paste(colnames(Z), collapse = "+"), "|| SUBJECT)")
# # unstructured random regression
# form <- paste0( "Y ~ Xf + (0+", paste(colnames(Z), collapse = "+"), "| SUBJECT)")
# for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
# mRR2b<-lmebreed(as.formula(form),
#            , data=DT)
# vc <- VarCorr(mRR2b); print(vc,comp=c("Variance"))


}
\keyword{datasets}
\seealso{The core function of the package \code{\link{lmebreed}} }
