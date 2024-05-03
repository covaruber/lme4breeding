\name{DT_cpdata}
\alias{DT_cpdata}
\alias{GT_cpdata}
\alias{MP_cpdata}
\docType{data}
\title{
Genotypic and Phenotypic data for a CP population
}
\description{
A CP population or F1 cross is the designation for a cross between 2 highly heterozygote individuals; i.e. humans, fruit crops, bredding populations in recurrent selection.

This dataset contains phenotpic data for 363 siblings for an F1 cross. These are averages over 2 environments evaluated for 4 traits; color, yield, fruit average weight, and firmness. The columns in the CPgeno file are the markers whereas the rows are the individuals. The CPpheno data frame contains the measurements for the 363 siblings, and as mentioned before are averages over 2 environments.

}
\usage{data("DT_cpdata")}
\format{
  The format is:
 chr "DT_cpdata"
}
\source{
This data was simulated for fruit breeding applications.
}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744

}
\examples{

# 
# ####=========================================####
# #### For CRAN time limitations most lines in the
# #### examples are silenced with one '#' mark,
# #### remove them and run the examples using
# #### command + shift + C |OR| control + shift + C
# ####=========================================####
data(DT_cpdata)
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# #### create the variance-covariance matrix
# A <- A.mat(GT) # additive relationship matrix
# A <- A + diag(1e-4, ncol(A), ncol(A))
# #### look at the data and fit the model
# head(DT)
# mix1 <- lmebreed(Yield~ (1|id) + (1|Rowf) + (1|Colf),
#                  relmat=list(id=Matrix::chol(A)),
#                  control = lmerControl(
#                    check.nobs.vs.nlev = "ignore",
#                    check.nobs.vs.rankZ = "ignore",
#                    check.nobs.vs.nRE="ignore"
#                  ),
#                  data=DT)
# vc <- VarCorr(mix1); print(vc,comp=c("Variance"))


}
\keyword{datasets}
\seealso{The core function of the package \code{\link{lmebreed}} }
