\name{DT_ige}
\alias{DT_ige}
\alias{A_ige}
\docType{data}
\title{
Data to fit indirect genetic effects.
}
\description{

This dataset contains phenotpic data for 98 individuals where they are measured with the purpose of identifying the effect of the neighbour in a focal individual. 

}
\usage{data("DT_ige")}
\format{
  The format is:
 chr "DT_ige"
}
\source{
This data was masked from a shared study.
}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\examples{

# 
# ####=========================================####
# #### For CRAN time limitations most lines in the 
# #### examples are silenced with one '#' mark, 
# #### remove them and run the examples
# ####=========================================####
# ####=========================================####
# #### EXAMPLES
# #### Different models with sommer
# ####=========================================####
 data(DT_ige)
# DT <- DT_ige
# # Indirect genetic effects model without covariance between DGE and IGE
# modIGE <- lmebreed(trait ~ block + (1|focal) + (1|neighbour),
#               data = DT)
# vc <- VarCorr(modIGE); print(vc,comp=c("Variance"))
# 
# # Add relationship matrices
# A_ige <- A_ige + diag(1e-4, ncol(A_ige), ncol(A_ige) )
# modIGE <- lmebreed(trait ~ block + (1|focal) + (1|neighbour),
#                    relmat = list(focal=Matrix::chol(A_ige),
#                                  neighbour=Matrix::chol(A_ige)),
#                    data = DT)
# vc <- VarCorr(modIGE); print(vc,comp=c("Variance"))
# 
# # Indirect genetic effects model with covariance between DGE and IGE using relationship matrices
# # Relationship matrix
# A_ige <- A_ige + diag(1e-4, ncol(A_ige), ncol(A_ige) )
# # Define 2 dummy variables to make a fake covariance
# # for two different random effects
# DT$fn <- DT$nn <- 1
# # Create the incidence matrix for the first random effect
# Zf <- Matrix::sparse.model.matrix( ~ focal-1, data=DT )
# colnames(Zf) <- gsub("focal","", colnames(Zf))
# Zf <- Zf[which(!is.na(DT$trait)),]
# # Create the incidence matrix for the second random effect
# Zn <- Matrix::sparse.model.matrix( ~ neighbour-1, data=DT )
# colnames(Zn) <- gsub("neighbour","", colnames(Zn))
# Zn <- Zn[which(!is.na(DT$trait)),]
# # Make inital values for incidence matrix but irrelevant
# # since these will be replaced by the addmat argument
# both <- (rep(colnames(Zf), nrow(DT)))[1:nrow(DT)]
# # Fit the model
# modIGE <- lmebreed(trait ~ block + (0+fn+nn|both),
#                    addmat = list(both=list(Zf,Zn)),
#                    relmat = list(both=Matrix::chol(A_ige)),
#                    data = DT)
# vc <- VarCorr(modIGE); print(vc,comp=c("Variance"))
# 
# blups <- ranef(modIGE)
# pairs(blups$both)
# cov2cor(vc$both)


}
\keyword{datasets}
\seealso{The core function of the package \code{\link{lmebreed}} }
