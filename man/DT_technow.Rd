\name{DT_technow}
\alias{DT_technow}
\alias{Ad_technow}
\alias{Af_technow}
\alias{Md_technow}
\alias{Mf_technow}
\docType{data}
\title{
Genotypic and Phenotypic data from single cross hybrids (Technow et al.,2014)
}
\description{
This dataset contains phenotpic data for 2 traits measured in 1254 single cross hybrids coming from the cross of Flint x Dent heterotic groups. In addition contains the genotipic data (35,478 markers) for each of the 123 Dent lines and 86 Flint lines. The purpose of this data is to demosntrate the prediction of unrealized crosses (9324 unrealized crosses, 1254 evaluated, total 10578 single crosses). We have added the additive relationship matrix (A) but can be easily obtained using the A.mat function on the marker data. Please if using this data for your own research cite Technow et al. (2014) publication (see References).

}
\usage{data("DT_technow")}
\format{
  The format is:
 chr "DT_technow"
}
\source{
This data was extracted from Technow et al. (2014).
}
\references{
If using this data for your own research please cite:

Technow et al. 2014. Genome properties and prospects of genomic predictions of hybrid performance in a Breeding program of maize. Genetics 197:1343-1355.

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\examples{

####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples using
#### command + shift + C |OR| control + shift + C
####=========================================####
data(DT_technow)
# DT <- DT_technow
# Md <- Md_technow
# Mf <- Mf_technow
# Md <- (Md*2) - 1
# Mf <- (Mf*2) - 1
# Ad <- A.mat(Md)
# Af <- A.mat(Mf)
# Ad <- Ad + diag(1e-4, ncol(Ad), ncol(Ad))
# Af <- Af + diag(1e-4, ncol(Af), ncol(Af))
# ###=========================================####
# ###=========================================####
# ans2 <- lmebreed(GY ~ (1|dent) + (1|flint),
#              relmat = list(dent=Matrix::chol(Ad), 
#                            flint=Matrix::chol(Af)),
#              data=DT)
# vc <- VarCorr(ans2); print(vc,comp=c("Variance"))
# ####################
# ### overlayed model
# ####################
# M <- rbind(Md,Mf)
# A <- A.mat(M)
# A <- A + diag(1e-4,ncol(A), ncol(A))
# Z <- with(DT, overlay(dent,flint) )
# # inital values for incidence matrix but irrelevant
# # since these will be replaced by admat argument
# fema <- (rep(colnames(Z), nrow(DT)))[1:nrow(DT)]
# #### model using overlay without relationship matrix
# ans2 <- lmebreed(GY ~ (1|fema),
#                  addmat = list(fema=Z),
#                  relmat = list(fema=Matrix::chol(A)),
#                  data=DT)
# vc <- VarCorr(ans2); print(vc,comp=c("Variance"))



}
\keyword{datasets}
\seealso{The core function of the package \code{\link{lmebreed}} }
