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

Giovanny Covarrubias-Pazaran & Johan Aparicio-Arce (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}
\examples{

data(DT_technow)
DT <- DT_technow
Md <- (Md_technow*2) - 1
Mf <- (Mf_technow*2) - 1
Ad <- A.mat(Md)
Af <- A.mat(Mf)
Ad <- Ad + diag(1e-4, ncol(Ad), ncol(Ad))
Af <- Af + diag(1e-4, ncol(Af), ncol(Af))
## simple model
ans2 <- lmebreed(GY ~ (1|dent) + (1|flint),
                 data=DT)
vc <- VarCorr(ans2); print(vc,comp=c("Variance"))

\donttest{

### with relationship matrices
ans2 <- lmebreed(GY ~ (1|dent) + (1|flint),
             relmat = list(dent=Matrix::chol(Ad),
                           flint=Matrix::chol(Af)),
             data=DT)
vc <- VarCorr(ans2); print(vc,comp=c("Variance"))

### overlayed model
M <- rbind(Md,Mf)
A <- A.mat(M)
A <- A + diag(1e-4,ncol(A), ncol(A))
Z <- with(DT, overlay(dent,flint) )
Z = Z[which(!is.na(DT$GY)),]
# inital values for incidence matrix but irrelevant
# since these will be replaced by admat argument
fema <- (rep(colnames(Z), nrow(DT)))[1:nrow(DT)]
#### model using overlay without relationship matrix
ans2 <- lmebreed(GY ~ (1|fema),
                 addmat = list(fema=Z),
                 relmat = list(fema=Matrix::chol(A)),
                 data=DT)
vc <- VarCorr(ans2); print(vc,comp=c("Variance"))

}

}
\keyword{datasets}
\seealso{The core function of the package \code{\link{lmebreed}} }
