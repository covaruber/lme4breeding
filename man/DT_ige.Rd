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

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}
\examples{

data(DT_ige)
DT <- DT_ige
# Indirect genetic effects model without covariance between DGE and IGE
modIGE <- lmebreed(trait ~ block + (1|focal) + (1|neighbour),
              data = DT)
vc <- VarCorr(modIGE); print(vc,comp=c("Variance"))


\donttest{

## Add relationship matrices
A_ige <- A_ige + diag(1e-4, ncol(A_ige), ncol(A_ige) )
modIGE <- lmebreed(trait ~ block + (1|focal) + (1|neighbour),
                   relmat = list(focal=A_ige,
                                 neighbour=A_ige),
                   data = DT)
vc <- VarCorr(modIGE); print(vc,comp=c("Variance"))

## Indirect genetic effects model with covariance between DGE and IGE using relationship matrices
## Relationship matrix
A_ige <- A_ige + diag(1e-4, ncol(A_ige), ncol(A_ige) )
## Define 2 dummy variables to make a fake covariance
## for two different random effects
DT$fn <- DT$nn <- 1
## Create the incidence matrix for the first random effect
Zf <- Matrix::sparse.model.matrix( ~ focal-1, data=DT )
colnames(Zf) <- gsub("focal","", colnames(Zf))
## Create the incidence matrix for the second random effect
Zn <- Matrix::sparse.model.matrix( ~ neighbour-1, data=DT )
colnames(Zn) <- gsub("neighbour","", colnames(Zn))
## Make inital values for incidence matrix but irrelevant
## since these will be replaced by the addmat argument
both <- (rep(colnames(Zf), nrow(DT)))[1:nrow(DT)]
## Fit the model
modIGE <- lmebreed(trait ~ block + (0+fn+nn|both),
                   addmat = list(both=list(Zf,Zn)),
                   relmat = list(both=A_ige),
                   data = DT)
vc <- VarCorr(modIGE); print(vc,comp=c("Variance"))
sigma(modIGE)^2 # error variance

blups <- ranef(modIGE)
PEV <- lapply(blups, function(x){attr(x, which="postVar")}) # take sqrt() for SEs
pairs(blups$both)
cov2cor(vc$both)

}


}
\keyword{datasets}
\seealso{The core function of the package \code{\link{lmebreed}} }
