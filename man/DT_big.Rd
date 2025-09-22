\name{DT_big}
\alias{DT_big}
\docType{data}
\title{
Big GxE Quantitative Genetic Simulation
}
\description{

A data frame and genetic markers for a population of 1000 individuals evaluated in 
25 environments and 2 replicates per environment giving a total of 50 K records
to show how to fit big GxE models with genomic information.

}
\usage{data("DT_big")}
\format{
  The format is:
 chr "DT_big"
}
\references{

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}
\examples{

data(DT_big)
DT = DT_big
M = apply(M_big,2,as.numeric)
rownames(M) <- rownames(M_big)

\donttest{

# compute the relationship matrix
MMT <- tcrossprod(M)
m <- sum(diag(MMT))/nrow(MMT)
MMT <- MMT/m
MMT <- MMT + diag(1e-05, ncol(MMT), ncol(MMT))

# let's fit the GxE model (diagonal first)
DT[,"envf_repf"] = paste(DT[,"envf"],DT[,"repf"],sep="_")
Z <- with(DT, smm(envf_repf))
for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
csdiagFormula <- paste0( "value ~ (1|envf) + (0+", paste(colnames(Z), collapse = "+"), "|| id)")

# fit the models
mixm2 <- lmebreed(formula=as.formula(csdiagFormula), data = DT, 
                  relmat = list(id = MMT), 
                  control = lmerControl( # how to control n iterations
                   calc.derivs = FALSE,
                   restart_edge = FALSE,
                   optCtrl = list(maxfun = 5000, maxeval = 5000)
                  ),
                  verbose=1L, rotation=TRUE
)

vc <- VarCorr(mixm2); print(vc,comp=c("Variance"))
ran2 = ranef(mixm2)
H0 <- ran2$id

}

}
\keyword{datasets}
\seealso{The core function of the package \code{\link{lmebreed}} }
