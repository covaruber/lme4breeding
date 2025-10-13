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

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}
\examples{

data(DT_cpdata)
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
## create the variance-covariance matrix
A <- A.mat(GT) # additive relationship matrix
A <- A + diag(1e-4, ncol(A), ncol(A))
## look at the data and fit the model
head(DT)

\donttest{

mix1 <- lmebreed(Yield~ (1|id) + (1|Rowf) + (1|Colf),
                 relmat=list(id=A),
                 data=DT)
vc <- VarCorr(mix1); print(vc,comp=c("Variance"))
sigma(mix1)^2 # error variance

# run one last iteration with imputed data
# to make sure you get predictions for every level
DT2 <- DT
DT2$Yield <- imputev(DT2$Yield)
mix2 <- update(mix1, suppressOpt = TRUE,
               start=getME(mix1, "theta"), 
               data=DT2)
predsMix2 <- ranef(mix2)
PEV <- lapply(predsMix2, function(x){attr(x, which="postVar")}) # take sqrt() for SEs

# if you don't want the imputed vector to have an effect in
# the predictions you can use the getMME function to use
# the extended model and get predictions without including the 
# imputed data (I know is a bit messy)
preds <- getMME(object=mix2, # extended model
                vc=VarCorr(mix1), # variance components
                recordsToUse = which(!is.na(DT$Yield)) # records to use for MME
                )
# now you could compare between both types of predictions, the last ones are in 
# theory the correct ones.
plot(preds$bu[2:364,], predsMix2$id[,1])

}

}
\keyword{datasets}
\seealso{The core function of the package \code{\link{lmebreed}} }
