\name{getMME}
\alias{getMME}

\title{
Build the mixed model equations from a lmebreed object
}
\description{
Uses the lmebreed object and builds the coefficient matrix (C) and returns its 
inverse and the solutions to the mixed model equations.
}
\usage{
getMME(object, vc=NULL, recordsToUse=NULL)
}

\arguments{
  \item{object}{
an object of class lmebreed.
  }
  \item{vc}{
optional variance components to force in the mixed model equations. This this to 
be the outpur of the \code{\link[lme4]{VarCorr}} function.
  }
  \item{recordsToUse}{
 a numeric vector of indices to say which records should be kept for forming the 
 mixed model equations and get solutions. This is particularly useful when we want
 to predict new individuals.
  }

}
\details{

Uses the lmebreed object and builds the coefficient matrix (C) and returns its inverse and the solutions to the mixed model equations. It is internally used by the ranef function when the user wants standard errors for the BLUPs.

}
\value{

\describe{
\item{$Ci}{inverse of the coefficient matrix.}
\item{$bu}{solutions to the mixed model equations}
}

}

\references{

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}

\examples{

data(DT_example)
DT <- DT_example
A <- A_example

ansMain <- lmebreed(Yield ~ Env + (1|Name),
                    relmat = list(Name = A ),
                    data=DT)
mme <- getMME(ansMain)

\donttest{

##########################################
## showing how to predict the individuals
## that didn't have records in the dataset
##########################################
data(DT_cpdata)
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
#### create the variance-covariance matrix
A <- A.mat(GT) # additive relationship matrix
A <- A + diag(1e-4, ncol(A), ncol(A))
#### look at the data and fit the model
head(DT)
mix1 <- lmebreed(Yield~ (1|id) + (1|Rowf) + (1|Colf),
                 relmat=list(id=A),
                 control = lmerControl(
                   check.nobs.vs.nlev = "ignore",
                   check.nobs.vs.rankZ = "ignore",
                   check.nobs.vs.nRE="ignore"
                 ),
                 data=DT)

vc <- VarCorr(mix1); print(vc,comp=c("Variance"))
# the new dataset includes more individuals
DT2 <- DT
DT2$Yield <- imputev(DT2$Yield)
mix1expanded <- lmebreed(Yield~ (1|id) + (1|Rowf) + (1|Colf),
                 relmat=list(id=A),
                 control = lmerControl(
                   check.nobs.vs.nlev = "ignore",
                   check.nobs.vs.rankZ = "ignore",
                   check.nobs.vs.nRE="ignore",
                   calc.derivs=TRUE,
                   optCtrl=list(maxeval=1)
                 ),
                 data=DT2)

vc <- VarCorr(mix1expanded); print(vc,comp=c("Variance"))
# predict the individuals that didn't have records in the dataset
res <- getMME(object=mix1expanded, vc=VarCorr(mix1), recordsToUse = which(!is.na(DT$Yield)) )


}


}
\seealso{\code{\link{lmebreed}}-- the core function of the package }
