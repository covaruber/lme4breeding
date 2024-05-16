\name{stackTrait}
\alias{stackTrait}
\title{Stacking traits in a dataset}
\description{
  \code{stackTrait} creates a dataset stacking traits in the long format to be used with the \code{\link{lmebreed}} solver for multi-trait models.
}
\usage{
  stackTrait(data, traits)
}

\arguments{
  \item{data}{a data frame with traits in wide format.}
  \item{traits}{variable names corresponding to the traits that should be in the long format.}
}
\value{
  \describe{
    
    \item{$res}{ a data frame with traits in long format.}
  }
}
\references{
  
Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}
\author{
  Giovanny Covarrubias-Pazaran
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
DT2 <- stackTrait(data=DT, traits = c("Yield","color"))
head(DT2$long)

\donttest{

mix1 <- lmebreed(valueS~ (0+trait|id),
                 relmat=list(id=A),
                 control = lmerControl(
                   check.nobs.vs.nlev = "ignore",
                   check.nobs.vs.rankZ = "ignore",
                   check.nobs.vs.nRE="ignore"
                 ),
                 data=DT2$long)
vc <- VarCorr(mix1); print(vc,comp=c("Variance"))
cov2cor(vc$id)

}

}
\seealso{The function for the \code{\link{lmebreed}} solver.}