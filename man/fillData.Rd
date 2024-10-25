\name{fillData}
\alias{fillData}
\title{Filling gaps for a dataset to balance}
\description{
  \code{fillData} creates a balanced dataset in the long format to be used with the \code{\link{lmebreed}} solver for multi-trait models.
}
\usage{
  fillData(data, toBalanceSplit=NULL, toBalanceFill=NULL)
}

\arguments{
  \item{data}{a data frame with traits in wide format.}
  \item{toBalanceSplit}{variable to split the dataset for balancing.}
  \item{toBalanceFill}{variable of the factor to balance across the levels of the toBalanceSplit variable.}
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

DT <- DT_example
A <- A_example
DT2 <- fillData(DT , toBalanceSplit="Env", toBalanceFill="Name" )
head(DT2)

\donttest{

## Compound simmetry (CS) model
ans2 <- lmebreed(Yield~Env+ (0+Env|Name),
                 relmat = list(Name=A),
                 rotation = TRUE,
                 data=DT2)
vc <- VarCorr(ans2); print(vc,comp=c("Variance"))

BLUP <- ranef(ans2, condVar=TRUE)$Name
SEs <- attr(BLUP, which="postVar")[,,]

}

}
\seealso{The function for the \code{\link{lmebreed}} solver.}