\name{balanceData}
\alias{balanceData}

\title{
Balance a dataset
}
\description{
Used to balance a dataset for given slope and intercept.
}
\usage{
balanceData(data, slope=NULL, intercept=NULL)
}

\arguments{
  \item{data}{
Dataset to be balanced.
}

  \item{slope}{
A character value indicating the column in the dataset that will serve as slope in the
random regression model.
}

\item{intercept}{
A character vector indicating the column(s) in the dataset that will serve as intercepts in the
random regression model. If more than one column name is provided the intercept value is assumed 
to be the concatenation of the variables provided.
}

}
\details{

It returns the same original dataset but balanced for the slopes given certain
intercept variables.

}
\value{

It returns the new balanced dataset.

}

\references{

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}

\examples{

DT = DT_big
M = apply(M_big,2,as.numeric)
rownames(M) <- rownames(M_big)

# unbalance the data intentionally to mimic a real scenario with missing records
nas <- sample(1:nrow(DT), round(nrow(DT)*.2))
DTu <- droplevels(DT[setdiff(1:nrow(DT), nas),])

# balance your data to be able to use the eigen decomposition
DTb <- balanceData(data=DTu, slope = "id", intercept = c("envf","repf"))
DTb$value <- imputev(x=DTb$value, by=DTb$env)

# now this dataset can be used with the argument 'rotation=TRUE' in lmebreed

}
\seealso{\code{\link{lmebreed}}-- the core function of the package }
