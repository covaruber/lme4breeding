\name{imputev}
\alias{imputev}
\title{Imputing a numeric or character vector}
\description{
This function is a very simple function to impute a numeric or character vector with the mean or median value of the vector.
}
\usage{
imputev(x, method="median",by=NULL)
}

\arguments{
  \item{x}{a numeric or character vector.}
  \item{method}{the method to choose between mean or median.}
  \item{by}{the level at which imputation will be done.}
}
\value{
\describe{

\item{$x}{a numeric or character vector imputed with the method selected.}

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
set.seed(1253)
x <- rnorm(100)
x[sample(1:100,10)] <- NA
imputev(x)
}
\seealso{The core function of the package \code{\link{lmeb}}}