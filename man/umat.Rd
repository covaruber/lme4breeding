\name{umat}
\alias{umat}
\title{Rotation matrix UDU' decomposition}
\description{
  `umat` creates the equivalent of a U matrix from the eigen decomposition of a relationship matrix of dimensions equal to the number of records equivalent to Lee and Van der Werf (2016).
}
\usage{
  umat(formula, relmat, data, addmat)
}

\arguments{
  \item{formula}{a formula expressing the factor to decompose.}
  \item{relmat}{an optional matrix of features explaining the levels of x. If not provided is assumed that the entire incidence matrix has been provided in x. But if provided, the decomposition occurs in the matrix M. }
  \item{data}{a dataset to be used for modeling.}
  \item{addmat}{additional matrices.}
}
\value{
  \describe{
    
    \item{$S3}{ A list with 3 elements:
    
    1) The U' matrix of dimensions n x n (eigen vectors), n being the number of records.
    
    2) The original U matrix of dimensions m x m (eigen vectors), m being the number of coefficients or levels in relmat.
    
    3) The D matrix of dimensions m x m (eigen values), m being the number of coefficients or levels in relmat.
    
    }
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
## create the variance-covariance matrix
A <- A.mat(GT) # additive relationship matrix
A <- A + diag(1e-4, ncol(A), ncol(A))
## look at the data and fit the model
head(DT)
xx <- umat(~id, relmat = list(id=A), data=DT)


}
\seealso{The core function of the package \code{\link{lmebreed}} }