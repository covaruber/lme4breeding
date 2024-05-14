\name{redmm}
\alias{redmm}
\title{Reduced Model Matrix}
\description{
  `redmm` reduces a model matrix by performing a singular value decomposition or Cholesky on an incidence matrix.
}
\usage{
  redmm(x, M = NULL, Lam=NULL, nPC=50, cholD=FALSE, returnLam=FALSE)
}

\arguments{
  \item{x}{as vector with values to form a model matrix or the complete incidence matrix itself for an effect of interest.}
  \item{M}{an optional matrix of features explaining the levels of x. If not provided is assumed that the entire incidence matrix has been provided in x. But if provided, the decomposition occurs in the matrix M. }
  \item{Lam}{a matrix of loadings in case is already available to avoid recomputing it.}
  \item{nPC}{number of principal components to keep from the matrix of loadings to form the model matrix.}
  \item{cholD}{should a Cholesky or a Singular value decomposition should be used. The default is the SVD.}
  \item{returnLam}{should the function return the loading matrix in addition to the incidence matrix. Default is FALSE.}
}
\value{
  \describe{
    
    \item{$S3}{ A list with 3 elements:
    
    1) The model matrix to be used in the mixed modeling.
    
    2) The reduced matrix of loadings (nPC columns).
    
    3) The full matrix of loadings.
    
    }
  }
}
\references{
  
 Giovanny Covarrubias-Pazaran & Johan Aparicio-Arce (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48. 
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
  
data(DT_technow)
DT <- DT_technow
Md <- Md_technow

M <- tcrossprod(Md)
Z = with(DT, redmm(x=dent, M=M, nPC=10))
custom <- (rep(colnames(Z), nrow(DT)))[1:nrow(DT)]
ans <- lmebreed(GY ~ (1|custom),
                addmat = list(custom=Z),
            data=DT)
vc <- VarCorr(ans); print(vc,comp=c("Variance"))

xx <- with(DT, redmm(x=dent, M=M, nPC=10, returnLam = TRUE))
u = tcrossprod(xx$Lam, t(as.matrix( ranef(ans)[[1]] ) ))


}
\seealso{The core function of the package \code{\link{lmebreed}} }