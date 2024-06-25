\name{rrm}
\alias{rrm}
\title{reduced rank matrix}
\description{
  \code{rrm} creates a reduced rank factor analytic matrix by selecting the n vectors of the L matrix of the Cholesky decomposition or the U vectors of the SVD decomposition (loadings or latent covariates) to create a new incidence matrix of latent covariates that can be used with the \code{\link{lmebreed}} solver to fit random regressions on the latent covariates.
}
\usage{
  rrm(x=NULL, H=NULL, nPC=2, returnGamma=FALSE, cholD=TRUE)
}

\arguments{
  \item{x}{vector of the dataset containing the variable to be used to form the incidence matrix.}
  \item{H}{two-way table of identifiers (rows; e.g., genotypes) by features (columns; e.g., environments) effects. Row names and column names are required. No missing data is allowed.}
  \item{nPC}{number of principal components to keep from the loadings matrix.}
  \item{returnGamma}{a TRUE/FALSE argument specifying if the function should return the matrix of loadings used to build the incidence matrix for the model. The default is FALSE so it returns only the incidence matrix.}
  \item{cholD}{a TRUE/FALSE argument specifying if the Cholesky decomposition should be calculated or the singular value decomposition should be used instead.}
  
}
\details{
This implementation of a version of the reduced rank factor analytic models uses the so-called principal component (PC) models (Meyer, 2009) which assumes specific effects (psi) are equal to 0. The model is as follows:

y = Xb + Zu + e

where the variance of u ~ MVN(0, Sigma)

Sigma = (Gamma_t Gamma) + Psi

\strong{Extended factor analytic model}:

y = Xb + Z(I Gamma)c + Zs + e = Xb + Z*c + Zs + e

where \code{y} is the response variable, \code{X} and \code{Z} are incidence matrices for fixed and random effects respectively, \code{I} is a diagonal matrix, \code{Gamma} are the factor loadings for \code{c} common factor scores, and \code{s} are the specific effects, \code{e} is the vector of residuals.   

\strong{Reduced rank model:}

y = Xb + Z(I Gamma)c + e = Xb + Z*c + e

which is equal to the one above but assumes specific effects = 0. 

\strong{The algorithm in rrm the following:}

1) uses a wide-format table of timevar (m columns) by idvar (q rows) named H to form the initial variance-covariance matrix (Sigma) which is calculated as Sigma = H'H of dimensions m x m (column dimensions, e.g., environments x environments). 

2) The Sigma matrix is then center and scaled. 

3) A Cholesky (L matrix) or SVD decomposition (U D V') is performed in the Sigma matrix.

4) n vectors from L (when Cholesky is used) or U sqrt(D) (when SVD is used) are kept to form Gamma. Gamma = L[,1:nPc] or Gamma = U[,1:nPC]. These are the so-called loadings (L for all loadings, Gamma for the subset of loadings).

4) Gamma is used to form a new incidence matrix as Z* = Z Gamma

5) This matrix is later used for the REML machinery to be used with the usc (unstructured) or smm (diagonal) structures to estimate variance components and factor scores. The resulting BLUPs from the mixed model are the optimized factor scores. Pretty much as a random regression over latent covariates.

This implementation does not update the loadings (latent covariates) during the REML process, only estimates the REML factor scores for fixed loadings. This is different to other software (e.g., asreml) where the loadings are updated during the REML process as well. 

BLUPs for genotypes in all locations can be recovered as:

u = Gamma * u_scores

The resulting loadings (Gamma) and factor scores can be thought as an equivalent to the classical factor analysis.

}
\value{
  \describe{
    
    \item{$Z}{ a incidence matrix Z* = Z Gamma which is the original incidence matrix for the timevar multiplied by the loadings.}
    
    \item{$Gamma}{ a matrix of loadings or latent covariates.}
        
    \item{$Sigma}{ the covariance matrix used to calculate Gamma.}
    
  }
}
\references{

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

  Meyer K (2009) Factor analytic models for genotype by environment type problems and structured covariance matrices. Genetics Selection Evolution, 41:21
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{

data(DT_h2)
DT <- DT_h2
DT=DT[with(DT, order(Env)), ]
head(DT)
indNames <- na.omit(unique(DT$Name))
A <- diag(length(indNames))
rownames(A) <- colnames(A) <- indNames

\donttest{
  
# fit diagonal model first to produce H matrix
Z <- with(DT, smm(Env))
diagFormula <- paste0( "y ~ Env + (0+", paste(colnames(Z), collapse = "+"), "|| Name)")
for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
ans1a <- lmebreed(as.formula(diagFormula),
                  relmat = list(Name = A ),
                  data=DT)
vc <- VarCorr(ans1a); print(vc,comp=c("Variance"))
H0 <- ranef(ans1a)$Name # GxE table

# reduced rank model
Z <- with(DT,  rrm(Env, H = H0, nPC = 3) ) 
Zd <- with(DT, smm(Env) )
faFormula <- paste0( "y ~ Env + (0+", paste(colnames(Z), collapse = "+"),
                     "| Name) + (0+",paste(colnames(Zd), collapse = "+"), "|| Name)")
for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
ansFA <- lmebreed(as.formula(faFormula),
                  relmat = list(Name = A ),
                  data=DT)
vc <- VarCorr(ansFA); print(vc,comp=c("Variance"))
  
}
\dontrun{

loadings=with(DT, rrm(Env, nPC = 3, H = H0, returnGamma = TRUE) )$Gamma
Gint <- loadings %*% vc$Name %*% t(loadings) ## replace * by matrix product
Gspec <- diag( unlist(lapply(vc[2:16], function(x){x[[1]]})) )
G <- Gint + Gspec
## lattice::levelplot(cov2cor(G))
colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
hv <- heatmap(cov2cor(G), col = colfunc(100), symm = TRUE)

u <- ranef(ansFA)$Name
uInter <- as.matrix(u[,1:3]) %*% t(as.matrix(loadings)) ## replace * by matrix product
uSpec <- as.matrix(u[,-c(1:3)])
u <- uSpec + uInter
  
}


}
\seealso{The \code{\link{lmebreed}} solver.}