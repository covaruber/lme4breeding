\name{A.mat}
\alias{A.mat}

\title{
Additive relationship matrix
}
\description{
Calculates the realized additive relationship matrix. 
}
\usage{
A.mat(X,min.MAF=NULL)
}

\arguments{
  \item{X}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} lines and \eqn{m} biallelic markers, 
coded as \{-1,0,1\}. Fractional (imputed) and missing values (NA) are allowed.
}
\item{min.MAF}{
Minimum minor allele frequency. The A matrix is not sensitive to rare alleles, so by default only monomorphic markers are removed.
}
}
\details{

For vanraden method: the marker matrix is centered by subtracting column means \eqn{M= X - ms} where ms is the coumn means. Then \eqn{A=M M'/c}, where \eqn{c = \sum_k{d_k}/k}, the mean value of the diagonal values of the \eqn{M M'} portion.

}
\value{
If return.imputed = FALSE, the \eqn{n \times n} additive relationship matrix is returned.

If return.imputed = TRUE, the function returns a list containing
\describe{
\item{$A}{the A matrix}
}
}

\references{

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}

\examples{

## random population of 200 lines with 1000 markers
X <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  X[i,] <- ifelse(runif(1000)<0.5,-1,1)
}

A <- A.mat(X)

 \donttest{
## take a look at the Genomic relationship matrix 
 colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
 hv <- heatmap(A[1:15,1:15], col = colfunc(100),Colv = "Rowv")
 str(hv)
 }

}
\seealso{\code{\link{lmebreed}} -- the core function of the package }
