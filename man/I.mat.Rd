\name{I.mat}
\alias{I.mat}

\title{
Identity relationship matrix
}
\description{
Shortcut to get an identity relationship matrix. 
}
\usage{
I.mat(x)
}

\arguments{
  \item{x}{
Vector of values to identify the unique values and form an iende.
}

}
\details{

Nothing but a shortcut to avoid few lines of code.

}
\value{
If return.imputed = FALSE, the \eqn{n \times n} identity matrix is returned.

\describe{
\item{$I}{the I matrix}
}
}

\references{

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}

\examples{

## random population of 200 lines with 1000 markers


}
\seealso{\code{\link{lmeb}} -- the core function of the package }
