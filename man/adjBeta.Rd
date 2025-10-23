\name{adjBeta}
\alias{adjBeta}
\title{Adjusting fixed effects for intercept}
\description{
This function is a very simple function to add the intercept to all fixed effects except for the first term.
}
\usage{
adjBeta(x)
}

\arguments{
  \item{x}{a numeric vector with fixed effects extracted by the fixef function.}
}
\value{
\describe{

\item{$x}{a numeric vector with the intercept added to all fixed effects except for the first term which corresponds to the intercept.}

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
data(DT_example)
DT <- DT_example
A <- A_example

ansMain <- lmeb(Yield ~ Env + (1|Name),
                    relmat = list(Name = Matrix::chol(A) ),
                    data=DT)

fixef(ansMain)
adjBeta(fixef(ansMain))
}
\seealso{The core function of the package \code{\link{lmeb}}}