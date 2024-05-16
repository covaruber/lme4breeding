\name{add.diallel.vars}
\alias{add.diallel.vars}
\title{add.diallel.vars}
\description{
  `add.diallel.vars` adds 4 columns to the provided diallel dataset. Specifically, the user provides a dataset with indicator variables for who is the male and female parent and the function returns the same dataset with 4 new dummy variables to allow the model fit of diallel models.
}
\usage{
  add.diallel.vars(df, par1="Par1", par2="Par2",sep.cross="-")
}

\arguments{
  \item{df}{a dataset with the two indicator variables for who is the male and female parent.}
  \item{par1}{the name of the column indicating who is the first parent (e.g. male).}
  \item{par2}{the name of the column indicating who is the second parent (e.g. female).}
  \item{sep.cross}{the character that should be used when creating the column for cross.id. A simple paste of the columns par1 and par2.}
}
\value{
  \describe{
    
    A new data set with the following 4 new dummy variables to allow the fit of complex diallel models:
    
    \item{is.cross}{returns a 0 if is a self and a 1 for a cross.}
    \item{is.self}{returns a 0 if is a cross and a 1 is is a self.}
    \item{cross.type}{returns a -1 for a direct cross, a 0 for a self and a 1 for a reciprocal cross.}
    \item{cross.id}{returns a column psting the par1 and par2 columns.}
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

data(DT_mohring)
DT <- DT_mohring
head(DT)
DT2 <- add.diallel.vars(DT,par1="Par1", par2="Par2")
head(DT2)
## see ?DT_mohring for an example on how to use the data to fit diallel models.

}
\seealso{The core function of the package \code{\link{lmebreed}} function and the \code{\link{DT_mohring}} example. }