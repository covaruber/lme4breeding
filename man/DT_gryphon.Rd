\name{DT_gryphon}
\alias{DT_gryphon}
\alias{A_gryphon}
\alias{P_gryphon}
\docType{data}
\title{
Gryphon data from the Journal of Animal Ecology
}
\description{
This is a dataset that was included in the Journal of animal ecology by Wilson et al. (2010; see references) to help users understand how to use mixed models with animal datasets with pedigree data. 

The dataset contains 3 elements:

gryphon; variables indicating the animal, the mother of the animal, sex of the animal, and two quantitative traits named 'BWT' and 'TARSUS'.

pedi; dataset with 2 columns indicating the sire and the dam of the animals contained in the gryphon dataset.

A; additive relationship matrix formed using the 'getA()' function used over the pedi dataframe.

}
\usage{data("DT_gryphon")}
\format{
  The format is:
 chr "DT_gryphon"
}
\source{
This data comes from the Journal of Animal Ecology. Please, if using this data cite Wilson et al. publication. If using our mixed model solver please cite Covarrubias' publication.
}
\references{

Wilson AJ, et al. (2010) An ecologist's guide to the animal model. Journal of Animal Ecology 79(1): 13-26.

Giovanny Covarrubias-Pazaran & Johan Aparicio-Arce (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48. 
 
}
\examples{

data(DT_gryphon)
DT <- DT_gryphon
A <- A_gryphon
P <- P_gryphon
#### look at the data
head(DT)

\donttest{

## fit the model with no fixed effects (intercept only)
mix1 <- lmebreed(BWT~ (1|ANIMAL),
                 relmat = list(ANIMAL=Matrix::chol(A)),
                 control = lmerControl(
                   check.nobs.vs.nlev = "ignore",
                   check.nobs.vs.rankZ = "ignore",
                   check.nobs.vs.nRE="ignore"
                 ),
                 data=DT)
vc <- VarCorr(mix1); print(vc,comp=c("Variance"))

}

}
\keyword{datasets}
\seealso{The core function of the package \code{\link{lmebreed}} }
