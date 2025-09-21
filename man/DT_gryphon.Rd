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

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

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
                 relmat = list(ANIMAL=A),
                 # how to control n iterations
                 # control = lmerControl( 
                 #  optCtrl = list(maxfun = 100, maxeval = 100)
                 #   ),
                 data=DT)
vc <- VarCorr(mix1); print(vc,comp=c("Variance"))
sigma(mix1)^2 # error variance
BLUP <- ranef(mix1, condVar=TRUE)
SEs <- attr(BLUP$ANIMAL, which="postVar")[,,]


### multi-trait model
traits <- c("BWT","TARSUS")
for(iTrait in traits){DT[,iTrait] <- scale(imputev(DT[,iTrait]))}
DTL <- reshape(DT[,c("ANIMAL", traits)], idvar = "ANIMAL", varying = traits,
               v.names = "value", direction = "long",
               timevar = "trait", times = traits )
DTL <- DTL[with(DTL, order(trait)), ]
head(DTL)

system.time(
  mix <- lmebreed(value ~ (0+trait|ANIMAL),
                  relmat = list(ANIMAL=A),
                  # rotation = TRUE, 
                  data=DTL)
)
vc <- VarCorr(mix); print(vc,comp=c("Variance"))
cov2cor(vc$ANIMAL)
sigma(mix)^2 # error variance


}

}
\keyword{datasets}
\seealso{The core function of the package \code{\link{lmebreed}} }
