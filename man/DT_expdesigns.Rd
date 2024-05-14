\name{DT_expdesigns}
\alias{DT_expdesigns}
\docType{data}
\title{
  Data for different experimental designs
}
\description{
  The following data is a list containing data frames for different type of experimental designs relevant in plant breeding:
  
  1) Augmented designs (2 examples)
  
  2) Incomplete block designs (1 example)
  
  3) Split plot design (2 examples)
  
  4) Latin square designs (1 example)
  
  5) North Carolina designs I,II and III
  
  How to fit each is shown at the Examples section. This may help you get introduced to experimental designs relevant to plant breeding. Good luck.
}

\format{
 Different based on the design.
}

\source{
  Datasets and more detail about them can be found in the agricolae package. Here we just show the datasets and how to analyze them using the \code{\link{lme4breeding}} package.
}

\references{

Giovanny Covarrubias-Pazaran & Johan Aparicio-Arce (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}

\examples{

data(DT_expdesigns)
DT <- DT_expdesigns
names(DT)
data1 <- DT$au1
head(data1)
## response variable: "yield"
## check indicator: "entryc" ('nc' for all unreplicated, but personal.name for checks)
## blocking factor: "block"
## treatments, personal names for replicated and non-replicated: "trt"
## check no check indicator: "new"
mix1 <- lmebreed(yield~entryc + (1|block)+(1|trt),
             data=data1)
vc <- VarCorr(mix1); print(vc,comp=c("Variance"))

}

