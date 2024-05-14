\name{DT_btdata}
\alias{DT_btdata}
\docType{data}
\title{
Blue Tit Data for a Quantitative Genetic Experiment
}
\description{

a data frame with 828 rows and 7 columns, with variables tarsus length (tarsus) and colour (back) measured on 828 individuals (animal). The mother of each is also recorded (dam) together with the foster nest (fosternest) in which the chicks were reared. The date on which the first egg in each nest hatched (hatchdate) is recorded together with the sex (sex) of the individuals.

}
\usage{data("DT_btdata")}
\format{
  The format is:
 chr "DT_btdata"
}
\references{

Giovanny Covarrubias-Pazaran & Johan Aparicio-Arce (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48. 

}
\examples{

data(DT_btdata)
DT <- DT_btdata
head(DT)

\donttest{

mix4 <- lmebreed(tarsus ~ sex + (1|dam) + (1|fosternest),
             data = DT)
vc <- VarCorr(mix4); print(vc,comp=c("Variance"))

}

}
\keyword{datasets}
\seealso{The core function of the package \code{\link{lmebreed}} }
