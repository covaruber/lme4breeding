\name{DT_polyploid}
\alias{DT_polyploid}
\alias{GT_polyploid}
\alias{MP_polyploid}
\docType{data}
\title{
Genotypic and Phenotypic data for a potato polyploid population
}
\description{
This dataset contains phenotpic data for 18 traits measured in 187 individuals from a potato diversity panel. In addition contains genotypic data for 221 individuals genotyped with 3522 SNP markers. Please if using this data for your own research make sure you cite Rosyara's (2015) publication (see References).

}
\usage{data("DT_polyploid")}
\format{
  The format is:
 chr "DT_polyploid"
}
\source{
This data was extracted from Rosyara (2016).
}
\references{
If using this data for your own research please cite:

Rosyara Umesh R., Walter S. De Jong, David S. Douches, Jeffrey B. Endelman. Software for genome-wide association studies in autopolyploids and its application to potato. The Plant Genome 2015.

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}
\examples{

data(DT_polyploid)
DT <- DT_polyploid
GT <- GT_polyploid
MP <- MP_polyploid
## convert markers to numeric format
numo <- atcg1234(data=GT[,1:100], ploidy=4);
numo$M[1:5,1:5];
numo$ref.allele[,1:5]

\donttest{

## plants with both genotypes and phenotypes
common <- intersect(DT$Name,rownames(numo$M))

## get the markers and phenotypes for such inds
marks <- numo$M[common,]; marks[1:5,1:5]
DT2 <- DT[match(common,DT$Name),];
DT2 <- as.data.frame(DT2)
DT2[1:5,]

## Additive relationship matrix, specify ploidy
A <- A.mat(marks)

## run as mixed model
A <- A + diag(1e-4,ncol(A),ncol(A))
ans <- lmebreed(tuber_shape~ (1|Name),
                relmat = list(Name=A),
                 # how to control n iterations
                 # control = lmerControl( 
                 #  optCtrl = list(maxfun = 100, maxeval = 100)
                 #   ),
                data=DT2)
vc <- VarCorr(ans); print(vc,comp=c("Variance"))
sigma(ans)^2 # error variance

}

}
\keyword{datasets}
\seealso{The core functions of the package \code{\link{lmebreed}} }
