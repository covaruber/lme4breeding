\name{lme4breeding-package}
\alias{lme4breeding}
\docType{package}
\author{
Giovanny Covarrubias-Pazaran
}
\title{

\strong{L}inear \strong{m}ixed \strong{e}quations \strong{4} \strong{B}reeding
\cr

}
\description{

lme4breeding is nice wrapper of the lme4 package that enables the use of especialized 
plant and animal breeding models that include relationship matrices among individuals 
(e.g., genomic relationship matrices) and complex covariance structures between 
factors (e.g., factor analytic structures) accelerated by the use of the eigen 
decomposition of relationship matrices. It uses all the lme4 machinery for linear
and non-linear models, for different response distributions opening a world of possibilities.

It took me several years to develop a package named sommer that allowed many of
the desired models. In May of 2024 I realized that few lines of code (exactly 100 lines)
would allow to tweak all the lme4 machinery to fit most plant and animal models 
popular today at a great speed enabled by the lme4 development team. 
I will not stop the development of the sommer package since it allows to fit
certain models at a greater speed than lme4breeding and other popular packages 
(in the p > n problem sommer still the best). The major advantage to use 
lme4breeding will be when you have balanced data for multiple traits or environments 
(not very likely but please look at the vignettes and make your own opinion) which 
will make this machinery extremely fast! if you have unbalanced data you may want 
to stick to the use of mmec until I discover how to adapt the eigen decomposition
to unbalanced data. I hope you enjoy it.

The \code{\link{lmebreed}} function is the core function of the package which is 
exactly the same function than lmer or glmer but with few added arguments 
\code{relmat} and \code{addmat} that allow the user to provide relationship matrices
and customized incidence matrices respectively. Also the argument \code{rotation} 
speeds up highly complex models. The lme4 machinery is designed to deal with a big 
number of records (r) since it works in the r > c problem and inverts a c x c matrix 
(being c the number of coefficients). There are \code{\link[lme4]{ranef}}, 
\code{\link[lme4]{fixef}}, \code{\link[lme4]{VarCorr}} functions to obtain 
variance-covariance components, BLUPs, BLUEs, residuals, fitted values, 
variances-covariances for fixed and random effects, etc.

}

\section{Functions for genetic analysis}{

The package provides kernels to the estimate additive (\code{\link{A.mat}})
relationship matrix for diploid and polyploid organisms. A good converter from 
letter code to numeric format is implemented in the function \code{\link{atcg1234}}, 
which supports higher ploidy levels than diploid. Additional functions for genetic 
analysis have been included such as build a genotypic hybrid marker matrix 
(\code{\link{build.HMM}}). If you need to use pedigree you need to convert your 
pedigree into a relationship matrix (use the `getA` function from the pedigreemm package).

}

\section{Functions for trial analysis}{

Recently, spatial modeling has been added added to lme4breeding using the two-dimensional
spline \code{\link{tps}} function.

}

\section{Keeping lme4breeding updated}{

The lme4breeding package is updated on CRAN every 4-months due to CRAN policies 
but you can find the latest source at https://github.com/covaruber/lme4breeding. 
This can be easily installed typing the following in the R console:

library(devtools)

install_github("covaruber/lme4breeding")

This is recommended if you reported a bug, was fixed and was immediately pushed 
to GitHub but not in CRAN until the next update.
}

\section{Tutorials}{

For tutorials on how to perform different analysis with lme4breeding please look 
at the vignettes by typing in the terminal:

\strong{vignette("lme4breeding.qg")}

\strong{vignette("lme4breeding.gxe")}

}

\section{Getting started}{
The package has been equiped with several datasets to learn how to use the 
lme4breeding package (and almost to learn all sort of quantitative genetic analysis): 

* \code{\link{DT_big}} simulated dataset containing 1K individuals in 50 environments to show how to fit big models.

* \code{\link{DT_btdata}} dataset contains an animal (birds) model.

* \code{\link{DT_cornhybrids}} dataset to perform genomic prediction in hybrid single crosses 

* \code{\link{DT_cpdata}} dataset to fit genomic prediction models within a biparental population coming from 2 highly heterozygous parents including additive, dominance and epistatic effects. 

* \code{\link{DT_fulldiallel}} dataset with examples to fit full diallel designs. 

* \code{\link{DT_gryphon}} data contains an example of an animal model including pedigree information.

* \code{\link{DT_halfdiallel}} dataset with examples to fit half diallel designs. 

* \code{\link{DT_h2}} to calculate heritability

* \code{\link{DT_ige}} dataset to show how to fit indirect genetic effect models.

* \code{\link{DT_legendre}} simulated dataset for random regression model.

* \code{\link{DT_mohring}} datasets with examples to fit half diallel designs.

* \code{\link{DT_polyploid}} to fit genomic prediction and GWAS analysis in polyploids. 

* \code{\link{DT_sleepstudy}} dataset to know how to translate lme4 models to lme4breeding models.

* \code{\link{DT_technow}} dataset to perform genomic prediction in hybrid single crosses 

* \code{\link{DT_wheat}} dataset to do genomic prediction in single crosses in species displaying only additive effects.

}

\section{Models Enabled}{

The machinery behind the scenes is lme4.

}

\section{Bug report and contact}{

If you have any questions or suggestions please post it in https://stackoverflow.com 
or https://stats.stackexchange.com

I'll be glad to help or answer any question. I have spent a valuable amount of
time developing this package. Please cite this package in your publication. 
Type 'citation("lme4breeding")' to know how to cite it.


}
\references{

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}
\keyword{ R package }
\examples{

data(DT_example)
DT <- DT_example
A <- A_example

ansMain <- lmebreed(Yield ~ Env + (1|Name),
                        relmat = list(Name = A ),
                        data=DT)
vc <- VarCorr(ansMain); print(vc,comp=c("Variance"))
BLUP <- ranef(ansMain, condVar=TRUE)
condVAR <- lapply(BLUP, function(x){attr(x, which="postVar")}) # take sqrt() for SEs


}