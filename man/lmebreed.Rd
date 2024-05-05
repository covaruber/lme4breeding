\name{lmebreed}
\alias{lmebreed}
\title{Fit mixed-effects models incorporating relationships}
\description{
  Fit linear or generalized linear mixed models incorporating the
  effects of relationships.
}
\usage{
lmebreed(formula, data, family = NULL, REML = TRUE,
           relmat = list(), addmat = list(), control = list(),
           start = NULL, verbose = FALSE, subset, weights,
           na.action, offset, contrasts = NULL, model = TRUE,
           x = TRUE, ...)
}
\arguments{
  \item{relmat}{an optional named list of relationship factors (relfactors) 
     or matrix decompositions (e.g., cholesky).
    The names of the elements must correspond to the names of grouping factors for
    random-effects terms in the \code{formula} argument.}
  \item{addmat}{an optional named list of customized incidence matrices.
    The names of the elements must correspond to the names of grouping factors for
    random-effects terms in the \code{formula} argument. Depending on the use-case
    the element in the list may be a single matrix or a list of matrices. Please see
    examples and vignettes to learn how to use it.}
  \item{formula}{as in \code{\link[lme4]{lmer}}}
  \item{data}{as in \code{\link[lme4]{lmer}}}
  \item{family}{as in \code{\link[lme4]{glmer}}}
  \item{REML}{as in \code{\link[lme4]{lmer}}}
  \item{control}{as in \code{\link[lme4]{lmer}}}
  \item{start}{as in \code{\link[lme4]{lmer}}}
  \item{verbose}{as in \code{\link[lme4]{lmer}}}
  \item{subset}{as in \code{\link[lme4]{lmer}}}
  \item{weights}{as in \code{\link[lme4]{lmer}}}
  \item{na.action}{as in \code{\link[lme4]{lmer}}}
  \item{offset}{as in \code{\link[lme4]{lmer}}}
  \item{contrasts}{as in \code{\link[lme4]{lmer}}}
  \item{model}{as in \code{\link[lme4]{lmer}}}
  \item{x}{as in \code{\link[lme4]{lmer}}}
  \item{\dots}{as in \code{\link[lme4]{lmer}}}
}
\details{
  All arguments to this function are the same as those to the function
  \code{\link[lme4]{lmer}} except \code{relmat} and \code{addmat} which must be 
  named lists.  Each name must correspond to the name of a grouping factor in a
  random-effects term in the \code{formula}.  The observed levels
  of that factor must be contained in the rownames and columnames of the relmat.  
  Each relmat is the (left) Cholesky factor of the relationship matrix restricted 
  to the observed levels and applied to the model matrix for that term. The incidence
  matrices in the addmat argument must match the dimensions of the final fit (pay 
  attention to missing data in responses).
  
  It is important to remember that when you use the \code{relmat} argument you are providing
  the square root of the relationship matrix and to recover the correct BLUPs for those effects
  you need to use the \code{\link{ranef}} function which internally multiple those BLUPs the 
  square root of the relationship matrix one more time to recover the correct BLUPs.
  
  
  \strong{Example Datasets}

The package has been equiped with several datasets to learn how to use the lme4breeding package:

* \code{\link{DT_halfdiallel}}, \code{\link{DT_fulldiallel}} and \code{\link{DT_mohring}} datasets have examples to fit half and full diallel designs.

* \code{\link{DT_h2}} to calculate heritability

* \code{\link{DT_cornhybrids}} and \code{\link{DT_technow}} datasets to perform genomic prediction in hybrid single crosses

* \code{\link{DT_wheat}} dataset to do genomic prediction in single crosses in species displaying only additive effects.

* \code{\link{DT_cpdata}} dataset to fit genomic prediction models within a biparental population coming from 2 highly heterozygous parents including additive, dominance and epistatic effects.

* \code{\link{DT_polyploid}} to fit genomic prediction and GWAS analysis in polyploids.

* \code{\link{DT_gryphon}} data contains an example of an animal model including pedigree information.

* \code{\link{DT_btdata}} dataset contains an animal (birds) model.

* \code{\link{DT_legendre}} simulated dataset for random regression model.

* \code{\link{DT_sleepstudy}} dataset to know how to translate lme4 models to sommer models.

* \code{\link{DT_ige}} dataset to show how to fit indirect genetic effect models.

}
\value{
  a \code{\linkS4class{lmebreed}} object.
}
\references{

Covarrubias-Pazaran G. Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 2016, 11(6): doi:10.1371/journal.pone.0156744

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.
   
}
\author{
  Giovanny Covarrubias-Pazaran
}
\examples{

data(DT_example)
DT <- DT_example
A <- A_example

ansMain <- lmebreed(Yield ~ Env + (1|Name),
                        relmat = list(Name = Matrix::chol(A) ),
                        data=DT)
vc <- VarCorr(ansMain); print(vc,comp=c("Variance"))

}
\keyword{models}
