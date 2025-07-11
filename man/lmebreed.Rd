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
           start = NULL, verbose = TRUE, subset, weights,
           na.action, offset, contrasts = NULL, model = TRUE,
           x = TRUE, dateWarning=TRUE, rotation=FALSE, 
           coefOutRotation=8, returnParams=FALSE,
           returnMod=FALSE, ...)
}
\arguments{
  \item{relmat}{an optional named list of relationship matrices (not the inverse). 
    Internally the Cholesky decomposition of those matrices will be computed.
    The names of the elements must correspond to the names of grouping factors for
    random-effects terms in the \code{formula} argument.}
  \item{addmat}{an optional named list of customized incidence matrices.
    The names of the elements must correspond to the names of grouping factors for
    random-effects terms in the \code{formula} argument. Depending on the use-case
    the element in the list may be a single matrix or a list of matrices. Please see
    examples and vignettes to learn how to use it.}
  \item{dateWarning}{an logical value indicating if you want to be warned when a new 
  version of lme4breeding is available on CRAN.}
  \item{rotation}{an logical value indicating if you want to compute the eigen decomposition 
  of the relationship matrix to rotate y and X and accelerate the computation. See details.}
  \item{coefOutRotation}{a numeric value denoting the inter-quantile outlier coefficient to
  be used in the rotation of the response when using the eigen decomposition to avoid 
  overshooting.}
    \item{returnParams}{a logical value indicating if you want to only get the incidence 
  matrices of the model.}
    \item{returnMod}{a logical value indicating if you want to force the model to use 
    your customized variance components without estimating them. The variance components
    should be provided in the \code{start} argument. If you want to provide the variance
    components from a previous model you can get the initial values by running:
    
    \code{ getME(mix1, 'sigma') }
    
    which returns a vector with theta values.
    
    }
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
  Each relmat is the relationship matrix restricted 
  to the observed levels and applied to the model matrix for that term. The incidence
  matrices in the addmat argument must match the dimensions of the final fit (pay 
  attention to missing data in responses).
  
  It is important to remember that when you use the \code{relmat} argument you are providing
  the square root of the relationship matrix and to recover the correct BLUPs for those effects
  you need to use the \code{\link[lme4]{ranef}} function which internally multiple those BLUPs the 
  square root of the relationship matrix one more time to recover the correct BLUPs.
  
  The argument \code{rotation} applies the eigen decomposition proposed by Lee and Van der Werf in 2016
  and makes the genetic evaluation totally sparse leading to incredible gains in speed compared 
  to the classical approach. Internally, the eigen decomposition UDU' is carried in the relationship 
  matrix. The U matrix is then taken to the n x n level (n being the number of records), and post-multiplied
  by a matrix of records presence (n x n) using the element-wise multiplication of two matrices (Schur product). 
  By default is not activated since this may not provide the exact same variance components than other software due to
  numerical reasons. If you would like to obtain the exact same variance components than other software you will
  have to keep \code{rotation=FALSE}. This will slow down considerably the speed. Normally when the rotation is 
  activated and variance components differ slightly with other software they will still provide highly similar 
  estimates at the speed of hundreds or thousands of times faster. Please consider this.
  
  Additional useful functions are; \code{\link{tps}} for spatial kernels, \code{\link{rrm}} 
  for reduced rank matrices, \code{\link{atcg1234}} for conversion of genetic markers, 
  \code{\link{overlay}} for overlay matrices, \code{\link{reshape}} for moving wide 
  format multi-trait datasets into long format, \code{\link{fillData}} for balancing 
  datasets for two variables.
  
  When using the optimizer argument inside the \code{\link[lme4]{lmerControl}} keep 
  in mind that the number of iterations is called differently depending on the optimizer. 
  For \code{\link[lme4]{Nelder_Mead}}, bobyqa and \code{\link[lme4]{nlminbwrap}} is
  called "maxfun" whereas for \code{\link[lme4]{nloptwrap}} is called "maxeval". 
  This should be passed inside a list in the \code{optCtrl} argument. For example:
  
\code{lmebreed(... , 
         control = lmerControl(
           optimizer="Nelder_Mead",
           optCtrl=list(maxfun=100)
         ), ...
        )}
        
  To predict values for unobserved levels you will need to impute the data and update
  your model with the new dataset and the initial starting values:
  
  \code{
    newModel <- update(oldModel, returnMod = TRUE,
                       start=getME(oldModel, 'sigma'), 
                       data=imputedData)
  }
  
  
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

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

Lee & Van der Werf (2016). MTG2: an efficient algorithm for multivariate linear mixed model analysis based on genomic information. Bioinformatics, 32(9), 1420-1422.
   
}
\author{
  Giovanny Covarrubias-Pazaran
}
\examples{

data(DT_example)
DT <- DT_example
A <- A_example

ansMain <- lmebreed(Yield ~ Env + (1|Name),
                        relmat = list(Name = A ),
                        data=DT)
vc <- VarCorr(ansMain); print(vc,comp=c("Variance"))

BLUP <- ranef(ansMain, condVar=TRUE)$Name
SEs <- attr(BLUP, which="postVar")[,,]

}
\keyword{models}
