\name{lmeb}
\alias{lmeb}
\alias{lmebreed}
\title{Fit breeding-related mixed-effects models}
\description{
  Fits linear or generalized linear mixed models incorporating known relationships 
  (e.g., genetic relationship matrices) and customized random effects (e.g., overlay models).
  Big-data genomic models can be fitted using the eigen decomposition proposed by Lee 
  and Van der Werf (2006).
}
\usage{
lmeb(formula, data, REML = TRUE, control = list(), start = NULL, 
         verbose = 1L, subset, weights, na.action, offset, contrasts = NULL, 
         calc.derivs = FALSE, nIters=100,
         # lmeb special params
         family = NULL, relmat = list(), addmat = list(), trace=1L,
         dateWarning=TRUE, rotation=FALSE,  rotationK=NULL, coefOutRotation=8, 
         returnFormula=FALSE, suppressOpt=FALSE, ...)
}
\arguments{

  \item{formula}{as in \code{\link[lme4]{lmer}}}
  
  \item{data}{as in \code{\link[lme4]{lmer}}}
  
  \item{REML}{as in \code{\link[lme4]{lmer}}}
  
  \item{control}{as in \code{\link[lme4]{lmer}}}
  
  \item{start}{as in \code{\link[lme4]{lmer}}}
  
  \item{verbose}{as in \code{\link[lme4]{lmer}}}
  
  \item{subset}{as in \code{\link[lme4]{lmer}}}
  
  \item{weights}{as in \code{\link[lme4]{lmer}}}
  
  \item{na.action}{as in \code{\link[lme4]{lmer}}}
  
  \item{offset}{as in \code{\link[lme4]{lmer}}}
  
  \item{contrasts}{as in \code{\link[lme4]{lmer}}}
  
  \item{calc.derivs}{as in \code{\link[lme4]{lmerControl}}. The intention of having
  this argument at the level of main argument is to set it as FALSE as default 
  to speed up computation. This means that in lme4breeding if you pass calc.derivs
  inside the lmerControl object it won't have any effect. You need to set it using the 
  argument in the lmeb function.}
  
  \item{nIters}{the number of iterations used by the \code{\link[lme4]{optimizeGlmer}} 
  and \code{\link[lme4]{optimizeLmer}} functions. Internally, it replaces the values of
  the \code{maxfun} and \code{maxeval} arguments in the \code{\link[lme4]{lmerControl}} 
  optCtrl passed.This means that in lme4breeding if you pass calc.derivs
  inside the lmerControl object it won't have any effect. You need to set it using the 
  argument in the lmeb function. }
  
  \item{family}{as in \code{\link[lme4]{glmer}}}
  
  \item{relmat}{an optional named list of relationship matrices between levels of a
    given random effect (not the inverse). 
    Internally the Cholesky decomposition of those matrices is computed to adjust the
    incidence matrices. The names of the elements in the list must correspond to 
    the names of slope factors for random-effects terms in the \code{formula} argument.}
    
  \item{addmat}{an optional named list of customized incidence matrices.
    The names of the elements must correspond to the names of grouping factors for
    random-effects terms in the \code{formula} argument. Depending on the use-case
    the element in the list may be a single matrix or a list of matrices. Please see
    examples and vignettes to learn how to use it.}
    
  \item{trace}{integer scalar. If > 0 verbose output is generated during the progress 
  of the model fit.}
    
  \item{dateWarning}{a logical value indicating if you want to be warned when a new 
  version of lme4breeding is available on CRAN. Default is TRUE.}
  
  \item{rotation}{a logical value indicating if you want to compute the eigen decomposition 
  of the relationship matrix to rotate the response vector y and the fixed effect matrix X 
  in order to accelerate the computation (Lee and Vander Werf, 2016). This argument requires the dataset to be balanced and without
  missing data for the slope variable, intercept variables, and the response involved in the model.
  Also make sure that your dataset is ordered by intercepts and slopes (e.g., \code{dt = dt[with(dt, order(intercept, slope)), ]}).
  See details to understand more about this argument and vignettes for examples.}
  
  \item{rotationK}{an integer value indicating the number of eigen vectors to compute
  when the rotation argument is set to TRUE. By default all eigen vectors are computed.}
  
  \item{coefOutRotation}{a numeric value denoting the inter-quantile outlier coefficient to
  be used in the rotation of the response when using the eigen decomposition to avoid 
  overshooting.}
  
  \item{returnFormula}{a logical value indicating if you want to only get the results from
  the \code{\link[lme4]{lFormula}} after the relmat and addmat arguments have been applied.
  This is normally just used for debugging.}
  
  \item{suppressOpt}{a logical value indicating if you want to force the model to use 
    your customized variance components without estimating them. It skips the 
    optimize(g)Lmer step to force the customized variance components. The variance components
    should be provided in the \code{start} argument. If you want to provide the variance
    components from a previous model you can get the initial values by running:
    
    \code{ getME(mix1, 'sigma') }
    
    which returns a vector with theta values.
    
    }
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
  
  When you use the \code{relmat} argument the square root of the relationship matrix
  will be computed internally. Therefore, to recover the correct BLUPs for those effects
  you need to use the \code{\link[lme4]{ranef}} function which internally multiple the
  obtained BLUPs by the square root of the relationship matrix one more time to recover the correct BLUPs.
  
  The argument \code{rotation} applies the eigen decomposition proposed by Lee and Van der Werf in 2016
  and makes the genetic evaluation totally sparse leading to incredible gains in speed compared 
  to the classical approach. Internally, the eigen decomposition UDU' is carried in the relationship 
  matrix. The U matrix is then taken to the n x n level (n being the number of records), and post-multiplied
  by a matrix of records presence (n x n) using the element-wise multiplication of two matrices (Schur product). 
  By default is not activated since this may not provide the exact same variance components than other software due to
  numerical reasons. If you would like to obtain the exact same variance components than other software you will
  have to keep \code{rotation=FALSE}. This will slow down considerably the speed. Normally when the rotation is 
  activated and variance components differ slightly with other software they will still provide extremely similar 
  estimates at the speed of hundreds or thousands of times faster. Please consider this.
  
  Additional useful functions are; \code{\link{tps}} for spatial kernels, \code{\link{rrm}} 
  for reduced rank matrices, \code{\link{atcg1234}} for conversion of genetic markers, 
  \code{\link{overlay}} for overlay matrices, \code{\link{reshape}} for moving wide 
  format multi-trait datasets into long format.
  
  Normally, using the optimizer argument inside the \code{\link[lme4]{lmerControl}} keep 
  in mind that the number of iterations is called differently depending on the optimizer. 
  For \code{\link[lme4]{Nelder_Mead}}, bobyqa and \code{\link[lme4]{nlminbwrap}} is
  called "maxfun" whereas for \code{\link[lme4]{nloptwrap}} is called "maxeval". 
  This should be passed inside a list in the \code{optCtrl} argument. For example:
  
\code{lmeb(... , 
         control = lmerControl(
           optimizer="Nelder_Mead",
           optCtrl=list(maxfun=100)
         ), ...
        )}
        
  But here, we have created the \code{nIters} argument to control the number of iterations
  in the optimizer. To predict values for unobserved levels you will need to impute the data and update
  your model with the new dataset and the initial starting values:
  
  \code{
    newModel <- update(oldModel, suppressOpt = TRUE,
                       start=getME(oldModel, 'sigma'), 
                       data=imputedData)
  }
  
  It is also worth mentioning that when the user does not provide the (g)lmerControl
  we take the freedom to specify it to avoid some common warning or error messages such
  as calculating the second derivatives for big models or not allowing to have a single 
  record per treatment level:
  
  \code{
  control <- (g)lmerControl(
        calc.derivs = FALSE,
        restart_edge = FALSE, 
        check.nobs.vs.nlev = "ignore", 
        check.nobs.vs.rankZ = "ignore", 
        check.nobs.vs.nRE="ignore" 
      )
  }
  
  This is important to notice because if the user decides to specify its own controls
  then the user should also consider some of these arguments that are the defaults in lmeb.
  
  \strong{Methods}
  
  Some important methods to keep in mind are:
  
  \code{\link[lme4]{fixef}}: allows you to extract fixed coefficients (BLUEs) 
  
  \code{\link[stats]{vcov}}: allows you to extract the standard errors of fixed effects
  
  \code{\link[lme4]{ranef}}: allows you to extract random coefficients (BLUPs) and their standard errors
  
  \code{\link{condVarRotated}}: allows you to extract the predicted error variance (PEV) matrix from a model
  
  \code{\link{mkMmeIndex}}: alllows you to extract the indices of the different fixed and random coefficients
  
  \code{\link{predict}}: allows you to extract fixed coefficients (BLUEs)
  
  \strong{Example Datasets}

The package has been equiped with several datasets to learn how to use the lme4breeding package:

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
\value{
  a \code{\linkS4class{lmeb}} object.
}
\references{

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

Lee & Van der Werf (2016). MTG2: an efficient algorithm for multivariate linear mixed model analysis based on genomic information. Bioinformatics, 32(9), 1420-1422.
   
}
\author{
  Giovanny Covarrubias-Pazaran
}
\examples{


nInds=300
nMarks=1000
nEnvs=2

#random population of nInds lines with nMarks markers
M <- matrix(rep(0,nInds*nMarks),nInds,nMarks)
for (i in 1:nInds) {
  M[i,] <- ifelse(runif(nMarks)<0.5,-1,1)
}
# compute the relationship matrix
MMT <- tcrossprod(M)
m <- sum(diag(MMT))/nrow(MMT)
MMT <- MMT/m
MMT <- MMT + diag(1e-05, ncol(MMT), ncol(MMT))
# get phenotypes
yl <- el <- il <- list()
for(iEnv in 1:nEnvs){
  #random phenotypes
  u <- rnorm(nMarks)
  g <- as.vector(crossprod(t(M),u))
  h2 <- 0.5  #heritability
  envEffect <- rnorm(1, mean = abs(rnorm(1)))
  yl[[iEnv]] <- g + rnorm(nInds,mean=0,sd=sqrt((1-h2)/h2*var(g))) + envEffect
  el[[iEnv]] <- rep(paste0("e",iEnv), length(g))
  il[[iEnv]] <-  paste0("i",1:length(g))
}
y <- unlist(yl)
env <- unlist(el)
ind <- unlist(il)
colnames(MMT) <- rownames(MMT) <- paste0("i",1:length(g))
# fit the mixed model with rotation
mix <- lmeb(y~env + (1|ind), relmat=list(ind=MMT), rotation = TRUE )
vc <- VarCorr(mix); print(vc,comp=c("Variance"))
# retrieve parameters
BLUP <- ranef(mix, condVar=TRUE)
condVAR <- lapply(BLUP, function(x){attr(x, which="postVar")}) # take sqrt() for SEs


}
\keyword{models}
