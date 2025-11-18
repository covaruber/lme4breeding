\name{predict.lmeb}
\alias{predict.lmeb}
\alias{predict.lmebreed}
\title{Predict form of a LMM fitted with lmeb/lmebreed}
\description{\code{predict} method for class \code{"lmeb"}.
}
\usage{
\method{predict}{lmeb}(object, hyperTable=NULL, classify, usePEV=FALSE, \dots)
}
\arguments{
 \item{object}{a mixed model of class \code{"lmeb"} fitted with \code{\link{lmeb}} }
 \item{hyperTable}{ a data frame with columns; variable, group, type, include,
     average. See the \code{\link{Dtable}} function to understand the format.
 
     A pure \strong{include} term means that the model matrices for that fixed or 
     random effect is filled with 1s for the positions where column names and row names match.

     An \strong{include and average} term means that the model matrices for that
     fixed or random effect is filled with 1/#1s in that row.

     A pure \strong{average} term alone means that all rows for such fixed or 
     random effect will be filled with 1/#levels in the effect.

     If a term is not considered \strong{include} or \strong{average} is then 
     totally ignored in the BLUP and SE calculation.

     The default behavior when the user doesn't provide the \code{hyperTable} is 
     to \strong{include and average} any fixed effect that is not part of classify. 
     Include any term making a perfect match with the \code{classify} argument 
     and \strong{include and average} any imperfect match with \code{classify} 
     argument (e.g., interactions).
 
 }
 \item{classify}{is a character value indicating which term we are computing
     the predictions for.}
  \item{usePEV}{is a logical value indicating whether we should use the conditional variance
  or the PEV for the computation of standard errors for the linear combination.}
 \item{\dots}{Further arguments to be passed.}
}
\details{
  This function allows to produce predictions specifying those variables that define the margins of the hypertable to be predicted (argument classify). Predictions are obtained for each combination of values of the specified variables that is present in the data set used to fit the model. See vignettes for more details. 
  
  For predicted values the pertinent design matrices X and Z together with BLUEs (b) and BLUPs (u) are multiplied and added together.
  
  \code{predicted.value = Db}
  
  \code{standard.error = D Cinv Dt}
  
 }
\value{

  \item{pvals}{the table of predictions according to the specified arguments.}
  
  \item{mapCondVar}{the map between the hypertable and the effects.}
  
  \item{hyperTable}{the table specifying the terms to include and terms to be averaged.}
  
  \item{b}{the fixed and random effects in a vector form.}
  
  \item{condVarMat}{the variance covariance for the predictions.}
  
  \item{D}{the model matrix for predictions as defined in Welham et al.(2004).}
  
  \item{classify}{the character value used to indicate which term we are computing the predictions for.}
  
  
  
 }
 \references{
  Welham, S., Cullis, B., Gogel, B., Gilmour, A., and Thompson, R. (2004). Prediction in linear mixed models. Australian and New Zealand Journal of Statistics, 46, 325 - 347.
 }
\author{Giovanny Covarrubias-Pazaran}
\seealso{\code{\link{predict}}, \code{\link{lmeb}}}
\keyword{models}
\examples{

data(DT_yatesoats, package="enhancer")
DT <- DT_yatesoats
m3 <- lmeb(Y ~ V + N + V:N +
                 (1|B) + (1|B:MP),
               data = DT)

#############################
## predict means for nitrogen
#############################
pp=predict(object=m3, classify="N")
pp$pvals
pp$hyperTable

#############################
## predict means for variety
#############################
pp=predict(object=m3, classify="V")
pp$pvals
pp$hyperTable
#############################
## predict means for nitrogen:variety
#############################
pp=predict(object=m3, classify="N:V")
pp$pvals
ht <- pp$hyperTable
ht[4,"include"]=1
ht[4,"average"]=0

pp2=predict(object=m3, classify="N:V", hyperTable=ht)
pp2$pvals

}
