\name{lmebreed-class}
\docType{class}
\alias{lmebreed-class}
\alias{fitted,lmebreed-method}
\alias{ranef,lmebreed-method}
\alias{residuals,lmebreed-method}
\alias{predict,lmebreed-method}
\title{Relationship-based mixed-effects model fits}
\description{A mixed-effects model fit by \code{\link{lmebreed}}.
  This class extends class \code{"\link[lme4:merMod-class]{merMod}"} class and includes one
  additional slot, \code{relfac}, which is a list of (left) Cholesky
  factors of the relationship matrices derived from 
  \code{"\linkS4class{lmebreed}"} objects.}
\section{Objects from the Class}{Objects are created by calls to the
  \code{\link{lmebreed}} function.}
\section{Slots}{
  \describe{
    \item{\code{relfac}:}{A list of relationship matrix factors.  All
      other slots are inherited from class \code{"\link[lme4:merMod-class]{merMod}"}.}
    \item{\code{udu}:}{A list of eigen decomposition elements.  All
      other slots are inherited from class \code{"\link[lme4:merMod-class]{merMod}"}.}
  }
}
\section{Extends}{
Class \code{"\link[lme4:merMod-class]{merMod}"}, directly.
}
\section{Methods}{
  \describe{
    \item{fitted}{\code{signature(object = "lmebreed")}: actually a
      non-method in that \code{fitted} doesn't apply to such objects
      because of the pre-whitening.}
    \item{ranef}{\code{signature(object = "lmebreed")}: back-transforms BLUPs and
      their conditional variances when models include
      the relationship between levels of random effects as returned for the object
      viewed as a \code{"\link[lme4:merMod-class]{merMod})"} object.}
    \item{residuals}{\code{signature(object = "lmebreed")}: also a
      non-method for the same reason as \code{fitted}}
    \item{predict}{\code{signature(object = "lmebreed")}: a function to compute 
     linear combinations of fixed and random effects as defined by Welham et al. 
     (2004). It has the following arguments:
     
     
     \code{predict(object, hyperTable=NULL, classify=NULL, ... )}
     
     \code{object}: is the model fit with the \code{\link{lmebreed}} function.
     
     \code{hyperTable}: a data frame with columns; variable, group, type, include,
     average. See the \code{\link{Dtable}} function to understand the format.
     
     A pure "include" term means that the model matrices for that fixed or random effect is filled with 1's for the positions where column names and row names match.

     An "include and average" term means that the model matrices for that fixed or random effect is filled with 1/#1's in that row.

     An "average" term alone means that all rows for such fixed or random effect will be filled with 1/#levels in the effect.

     If a term is not considered "include" or "average" is then totally ignored in the BLUP and SE calculation.

     The default behavior when the user doesn't provide the hyperTable is to include and average any fixed effect that is not part of classify. Include any term making a perfect match with the classify argument and include and average any imperfect match with classify argument (e.g., interactions).
     
     \code{classify}: is a character value indicating which term we are computing
     the predictions for.


    }
    
  }
}
\seealso{\code{\link{lmebreed}}
}
\examples{
showClass("lmebreed")

data(DT_example)
DT <- DT_example
A <- A_example

## Compound simmetry (CS) model
ans1 <- lmebreed(Yield~Env + (1|Name) + (1|Env:Name),
                 data=DT)
                 
pp <- predict(ans1, classify="Name")
head(pp$pvals)

image(pp$D)

pp$hyperTable

}
\keyword{classes}
