\name{lmeb-class}
\docType{class}
\alias{lmeb-class}
\alias{fitted,lmeb-method}
\alias{ranef,lmeb-method}
\alias{residuals,lmeb-method}
\title{Relationship-based mixed-effects model fits}
\description{A mixed-effects model fit by \code{\link{lmeb}}.
  This class extends class \code{"\link[lme4:merMod-class]{merMod}"} class and includes one
  additional slot, \code{relfac}, which is a list of (left) Cholesky
  factors of the relationship matrices derived from 
  \code{"\linkS4class{lmeb}"} objects.}
\section{Objects from the Class}{Objects are created by calls to the
  \code{\link{lmeb}} function.}
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
    \item{fitted}{\code{signature(object = "lmeb")}: actually a
      non-method in that \code{fitted} doesn't apply to such objects
      because of the pre-whitening.}
    \item{ranef}{\code{signature(object = "lmeb")}: back-transforms BLUPs and
      their conditional variances when models include
      the relationship between levels of random effects as returned for the object
      viewed as a \code{"\link[lme4:merMod-class]{merMod})"} object.}
    \item{residuals}{\code{signature(object = "lmeb")}: also a
      non-method for the same reason as \code{fitted}}
    
  }
}
\seealso{\code{\link{lmeb}}
}
\examples{
showClass("lmeb")

data(DT_example)
DT <- DT_example
A <- A_example

## Compound simmetry (CS) model
ans1 <- lmeb(Yield~Env + (1|Name) + (1|Env:Name),
                 data=DT)
                 
# fitted(ans1)
# residuals(ans1)
# rr <- ranef(ans1)

}
\keyword{classes}
