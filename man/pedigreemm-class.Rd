\name{pedigreemm-class}
\docType{class}
\alias{pedigreemm-class}
\alias{fitted,pedigreemm-method}
\alias{ranef,pedigreemm-method}
\alias{residuals,pedigreemm-method}
\title{Pedigree-based mixed-effects model fits}
\description{A mixed-effects model fit by \code{\link{pedigreemm}}.
  This class extends class \code{"\link[lme4:merMod-class]{merMod}"} class and includes one
  additional slot, \code{relfac}, which is a list of (left) Cholesky
  factors of the relationship matrices derived from 
  \code{"\linkS4class{pedigree}"} objects.}
\section{Objects from the Class}{Objects are created by calls to the
  \code{\link{pedigreemm}} function.}
\section{Slots}{
  \describe{
    \item{\code{relfac}:}{A list of relationship matrix factors.  All
      other slots are inherited from class \code{"\link[lme4:merMod-class]{merMod}"}.}
  }
}
\section{Extends}{
Class \code{"\link[lme4:merMod-class]{merMod}"}, directly.
}
\section{Methods}{
  \describe{
    \item{fitted}{\code{signature(object = "pedigreemm")}: actually a
      non-method in that \code{fitted} doesn't apply to such objects
      because of the pre-whitening.}
    \item{ranef}{\code{signature(object = "pedigreemm")}: incorporates
      the pedigree into the random effects as returned for the object
      viewed as a \code{"\link[lme4:merMod-class]{merMod})"} object.}
    \item{residuals}{\code{signature(object = "pedigreemm")}: also a
      non-method for the same reason as \code{fitted}}
  }
}
\seealso{\code{\link{pedigreemm}}
}
\examples{
showClass("pedigreemm")
}
\keyword{classes}
