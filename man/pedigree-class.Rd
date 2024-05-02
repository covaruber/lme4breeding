\name{pedigree-class}
\docType{class}
\alias{pedigree-class}
\alias{pedinbred-class}
%\alias{chol,pedigree-method}
\alias{coerce,pedigree,sparseMatrix-method}
\alias{coerce,pedigree,data.frame-method}
\alias{head,pedigree-method}
\alias{show,pedigree-method}
\alias{tail,pedigree-method}
%
\title{Class "pedigree"}
\description{
  Objects of class \code{"pedigree"} represent a set of
  individuals that can have two parents including their parent-child
  relations.  The terminology has been taken from cattle breeding.
  The \code{"pedinbred"} class is an extension of the pedigree class
  with an additional slot of the inbreeding coefficients.
}
\section{Objects from the Class}{
  Objects in the \code{"pedigree"} class can be created by calls of the
  form \code{new("pedigree", ...)}, or more conveniently,
  \code{\link{pedigree}(sire= ., dam = ., label =.)}.

  Objects of the \code{"pedinbred"} class are created by coercing a
  pedigree to class \code{"pedinbred"}.
}
\section{Slots}{
  \describe{
    \item{\code{sire}:}{integer vector (with some \code{NA} entries),
      denoting a \emph{previous} entry in the pedigree corresponding to
      the current entry's \dQuote{father}.}
    \item{\code{dam}:}{similarly as \code{sire} for the \dQuote{mother}
      of each entry.}
    \item{\code{label}:}{a \code{"character"} vector of the same length
      as \code{sire} and \code{dam} giving a unique ID for the
      corresponding entry.}
    \item{\code{F}:}{(class \code{"pedinbred"} only) a numeric vector of
      inbreeding coefficients.}
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "pedigree", to = "sparseMatrix")}:
      returns a sparse, unit lower-triangular matrix which is the inverse of the
      "L" part of the "LDL'" form of the Cholesky factorization of the
      relationship matrix.  All non-zero elements below the diagonal
      are -0.5.}
    \item{coerce}{\code{signature(from = "pedigree", to = "data.frame")}: ... }
    \item{head}{\code{signature(x = "pedigree")}: ... }
    \item{show}{\code{signature(object = "pedigree")}: ... }
    \item{tail}{\code{signature(x = "pedigree")}: ... }
  }
}
\references{
  R. A. Mrode, \emph{Linear Models for the Prediction of Animal Breeding
    Values, 2nd ed}, CABI Publishing, 2005.

  2010. A.I. Vazquez, D.M. Bates, G.J.M. Rosa, D. Gianola and K.A. Weigel.
  Technical Note: An R package for fitting generalized linear mixed models
  in animal breeding. Journal of Animal Science, 88:497-504.

}
\seealso{
  \code{\link{pedigree}}, \code{\link{inbreeding}}
}
\examples{
## Rather use, pedigree()!   The following is "raw code":
## Example from chapter 2 of Mrode (2005)
p1 <- new("pedigree",
          sire = as.integer(c(NA,NA,1, 1,4,5)),
          dam  = as.integer(c(NA,NA,2,NA,3,2)),
          label = as.character(1:6))
p1
(dtc <- as(p1, "sparseMatrix"))   # T-inverse in Mrode's notation
solve(dtc)
inbreeding(p1)
}
\keyword{classes}
