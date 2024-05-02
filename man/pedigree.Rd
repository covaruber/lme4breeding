\name{pedigree}
\alias{pedigree}
\title{Pedigree Constructor}
\description{
  Construct an object of class \code{"pedigree"}, more
  conveniently than by \code{new("pedigree", ....)}.
}
\usage{
pedigree(sire, dam, label)
}
\arguments{
    \item{sire}{numeric vector (with some \code{NA} entries) of integer IDs,
      denoting a \emph{previous} entry in the pedigree corresponding to
      the current entry's \dQuote{father}.}
    \item{dam}{similarly as \code{sire} for the \dQuote{mother} of each entry.}
    \item{label}{a vector coercable to \code{"character"} of the same
      length as \code{sire} and \code{dam} giving a unique ID for the
      corresponding entry.}
}
\references{
   2010. A.I. Vazquez, D.M. Bates, G.J.M. Rosa, D. Gianola and K.A. Weigel.
   Technical Note: An R package for fitting generalized linear mixed models
   in animal breeding. Journal of Animal Science, 88:497-504.
}
\value{
  an object of formal class \code{"\linkS4class{pedigree}"}.
}
\seealso{the \code{\linkS4class{pedigree}} class.}
\examples{
example("pedigree-class") ## 'p1' pedigree object `the hard way'

ped <- pedigree(sire = c(NA,NA,1, 1,4,5),
                dam  = c(NA,NA,2,NA,3,2), label= 1:6)
## note that 'label' is coerced to character automatically
ped
stopifnot(identical(ped, p1))
}
\keyword{misc}

