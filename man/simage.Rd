\name{simage}
\alias{simage}
\title{Image of sparsity between two variables}
\description{\code{image} for sparsity.
}
\usage{
simage(data, Var1=NULL, Var2=NULL, \dots)
}
\arguments{
 \item{data}{model data of class \code{"data.frame"}}
 \item{Var1}{variable to set in the x axis}
 \item{Var2}{variable to set in the y axis}
 \item{\dots}{Further arguments to be passed to the image function.}
}
\value{
  vector of image
 }
\author{Giovanny Covarrubias}
\examples{

# row x column combinations
data(DT_cpdata)
DT <- DT_cpdata
head(DT)
simage(data=DT, Var1 = "Rowf", Var2 = "Colf")

# dent x flint combinations
data(DT_technow)
DT <- DT_technow
head(DT)
simage(data=DT, Var1 = "dent", Var2 = "flint")

# male x female combinations
data(DT_halfdiallel)
DT <- DT_halfdiallel
head(DT)
simage(data=DT, Var1 = "male", Var2 = "female")

}
\seealso{\code{\link{image}}, \code{\link{lmeb}}}
\keyword{models}

