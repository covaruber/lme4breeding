\name{pedigreemm}
\alias{pedigreemm}
\title{Fit mixed-effects models incorporating pedigrees}
\description{
  Fit linear or generalized linear mixed models incorporating the
  effects of a pedigree.
}
\usage{
pedigreemm(formula, data, family = NULL, REML = TRUE,
           pedigree = list(), control = list(),
           start = NULL, verbose = FALSE, subset, weights,
           na.action, offset, contrasts = NULL, model = TRUE,
           x = TRUE, ...)
}
\arguments{
  \item{pedigree}{a named list of relationship factors (relfactors) 
     or matrix decompositions (e.g., cholesky).
    The names must correspond to the names of grouping factors for
    random-effects terms in the \code{formula} argument.}
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
  \code{\link[lme4]{lmer}} except \code{pedigree} which must be a named list of
  \code{\linkS4class{pedigree}} objects.  Each name (frequently there is
  only one) must correspond to the name of a grouping factor in a
  random-effects term in the \code{formula}.  The observed levels
  of that factor must be contained in the pedigree.  For each pedigree
  the (left) Cholesky factor of the
  relationship matrix restricted to the observed levels is calculated
  using \code{\link{relfactor}} and applied to the model matrix for that
  term.
}
\value{
  a \code{\linkS4class{pedigreemm}} object.
}
\references{
   2010. A.I. Vazquez, D.M. Bates, G.J.M. Rosa, D. Gianola and K.A. Weigel. 
   Technical Note: An R package for fitting generalized linear mixed models 
   in animal breeding. Journal of Animal Science, 88:497-504.
}
\seealso{\code{\linkS4class{pedigreemm}}, \code{\linkS4class{pedigree}},
  \code{\link{relfactor}}.}
\examples{

  p1 <- new("pedigree",
           sire = as.integer(c(NA,NA,1, 1,4,5)),
           dam  = as.integer(c(NA,NA,2,NA,3,2)),
           label = as.character(1:6))
  A<-getA(p1)
  cholA<-chol(A)  
  varU<-0.4; varE<-0.6; rep<-20
  n<-rep*6
  set.seed(108)
  bStar<- rnorm(6, sd=sqrt(varU))
  b<-crossprod(as.matrix(cholA),bStar)
  ID <- rep(1:6, each=rep)
  e0<-rnorm(n, sd=sqrt(varE))
  y<-b[ID]+e0
  fm1 <- pedigreemm(y ~ (1|ID) , pedigree = list(ID =  relfactor(p1, unique(ID) )) )
  table(y01<-ifelse(y<1.3,0,1))
  fm2 <- pedigreemm(y01 ~ (1|ID) , pedigree = list(ID =  relfactor(p1, unique(ID) )), family = 'binomial')
}
\keyword{models}
