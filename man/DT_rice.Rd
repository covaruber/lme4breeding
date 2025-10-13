\name{DT_rice}
\alias{DT_rice}
\alias{GT_rice}
\alias{GTn_rice}
\title{Rice lines dataset}
\description{

Information from a collection of 413 rice lines.  The DT_rice data set is from 
Rice Diversity Org. Program. The lines are genotyped with 36,901 SNP markers and phenotyped for more than 30 traits. This data set was included in the package to play with it. If using it for your research make sure you cite the original publication from Zhao et al.(2011).

}

\usage{
  data(DT_rice)
}

\format{
 RicePheno contains the phenotypes
 RiceGeno contains genotypes letter code
 RiceGenoN contains the genotypes in numerical code using atcg1234 converter function
}

\source{
  Rice Diversity Organization
  http://www.ricediversity.org/data/index.cfm.
}

\references{

Keyan Zhao, Chih-Wei Tung, Georgia C. Eizenga, Mark H. Wright, M. Liakat Ali, Adam H. Price, Gareth J. Norton, M. Rafiqul Islam, Andy Reynolds, Jason Mezey, Anna M. McClung, Carlos D. Bustamante & Susan R. McCouch (2011). Genome-wide association mapping reveals a rich genetic architecture of complex traits in Oryza sativa. Nat Comm 2:467 DOI: 10.1038/ncomms1467, Published Online 13 Sep 2011.

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}
\examples{

data(DT_rice)
DT <- DT_rice
GT <- GT_rice
GTn <- GTn_rice
head(DT)
M <- atcg1234(GT)

\donttest{

### univariate model
A <- A.mat(M$M)
A <- A + diag(1e-4, ncol(A), ncol(A))
mix <- lmebreed(Protein.content ~ (1|geno),
                relmat = list(geno=A),
                data=DT)
vc <- VarCorr(mix); print(vc,comp=c("Variance"))
sigma(mix)^2 # error variance

### multi-trait model
traits <- c("Flowering.time.at.Arkansas" ,"Seed.volume", "Protein.content")
for(iTrait in traits){DT[,iTrait] <- scale(DT[,iTrait])}
DTL <- reshape(DT[,c("geno", traits)], idvar = "geno", varying = traits,
        v.names = "value", direction = "long",
        timevar = "trait", times = traits )
DTL <- DTL[with(DTL, order(trait)), ]
head(DTL)

system.time(
  mix <- lmebreed(value ~ (0+trait|geno),
                  relmat = list(geno=A),
                  rotation = TRUE,
                  data=DTL)
)
vc <- VarCorr(mix); print(vc,comp=c("Variance"))
vc$geno
sigma(mix)^2 # error variance

BLUP <- ranef(mix, condVar=TRUE)
PEV <- lapply(BLUP, function(x){attr(x, which="postVar")}) # take sqrt() for SEs

}


}
\keyword{datasets}
\seealso{The core function of the package \code{\link{lmebreed}} }
