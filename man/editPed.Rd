\name{editPed}
\alias{editPed}
\title{Complete and Order a Pedigree}
\usage{editPed(sire, dam, label, verbose)}
\description{This function helps to prepare a pedigree to generate a pedigree object}
\details{
  The function takes a vector of sires, another for dams and a final one for subjects
  all of the same length, convert them to character. If there are dams or sires not
  declared as subjects the function generates them. Finally, it orders the pedigree. 
  The output can be used to build a pedigree object \code{ped}
}
%\references{}
\value{A data frame with strings as characters. All subjects are in the label column, 
       and all subjects will appear in this column before appering as sires or dams.}
\arguments{
    \item{sire}{a vector (with some \code{NA} entries) with the father IDs}
    \item{dam}{similarly as \code{sire} for the \dQuote{mother} of each entry. 
               The vector must be of the same length than the one for the sire}
    \item{label}{a vector with the subjects id. Giving a unique ID for the corresponding 
      entry. The length as \code{sire} and \code{dam} should be the same}
    \item{verbose}{logical entry inquiring whether to print line that the program is evaluating.
                    The default is FALSE.} 
}
\examples{
  #(1)
  pede<-data.frame(sire=as.character(c(NA,NA,NA,NA,NA,1,3,5,6,4,8,1,10,8)),
             dam= as.character(c(NA,NA,NA,NA,NA,2,2,NA,7,7,NA,9,9,13)),
             label=as.character(1:14))
   #scrambled original pedigree:
   (pede<- pede[sample(replace=FALSE, 1:14),]  )
   (pede<- editPed(sire=pede$sire, dam= pede$dam, label=pede$label)) 
   ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))
#################################################################################################
   #(2) With missing labels
   pede<-data.frame(sire=as.character(c(NA,1,3,5,6,4,8,1,10,8)),
             dam= as.character(c(NA,2,2,NA,7,7,NA,9,9,13)),
             label=as.character(5:14))
   #scrambled original pedigree:
   (pede<- pede[sample(replace=FALSE, 1:10),]  )
   (pede<- editPed(sire=pede$sire, dam= pede$dam, label=pede$label)) 
   ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))

#################################################################################################
   #(2) A larger pedigree
 #Useing pedCows pedigree
# str(pedCows)
# pede<-data.frame(id=pedCows@label, sire=pedCows@sire, dam=pedCows@dam)
# pede<-pede[sample(1:nrow(pede),replace=FALSE),]
# pede<- editPed(sire=pede$sire, dam=pede$dam, label=pede$id)
# ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))
}
\keyword{array}
\keyword{algebra}
