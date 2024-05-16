#### "relmat" class methods


lmebreed <-
  function(formula, data, family = NULL, REML = TRUE, relmat = list(), 
           addmat=list(), 
           # initDerivs=NULL, initCov=NULL,
           control = list(), start = NULL, verbose = FALSE, 
           subset, weights, na.action, offset, contrasts = NULL,
           model = TRUE, x = TRUE, dateWarning=TRUE, returnParams=FALSE, ...)
  {
    my.date <- "2024-09-01" # expiry date
    your.date <- Sys.Date()
    ## if your month is greater than my month you are outdated
    if(dateWarning){
      if (your.date > my.date) {
        cat("Version out of date. Please update lme4breeding to the newest version using:\ninstall.packages('lme4breeding') in a new session\n Use the 'dateWarning' argument to disable the warning message.")
      }
    }
    gaus <- FALSE
    if (is.null(family)) {
      gaus <- TRUE
    } else {
      ## copied from glm()
      if (is.character(family)) {
        family <- get(family, mode = "function", envir = parent.frame())
      }
      if (is.function(family)){
        family <- family()
      } 
      if (!inherits(family, "family")) stop("unknown family type")
      gaus <- family$family == "gaussian" && family$link == "identity"
    }
    mc <- match.call()
    lmerc <- mc                         # create a call to lmer 
    # lmerc$formula <- formula; lmerc$data <- data; lmerc$control <- control
    lmerc[[1]] <- if (gaus){as.name("lmer")}else{as.name("glmer")} 
    lmerc$relmat <- NULL
    lmerc$addmat <- NULL
    lmerc$dateWarning <- NULL
    lmerc$returnParams <- NULL
    if (!gaus) {lmerc$REML <- NULL}
    if (!length(relmat) & !length(addmat))  {
      return(eval.parent(lmerc))
    }            # call [g]lmer instead
    stopifnot(is.list(relmat),        # check the relmat argument
              length(names(relmat)) == length(relmat),
              all( sapply(relmat, inherits, what = c("relmat","matrix","dtCMatrix","ddiMatrix"))  ))
    
    lmf <- eval(lmerc, parent.frame()) 
    relfac <- relmat          # copy the relmat list for relfactor
    pnms <- names(relmat)
    pnms2 <- names(addmat)
    pp <- lmf@pp
    resp <- lmf@resp
    fl <- lmf@flist
    stopifnot(all(pnms %in% names(fl)))
    asgn <- attr(fl, "assign")
    Zt <- pp$Zt # Matrix::image(Zt)  Matrix::image(as(addmat[[1]], Class="dgCMatrix"))
    ##############################
    ## replace additional matrices
    for (i in seq_along(addmat)) {
      tn0 <- which(match(pnms2[i], names(fl)) == asgn)
      for(j in 1:length(tn0)){ # diagonal and unstructured models require to multiple all matrices by the same relfactor
        ind <- (lmf@Gp)[tn0[j]:(tn0[j]+1L)]
        rowsi <- (ind[1]+1L):ind[2]
        covariate <- unlist(lmf@cnms[tn0])[j]
        if(covariate %in% names(data)){ # if is a random regression
          covariateZ <- Matrix::sparse.model.matrix(as.formula(paste("~",covariate,"-1")), data=lmf@frame)
          if(is.list(addmat[[i]])){ # user has different matrices for the same effect (e.g., indirect genetic effects)
            provZ <- addmat[[i]][[j]]
          }else{ # user has a single matrix for a given effect
            provZ <- addmat[[i]]
          }
          covariateZ <- covariateZ %*% Matrix::Matrix(1, nrow=1, ncol=ncol(provZ)) # expand the covariate
          provZt <- t(provZ*covariateZ)
          Zt[rowsi,] <- provZt[rownames(Zt[rowsi,]),]
        }else{ # is an intercept
          provZt <- t(addmat[[i]])
          Zt[rowsi,] <- provZt[rownames(Zt[rowsi,]),]
        }
      }
    }
    #############################
    ## use the relfactors
    for (i in seq_along(relmat)) {
      tn <- which(match(pnms[i], names(fl)) == asgn)
      for(j in 1:length(tn)){ # diagonal and unstructured models require to multiple all matrices by the same relfactor
        ind <- (lmf@Gp)[tn[j]:(tn[j]+1L)]
        rowsi <- (ind[1]+1L):ind[2]
        pick <- intersect( rownames(Zt), rownames(relfac[[i]])  )
        toAdd <- setdiff( rownames(relfac[[i]]), rownames(Zt) )
        if(length(pick)==0){stop(paste("The names on your relmat does not coincide with the names in your factor",pnms[i]))}
        provRelFac <- relfac[[i]][pick,pick]
        if(nrow(Zt[rowsi,]) == nrow(provRelFac)){ # regular model
          Zt[rowsi,] <- provRelFac %*% Zt[rowsi,]
        }else{ # unstructured model
          mm <- Matrix::Diagonal( length(lmf@cnms[[pnms[i]]]) )
          Zt[rowsi,] <- Matrix::kronecker(provRelFac, mm) %*% Zt[rowsi,]
        }
      }
    }
    reTrms <- list(Zt=Zt,theta=lmf@theta,Lambdat=pp$Lambdat,Lind=pp$Lind,
                   lower=lmf@lower,flist=lmf@flist,cnms=lmf@cnms, Gp=lmf@Gp)
    dfl <- list(fr=lmf@frame, X=pp$X, reTrms=reTrms, start=lmf@theta)
    if(returnParams){
      return(dfl)
    }else{
      if (gaus) {
        dfl$REML = resp$REML > 0L
        devfun <- do.call(mkLmerDevfun,dfl)
        opt <- optimizeLmer(devfun, optimizer="Nelder_Mead",...)
      } else {
        dfl$family <- family
        devfun <- do.call(mkGlmerDevfun,dfl)
        opt <- optimizeGlmer(devfun, optimizer="Nelder_Mead",...)
      }
      mm <- mkMerMod(environment(devfun), opt, reTrms, lmf@frame, mc)
      cls <- if (gaus){"lmerlmebreed"}else{"glmerlmebreed"} 
      ans <- do.call(new, list(Class=cls, relfac=relfac,
                               frame=mm@frame, flist=mm@flist, cnms=mm@cnms, Gp=mm@Gp,
                               theta=mm@theta, beta=mm@beta,u=mm@u,lower=mm@lower,
                               devcomp=mm@devcomp, pp=mm@pp,resp=mm@resp,optinfo=mm@optinfo))
      ans@call <- evalq(mc)
      return(ans)
    }
    
  }

setMethod("ranef", signature(object = "lmebreed"),
          function(object, condVar = FALSE, drop = FALSE, whichel = names(ans), relmat = TRUE, ...)
          {
            ans <- lme4::ranef(object, condVar, drop = FALSE) # as(object, "merMod")
            ans <- ans[whichel]
            if (relmat) {
              rf <- object@relfac
              for (nm in names(rf)) { # nm <- names(rf)[1]
                dm <- data.matrix(ans[[nm]])
                cn <- colnames(dm)
                rn <- rownames(dm)
                dm <- t(as.matrix( t(dm) %*% rf[[nm]][rn,rn] )) # rotate BLUPs
                colnames(dm) <- cn
                rownames(dm) <- rn
                for(kCol in 1:ncol(ans[[nm]])){
                  ans[[nm]][[kCol]] <- as.numeric(dm[,kCol])# data.frame(dm[[kCol]], check.names = FALSE)
                }
                # replace postVar if condVar=TRUE
                if (condVar){
                  # if(simpleCondVar){
                  #   tn <- which(match(nm, names(object@flist)) == attr(object@flist, "assign") )
                  #   for(j in 1:length(tn)){
                  #     ind <- (object@Gp)[tn[j]:(tn[j]+1L)]
                  #     rowsi <- (ind[1]+1L):ind[2]
                  #     if(is.list(attr(ans[[nm]], which="postVar"))){
                  #       prov <- t(rf[[nm]][rn,rn]) %*% diag(attr(ans[[nm]], which="postVar")[[j]][,,]) %*% (rf[[nm]][rn,rn])
                  #       attr(ans[[nm]], which="postVar")[[j]][,,] <- diag(prov)
                  #     }else{
                  #       prov <- t(rf[[nm]][rn,rn]) %*% diag(attr(ans[[nm]], which="postVar")[,,]) %*% (rf[[nm]][rn,rn])
                  #       attr(ans[[nm]], which="postVar")[,,] <- prov
                  #     }
                  #   }
                  # }else{
                    Ci <- getMME(object)$Ci
                    tn <- which(match(nm, names(object@flist)) == attr(object@flist, "assign") )
                    for(j in 1:length(tn)){
                      ind <- (object@Gp)[tn[j]:(tn[j]+1L)]
                      rowsi <- (ind[1]+1L):ind[2]
                      if(is.list(attr(ans[[nm]], which="postVar"))){
                        attr(ans[[nm]], which="postVar")[[j]][,,] <- diag(Ci[rowsi,rowsi])
                      }else{
                        attr(ans[[nm]], which="postVar")[,,] <- diag(Ci[rowsi,rowsi])
                      }
                    }
                  # }
                }
              }
            }
            return(ans)
          })


setMethod("fitted", signature(object = "lmebreed"),
          function(object, ...) {
            stop("fitted() not applicable to lmebreed objects")
          })


setMethod("residuals", signature(object = "lmebreed"),
          function(object, ...) {
            stop("residuals() not applicable to lmebreed objects")
          })

# setMethod("predict", signature(object = "lmebreed"),
#           function(object, Dtable=NULL, D, ...){
#             if(is.character(D)){classify <- D}else{classify="id"} # save a copy before D is overwriten
#             # complete the Dtable withnumber of effects in each term
#             xEffectN <- lapply(object$partitionsX, as.vector)
#             nz <- unlist(lapply(object$uList,function(x){nrow(x)*ncol(x)}))
#             # add a value but if there's no intercept consider it
#             lenEff <- length(xEffectN)
#             toAdd <- xEffectN[[lenEff]];
#             if(length(toAdd) > 0){ # there's a value in the last element of xEffectN
#               add <- max(toAdd) + 1 # this is our last column of fixed effects
#             }else{ # there's not a value in the last element of xEffectN
#               if(lenEff > 1){
#                 toAdd2 <- xEffectN[[lenEff-1]]
#                 add <- max(toAdd2) + 1
#               }
#             }
#             zEffectsN <- list()
#             for(i in 1:length(nz)){
#               end= add + nz[i] - 1
#               zEffectsN[[i]] <- add:end
#               add = end + 1
#             }
#             names(zEffectsN) <- names(nz)
#             effectsN = c(xEffectN,zEffectsN)
#             # fill the Dt table for rules
#             if(is.null(Dtable) & is.character(D) ){ # if user didn't provide the Dtable but D is character
#               Dtable <- object$Dtable # we extract it from the model
#               termsInDtable <- apply(data.frame(Dtable$term),1,function(xx){all.names(as.formula(paste0("~",xx)))})
#               termsInDtable <- lapply(termsInDtable, function(x){intersect(x,colnames(object$data))})
#               termsInDtable <- lapply(termsInDtable, function(x){return(unique(c(x, paste(x, collapse=":"))))})
#               # term identified
#               termsInDtableN <- unlist(lapply(termsInDtable,length))
#               pickTerm <- which( unlist(lapply(termsInDtable, function(xxx){ifelse(length(which(xxx == D)) > 0, 1, 0)})) > 0)
#               if(length(pickTerm) == 0){
#                 isInDf <- which(colnames(object$data) %in% D)
#                 if(length(isInDf) > 0){ # is in data frame but not in model
#                   stop(paste("Predict:",classify,"not in the model but present in the original dataset. You may need to provide the
#                    Dtable argument to know how to predict", classify), call. = FALSE)
#                 }else{
#                   stop(paste("Predict:",classify,"not in the model and not present in the original dataset. Please correct D."), call. = FALSE)
#                 }
#               }
#               # check if the term to predict is fixed or random
#               pickTermIsFixed = ifelse("fixed" %in% Dtable[pickTerm,"type"], TRUE,FALSE)
#               ## 1) when we predict a random effect, fixed effects are purely "average"
#               if(!pickTermIsFixed){Dtable[which(Dtable$type %in% "fixed"),"average"]=TRUE; Dtable[pickTerm,"include"]=TRUE}
#               ## 2) when we predict a fixed effect, random effects are ignored and the fixed effect is purely "include"
#               if(pickTermIsFixed){Dtable[pickTerm,"include"]=TRUE}
#               ## 3) for predicting a random effect, the interactions are ignored and only main effect is "included", then we follow 1)
#               ## 4) for a model with pure interaction trying to predict a main effect of the interaction we "include" and "average" the interaction and follow 1)
#               if(length(pickTerm) == 1){ # only one effect identified
#                 if(termsInDtableN[pickTerm] > 1){ # we are in #4 (is a pure interaction model)
#                   Dtable[pickTerm,"average"]=TRUE
#                 }
#               }else{# more than 1, there's main effect and interactions, situation #3
#                 main <- which(termsInDtableN[pickTerm] == min(termsInDtableN[pickTerm])[1])
#                 Dtable[pickTerm,"include"]=FALSE;  Dtable[pickTerm,"average"]=FALSE # reset
#                 Dtable[pickTerm[main],"include"]=TRUE
#               }
#             }
#             ## if user has provided D as a classify then we create the D matrix
#             if(is.character(D)){
#               # create model matrices to form D
#               P <- sparse.model.matrix(as.formula(paste0("~",D,"-1")), data=object$data)
#               colnames(P) <- gsub(D,"",colnames(P))
#               tP <- t(P)
#               W <- object$W
#               D = tP %*% W
#               colnames(D) <- c(rownames(object$b),rownames(object$u))
#               rd <- rownames(D)
#               cd <- colnames(D)
#               for(jRow in 1:nrow(D)){ # for each effect add 1's where missing
#                 myMatch <- which(cd == rd[jRow])
#                 if(length(myMatch) > 0){D[jRow,myMatch]=1}
#               }
#               # apply rules in Dtable
#               for(iRow in 1:nrow(Dtable)){
#                 w <- effectsN[[iRow]]
#                 # include/exclude rule
#                 if(Dtable[iRow,"include"]){ # set to 1
#                   subD <- D[,w,drop=FALSE]
#                   subD[which(subD > 0, arr.ind = TRUE)] = 1
#                   D[,w] <- subD
#                   # average rule
#                   if(Dtable[iRow,"average"]){ # set to 1
#                     # average the include set
#                     for(o in 1:nrow(subD)){
#                       v <- which(subD[o,] > 0);  subD[o,v] <- subD[o,v]/length(v)
#                     }
#                     D[,w] <- subD
#                   }
#                 }else{ # set to zero
#                   if(Dtable[iRow,"average"]){ # set to 1
#                     subD <- D[,w,drop=FALSE] + 1
#                     subD <- subD/subD
#                     subD[which(subD > 0, arr.ind = TRUE)] = subD[which(subD > 0, arr.ind = TRUE)]/ncol(subD)
#                     D[,w] <- subD
#                   }else{
#                     D[,w] <- D[,w] * 0
#                   }
#                 }
#               }
#               interceptColumn <- unique(c(grep("Intercept",rownames(object$b) ))) # ,which(rownames(object$b)=="1")
#               if(length(interceptColumn) > 0){D[,interceptColumn] = 1}
#             }else{ }# user has provided D as a matrix to do direct multiplication
#             ## calculate predictions and standard errors
#             bu <- object$bu
#             predicted.value <- D %*% bu
#             vcov <- D %*% object$Ci %*% t(D)
#             std.error <- sqrt(diag(vcov))
#             pvals <- data.frame(id=rownames(D),predicted.value=predicted.value[,1], std.error=std.error)
#             if(is.character(classify)){colnames(pvals)[1] <- classify}
#             return(list(pvals=pvals,D=D,vcov=vcov, Dtable=Dtable))
#           }
# )
