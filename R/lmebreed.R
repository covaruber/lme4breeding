#### "relmat" class methods
lmebreed <-  function(formula, data, REML = TRUE, control = list(), start = NULL, 
                      verbose = 1L, subset, weights, na.action, offset, contrasts = NULL,
                      calc.derivs=FALSE,
                      # new params
                      family = NULL, relmat = list(),  addmat=list(), trace=1L,
                      dateWarning=TRUE, rotation=FALSE, rotationK=NULL, coefOutRotation=8, 
                      returnParams=FALSE, returnMod=FALSE, ...)
{
  my.date <- "2026-01-01" # expiry date
  your.date <- Sys.Date()
  ## if your month is greater than my month you are outdated
  if(dateWarning){
    if (your.date > my.date) {
      warning("Version out of date. Please update lme4breeding to the newest version using:\ninstall.packages('lme4breeding') in a new session\n Use the 'dateWarning' argument to disable the warning message.")
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
    if (!inherits(family, "family")) {stop("unknown family type")}
    gaus <- family$family == "gaussian" && family$link == "identity"
  }
  mc <- match.call()
  lmerc <- mc                         # create a call to lmer 
  ## >>>>>>>>>>>> create control if user didn't specify it for the 3 things we want to force
  if(length(control)==0){
    if(gaus){
      control <- lmerControl(
        calc.derivs = FALSE,
        restart_edge = FALSE,
        check.nobs.vs.nlev = "ignore",
        check.nobs.vs.rankZ = "ignore",
        check.nobs.vs.nRE="ignore"
      )
    }else{
      control <- glmerControl(
        calc.derivs = FALSE,
        restart_edge = FALSE,
        check.nobs.vs.nlev = "ignore",
        check.nobs.vs.rankZ = "ignore",
        check.nobs.vs.nRE="ignore"
      )
    }
  }else{ # if user provides a control force this 3 controls
    control$checkControl$check.nobs.vs.nlev = "ignore"
    control$checkControl$check.nobs.vs.rankZ = "ignore"
    control$checkControl$check.nobs.vs.nRE="ignore"
    control$calc.derivs = calc.derivs
  }
  # lmerc$formula <- formula; lmerc$data <- data; 
  # lmerc$control <- control # only if we are using it 
  ## silence additional parameters from lme4breeding that don't apply to lmer
  lmerc[[1]] <- if (gaus){as.name("lmer")}else{as.name("glmer")} 
  lmerc$family <- family
  lmerc$control <- control
  lmerc$start <- start
  if (!gaus) {lmerc$REML <- NULL}
  ## if there are no relmats or additional matrices just return th regular lmer model
  if (!length(relmat) & !length(addmat))  {
    lmerc$relmat <- NULL # remove relmat from the match call to avoid errors when evaluating the call
    lmerc$addmat <- NULL # remove relmat from the match call
    lmerc$trace <- NULL # remove relmat from the match call
    lmerc$dateWarning=NULL; lmerc$rotation=NULL; lmerc$rotationK=NULL
    lmerc$coefOutRotation=NULL; lmerc$returnParams=NULL; lmerc$returnMod=NULL
    mm <- eval.parent(lmerc)
    cls <- if (gaus){"lmerlmebreed"}else{"glmerlmebreed"} 
    # put it in a lmebreed object
    ans <- do.call(new, list(Class=cls, relfac=list(), udu=list(), #goodRecords=goodRecords,
                             frame=mm@frame, flist=mm@flist, cnms=mm@cnms, Gp=mm@Gp,
                             theta=mm@theta, beta=mm@beta,u=mm@u,lower=mm@lower,
                             devcomp=mm@devcomp, pp=mm@pp,resp=mm@resp,optinfo=mm@optinfo))
    ans@call <- evalq(mc)
    return(ans)
  }            # call [g]lmer instead
  
  ## DO TRANSFORMATION BEFORE EVALUATING THE CALL
  '%!in%' <- function(x,y)!('%in%'(x,y)) 
  response <- all.vars(formula)[1]
  # control to ignore relmats if there's no match with formula vars
  if(length(relmat) > 0){
    if( length(intersect(names(relmat), all.vars(formula) )) == 0){
      relmat <- list()
    }
  }
  # control to ignore addmats if there's no match with formula vars
  if(length(addmat) > 0){
    if( length(intersect(names(addmat), all.vars(formula) )) == 0){
      addmat <- list()
    }
  }
  if( response %!in% colnames(data) ){stop("Response selected in your formula is not part of the dataset provided.", call. = FALSE)}
  goodRecords <- which(!is.na(data[,response]))
  udu <- list()
  # >>>>>>>>> get cholesky factor (and new response if rotation)
  if(length(relmat) > 0){ 
    if(rotation){ # if UDU decomposition + Cholesky is requested
      if(length(relmat) > 1){warning("Rotation is only reported to be accurate with one relationship matrix. Make sure you are using the same relationship matrix for the different random effects for the rotation approach.", call. = FALSE)}
      for(iRel in 1:length(relmat)){
        idsOrdered <- as.character(unique(data[,names(relmat)[iRel]])) # when we rotate we need to have relmat already ordered before creating the matrices
        relmat[[iRel]] = relmat[[iRel]][ idsOrdered , idsOrdered ]
      }
      # only the first relmat will be used so if more, the rotation will only work if it is the same relmat in the next random effects
      udu <- umat(formula=as.formula(paste("~", paste(names(relmat), collapse = "+"))), relmat = relmat, 
                  data=data, addmat = addmat, k=rotationK)
      # if rotation we impute the response
      data[,response] <- imputev(x=data[,response],method="median")#, by=data[udu$effect])
      newValues <- udu$Utn %*% Matrix::Matrix(data[,response])
      newValues <- newValues[,1]
      outlier <- grDevices::boxplot.stats(x=newValues,coef=coefOutRotation )$out
      if(length(outlier) > 0){newValues[which(newValues %in% outlier)] = mean(newValues[which(newValues %!in% outlier)])}
      data[,response] <- newValues
      if(trace){message(magenta("* Rotation of response finished."))}
      for(iD in names(udu$D)){
        relmat[[iD]] <- Matrix::chol(udu$D[[iD]])
      }
      udu$newValues <- newValues
      lmerc$data <- data
      if(trace){message(magenta("* Cholesky decomposition finished."))}
    }else{ # classical approach, just cholesky
      for (i in seq_along(relmat)) {
        relmat[[i]] <- Matrix::chol(relmat[[i]])
      }
      if(trace){message(magenta("* Cholesky decomposition finished."))}
    }
  }
  
  stopifnot(is.list(relmat),        # check the relmat argument
            length(names(relmat)) == length(relmat),
            all( sapply(relmat, inherits, what = c("relmat","matrix","Matrix"))  ))
  
  # >>>>>>>>> use the match call to parse the data and formula
  lmerc[[1]] <- if (gaus){as.name("lFormula")}else{as.name("glFormula")} # change from model to lFormula
  lmerc$relmat <- NULL # remove relmat from the match call
  lmerc$addmat <- NULL # remove relmat from the match call
  lmerc$dateWarning <- NULL # remove relmat from the match call
  lmerc$rotation <- NULL # remove relmat from the match call
  lmerc$rotationK <- NULL # remove relmat from the match call
  lmerc$coefOutRotation <- NULL # remove relmat from the match call
  lmerc$returnParams <- NULL # remove relmat from the match call
  lmerc$returnMod <- NULL # remove relmat from the match call
  
  suppressWarnings( lmod <- eval.parent(lmerc) , classes = "warning") # necesary objects from lFormula
  
  # >>>>>>>>> extract relevant matrices from the lFormula object
  relfac <- relmat          # copy te relmat list for relfactor
  pnms <- names(relmat)
  pnms2 <- names(addmat)
  fl <- lmod$reTrms$flist
  stopifnot(all(pnms %in% names(fl)))
  asgn <- attr(fl, "assign")
  Zt <- lmod$reTrms$Zt
  ##############################
  ## transform X if rotation is needed
  if(length(relmat) > 0){
    if(rotation){
      if(ncol(udu$Utn) != nrow(lmod$X)){stop("Rotation approach requires your dataset to be balanced and imputed.")}
      lmod$X <- (udu$Utn %*% lmod$X) * lmod$X
      if(trace){message(magenta("* Rotation applied to the X matrix."))}
      udu$Utn <- NULL # avoid storing a big matrix after the multiplication
    }
  }
  ##############################
  
  # >>>>>>>>>> apply addmat (additional matrices)
  for (i in seq_along(addmat)) {
    tn0 <- which(match(pnms2[i], names(fl)) == asgn)
    for(j in 1:length(tn0)){ # diagonal and unstructured models require to multiple all matrices by the same relfactor
      ind <- (lmod$reTrms$Gp)[tn0[j]:(tn0[j]+1L)]
      rowsi <- (ind[1]+1L):ind[2]
      covariate <- unlist(lmod$reTrms$cnms[tn0])[j]
      if(covariate %in% names(data)){ # if is a random regression
        covariateZ <- Matrix::sparse.model.matrix(as.formula(paste("~",covariate,"-1")), data=lmod$fr)
        if(is.list(addmat[[i]])){ # user has different matrices for the same effect (e.g., indirect genetic effects)
          provZ <- addmat[[i]][[j]][goodRecords,]
        }else{ # user has a single matrix for a given effect
          provZ <- addmat[[i]][goodRecords,]
        }
        covariateZ <- covariateZ %*% Matrix::Matrix(1, nrow=1, ncol=ncol(provZ)) # expand the covariate
        provZt <- t(provZ*covariateZ)
        Zt[rowsi,] <- provZt[rownames(Zt[rowsi,]),]
      }else{ # is an intercept
        provZt <- t(addmat[[i]][goodRecords,])
        Zt[rowsi,] <- provZt[rownames(Zt[rowsi,]),]
      }
    }
    if(trace){message(magenta("* Additional matrices (addmat) added."))}
  }
  
  # >>>>>>>>> time to apply the relmat
  for (i in seq_along(relmat)) { # for each relationship matrix
    tn <- which(match(pnms[i], names(fl)) == asgn) # match relmat names with random effects names
    for(j in 1:length(tn)){ # for each random effect matching this relationship matrix (diagonal and unstructured models require to multiple all incidence matrices by the same relfactor)
      ind <- (lmod$reTrms$Gp)[tn[j]:(tn[j]+1L)] # which columns match this random effect
      rowsi <- (ind[1]+1L):ind[2] # first to last column from Z
      colnamesRelFac <- colnames(relfac[[i]])
      
      if( mean(table(colnamesRelFac)) > 1 ){  # is this complex because we may have a relationship matrix with repeated names
        toBeRemoved <- character()
        namesProvRelFac <- character() 
        foundV <- numeric()
        for(p in which( rownames(Zt) %in% rownames(relfac[[i]]) ) ){ # p=1
          found <- which(colnamesRelFac %in% rownames(Zt)[p])
          found <- setdiff(found, toBeRemoved)[1]
          toBeRemoved <- c(toBeRemoved, found[1])
          if(!is.na(found)){
            foundV <- c(foundV,found)
            namesProvRelFac <- c(namesProvRelFac, colnamesRelFac[found] )
          }
        }
        provRelFac <- relfac[[i]][foundV,foundV] 
        colnames(provRelFac) <- rownames(provRelFac) <- namesProvRelFac
        relfac[[i]] <- provRelFac
      }else{
        pick <- intersect( rownames(Zt), rownames(relfac[[i]])  ) # match names in relmat and Z matrix
        if(length(pick)==0){stop(paste("The names on your relmat does not coincide with the names in your factor",pnms[i],". Maybe you didn't code it as factor?"))}
        provRelFac <- relfac[[i]][pick,pick] # only pick portion of relmat that coincides with Z
      }
      if(nrow(Zt[rowsi,]) == nrow(provRelFac)){ # regular model (single random intercept)
        provRelFac <- as(as(as( provRelFac,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
        ZtL <- list() # we have to do this because filling by rows a Column-oriented matrix is extremely slow so it is faster to cut and paste
        if(min(rowsi) > 1){ZtL[[1]] <- Zt[1:(min(rowsi)-1),]}
        ZtL[[2]] <- provRelFac %*% Zt[rowsi,] 
        if(max(rowsi) < nrow(Zt)){ZtL[[3]] <- Zt[(max(rowsi)+1):nrow(Zt),]}
        Zt <- do.call(rbind, ZtL)
      }else{ # complex model (multiple random intercepts)
        mm <- Matrix::Diagonal( length(lmod$reTrms$cnms[[pnms[i]]]) )
        if(length(rowsi) != ncol(provRelFac)*ncol(mm) ){stop(paste("Relationship matrix dimensions of ",pnms[i],"do not conform with the random effect, please review."), call. = FALSE)}
        provRelFac <- as(as(as( provRelFac,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
        ZtL <- list()
        if(min(rowsi) > 1){ZtL[[1]] <- Zt[1:(min(rowsi)-1),]}
        ZtL[[2]] <- Matrix::kronecker(provRelFac, mm, make.dimnames = TRUE) %*% Zt[rowsi,]
        rownames(ZtL[[2]]) <- rownames(Zt[rowsi,])
        if(max(rowsi) < nrow(Zt)){ZtL[[3]] <- Zt[(max(rowsi)+1):nrow(Zt),]}
        Zt <- do.call(rbind, ZtL)
      }
    }
  }
  
  if(trace){message(magenta("* Relfactors (relmat) applied to Z"))}
  
  reTrms <- list(Zt=Zt,theta=if(is.null(start)){lmod$reTrms$theta}else{start},Lambdat=lmod$reTrms$Lambdat,Lind=lmod$reTrms$Lind,
                 lower=lmod$reTrms$lower,flist=lmod$reTrms$flist,cnms=lmod$reTrms$cnms, Gp=lmod$reTrms$Gp)
  lmod <- list(fr=lmod$fr, X=lmod$X, reTrms=reTrms, formula=formula, verbose=verbose,
               start=if(is.null(start)){lmod$reTrms$theta}else{start},
               control=control)
  # print(str(lmod))
  if(length(control) == 0){
    if(gaus){ # if user calls a gaussian response family
      control <- lmerControl()
      lmod$control <- control
    }else{ # any other family response
      control <- glmerControl()
      control$optimizer <- control$optimizer[1]
      lmod$control <- control
    }
  }
  if(returnParams){ # if user only wants the incidence matrices
    return(lmod)
  }else{
    if(trace){message(magenta("* Optimizing ..."))}
    if (gaus) { # gaussian distribution
      lmod$REML = REML # TRUE# resp$REML > 0L
      suppressWarnings( devfun <- do.call(mkLmerDevfun, lmod ), classes = "warning") # creates a deviance function
      if(returnMod){ # user wants to force variance components without fitting a model
        opt <- list(par = start, fval = devfun(start), feval = 1, conv = 0)
      }else{ # # user wants to optimize the varcomp optimizer
        suppressWarnings( opt <- optimizeLmer(devfun, optimizer = lmod$control$optimizer, 
                                              control = lmod$control$optCtrl, 
                                              verbose=lmod$verbose, 
                                              calc.derivs=lmod$control$calc.derivs,
                                              restart_edge=lmod$control$restart_edge,
                                              boundary.tol=lmod$control$boundary.tol,
                                              use.last.params=lmod$control$use.last.params, ...)   , classes = "warning") # need to pass control 
      } 
    } else { # exponential family of distributions
      lmod$family <- family
      suppressWarnings( devfun <- do.call(mkGlmerDevfun,lmod) , classes = "warning") # creates a deviance function
      if(returnMod){ # user wants to force variance components without optimizing
        opt <- list(par = start, fval = devfun(start), feval = 1, conv = 0)
      }else{ # user wants to optimize the varcomp optimizer
        suppressWarnings( opt <- optimizeGlmer(devfun, optimizer = lmod$control$optimizer, 
                                               control = lmod$control$optCtrl, 
                                               verbose=lmod$verbose, 
                                               calc.derivs=lmod$control$calc.derivs,
                                               restart_edge=lmod$control$restart_edge,
                                               boundary.tol=lmod$control$boundary.tol,
                                               use.last.params=lmod$control$use.last.params,  ...)  ) # need to pass control 
      } 
    }
    if(trace){message(magenta("* Done!!"))}
    # make results in a mkMerMod object format
    suppressWarnings( mm <- mkMerMod(environment(devfun), opt, lmod$reTrms, lmod$fr, mc), classes = "warning" )
    cls <- if (gaus){"lmerlmebreed"}else{"glmerlmebreed"} 
    ans <- do.call(new, list(Class=cls, relfac=relfac, udu=udu, #goodRecords=goodRecords,
                             frame=mm@frame, flist=mm@flist, cnms=mm@cnms, Gp=mm@Gp,
                             theta=mm@theta, beta=mm@beta,u=mm@u,lower=mm@lower,
                             devcomp=mm@devcomp, pp=mm@pp,resp=mm@resp,optinfo=mm@optinfo))
    ans@call <- evalq(mc)
    return(ans)
  }
  
}

setMethod("ranef", signature(object = "lmebreed"),
          function(object, condVar = TRUE, drop = FALSE, whichel = names(ans), includePEV=TRUE, ...)  {
            # print("new")
            relmat <- ifelse(length(object@relfac) > 0, TRUE, FALSE)
            ans <- lme4::ranef(object, condVar, drop = FALSE) # as(object, "merMod")
            ans <- ans[whichel]
            if(condVar){
              namesCi <- mkMmeIndex(object) # rbind(namesBlue, namesBlup)
              Ci <- getCi(object)
            }
            if (relmat) { # transform back when relfac was used
              rf <- object@relfac
              for (nm in names(rf)) { # for each random effect # nm <- names(rf)[1]
                dm <- data.matrix(ans[[nm]])
                cn <- colnames(dm)
                rn <- rownames(dm)
                message(magenta(paste("Rotating back BLUPs by transpose of Cholesky (u=L'u*) for:", nm)))
                dm <- t(as.matrix( t(dm) %*% rf[[nm]][rn,rn] )) # rotate BLUPs by the relfactor
                # rotate one more time if a UDU rotation was used in the response
                if(length(object@udu) > 0){ # rotation was used
                  if( nm %in% names(object@udu$U) ){ # this only happens if there was a single relmat
                    dm <- object@udu$U[[nm]][rn,rn] %*% dm
                  }
                }
                colnames(dm) <- cn
                rownames(dm) <- rn
                for(kCol in 1:ncol(ans[[nm]])){ # put in a nice shape
                  ans[[nm]][[kCol]] <- as.numeric(dm[,kCol])# data.frame(dm[[kCol]], check.names = FALSE)
                }
                # replace postVar if condVar=TRUE
                if (condVar){
                  namesCiNm <- namesCi[which(namesCi$group == nm),]
                  intercepts <- unique(namesCiNm$variable)
                  for(j in 1:length(intercepts)){ # iInter = 1 # intercepts[1]
                    iInter <- intercepts[j]
                    v <- namesCiNm[which(namesCiNm$variable == iInter), "index"]
                    if(is.list(attr(ans[[nm]], which="postVar"))){ # diagonal model
                      attr(ans[[nm]], which="postVar")[[j]][,,] <- diag(Ci)[v]
                    }else{ # unstructured model # fill the diagonal element corresponding to the intercept level
                      attr(ans[[nm]], which="postVar")[j,j,] <- diag(Ci)[v]
                    }
                  }
                  
                }
              }
              
            }
            # before returning store the Ci and its names
            if(all(c(includePEV, condVar))){ # if both are TRUE add the PEV
              attr(ans, which="PEV") = Ci
              attr(ans, which="namesCi") = namesCi
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

setMethod("predict", signature(object = "lmebreed"),
          function(object, hyperTable=NULL, classify=NULL, ...)  {
            
            if(is.null(classify)){
              stop("Please provide the classify argument to build the D matrix.", call. = FALSE )
            }
            '%!in%' <- function(x,y)!('%in%'(x,y))
            if(is.null(hyperTable)){
              message(magenta("hyperTable argument not provided. Building a hyper table based on the classify argument. Please check the output to ensure the different effects have been grouped and average as you expected."))
              hyperTable <- Dtable(object)
              hyperTable$include[which(hyperTable$group %in% classify)]=1
              hyperTable$include[which(hyperTable$group %!in% classify & hyperTable$type == "fixed")]=1
              hyperTable$average[which(hyperTable$group %!in% classify & hyperTable$type == "fixed")]=1
            }
            # get all information from the model
            BLUP <- ranef(object, condVar=TRUE)
            # get D table information
            namesCi <- attr(BLUP, which="namesCi")
            # get inverse of coefficient matrix
            Ci <- attr(BLUP, which="PEV")
            # get coefficients
            b <- c( 
              as.vector( fixef(object) ),
              unlist( lapply(BLUP, function(x){
                unlist(lapply(x, function(y){as.vector(y)}), use.names = FALSE )
              }), use.names = FALSE )
            )
            # create the D matrix of linear combination
            
            classifys <- unique(c( all.vars(as.formula(paste("~",classify,"-1"))), classify))
            Zd <- Matrix::sparse.model.matrix(as.formula(paste("~",classify,"-1")), data=object@frame ) 
            levsOr <- colnames(Zd)
            for(iClassify in classifys){
              colnames(Zd) <- gsub(iClassify,"",colnames(Zd))
            }
            levs <- colnames(Zd)
            D <- Matrix::Matrix(0, ncol=length(b), nrow=length(levs))
            rownames(D) <- levs
            # now add the rules specified in the hyperTable
            for(iRow in 1:nrow(hyperTable)){ # iRow=1
              iVar <- hyperTable[iRow,"variable"]
              iGroup <- hyperTable[iRow,"group"]
              if(hyperTable[iRow,"include"]>0){
                # if(hyperTable[iRow,"group"] %in% classifys){ # if this is a classify
                  for (jRow in 1:nrow(D)) { # jRow=3
                    v <- which(namesCi[,"variable"]==iVar & namesCi[,"group"]==iGroup )
                    w <- which(namesCi[,"level"] %in%
                                 c(
                                   rownames(D)[jRow], 
                                   levsOr[jRow],
                                   strsplit(rownames(D)[jRow], split = ":")[[1]],
                                   strsplit(levsOr[jRow], split = ":")[[1]],
                                   "(Intercept)"
                                 )
                    )
                    myMatch <- intersect(v,w)
                    # print(myMatch)
                    if (length(myMatch) > 0) {
                      D[jRow, myMatch] = 1
                    }
                  }
                # }else{ # this group is not part of classify
                #   v <- which(namesCi[,"variable"]==iVar & namesCi[,"group"]==iGroup )
                #   D[,v]=1
                # }
              }
              if(hyperTable[iRow,"average"]>0){
                D[,v]= D[,v]/length(v)
              }
              
            }
            # compute the predicted values and std errors
            predicted.value <- D %*% b
            vcov <- D %*% Ci %*% t(D)
            std.error <- sqrt(diag(vcov))
            pvals <- data.frame(id = rownames(D),
                                predicted.value = predicted.value[,1], 
                                std.error = std.error)
            # compile results
            ans <- list(pvals=pvals, b=b, Ci=Ci, D=D, hyperTable=hyperTable, classify=classify )
            return(ans)
          })



