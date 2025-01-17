#### "relmat" class methods
lmebreed <-
  function(formula, data, family = NULL, REML = TRUE, relmat = list(), 
           addmat=list(), 
           control = list(), start = NULL, verbose = TRUE, 
           subset, weights, na.action, offset, contrasts = NULL,
           model = TRUE, x = TRUE, dateWarning=TRUE, returnParams=FALSE, 
           rotation=FALSE, coefOutRotation=8, ...)
  {
    my.date <- "2025-03-01" # expiry date
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
    # lmerc$formula <- formula; lmerc$data <- data; 
    # lmerc$control <- control # only if we are using it 
    ## silence additional parameters from lme4breeding that don't apply to lmer
    lmerc[[1]] <- if (gaus){as.name("lmer")}else{as.name("glmer")} 
    lmerc$relmat <- NULL
    lmerc$addmat <- NULL
    lmerc$dateWarning <- NULL
    lmerc$returnParams <- NULL
    lmerc$rotation <- NULL 
    lmerc$coefOutRotation <- NULL
    lmerc$family <- family
    ##
    if (!gaus) {lmerc$REML <- NULL}
    ## if there are no relmats or additional matrices just return th regular lmer model
    if (!length(relmat) & !length(addmat))  {
      return(eval.parent(lmerc))
    }            # call [g]lmer instead
    stopifnot(is.list(relmat),        # check the relmat argument
              length(names(relmat)) == length(relmat),
              all( sapply(relmat, inherits, what = c("relmat","matrix","Matrix"))  ))
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
    if(length(relmat) > 0){ # get cholesky factor
      if(rotation){ # if UDU decomposition + Cholesky is requested
        if(length(relmat) > 1){warning("Rotation is only reported to be accurate with one relationship matrix.", call. = FALSE)}
        udu <- umat(formula=as.formula(paste("~", paste(names(relmat), collapse = "+"))), relmat = relmat, data=data, addmat = addmat)
        # if rotation we impute
        data[,response] <- imputev(x=data[,response],method="median",by=udu$record)
        goodRecords <- 1:nrow(data)          
        newValues <- (udu$Utn[goodRecords,goodRecords] %*% data[goodRecords,response])[,1]
        # newValues <- (udu$Utn[,goodRecords] %*% data[goodRecords,response])[,1]
        outlier <- grDevices::boxplot.stats(x=newValues,coef=coefOutRotation )$out
        if(length(outlier) > 0){newValues[which(newValues %in% outlier)] = mean(newValues[which(newValues %!in% outlier)])}
        data[goodRecords,response] <- newValues
        # data[,response] <- newValues
        lmerc$data <- data
        if(verbose){message("* Rotation of response finished.")}
        for(iD in names(udu$D)){
          relmat[[iD]] <- Matrix::chol(udu$D[[iD]])
        }
        if(verbose){message("* Cholesky decomposition finished.")}
      }else{ # classical approach, just cholesky
        for (i in seq_along(relmat)) {
          relmat[[i]] <- Matrix::chol(relmat[[i]])
        }
        if(verbose){message("* Cholesky decomposition finished.")}
      }
    }
    ## END OF TRANSFORMATION
    # lmod <- do.call(lFormula, as.list(lmerc)) # basic elements to modify for our purposes
    suppressWarnings(lmod <- do.call(lFormula, as.list(lmerc)), classes = "warning")
    relfac <- relmat          # copy te relmat list for relfactor
    pnms <- names(relmat)
    pnms2 <- names(addmat)
    # resp <- lmf@resp
    fl <- lmod$reTrms$flist
    stopifnot(all(pnms %in% names(fl)))
    asgn <- attr(fl, "assign")
    Zt <- lmod$reTrms$Zt # Matrix::image(Zt)  Matrix::image(as(addmat[[1]], Class="dgCMatrix"))
    ##############################
    ## transform X if rotation is needed
    if(length(relmat) > 0){
      if(rotation){
        # toZero <- which(lmod$X == 0, arr.ind = TRUE)
        # lmod$X <- udu$Utn[goodRecords,goodRecords] %*% lmod$X
        # lmod$X[toZero] =0
        # lmod$X <- udu$Utn %*% lmod$X
        if(verbose){message("* Rotation applied to the X matrix.")}
      }
    }
    ##############################
    ## replace additional matrices
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
      if(verbose){message("* Additional matrices (addmat) added.")}
    }
    #############################
    ## use the relfactors
    for (i in seq_along(relmat)) {
      tn <- which(match(pnms[i], names(fl)) == asgn)
      for(j in 1:length(tn)){ # diagonal and unstructured models require to multiple all matrices by the same relfactor
        ind <- (lmod$reTrms$Gp)[tn[j]:(tn[j]+1L)]
        rowsi <- (ind[1]+1L):ind[2]
        pick <- intersect( rownames(Zt), rownames(relfac[[i]])  )
        toAdd <- setdiff( rownames(relfac[[i]]), rownames(Zt) )
        if(length(pick)==0){stop(paste("The names on your relmat does not coincide with the names in your factor",pnms[i],". Maybe you didn't code it as factor."))}
        provRelFac <- relfac[[i]][pick,pick]
        if(nrow(Zt[rowsi,]) == nrow(provRelFac)){ # regular model
          Zt[rowsi,] <- provRelFac %*% Zt[rowsi,]
          # Zt[rowsi,] <- t( t(udu$Utn[goodRecords,goodRecords]) %*% t(provRelFac %*% Zt[rowsi,] ) )
        }else{ # unstructured model
          mm <- Matrix::Diagonal( length(lmod$reTrms$cnms[[pnms[i]]]) )
          Zt[rowsi,] <- Matrix::kronecker(provRelFac, mm) %*% Zt[rowsi,]
          # Zt[rowsi,] <- t( t(udu$Utn[goodRecords,goodRecords]) %*%  t(Matrix::kronecker(provRelFac, mm) %*% Zt[rowsi,]) )
        }
      }
    }
    if(verbose){message("* Relfactors (relmat) applied to Z")}
    reTrms <- list(Zt=Zt,theta=lmod$reTrms$theta,Lambdat=lmod$reTrms$Lambdat,Lind=lmod$reTrms$Lind,
                   lower=lmod$reTrms$lower,flist=lmod$reTrms$flist,cnms=lmod$reTrms$cnms, Gp=lmod$reTrms$Gp)
    dfl <- list(fr=lmod$fr, X=lmod$X, reTrms=reTrms, start=lmod$reTrms$theta, control=control)
    if(length(control) == 0){
      if(gaus){ # if user calls a gaussian response family
        control <- lmerControl()
        dfl$control <- control
      }else{ # any other family response
        control <- glmerControl()
        control$optimizer <- control$optimizer[1]
        dfl$control <- control
      }
    }
    if(returnParams){
      return(dfl)
    }else{
      if(verbose){message("* Optimizing ...")}
      if (gaus) {
        dfl$REML = REML # TRUE# resp$REML > 0L
        # print(dfl$REML)
        devfun <- do.call(mkLmerDevfun, dfl )
        opt <- optimizeLmer(devfun, optimizer = dfl$control$optimizer, control = dfl$control$optCtrl, ...) # need to pass control 
      } else {
        dfl$family <- family
        devfun <- do.call(mkGlmerDevfun,dfl)
        opt <- optimizeGlmer(devfun, optimizer = dfl$control$optimizer, control = dfl$control$optCtrl,  ...) # need to pass control 
      }
      if(verbose){message("* Done!!")}
      # make results a mkMerMod object
      mm <- mkMerMod(environment(devfun), opt, reTrms, lmod$fr, mc)
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
          function(object, condVar = FALSE, drop = FALSE, whichel = names(ans), ...)
          {
            relmat <- ifelse(length(object@relfac) > 0, TRUE, FALSE)
            ans <- lme4::ranef(object, condVar, drop = FALSE) # as(object, "merMod")
            ans <- ans[whichel]
            if (relmat) { # transform back when relfac was used
              rf <- object@relfac
              for (nm in names(rf)) { # nm <- names(rf)[1]
                dm <- data.matrix(ans[[nm]])
                cn <- colnames(dm)
                rn <- rownames(dm)
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
                  X <- getME(object, "X") 
                  Ci <- getMME(object)$Ci
                  tn <- which(match(nm, names(object@flist)) == attr(object@flist, "assign") )
                  for(j in 1:length(tn)){ # j=1
                    ind <- (object@Gp)[tn[j]:(tn[j]+1L)]
                    rowsi <- ( (ind[1]+1L):ind[2] ) + ncol(X)
                    CiSub <- Ci[rowsi,rowsi]
                    if(length(object@udu) > 0){ # rotation was used
                      if( nm %in% names(object@udu$U) ){ # this only happens if there was a single relmat
                        CiSub <- (object@udu$U[[nm]][rn,rn]) %*% CiSub[rn,rn] %*% t(object@udu$U[[nm]][rn,rn])
                      }
                    }
                    if(is.list(attr(ans[[nm]], which="postVar"))){ # unstructured model
                      attr(ans[[nm]], which="postVar")[[j]][,,] <- diag(CiSub)
                    }else{ # simple model
                      attr(ans[[nm]], which="postVar")[,,] <- diag(CiSub)
                    }
                  }
                  
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
