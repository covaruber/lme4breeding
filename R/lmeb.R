#### "relmat" class methods
lmebreed <-  lmeb <- function(formula, data, REML = TRUE, control = list(), start = NULL, 
                      verbose = 1L, subset, weights, na.action, offset, contrasts = NULL,
                      calc.derivs=FALSE, nIters=100,
                      # new params
                      family = NULL, relmat = list(),  addmat=list(), trace=1L,
                      dateWarning=TRUE, rotation=FALSE, rotationK=NULL, coefOutRotation=8, 
                      returnFormula=FALSE, suppressOpt=FALSE, ...)
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
  lmerc <- mc  # create a call to lmer (we need it twice)
  
  ## >>>>>>>>>>>>
  ## >>>>>>>>>>>> add addmat variables
  if( any( names(addmat) %in% all.vars(formula) ) ){ # at least one addmat
    for(iAddMat in names(addmat)){ # iAddMat = names(addmat)[1]
      if(!missing(data)){
        if(is.list(addmat[[iAddMat]])){ # indirect genetic effects so addmat is a list nested in another list
          Zaddmat <- addmat[[iAddMat]][[1]]
        }else{Zaddmat <- addmat[[iAddMat]]} # simple list of addmats
        nTimes <- nrow(data)/ncol(Zaddmat) + 2
        data[,iAddMat] <- (rep(colnames(Zaddmat), nTimes ))[1:nrow(data)]
        Zaddmat <- NULL
      }else{
        checkExistIAddMat <- exists(iAddMat)
        if(!checkExistIAddMat){ # if doesn't exist the variable in the environment create
          if(is.list(addmat[[iAddMat]])){ # indirect genetic effects so addmat is a list nested in another list
            Zaddmat <- addmat[[iAddMat]][[1]]
          }else{Zaddmat <- addmat[[iAddMat]]} # simple list of addmats
          respy <- get(all.vars(formula)[1]) # get response variable
          nTimes <- length(respy)/ncol(Zaddmat) + 2 # how many times to repeate a vector
          newVariable <- (rep(colnames(Zaddmat), nTimes ))[1:length(respy)]
          assign(iAddMat, newVariable); Zaddmat <- NULL
        }else{stop(paste("Your variable",iAddMat,"does not exist in the environment and you have not provided a data argument."), call. = FALSE)}
      }
    }
  }
  ## >>>>>>>>>>>> 
  ## >>>>>>>>>>>> create a new formula (lme4 cannot interpret properly || )
  '%!in%' <- function(x,y)!('%in%'(x,y))
  if(!missing(data)){
    classDT <- unlist(lapply(data, class))
    data$unitsR <- 1:nrow(data)
  }else{
    vars <- all.vars(formula)
    unitsR <- 1:length(get(vars[1]))
    classDT <- unlist(lapply(as.list(vars), function(x){class(get(x))})); 
    names(classDT) <- vars
  }
  j <- paste(deparse(formula), collapse="")
  match0 <- gregexpr("\\(([^()]|\\([^()]*\\))*\\)",j,perl=TRUE)
  extracted <- regmatches(j,match0)[[1]]
  randomTerms <- gsub("^\\(|\\)$","", extracted)
  for(h in 1:3){randomTerms <- gsub("\\)","", sub(".*?\\(","",randomTerms) )}
  randomTerms0 <- randomTerms
  for(i in 1:length(randomTerms)){ # i=1 # for each random term
    ithRandomTerm <- strsplit(randomTerms[i], split="[|]")[[1]]
    intercept <- ithRandomTerm[1]
    intercept <- strsplit(intercept,"[+]")[[1]]
    intercept <- gsub(" ","",intercept)
    slope <- ithRandomTerm[-c(1)]
    slope <- gsub(" ","",slope)
    for(k in 1:length(intercept)){ # k=2 # for each intercept int1+int2+int3 | slope
      interceptK <- intercept[k]
      if(!missing(data)){
        checkExistInterK <- length(which(interceptK %in% names(classDT))) > 0
      }else{
        checkExistInterK <- exists(interceptK)
      }
      if(checkExistInterK){ # if the kth intercept is part of the model.frame or is in the environment
        if( classDT[interceptK] %in% c("factor","character") ){ # if is a character of factor get levels and add
          
          if(!missing(data)){ # we can add the dummy variable to the data
            variableForInterK <- gsub("[^a-zA-Z0-9._]", "", data[,interceptK]) # we remove special characters from intercept variable except dots or underscores
            # variableForInterK <- gsub("[+-]","",data[,interceptK])
            Z <- smm(variableForInterK) # add new dummy columns
            for(l in 1:ncol(Z)){data[,colnames(Z)[l]] <- Z[,l]}
          }else{ # we have to create the variables and put them in the environment
            # variableForInterK <- get(interceptK)
            variableForInterK <- gsub("[^a-zA-Z0-9._]", "", get(interceptK)) # we remove special characters from intercept variable except dots or underscores
            Z <- smm(variableForInterK)
            for(l in 1:ncol(Z)){assign(colnames(Z)[l], Z[,l] )}
          }
          levsIntercept <- sort(as.character(unique(variableForInterK)))
          if("unitsR" %in% slope){ # if this is a random effect for unitsR
            unitsR <- NA
            for(m in levsIntercept){ # m = levsIntercept[1]
              sam <- which(variableForInterK == m)
              unitsR[sam] <- 1:length(sam)
            }; unitsR <- as.factor(unitsR)
            if(!missing(data)){data$unitsR <- unitsR}
          } # "IHYB23Rattray"
        }else{
          levsIntercept <- interceptK
          if("unitsR" %in% slope){
            unitsR <- as.factor(1:length(variableForInterK))
            if(!missing(data)){data$unitsR <- unitsR}
          }
        } # if kth intercept is a numeric covariate already
      }else{ # if is not part of the model frame (e.g., 0 or 1)
        levsIntercept <- interceptK
        if(("unitsR" %in% slope) & (levsIntercept %!in% c("0","1"))){
          unitsR <- as.factor(1:length(variableForInterK))
          if(!missing(data)){data$unitsR <- unitsR}
        }
      }
      randomTerms[i] <- gsub( interceptK, paste(levsIntercept,collapse = " + ") , randomTerms[i] )
    }
  }
  newJ <- j # copy the formula
  for(h in 1:length(randomTerms)){ # replace with new terms
    newJ <- gsub(randomTerms0[h],randomTerms[h],newJ, fixed = TRUE)
  }
  lmerc$formula <- as.formula(newJ) # save the new formula in the call
  if(!missing(data)){lmerc$data <- data} # save the new dataset in the call
  
  ## >>>>>>>>>>>>
  ## >>>>>>>>>>>> create control if user didn't specify it for the 3 things we want to force
  if(length(control)==0){
    if(gaus){
      control <- lmerControl(
        calc.derivs = FALSE,
        restart_edge = FALSE,
        check.nobs.vs.nlev = "ignore",
        check.nobs.vs.rankZ = "ignore",
        check.nobs.vs.nRE="ignore",
        optCtrl = list(maxfun = nIters, maxeval = nIters)
      )
    }else{
      control <- glmerControl(
        calc.derivs = FALSE,
        restart_edge = FALSE,
        check.nobs.vs.nlev = "ignore",
        check.nobs.vs.rankZ = "ignore",
        check.nobs.vs.nRE="ignore",
        optCtrl = list(maxfun = nIters, maxeval = nIters)
      )
    }
  }else{ # if user provides a control force these controls
    control$checkControl$check.nobs.vs.nlev = "ignore"
    control$checkControl$check.nobs.vs.rankZ = "ignore"
    control$checkControl$check.nobs.vs.nRE="ignore"
    control$calc.derivs = calc.derivs
    control$optCtrl$maxfun=nIters
    control$optCtrl$maxeval=nIters
  }
  # lmerc$formula <- formula; lmerc$data <- data; 
  # lmerc$control <- control # only if we are using it 
  ## silence additional parameters from lme4breeding that don't apply to lmer
  lmerc[[1]] <- if (gaus){as.name("lmer")}else{as.name("glmer")} 
  lmerc$family <- family
  lmerc$control <- control
  lmerc$start <- start
  lmerc$verbose <- verbose # remove relmat from the match call
  if (!gaus) {lmerc$REML <- NULL}
  ## if there are no relmats or additional matrices just return th regular lmer model
  if (!length(relmat) & !length(addmat))  {
    lmerc$nIters <- NULL # remove relmat
    lmerc$relmat <- NULL # remove relmat from the match call to avoid errors when evaluating the call
    lmerc$addmat <- NULL # remove relmat from the match call
    lmerc$trace <- NULL # remove relmat from the match call
    lmerc$dateWarning=NULL; lmerc$rotation=NULL; lmerc$rotationK=NULL
    lmerc$coefOutRotation=NULL; lmerc$returnFormula=NULL; lmerc$suppressOpt=NULL
    if(returnFormula){
      lmerc[[1]] <- if (gaus){as.name("lFormula")}else{as.name("glFormula")} 
      suppressWarnings( mm <- eval.parent(lmerc), classes = "warning")
      return(mm)
    }else{
      suppressWarnings( mm <- eval.parent(lmerc), classes = "warning")
      cls <- if (gaus){"lmerMod"}else{"glmerMod"} 
      # put it in a lmeb object
      ans <- do.call(new, list(Class=cls, relfac=list(), udu=list(), 
                               frame=mm@frame, flist=mm@flist, cnms=mm@cnms, Gp=mm@Gp,
                               theta=mm@theta, beta=mm@beta,u=mm@u,lower=mm@lower,
                               devcomp=mm@devcomp, pp=mm@pp,resp=mm@resp,optinfo=mm@optinfo))
      ans@call <- evalq(mc)
      return(ans)
    }
  }            # call [g]lmer instead
  
  stopifnot(is.list(relmat),        # check the relmat argument
            length(names(relmat)) == length(relmat),
            all( sapply(relmat, inherits, what = c("relmat","matrix","Matrix"))  ))
  
  # >>>>>>>>> use the match call to parse the data and formula
  lmerc[[1]] <- if (gaus){as.name("lFormula")}else{as.name("glFormula")} # change from model to lFormula
  lmerc$nIters <- NULL # remove nIters
  lmerc$relmat <- NULL # remove relmat from the match call
  lmerc$addmat <- NULL # remove relmat from the match call
  lmerc$dateWarning <- NULL # remove relmat from the match call
  lmerc$rotation <- NULL # remove relmat from the match call
  lmerc$rotationK <- NULL # remove relmat from the match call
  lmerc$coefOutRotation <- NULL # remove relmat from the match call
  lmerc$returnFormula <- NULL # remove relmat from the match call
  lmerc$suppressOpt <- NULL # remove relmat from the match call
  
  response <- all.vars(formula)[1]
  if(rotation){ # if user want rotation we need to impute in advance
    if(!missing(data)){
      lmerc$data <- data
      missed <- which(is.na(data[,response]))
      if(length(missed)>0){
        if(trace){message(magenta("* Response imputed for rotation."))}
      }
      lmerc$data[,response] <- imputev(data[,response])
    }
  }
  suppressWarnings( lmod <- eval.parent(lmerc) , classes = "warning") # necesary objects from lFormula
  # return(lmod)
  
  ## DO ROTATION OF RESPONSE AND RELMATS IF REQUIRED (lmod$fr[,response])
  '%!in%' <- function(x,y)!('%in%'(x,y)) 
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
  if( response %!in% colnames(lmod$fr) ){stop("Response selected in your formula is not part of the dataset provided.", call. = FALSE)}
  # goodRecords <- which(!is.na(lmod$fr[,response]))
  udu <- list()
  # >>>>>>>>> get cholesky factor (and new response if rotation)
  if(length(relmat) > 0){ 
    if(rotation){ # if UDU decomposition + Cholesky is requested
      if(length(relmat) > 1){warning("Rotation is only reported to be accurate with one relationship matrix. Make sure you are using the same relationship matrix for the different random effects for the rotation approach.", call. = FALSE)}
      for(iRel in 1:length(relmat)){
        idsOrdered <- as.character(unique(lmod$fr[,names(relmat)[iRel]])) # when we rotate we need to have relmat already ordered before creating the matrices
        relmat[[iRel]] = relmat[[iRel]][ idsOrdered , idsOrdered ]
      }
      if(trace){message(magenta("* Rotation of response step."))}
      # only the first relmat will be used so if more, the rotation will only work if it is the same relmat in the next random effects
      udu <- umat(formula=as.formula(paste("~", paste(names(relmat), collapse = "+"))), relmat = relmat, 
                  data=lmod$fr, addmat = addmat, k=rotationK)
      # if rotation we impute the response
      lmod$fr[,response] <- imputev(x=lmod$fr[,response],method="median")#, by=data[udu$effect])
      
      Ut <- t(udu$U[[1]]) # extract Ut
      U <- udu$U[[1]] # extract U
      sx0 <- seq(1,length(lmod$fr[,response]), ncol(Ut) ) # start of each piece
      ex0 <- seq(ncol(Ut),length(lmod$fr[,response]), ncol(Ut) ) # end of each piece
      newValues <- Matrix::Matrix(lmod$fr[,response])
      for(iPiece in 1:length(sx0)){ # rotate response piece by piece
        newValues[sx0[iPiece]:ex0[iPiece],] <- Ut %*% newValues[sx0[iPiece]:ex0[iPiece],]
      }
      # newValues <- udu$Utn %*% Matrix::Matrix(lmod$fr[,response])
      newValues <- newValues[,1]
      
      outlier <- grDevices::boxplot.stats(x=newValues,coef=coefOutRotation )$out
      if(length(outlier) > 0){newValues[which(newValues %in% outlier)] = mean(newValues[which(newValues %!in% outlier)])}
      lmod$fr[,response] <- newValues
      if(trace){message(magenta("* Cholesky of relmats step."))}
      for(iD in names(udu$D)){
        relmat[[iD]] <- Matrix::chol(udu$D[[iD]])
      }
      udu$newValues <- newValues
      # lmerc$data <- data
    }else{ # classical approach, just cholesky
      if(trace){message(magenta("* Cholesky of relmats step."))}
      for (i in seq_along(relmat)) {
        relmat[[i]] <- Matrix::chol(relmat[[i]])
      }
    }
  }
  
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
      for(iPiece in 1:length(sx0)){ # rotate X piece by piece
        lmod$X[sx0[iPiece]:ex0[iPiece],] <- Ut %*% lmod$X[sx0[iPiece]:ex0[iPiece],]
      }
      # lmod$X <- udu$Utn %*% lmod$X # previous approach
      if(trace){message(magenta("* Rotation applied to the X matrix."))}
    }
  }
  ##############################
 
  # >>>>>>>>>> apply addmat (additional matrices)
  for (i in seq_along(addmat)) {
    if(trace){message(magenta("* Merging additional matrices."))}
    if(!missing(data)){
      goodRecords <- which(!is.na(data[,response]))
    }else{
      provResponse <- get(response)
      goodRecords <- which(!is.na(provResponse))
    } # this second option may still not work
    tn0 <- which(match(pnms2[i], names(fl)) == asgn)
    for(j in 1:length(tn0)){ # diagonal and unstructured models require to multiple all matrices by the same relfactor
      ind <- (lmod$reTrms$Gp)[tn0[j]:(tn0[j]+1L)]
      rowsi <- (ind[1]+1L):ind[2]
      covariate <- unlist(lmod$reTrms$cnms[tn0])[j]
      if(covariate %in% colnames(lmod$fr)){ # if is a random regression
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
  }
  
  # >>>>>>>>> time to apply the relmat
  namR <- unique(names(lmod$reTrms$cnms))
  if(any(namR %in% names(relmat) )){
    if(trace){message(magenta("* Postmultiplying LZ' step."))}
  }
  for (i in seq_along(namR)) { # for each random effect readjust # Zt i=2
    tn <- which(match(namR[i], names(fl)) == asgn) # match relmat names with random effects names
    # unstructured tn=1, other models tn=n.levels.intercept
    for(j in 1:length(tn)){ # for each intercept matching this relationship matrix (diagonal and unstructured models require to multiple all incidence matrices by the same relfactor)
      ind <- (lmod$reTrms$Gp)[tn[j]:(tn[j]+1L)] # which columns match this random effect
      rowsi <- (ind[1]+1L):ind[2] # first to last column from Z
      # reorder relfac 
      if(namR[i] %in% names(relmat) ){ # this random effect has a relationship matrix so we need to adjust
        colnamesRelFac <- colnames(relfac[[namR[i]]])
        if( mean(table(colnamesRelFac)) > 1 ){  # is this complex because we may have a relationship matrix with repeated names
          toBeRemoved <- character()
          namesProvRelFac <- character() 
          foundV <- numeric()
          for(p in which( rownames(Zt) %in% rownames(relfac[[namR[i]]]) ) ){ # p=1
            found <- which(colnamesRelFac %in% rownames(Zt)[p])
            found <- setdiff(found, toBeRemoved)[1]
            toBeRemoved <- c(toBeRemoved, found[1])
            if(!is.na(found)){
              foundV <- c(foundV,found)
              namesProvRelFac <- c(namesProvRelFac, colnamesRelFac[found] )
            }
          }
          provRelFac <- relfac[[namR[i]]][foundV,foundV] 
          colnames(provRelFac) <- rownames(provRelFac) <- namesProvRelFac
          relfac[[namR[i]]] <- provRelFac
        }else{
          pick <- intersect( rownames(Zt), rownames(relfac[[namR[i]]])  ) # match names in relmat and Z matrix
          if(length(pick)==0){stop(paste("The names on your relmat does not coincide with the names in your factor",pnms[i],". Maybe you didn't code it as factor?"))}
          provRelFac <- relfac[[namR[i]]][pick,pick] # only pick portion of relmat that coincides with Z
        }
      }
      # multiply by the provRelFac or by the Utn matrix
      if( length(lmod$reTrms$cnms[[j]]) == 1 ){ # regular model (intercept || slope) OR (1 | slope )
        
        ZtL <- list() # we have to do this because filling by rows a Column-oriented matrix is extremely slow so it is faster to cut and paste
        
        if(namR[i] %in% names(relmat) ){ # if random effect has a relationship matrix
          # left part
          if(min(rowsi) > 1){ZtL[[1]] <- Zt[1:(min(rowsi)-1),]}
          # central part
          provRelFac <- as(as(as( provRelFac,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
          ZtL[[2]] <- provRelFac %*% Zt[rowsi,] 
          # right part
          if(max(rowsi) < nrow(Zt)){ZtL[[3]] <- Zt[(max(rowsi)+1):nrow(Zt),]}
          Zt <- do.call(rbind, ZtL) # bind all
        }else{ # random effect does not have a relationship matrix
          if(rotation){  # if Lee and Vander Werf 2016 rotation is requested a non-relationship random effect needs to be adjusted by U
            # left part
            if(min(rowsi) > 1){ZtL[[1]] <- Zt[1:(min(rowsi)-1),]}
            # central part
            for(iPiece in 1:length(sx0)){ # rotate Z piece by piece
              Zt[rowsi,sx0[iPiece]:ex0[iPiece]] <- Zt[rowsi,sx0[iPiece]:ex0[iPiece]] %*% U
            }; ZtL[[2]] <- Zt[rowsi,]
            # ZtL[[2]] <- Zt[rowsi,] %*% t(udu$Utn)
            if(trace){message(magenta("* Rotation applied to other Z matrices."))}
            # right part
            if(max(rowsi) < nrow(Zt)){ZtL[[3]] <- Zt[(max(rowsi)+1):nrow(Zt),]}
            Zt <- do.call(rbind, ZtL) # bind all
          }
        }
        
      }else{ # complex model (intercept | slope)
        mm <- Matrix::Diagonal( length(lmod$reTrms$cnms[[j]]) )
        ZtL <- list()
        if(namR[i] %in% names(relmat) ){ # if random effect has a relmat
          # left part
          if(min(rowsi) > 1){ZtL[[1]] <- Zt[1:(min(rowsi)-1),]}
          # central part
          if(length(rowsi) != ncol(provRelFac)*ncol(mm) ){stop(paste("Relationship matrix dimensions of ",pnms[i],"do not conform with the random effect, please review."), call. = FALSE)}
          provRelFac <- as(as(as( provRelFac,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
          ZtL[[2]] <- Matrix::kronecker(provRelFac, mm, make.dimnames = TRUE) %*% Zt[rowsi,]
          rownames(ZtL[[2]]) <- rownames(Zt[rowsi,])
          # right part
          if(max(rowsi) < nrow(Zt)){ZtL[[3]] <- Zt[(max(rowsi)+1):nrow(Zt),]}
          Zt <- do.call(rbind, ZtL) # bind all
        }else{ # if random effect has no relmat
          if(rotation){  # if Lee and Vander Werf 2016 rotation is requested a non-relationship random effect needs to be adjusted by U
            # left part
            if(min(rowsi) > 1){ZtL[[1]] <- Zt[1:(min(rowsi)-1),]}
            # central part
            for(iPiece in 1:length(sx0)){ # rotate Z piece by piece
              Zt[rowsi,sx0[iPiece]:ex0[iPiece]] <- Zt[rowsi,sx0[iPiece]:ex0[iPiece]] %*% U
            }; ZtL[[2]] <- Zt[rowsi,]
            # ZtL[[2]] <- Zt[rowsi,] %*% t(udu$Utn)
            if(trace){message(magenta("* Rotation applied to the Z matrices."))}
            # right part
            if(max(rowsi) < nrow(Zt)){ZtL[[3]] <- Zt[(max(rowsi)+1):nrow(Zt),]}
            Zt <- do.call(rbind, ZtL) # bind all
          }
        } # end of: if random effect has no relmat
        
      }# end of : type of random effect
      
    } # enf of for each intercept
    provRelFac <- NULL
  } # end of for each random effect
  
  reTrms <- list(Zt=Zt,theta=if(is.null(start)){lmod$reTrms$theta}else{start},Lambdat=lmod$reTrms$Lambdat,Lind=lmod$reTrms$Lind,
                 lower=lmod$reTrms$lower,flist=lmod$reTrms$flist,cnms=lmod$reTrms$cnms, Gp=lmod$reTrms$Gp)
  lmod <- list(fr=lmod$fr, X=lmod$X, reTrms=reTrms, formula=formula, verbose=verbose,
               start=if(is.null(start)){lmod$reTrms$theta}else{start},
               control=control)
  
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
  if(returnFormula){ # if user only wants the incidence matrices
    return(lmod)
  }else{
    if(trace){message(magenta("* Optimization step ..."))}
    if (gaus) { # gaussian distribution
      lmod$REML = REML # TRUE# resp$REML > 0L
      suppressWarnings( devfun <- do.call(mkLmerDevfun, lmod ), classes = "warning") # creates a deviance function
      if(suppressOpt){ # user wants to force variance components without fitting a model
        opt <- list(par = start, fval = devfun(start), feval = 1, conv = 0)
      }else{ # # user wants to optimize the varcomp optimizer
        suppressWarnings( opt <- optimizeLmer(devfun, optimizer = lmod$control$optimizer, 
                                              control = lmod$control$optCtrl, 
                                              verbose=lmod$verbose, 
                                              calc.derivs=lmod$control$calc.derivs,
                                              restart_edge=lmod$control$restart_edge,
                                              boundary.tol=lmod$control$boundary.tol,
                                              use.last.params=lmod$control$use.last.params)   , classes = "warning") # need to pass control 
      } 
    } else { # exponential family of distributions
      lmod$family <- family
      suppressWarnings( devfun <- do.call(mkGlmerDevfun,lmod) , classes = "warning") # creates a deviance function
      if(suppressOpt){ # user wants to force variance components without optimizing
        opt <- list(par = start, fval = devfun(start), feval = 1, conv = 0)
      }else{ # user wants to optimize the varcomp optimizer
        suppressWarnings( opt <- optimizeGlmer(devfun, optimizer = lmod$control$optimizer[1], # only first optimizer used 
                                               control = lmod$control$optCtrl, 
                                               verbose=lmod$verbose, 
                                               calc.derivs=lmod$control$calc.derivs,
                                               restart_edge=lmod$control$restart_edge,
                                               boundary.tol=lmod$control$boundary.tol,
                                               use.last.params=lmod$control$use.last.params)  ) # need to pass control 
      } 
    }
    if(trace){message(magenta("* Done!!"))}
    # make results in a mkMerMod object format
    suppressWarnings( mm <- mkMerMod(environment(devfun), opt, lmod$reTrms, lmod$fr, mc), classes = "warning" )
    cls <- if (gaus){c("lmerMod")}else{c("glmerMod")} 
    ans <- do.call(new, list(Class=cls, relfac=relfac, udu=udu, #goodRecords=goodRecords,
                             frame=mm@frame, flist=mm@flist, cnms=mm@cnms, Gp=mm@Gp,
                             theta=mm@theta, beta=mm@beta,u=mm@u,lower=mm@lower,
                             devcomp=mm@devcomp, pp=mm@pp,resp=mm@resp,optinfo=mm@optinfo))
    ans@call <- evalq(mc)
    return(ans)
  }
  
}

setMethod("ranef", signature(object = "lmeb"),
          function(object, condVar = TRUE, drop = FALSE, whichel = names(ans), includeCVM=TRUE, ...)  {
            # print("new")
            relmat <- ifelse(length(object@relfac) > 0, TRUE, FALSE)
            if(relmat){rf <- object@relfac}
            ans <- lme4::ranef(object, condVar=FALSE, drop = FALSE) # extracts condVar 1st time
            ans <- ans[whichel]
            if(condVar){
              mapCondVar <- mkMmeIndex(object) # rbind(namesBlue, namesBlup)
              condVarMat <- condVarRotated(object) # internally we're extracting condVar 2nd time
            }
            for (nm in names(object@flist)) { # for each random effect # nm <- names(rf)[1]
              dm <- data.matrix(ans[[nm]])
              cn <- colnames(dm)
              rn <- rownames(dm)
              
              if (nm %in% names(object@relfac) ) { # transform back when relfac was used for this random effect
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
              }
              if (condVar){ # if conditional variance was requested put it in our desired shape
                mapCondVarNm <- mapCondVar[which(mapCondVar$group == nm),]
                intercepts <- unique(mapCondVarNm$variable)
                postVarNm <- matrix(NA, ncol=length(intercepts), nrow=nrow(mapCondVarNm)/length(intercepts) )
                for(j in 1:length(intercepts)){ # iInter = 1 # intercepts[1]
                  iInter <- intercepts[j]
                  v <- mapCondVarNm[which(mapCondVarNm$variable == iInter), "index"]
                  postVarNm[,j] <- diag(condVarMat)[v]
                }
                rownames(postVarNm) <- rownames(dm)
                colnames(postVarNm) <- colnames(dm)
                attr(ans[[nm]], which="postVar") <- postVarNm
              }
              
            }
            # before returning store the condVarMat and its names
            if(all(c(includeCVM, condVar))){ # if both are TRUE add the PEV
              attr(ans, which="condVarMat") = condVarMat
              attr(ans, which="mapCondVar") = mapCondVar
            }
            return(ans)
          })

setMethod("fitted", signature(object = "lmeb"),
          function(object, ...) {
            W <- do.call(cbind, getME(object = object, c("X","Z")) )
            b <- rbind( as.matrix(fixef(object)),
                        getME(object = object, c("b")) )
            y.hat <- W%*%b
            return(y.hat)
          })


setMethod("residuals", signature(object = "lmeb"),
          function(object, ...) {
            getME(object, "y") - fitted(object)
          })




