##################################################################################################
#Startup function
#this function is executed once the library is loaded
.onAttach = function(library, pkg)
{
  Rv = R.Version()
  if(!exists("getRversion", baseenv()) || (getRversion() < "3.5.0"))
    stop("This package requires R 3.5.0 or later")
  if(interactive()) {
    packageStartupMessage(green(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(green(paste("[] Linear Mixed Equations 4 Breeding (lme4breeding) 1.1.0 (2025-12) []",sep="")),appendLF=TRUE)
    packageStartupMessage(paste0(green("[] Author: Giovanny Covarrubias-Pazaran",paste0(bgGreen(white(" ")), bgWhite(magenta("M")), bgRed(white(" ")),"  ", bgRed(bold(yellow(" (") )),bgRed(bold(white("W"))), bgRed(bold(yellow(") "))) ) ,"                 []")),appendLF=TRUE)
    packageStartupMessage(green("[] Special thanks to the lme4 dev team (Bolker, Bates, et al.)      []"),appendLF=TRUE)
    packageStartupMessage(green("[] Type 'vignette('lmebreed.gxe')' for a short tutorial             []"),appendLF=TRUE)
    packageStartupMessage(green("[] Type 'citation('lme4breeding')' to know how to cite it           []"),appendLF=TRUE)
    packageStartupMessage(green(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(green("lme4breeding is updated on CRAN every 3-months due to CRAN policies"),appendLF=TRUE)
    packageStartupMessage(green("Source code is available at https://github.com/covaruber/lme4breeding"),appendLF=TRUE)
  }
  invisible()
}

###
umat <- function(formula, relmat, data, addmat, k=NULL){
  
  if(missing(data)){stop("Please provide the dataset where we can extract find the factor in formula.")}
  if(missing(relmat)){stop("Please provide the relationship matrix where we will apply the eigen decomposition.")}
  if(missing(formula)){stop("Please provide the formula with the factor to do the decomposition on.")}
  if(!inherits(formula, "formula")){stop("Please provide the formula as.formula().")}
  if(!is.list(relmat)){stop("relmat argument should be a named list of relationship matrices", call. = FALSE)}
  
  idProvided <- all.vars(formula)
  if(length(idProvided) > 1){message(magenta("Only one dense relationship matrix can be eigen decomposed. You have provided more than one relmat. Only the first one will be used for the rotation."))}
  
  # build the U nxn matrix
  Ul <- Dl <- Zu <- nLev <- list()
  nLev <- numeric()
  for(iProv in idProvided){ # only one is actually accepted for now according to theory
    nTraits <- max(table(data[,idProvided]))
    levsInA = unique(data[,idProvided])
    if(iProv %in% colnames(data)){
      Z <- sparse.model.matrix(as.formula(paste("~",iProv,"-1")), data=data)
      colnames(Z) <- gsub(iProv,"", colnames(Z))
    }else{
      if(iProv %in% names(addmat)){
        if(is.list(addmat[[iProv]])){ # indirect genetic effects
          Z <- Reduce("+", addmat[[iProv]])
        }else{ # single model
          Z <- addmat[[iProv]]
        }
      }else{
        stop(paste("Your term", iProv, "is neither in the dataset nor in addmat, please correct."), call. = FALSE)
      }
    }
    if(is.null(k)){k<- nrow(relmat[[iProv]])}
    suppressWarnings( UD <- RSpectra::eigs_sym(relmat[[iProv]], k=k, which = "LM"), classes = "warning")
    difference <- ncol(relmat[[iProv]]) - k
    if(difference > 0){
      UD$values <- c(UD$values, rep(1e-6, difference))
      UD$vectors <- cbind(UD$vectors, matrix(0, nrow=nrow(UD$vectors), ncol=difference))
    }
    U<-UD$vectors
    D<-Matrix::Diagonal(x=UD$values)# This will be our new 'relationship-matrix'
    rownames(D) <- colnames(D) <- rownames(relmat[[iProv]])
    rownames(U) <- colnames(U) <- rownames(relmat[[iProv]])
    common <- intersect(colnames(U), colnames(Z))
    Ul[[iProv]]<- U[common,common]
    Dl[[iProv]]<-D[common,common]# This will be our new 'relationship-matrix'
    # Utn <- Matrix::kronecker(Matrix::Diagonal(n=nTraits),t(U[common,common]))
  }
  # UnList <- list()
  # for(iel in 1:1){ # we only use the first relationship matrix for the rotation # length(Ul)
  #   UnList[[iel]] <- Matrix::kronecker(Matrix::Diagonal(n=nLev[[iel]]),t(Ul[[iel]]))
  # }
  # Utn <- Reduce("+",UnList)
  
  if( (nrow(U)*nTraits) != nrow(data)){
    stop("The eigen decomposition only works when the rotated effect is balanced (equal number of reps for each level). \n Please ensure you fill the dataset to make it balanced for the \n 'relmat' terms or set 'rotation' to FALSE.", call. = FALSE)
  }
  
  return(list(D=Dl, U=Ul, # RRt=ZrZrt, 
              effect=idProvided, 
              record=idProvided # data$recordF
  ))
}

umat2 <- function(formulaU, relmat, data, addmat, k=NULL, trace=TRUE){
  
  if(missing(data)){stop("Please provide the dataset where we can extract find the factor in formulaU.")}
  if(missing(relmat)){stop("Please provide the relationship matrix where we will apply the eigen decomposition.")}
  if(missing(formulaU)){stop("Please provide the formulaU with the factor to do the decomposition on.")}
  if(!inherits(formulaU, "formula")){stop("Please provide the formulaU as.formulaU().")}
  if(!is.list(relmat)){stop("relmat argument should be a named list of relationship matrices", call. = FALSE)}
  
  idProvided <- all.vars(formulaU)
  if(length(idProvided) > 1){message(magenta("*** More than one relationship matrix present. VCs will be biased."))}
  
  # build the U nxn matrix
  Ul <- Dl <-  Utnl <- list()
  counter = 1
  for(iProv in idProvided){ # iProv = idProvided[1]
    if(trace){message(magenta(paste("** Rotating", iProv,"effect")))}
    # get number of replicates
    if(iProv %in% colnames(data)){
      Z <- sparse.model.matrix(as.formula(paste("~",iProv,"-1")), data=data)
      colnames(Z) <- gsub(iProv,"", colnames(Z))
    }else{
      if(iProv %in% names(addmat)){
        if(is.list(addmat[[iProv]])){ # indirect genetic effects
          Z <- Reduce("+", addmat[[iProv]])
        }else{ # single model
          Z <- addmat[[iProv]]
        }
      }else{
        stop(paste("Your term", iProv, "is neither in the dataset nor in addmat, please correct."), call. = FALSE)
      }
    }
    n <- rep(1,nrow(Z))
    Z2 <- apply(Z,2,cumsum)
    tabIprov <- table(data[,idProvided])
    if(var(tabIprov) > 0){message(magenta(paste("*** Levels in relmat for",iProv,"are unbalanced so VCs will be biased.")))}
    for(j in 1: max(tabIprov)){ # j=2
      where <- apply(Z2,2,function(x){
        v <- which(x==j)
        if(length(v)>0){return(min(v))}else{return(NA)}
      })
      start <- min(where, na.rm = TRUE)
      end <- max(where, na.rm = TRUE)
      n[start:end] <-  n[start:end]*j
    }
    # calculate U and D
    if(is.null(k)){k<- nrow(relmat[[iProv]])}
    suppressWarnings( UD <- RSpectra::eigs_sym(relmat[[iProv]], k=k, which = "LM"), classes = "warning")
    difference <- ncol(relmat[[iProv]]) - k
    if(difference > 0){
      UD$values <- c(UD$values, rep(1e-6, difference))
      UD$vectors <- cbind(UD$vectors, matrix(0, nrow=nrow(UD$vectors), ncol=difference))
    }
    U<-UD$vectors
    D<-Matrix::Diagonal(x=UD$values)# This will be our new 'relationship-matrix'
    rownames(D) <- colnames(D) <- rownames(relmat[[iProv]])
    rownames(U) <- colnames(U) <- rownames(relmat[[iProv]])
    common <- intersect(colnames(U), colnames(Z))
    Ul[[iProv]]<- U[common,common]
    Dl[[iProv]]<-D[common,common]# This will be our new 'relationship-matrix'
    
    Ut <- as(as(as(  t(U),  "dMatrix"), "generalMatrix"), "CsparseMatrix") 
    Utl <- list()
    for(j in 1: max(tabIprov)){ # j=1
      v <- which(n == j)
      w <- as.character(data[v,idProvided])
      Utl[[j]] <- Ut[w,w]#/(j^2)
    }
    if(counter == 1){
      Utn  <- do.call(Matrix::bdiag,Utl)
    }else{
      Utn  <- Utn + do.call(Matrix::bdiag,Utl)
    }
    counter <- counter + 1
    # we still need addmat version
  }
  return(list(Utn=Utn, D=Dl, U=Ul, # RRt=ZrZrt, 
              effect=idProvided, 
              record=idProvided # data$recordF
  ))
}

balanceData <- function(data, slope=NULL, intercept=NULL){
  if(is.null(slope)){stop("Please provide the column name corresponsing to the slope (e.g., treatments).", call. = FALSE)}
  if(is.null(intercept)){stop("Please provide the column name corresponsing to the slope (e.g., treatments).", call. = FALSE)}
  slopeLevs = unique(data[,slope])
  data[,paste0(intercept, collapse = "_")] <- apply(data[,intercept],1,function(x){paste0(x, collapse = "_")})
  interLevs = unique(data[,paste0(intercept, collapse = "_")])
  balanced = expand.grid(slopeLevs, interLevs)
  colnames(balanced) <- c(slope,paste0(intercept, collapse = "_"))
  balanced = merge(balanced, data, by=c(paste0(intercept, collapse = "_"), slope), 
                   all.x = TRUE)
  balanced = balanced[ order(balanced[,paste0(intercept, collapse = "_")], balanced[,slope]), ]
  return(balanced)
}

###
adjBeta <- function(x){
  if(length(x) > 1){
    x[2:length(x)] <- x[2:length(x)] + x[1]
  }
  return(x)
}
###


####
getMME <- function(object, vc=NULL, recordsToUse=NULL){
  if(is.null(vc)){
    vc <- VarCorr(object) #object %>% VarCorr %>% as_tibble # extract estimated variance components (vc)
  }
  # R = varcov-matrix for error term
  n    <- length(object@resp$y) # object %>% summary %>% pluck(residuals) %>% length # numer of observations
  vc_e <- attr(VarCorr(object), "sc")^2
  # vc_e <- vc %>% filter(grp=="Residual") %>% pull(vcov)     # error vc
  Ri    <- Matrix::Diagonal(n)*(1/vc_e)                                      # R matrix = I_n * vc_e
  
  # Design Matrices
  X <- getME(object, "X") #%>% as.matrix # Design matrix fixed effects
  Z <- getME(object, "Z") #%>% as.matrix # Design matrix random effects
  y <- getME(object, "y") #%>% as.matrix # Design matrix random effects
  
  if(!is.null(recordsToUse)){
    X <- X[recordsToUse,, drop=FALSE]
    Z <- Z[recordsToUse,, drop=FALSE]
    y <- y[recordsToUse]
    Ri <- Ri[recordsToUse,recordsToUse]
  }
  # Mixed Model Equation (HENDERSON 1986; SEARLE et al. 2006)
  C11 <- t(X) %*% Ri %*% X
  C12 <- t(X) %*% Ri %*% Z
  C21 <- t(Z) %*% Ri %*% X
  C22 <- t(Z) %*% Ri %*% Z #+ solve(G) 
  
  C <- rbind(cbind(C11, C12),  
             cbind(C21, C22)) #%>% as.matrix # Combine components into one matrix C
  
  RHS1 <- t(X) %*% Ri %*% y
  RHS2 <- t(Z) %*% Ri %*% y
  RHS <- rbind(RHS1, RHS2)
  # G = varcov-matrx for all random effects
  # subset of G regarding genotypic effects
  # is a diagonal matrix because Z has the relationship incorporated
  fl <- object@flist
  asgn <- attr(fl, "assign")
  pnms <- names(object@flist)# names(object@relfac)
  Gi <- Matrix::Matrix(0, nrow = nrow(C), ncol=ncol(C))
  # Zt <- getME(object, "Zt")
  # vc <- VarCorr(object);
  for(iFac in pnms){ # iFac = pnms[1]
    tn <- which(match(iFac, names(fl)) == asgn)
    
    ####
    vcov <- do.call(Matrix::bdiag, vc[tn]) # var covar matrix
    vcov <- vcov + diag(1e-6,ncol(vcov),ncol(vcov))
    LLt <- Matrix::Diagonal( length(unique(object@flist[[iFac]])) ) # digonal for unique number of levels
    rowsi <- list()
    for(j in 1:length(tn)){ # j=1
      ind <- (object@Gp)[tn[j]:(tn[j]+1L)]
      rowsi[[j]] <- ((ind[1]+1L):ind[2])+1
    }
    Gi[unlist(rowsi),unlist(rowsi)] <- kronecker( LLt , solve( Matrix::nearPD( vcov )$mat ) )
    ##
    # for(j in 1:length(tn)){ # j=1
    #   ind <- (object@Gp)[tn[j]:(tn[j]+1L)]
    #   rowsi <- ((ind[1]+1L):ind[2])+1
    #   LLt <- Matrix::Diagonal( length(unique(object@flist[[iFac]])) )
    #   if(length(diag(vc[[iFac]])) > 0){
    #     Gi[rowsi,rowsi] <- kronecker( LLt , solve( Matrix::nearPD( vc[[iFac]] )$mat ) )
    #   }else{
    #     Gi[rowsi,rowsi] <- kronecker( LLt ,  vc[[iFac]]  )
    #   }
    # }
  }
  # incomplete block part of G matrix = I * vc.b
  
  C <- C + Gi + diag(1e-4, ncol(C), ncol(C))
  
  # Mixed Model Equation Solutions 
  C_inv <- solve(C)# %>% solve   
  rownames(C_inv) <- colnames(C_inv) <- c(colnames(X), colnames(Z))
  bu <- C_inv%*% RHS
  rownames(bu) <- rownames(C_inv)
  
  # when relfac is present save them in block diagonal and multiple bu and C by it
  if(length(object@relfac) > 0){
    ROT <- Matrix::Diagonal(n=nrow(C))#Matrix(0, nrow = nrow(C), ncol=ncol(C))
    for(iFac in pnms){ # iFac = pnms[1]
      tn <- which(match(iFac, names(fl)) == asgn)
      for(j in 1:length(tn)){ # j=1
        ind <- (object@Gp)[tn[j]:(tn[j]+1L)]
        rowsi <- ( (ind[1]+1L):ind[2] ) + ncol(X)
        if( iFac %in% names(object@relfac) ){
          pick <- rownames(C)[rowsi] # intersect( colnames(C), rownames(  object@relfac[[iFac]] )  )
          ROT[rowsi,rowsi] <- object@relfac[[iFac]][pick,pick]
        }
      }
    }
    # rotate blups and Ci matrix
    rn <- rownames(C_inv)
    buROT <- t(as.matrix( t( bu ) %*% ROT ))
    C_invROT <- t(ROT) %*% C_inv %*% (ROT)
    rownames(buROT) <- rn
    colnames(C_invROT) <- rownames(C_invROT) <- rn
    return(list(Ci=C_invROT, bu=buROT, RHS=RHS))
  }else{
    return(list(Ci=C_inv, bu=bu, RHS=RHS))
  }
  
  
}

smm <- function(x){
  if(is.matrix(x)){
    dummy <- x
    mm <- diag(1,ncol(x))
  }else{
    if(!is.character(x) & !is.factor(x)){
      namess <- as.character(substitute(list(x)))[-1L]
      dummy <- as(as(as( Matrix(x,ncol=1),  "dMatrix"), "generalMatrix"), "CsparseMatrix") 
      colnames(dummy) <- namess
      mm <- diag(ncol(dummy));
    }else{
      dummy <- x
      levs <- na.omit(unique(dummy))
      if(length(levs) > 1){
        dummy  <- Matrix::sparse.model.matrix(~dummy-1,na.action = na.pass)
        colnames(dummy) <- gsub("dummy","",colnames(dummy))
      }else{
        vv <- which(!is.na(dummy));
        dummy <- matrix(0,nrow=length(dummy))
        dummy[vv,] <- 1; colnames(dummy) <- levs
      }
      mm <- diag(1,ncol(dummy))
    }
  }
  return(dummy)
}

condVarFun <- function (object, scaled = TRUE) {
  Lamt <- getME(object, "Lambdat")
  L <- getME(object, "L")
  LL <- solve(L, Lamt, system = "A")
  cc <- crossprod(Lamt, LL)
  if (scaled) 
    cc <- sigma(object)^2 * cc
  cc
}

condVarRotated <- function(object){
  vcovBlue <- vcov(object )
  condVarBlup <- condVarFun(object)
  condVarFull <- Matrix::bdiag( vcovBlue, condVarBlup )
  mapCondVar <- mkMmeIndex(object)# rbind(namesBlue, namesBlup)
  condVarFull@Dimnames[[1]] <- mapCondVar$level
  condVarFull@Dimnames[[2]] <- mapCondVar$variable
  # rotate back cholesky (triangular dense matrix)
  relmat <- ifelse(length(object@relfac) > 0, TRUE, FALSE) # control to know if we should rotate
  if(relmat){
    message(magenta(paste("Rotating back conditional variance using Cholesky factors")))
    groups <- unique(mapCondVar[,c("variable","group")])
    Ll <- vl <- list()
    for(iGroup in 1:nrow(groups)){ # iGroup = 2
      v <- which( ( mapCondVar[,"variable"] == groups[iGroup,"variable"] ) & ( mapCondVar[,"group"] == groups[iGroup,"group"] ) )
      vl[[iGroup]] <- v
      if(groups[iGroup,"group"] %in% names( object@relfac )){ # extract relfac
        # L[v,v] <- as(as(as( object@relfac[[ groups[iGroup,"group"] ]][ mapCondVar[v,"level"], mapCondVar[v,"level"] ] ,  "dMatrix"), "generalMatrix"), "CsparseMatrix") 
        Ll[[iGroup]] <- as(as(as( object@relfac[[ groups[iGroup,"group"] ]][ mapCondVar[v,"level"], mapCondVar[v,"level"] ] ,  "dMatrix"), "generalMatrix"), "CsparseMatrix") 
      }else{ # build a digonal
        # L[v,v] <- Matrix::Diagonal(n=length(v))
        Ll[[iGroup]] <- Matrix::Diagonal(n=length(v))
      }
    }
    L <- do.call(bdiag, Ll)
    L <- L[unlist(vl), unlist(vl)]
    condVarFull <- crossprod(L, condVarFull %*% L)
    L <- NULL
    condVarFull@Dimnames[[1]] <- mapCondVar$level
    condVarFull@Dimnames[[2]] <- mapCondVar$variable
  }
  # rotate back eigen (dense eigen vectors)
  eigmat <- ifelse(length(object@udu) > 0, TRUE, FALSE) # control to know if we should rotate
  if(eigmat){
    message(magenta(paste("Rotating back conditional variance using Eigen factors")))
    groups <- unique(mapCondVar[,c("variable","group")])
    Ul <- vl <- list()
    for(iGroup in 1:nrow(groups)){ # iGroup = 2
      v <- which( ( mapCondVar[,"variable"] == groups[iGroup,"variable"] ) & ( mapCondVar[,"group"] == groups[iGroup,"group"] ) )
      vl[[iGroup]] <- v
      if(groups[iGroup,"group"] %in% names( object@udu$U ) ){ # extract relfac
        # U[v,v] <- as(as(as( object@udu$U[[ groups[iGroup,"group"] ]][ mapCondVar[v,"level"], mapCondVar[v,"level"] ] ,  "dMatrix"), "generalMatrix"), "CsparseMatrix") 
        Ul[[iGroup]] <- as(as(as( object@udu$U[[ groups[iGroup,"group"] ]][ mapCondVar[v,"level"], mapCondVar[v,"level"] ] ,  "dMatrix"), "generalMatrix"), "CsparseMatrix") 
      }else{ # build a digonal
        # U[v,v] <- Matrix::Diagonal(n=length(v))
        Ul[[iGroup]] <- Matrix::Diagonal(n=length(v))
      }
    }
    U <- do.call(bdiag, Ul)
    U <- U[unlist(vl), unlist(vl)]
    condVarFull <- tcrossprod(U%*%condVarFull,U)
    U <- NULL
    condVarFull@Dimnames[[1]] <- mapCondVar$level
    condVarFull@Dimnames[[2]] <- mapCondVar$variable
  }
  return(condVarFull)
}

mkMmeIndex <- function(object) {
  # get information about the
  # dimensions of the random
  # effects
  rp <- rePos$new(object)
  # list with the levels for
  # each grouping factor (one
  # character-valued list
  # element per factor)
  levs <- lapply(rp$flist, levels)
  # names of the variables
  # associated with each random
  # effect
  variableNames <- groupl <- levsl <- list(); counter=1
  for(i in 1:length(rp$cnms)){
    variableNames[[counter]] <- rep(rp$cnms[[i]], rp$nlevs[ attr(rp$flist, "assign")[i] ] )
    groupl[[counter]] <- rep( rep( names(rp$cnms)[i],  length(rp$cnms[[i]]) ), rp$nlevs[ attr(rp$flist, "assign")[i] ] )
    levsl[[counter]] <- sort( rep(levs[[attr(rp$flist, "assign")[i]] ], length(rp$cnms[[i]])) )
    counter=counter+1
  }
  variableNames <- unlist(variableNames, use.names = FALSE)
  groupl <- unlist(groupl, use.names = FALSE)
  level <- unlist(levsl, use.names = FALSE) # unlist(levs[attr(rp$flist, "assign")], use.names = FALSE)
  namesBlup <- data.frame(index=1:length(level), level, variable=variableNames, group=groupl, type="random")
  
  vcovBlue <- vcov(object )
  namesBlue <- data.frame(index=1:ncol(vcovBlue), level=colnames(vcovBlue),
                          variable="(Intercept)", group= colnames(vcovBlue), type="fixed" )
  fixedTerms <- terms(object)
  fixedTerms <- attr(fixedTerms,"term.labels")
  for(iF in fixedTerms){ # iF = fixedTerms[1]
    Xi <- sparse.model.matrix(as.formula(paste("~",iF,"-1")), data=object@frame)
    namesBlue$group[which(namesBlue$level %in% colnames(Xi))] <- iF
  }
  
  namesBlup$index <- namesBlup$index + nrow(namesBlue)
  res <- rbind(namesBlue, namesBlup)
  
  # variableNames <-
  #   mapply(rep, rp$cnms,
  #          times = rp$nlevs[attr(rp$flist, "assign")],
  #          SIMPLIFY = FALSE) %>%
  #   unlist(use.names = FALSE)
  # construct the output data
  # frame
  # xx <- rep %>%
  #   mapply(levs[attr(rp$flist, "assign")], # levels associated with each RE term
  #          each = rp$ncols,                # num vars associated with each RE term
  #          SIMPLIFY = FALSE) %>%
  #   melt() %>%
  #   setNames(c("level", "group")) %>%
  #   mutate(variable = variableNames) %>%
  #   mutate(index = row_number()) %>%
  #   select(index, level, variable, group)
  return(res)
}

Dtable <- function(object){
  
  mapCondVar <- mkMmeIndex(object) # rbind(namesBlue, namesBlup)
  
  tab <- unique(mapCondVar[,c("variable","group","type")])
  tab[,"include"]=0
  tab[,"average"]=0
  return(tab)
  
}
