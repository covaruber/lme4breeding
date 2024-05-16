##################################################################################################
#Startup function
#this function is executed once the library is loaded
.onAttach = function(library, pkg)
{
  Rv = R.Version()
  if(!exists("getRversion", baseenv()) || (getRversion() < "3.5.0"))
    stop("This package requires R 3.5.0 or later")
  if(interactive()) {
    packageStartupMessage(blue(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(blue(paste("[] Linear Mixed Equations 4 Breeding (lme4breeding) 1.0.1 (2024-06) []",sep="")),appendLF=TRUE)
    packageStartupMessage(paste0(blue("[] Author: Giovanny Covarrubias-Pazaran",paste0(bgGreen(white(" ")), bgWhite(magenta("*")), bgRed(white(" "))),"                        []")),appendLF=TRUE)
    packageStartupMessage(blue("[] Dedicated to the University of Chapingo and UW-Madison           []"),appendLF=TRUE)
    packageStartupMessage(blue("[] Type 'vignette('lmebreed.gxe')' for a short tutorial             []"),appendLF=TRUE)
    packageStartupMessage(blue("[] Type 'citation('lme4breeding')' to know how to cite it           []"),appendLF=TRUE)
    packageStartupMessage(blue(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(blue("lme4breeding is updated on CRAN every 4-months due to CRAN policies"),appendLF=TRUE)
    packageStartupMessage(blue("Source code is available at https://github.com/covaruber/lme4breeding"),appendLF=TRUE)
  }
  invisible()
}

###
umat <- function(formula, relmat, data){
  
  if(missing(data)){stop("Please provide the dataset where we can extract find the factor in formula.")}
  if(missing(relmat)){stop("Please provide the relationship matrix where we will apply the eigen decomposition.")}
  if(missing(formula)){stop("Please provide the formula with the factor to do the decomposition on.")}
  if(!inherits(formula, "formula")){stop("Please provide the formula as.formula().")}
  
  idProvided <- all.vars(formula)
  if(length(idProvided) > 1){stop("Only one factor can be provided in the formula.")}
  data$record <- NA
  ids <- as.character(na.omit(unique(data[,idProvided])))
  for(iId in ids){ # iId<- ids[1]
    found <- which(data[,idProvided] == iId)
    data[found,"record"] <- 1:length(found)
  }
  data$recordF <- as.factor(data$record)
  nLev <- length(levels(data$recordF))
  if(nLev > 1){
    Zr <- sparse.model.matrix(~recordF-1, data=data)
  }else{
    Zr <- Matrix::Matrix(1, ncol=1, nrow=nrow(data))
  }
  
  # dim(Zr)
  part0 <- Zr %*% t(Zr)
  # Matrix::image(part0)
  Z <- sparse.model.matrix(as.formula(paste("~",idProvided,"-1")), data=data)
  colnames(Z) <- gsub(idProvided,"", colnames(Z))
  UD <- eigen(relmat)
  U<-UD$vectors
  D<-Matrix::Diagonal(x=UD$values)# This will be our new 'relationship-matrix'
  rownames(D) <- colnames(D) <- rownames(relmat)
  rownames(U) <- colnames(U) <- rownames(relmat)
  
  # part1 <- Z%*%U[colnames(Z),colnames(Z)]%*%t(Z)
  part1 <- U[as.character(data[,idProvided]), as.character(data[, idProvided])]
  W0 <- part0 * part1
  return(list(Utn=t(W0), D=D, U=U, RRt=part0))
}
###
adjBeta <- function(x){
  if(length(x) > 1){
    x[2:length(x)] <- x[2:length(x)] + x[1]
  }
  return(x)
}
###
leg <- function(x,n=1,u=-1,v=1, intercept=TRUE, intercept1=FALSE){
  
  init0 <- as.character(substitute(list(x)))[-1L]
  
  if(system.file(package = "orthopolynom") == ""){
    stop("Please install the orthopolynom package to use this function.",call. = FALSE)
  }
  requireNamespace("orthopolynom",quietly=TRUE)
  (leg4coef <- orthopolynom::legendre.polynomials(n=n, normalized=TRUE))
  leg4 <- as.matrix(as.data.frame(orthopolynom::polynomial.values(polynomials=leg4coef,
                                                                  x=orthopolynom::scaleX(x, u=u, v=v))))
  colnames(leg4) <- paste("leg",0:(ncol(leg4)-1),sep="")
  if(!intercept){
    leg4 <- leg4[, 2:ncol(leg4), drop = FALSE]
  }
  if(intercept1){
    leg4 <- leg4*sqrt(2)
    # leg4[,1] <- leg4[,1]*sqrt(2)
  }
  attr(leg4,"variables") <- c(init0)
  return(leg4)
}

####
getMME <- function(object, vc=NULL, recordsToKeep=NULL){
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
  
  if(!is.null(recordsToKeep)){
    X <- X[recordsToKeep,, drop=FALSE]
    Z <- Z[recordsToKeep,, drop=FALSE]
    y <- y[recordsToKeep]
    Ri <- Ri[recordsToKeep,recordsToKeep]
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
  vc <- VarCorr(object);
  for(iFac in pnms){ # iFac = pnms[1]
    tn <- which(match(iFac, names(fl)) == asgn)
    for(j in 1:length(tn)){
      ind <- (object@Gp)[tn[j]:(tn[j]+1L)]
      rowsi <- (ind[1]+1L):ind[2]
      LLt <- Matrix::Diagonal(length(rowsi))
      if(sum(diag(vc[[iFac]])) > 0){
        Gi[rowsi,rowsi] <- kronecker( LLt , solve( vc[[iFac]] ) )
      }else{
        Gi[rowsi,rowsi] <- kronecker( LLt ,  vc[[iFac]]  )
      }
    }
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
    ROT <- Matrix::Diagonal(nrow(C))#Matrix(0, nrow = nrow(C), ncol=ncol(C))
    for(iFac in pnms){ # iFac = pnms[1]
      tn <- which(match(iFac, names(fl)) == asgn)
      for(j in 1:length(tn)){
        ind <- (object@Gp)[tn[j]:(tn[j]+1L)]
        rowsi <- (ind[1]+1L):ind[2] + ncol(X)
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
    return(list(Ci=C_invROT, bu=buROT))
  }else{
    return(list(Ci=C_inv, bu=bu))
  }
  
  
}

stackTrait <- function(data, traits){
  
  '%!in%' <- function(x,y)!('%in%'(x,y)) 
  dataScaled <- data
  # remove traits not present in the dataset
  traits <- intersect(traits, colnames(data) )
  # check that traits are numeric
  columnTypesTraits <- unlist(lapply(data[traits],class)) 
  badones <- which(columnTypesTraits %!in% c("integer","numeric") )
  if(length(badones) > 0){stop("some of your selected traits are not numeric. Please correct the traits provided.", call. = FALSE)}
  # identify possible idvars
  idvars <- setdiff(colnames(data), traits)
  for(iTrait in traits){
    dataScaled[,iTrait] <- scale(dataScaled[,iTrait])
  }
  columnTypes <- unlist(lapply(data[idvars],class)) 
  columnTypes <- columnTypes[which(columnTypes %in% c("factor","character","integer"))]
  idvars <- intersect(idvars,names(columnTypes))
  data2 <- reshape(data[,c(idvars, traits)], idvar = idvars, varying = traits,
                   timevar = "trait",
                   times = traits,v.names = "value", direction = "long")
  data2Scaled <- reshape(dataScaled[,c(idvars, traits)], idvar = idvars, varying = traits,
                         timevar = "trait",
                         times = traits,v.names = "value", direction = "long")
  data2 <- as.data.frame(data2)
  data2$valueS <- as.vector(unlist(data2Scaled$value))
  rownames(data2) <- NULL
  varG <- cov(data[,traits], use="pairwise.complete.obs")
  # varG <- apply(data[,traits],2,var, na.rm=TRUE) 
  mu <- apply(data[,traits],2,mean, na.rm=TRUE) 
  return(list(long=data2, varG=varG, mu=mu))
}

add.diallel.vars <- function(df, par1="Par1", par2="Par2",sep.cross="-"){
  # Dummy variables for selfs, crosses, combinations
  df[,"is.cross"] <- ifelse(df[,par1] == df[,par2], 0, 1)
  df[,"is.self"] <- ifelse(df[,par1] == df[,par2], 1, 0)
  df[,"cross.type"] <- ifelse(as.character(df[,par1]) < as.character(df[,par2]), -1,
                              ifelse(as.character(df[,par1]) == as.character(df[,par2]), 0, 1))
  # Dummy variable for the combinations, ignoring the reciprocals
  df[,"cross.id"]<-factor(ifelse(as.character(df[,par1]) <= as.character(df[,par2]),
                                 paste(df[,par1], df[,par2], sep =sep.cross),
                                 paste(df[,par2], df[,par1], sep =sep.cross)) )
  return(df)
}

overlay<- function (..., rlist = NULL, prefix = NULL, sparse=FALSE){
  init <- list(...) # init <- list(DT$femalef,DT$malef)
  ## keep track of factor variables
  myTypes <- unlist(lapply(init,class))
  init0 <- init
  ##
  init <- lapply(init, as.character)
  namesInit <- as.character(substitute(list(...)))[-1L] # names <- c("femalef","malef")
  dat <- as.data.frame(do.call(cbind, init))
  dat <- as.data.frame(dat)
  ## bring back the levels
  for(j in 1:length(myTypes)){
    if(myTypes[j]=="factor"){
      levels(dat[,j]) <- c(levels(dat[,j]),setdiff(levels(init0[[j]]),levels(dat[,j]) ))
    }
  }
  ##
  if (is.null(dim(dat))) {
    stop("Please provide a data frame to the overlay function, not a vector.\\n",
         call. = FALSE)
  }
  if (is.null(rlist)) {
    rlist <- as.list(rep(1, dim(dat)[2]))
  }
  ss1 <- colnames(dat)
  dat2 <- as.data.frame(dat[, ss1])
  head(dat2)
  colnames(dat2) <- ss1
  femlist <- list()
  S1list <- list()
  for (i in 1:length(ss1)) {
    femlist[[i]] <- ss1[i]
    dat2[, femlist[[i]]] <- as.factor(dat2[, femlist[[i]]])
    if(sparse){
      S1 <- Matrix::sparse.model.matrix(as.formula(paste("~", femlist[[i]],
                                                         "-1")), dat2)
    }else{
      S1 <- model.matrix(as.formula(paste("~", femlist[[i]],
                                          "-1")), dat2)
    }
    colnames(S1) <- gsub(femlist[[i]], "", colnames(S1))
    S1list[[i]] <- S1
  }
  levo <- sort(unique(unlist(lapply(S1list, function(x) {
    colnames(x)
  }))))
  if(sparse){
    S3 <- Matrix(0, nrow = dim(dat2)[1], ncol = length(levo))
  }else{
    S3 <- matrix(0, nrow = dim(dat2)[1], ncol = length(levo))
  }
  
  rownames(S3) <- rownames(dat2)
  colnames(S3) <- levo
  for (i in 1:length(S1list)) {
    if (i == 1) {
      S3[rownames(S1list[[i]]), colnames(S1list[[i]])] <- S1list[[i]] *
        rlist[[i]]
    }
    else {
      S3[rownames(S1list[[i]]), colnames(S1list[[i]])] <- S3[rownames(S1list[[i]]),
                                                             colnames(S1list[[i]])] + (S1list[[i]][rownames(S1list[[i]]),
                                                                                                   colnames(S1list[[i]])] * rlist[[i]])
    }
  }
  if (!is.null(prefix)) {
    colnames(S3) <- paste(prefix, colnames(S3), sep = "")
  }
  attr(S3,"variables") <- namesInit
  return(S3)
}


## VS structures for lmebreed
redmm <- function (x, M = NULL, Lam=NULL, nPC=50, cholD=FALSE, returnLam=FALSE) {
  
  if(system.file(package = "RSpectra") == ""){
    stop("Please install the RSpectra package to use the redmm() function.",call. = FALSE)
  }else{
    requireNamespace("RSpectra",quietly=TRUE)
  }
  
  if(is.null(M)){
    # stop("M cannot be NULL. We need a matrix of features that defines the levels of x")
    smd <- RSpectra::svds(x, k=nPC, which = "LM")
    if(is.null(Lam)){
      Lam0 <- smd$u
      Lam = Lam0[,1:min(c(nPC,ncol(x))), drop=FALSE]
      rownames(Lam) <- rownames(x)
      colnames(Lam) <- paste0("nPC",1:nPC)
    }else{
      Lam0=Lam
      Lam = Lam0[,1:min(c(nPC,ncol(M))), drop=FALSE]
      rownames(Lam) <- rownames(M)
      colnames(Lam) <- paste0("nPC",1:nPC)
    }
    Zstar <- Lam
  }else{
    
    if (inherits(x, "dgCMatrix") | inherits(x, "matrix")) {
      notPresentInM <- setdiff(colnames(Z),rownames(M))
      notPresentInZ <- setdiff(rownames(M),colnames(x))
    }else{
      notPresentInM <- setdiff(unique(x),rownames(M))
      notPresentInZ <- setdiff(rownames(M),unique(x))
    }
    if(is.null(Lam)){ # user didn't provide a Lambda matrix
      if(nPC == 0){ # user wants to use the full marker matrix
        Lam <- Lam0 <- M
      }else{ # user wants to use the PCA method
        nPC <- min(c(nPC, ncol(M)))
        if(cholD){
          smd <- try(chol(M) , silent = TRUE)
          if(inherits(smd, "try-error")){smd <- try(chol((M+diag(1e-5,nrow(M),nrow(M))) ) , silent = TRUE)}
          Lam0 = t(smd)
        }else{
          smd <- RSpectra::svds(M, k=nPC, which = "LM")
          Lam0 <- smd$u
        }
        Lam = Lam0[,1:min(c(nPC,ncol(M))), drop=FALSE]
        rownames(Lam) <- rownames(M)
        colnames(Lam) <- paste0("nPC",1:nPC)
      }
    }else{ # user provided it's own Lambda matrix
      Lam0=Lam
      Lam = Lam0[,1:min(c(nPC,ncol(M))), drop=FALSE]
      rownames(Lam) <- rownames(M)
      colnames(Lam) <- paste0("nPC",1:nPC)
    }
  }
  if (inherits(x, "dgCMatrix") | inherits(x, "matrix")) {
    Z <- x
  }else{
    if (!is.character(x) & !is.factor(x)) {
      namess <- as.character(substitute(list(x)))[-1L]
      Z <- Matrix(x, ncol = 1)
      colnames(Z) <- namess
    }else {
      dummy <- x
      levs <- na.omit(unique(dummy))
      if (length(levs) > 1) {
        Z <- Matrix::sparse.model.matrix(~dummy - 1, na.action = na.pass)
        colnames(Z) <- gsub("dummy", "", colnames(Z))
      } else {
        vv <- which(!is.na(dummy))
        Z <- Matrix(0, nrow = length(dummy))
        Z[vv, ] <- 1
        colnames(Z) <- levs
      }
    }
  }
  
  if(is.null(M)){
    Zstar <- Lam
  }else{
    Zstar <- as.matrix(Z %*% Lam[colnames(Z),])
  }
  
  if(returnLam){
    return(list(Z = Zstar, Lam=Lam, Lam0=Lam0)) 
  }else{return(Zstar)}
  
}

imputev <- function (x, method = "median") {
  if (is.numeric(x)) {
    if (method == "mean") {
      x[which(is.na(x))] <- mean(x, na.rm = TRUE)
    }
    else if (method == "median") {
      x[which(is.na(x))] <- median(x, na.rm = TRUE)
    }
    else {
      x[which(is.na(x))] <- mean(x, na.rm = TRUE)
    }
  }
  else {
    if (method == "mean") {
      stop("Method 'mean' is not available for non-numeric vectors.", 
           call. = FALSE)
    }
    else if (method == "median") {
      tt <- table(x)
      x[which(is.na(x))] <- names(tt)[which(tt == max(tt))]
    }
    else {
      x[which(is.na(x))] <- names(tt)[which(tt == max(tt))]
    }
  }
  return(x)
}

rrc <- function(x=NULL, H=NULL, nPC=2, returnGamma=FALSE, cholD=TRUE){
  if(is.null(x) ){stop("Please provide the x argument.", call. = FALSE)}
  if(is.null(H) ){stop("Please provide the x argument.", call. = FALSE)}
  # these are called PC models by Meyer 2009, GSE. This is a reduced rank implementation
  # we produce loadings, the Z*L so we can use it to estimate factor scores in mmec()
  
  Y <- apply(H,2, imputev)
  Sigma <- cov(scale(Y, scale = TRUE, center = TRUE)) # surrogate of unstructured matrix to start with
  Sigma <- as.matrix(nearPD(Sigma)$mat)
  # GE <- as.data.frame(t(scale( t(scale(Y, center=T,scale=F)), center=T, scale=F)))  # sum(GE^2)
  if(cholD){
    ## OPTION 2. USING CHOLESKY
    Gamma <- t(chol(Sigma)); # LOADINGS  # same GE=LL' from cholesky  plot(unlist(Gamma%*%t(Gamma)), unlist(GE))
  }else{
    ## OPTION 1. USING SVD
    U <- svd(Sigma)$u;  # V <- svd(GE)$v
    D <- diag(svd(Sigma)$d)
    Gamma <- U %*% sqrt(D); # LOADINGS
    rownames(Gamma) <- colnames(Gamma) <- rownames(Sigma)
  }
  colnamesGamma <- colnames(Gamma)
  rownamesGamma <- rownames(Gamma)
  Gamma <- Gamma[,1:nPC, drop=FALSE]; 
  colnames(Gamma) <- colnamesGamma[1:nPC]
  rownames(Gamma) <- rownamesGamma
  ##
  rownames(Gamma) <- gsub("v.names_","",rownames(Gamma))#rownames(GE)#levels(dataset$Genotype);  # rownames(Se) <- colnames(GE)#levels(dataset$Environment)
  colnames(Gamma) <- paste("PC", 1:ncol(Gamma), sep =""); # 
  ######### GEreduced = Sg %*% t(Se) 
  # if we want to merge with PCs for environments
  dtx <- data.frame(timevar=x)
  dtx$index <- 1:nrow(dtx)
  Z <- Matrix::sparse.model.matrix(~timevar -1, na.action = na.pass, data=dtx)
  colnames(Z) <- gsub("timevar","",colnames(Z))
  Zstar <- Z%*%Gamma[colnames(Z),] # we multiple original Z by the LOADINGS
  Zstar <- as.matrix(Zstar)
  rownames(Z) <- NULL
  
  if(returnGamma){
    return(list(Gamma=Gamma, H=H, Sigma=Sigma))
  }else{
    return(Zstar)
  }
}

dsc <- function(x, thetaC=NULL, theta=NULL){
  if(is.matrix(x)){
    dummy <- x
    mm <- diag(1,ncol(x))
  }else{
    if(!is.character(x) & !is.factor(x)){
      namess <- as.character(substitute(list(x)))[-1L]
      dummy <- as(Matrix(x,ncol=1), Class = "dgCMatrix"); colnames(dummy) <- namess
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
  colnames(mm) <- rownames(mm) <- colnames(dummy)
  bnmm <- mm*0.15
  if(nrow(bnmm) > 1){
    bnmm[upper.tri(bnmm)]=bnmm[upper.tri(bnmm)]/2
  }
  if(!is.null(thetaC)){
    mm <- thetaC
    colnames(mm) <- rownames(mm) <- colnames(dummy)
  }
  if(!is.null(theta)){
    bnmm <- theta
    colnames(bnmm) <- rownames(bnmm) <- colnames(dummy)
  }
  mm[lower.tri(mm)]=0
  return(list(Z=dummy,thetaC=mm, theta=bnmm))
}

atcg1234 <- function(data, ploidy=2, format="ATCG", maf=0, multi=TRUE, silent=FALSE, 
                     by.allele=FALSE, imp=TRUE, ref.alleles=NULL){
  
  impute.mode <- function(x) {
    ix <- which(is.na(x))
    if (length(ix) > 0) {
      x[ix] <- as.integer(names(which.max(table(x))))
    }
    return(x)
  }
  ##### START GBS.TO.BISNP DATA ######
  gbs.to.bisnp <- function(x) {
    y <- rep(NA,length(x))
    y[which(x=="A")] <- "AA"
    y[which(x=="T")] <- "TT"
    y[which(x=="C")] <- "CC"
    y[which(x=="G")] <- "GG"
    y[which(x=="R")] <- "AG"
    y[which(x=="Y")] <- "CT"
    y[which(x=="S")] <- "CG"
    y[which(x=="W")] <- "AT"
    y[which(x=="K")] <- "GT"
    y[which(x=="M")] <- "AC"
    y[which(x=="+")] <- "++"
    y[which(x=="0")] <- "NN"
    y[which(x=="-")] <- "--"
    y[which(x=="N")] <- NA
    return(y)
  }
  ##### END GBS.TO.BISNP DATA ######
  imputeSNP <- function(data){
    #######
    data2 <- apply(data,2,function(x){
      areNA <- which(is.na(x))
      if(length(areNA)>0){
        pos.all <- table(data[,1])
        totake <- names(pos.all)[which(pos.all == max(pos.all))]
        x[areNA] <- totake
      }
      return(x)
    })
    #######
    return(data2)
  }
  #### apply with progress bar ######
  apply_pb <- function(X, MARGIN, FUN, ...){
    env <- environment()
    pb_Total <- sum(dim(X)[MARGIN])
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total,
                         style = 3)
    
    wrapper <- function(...)
    {
      curVal <- get("counter", envir = env)
      assign("counter", curVal +1 ,envir= env)
      setTxtProgressBar(get("pb", envir= env),
                        curVal +1)
      FUN(...)
    }
    res <- apply(X, MARGIN, wrapper, ...)
    close(pb)
    res
  }
  ###### zero.one function
  zero.one <- function(da){
    # this function takes a matrix of markers in biallelic format and returns a matrix of
    # presense/absense of alleles
    mar.nam <- colnames(da)#unique(gsub("\\.\\d","", names(da))) # find a dot and a number after the dot
    mat.list <- list(NA) # list of matrices for each marker
    wi=0 # counter
    if(!silent){
      count <- 0
      tot <- length(mar.nam)
      pb <- txtProgressBar(style = 3)
      setTxtProgressBar(pb, 0)
    }
    for(i in 1:length(mar.nam)){ # for each marker
      wi=wi+1
      if(!silent){
        count <- count + 1
      }
      
      v <- which(colnames(da)==mar.nam[i])#grep(mar.nam[i], colnames(da))
      
      if(length(v)==0){
        qqqqq <- grep(mar.nam[i-1],names(da))
        qqqqq2 <- names(da)[qqqqq[length(qqqqq)] + 1]
        
        stop(paste("Marker",qqqqq2,"has a problem"), call.=FALSE)
      }else if(length(v) == 1){ # for markers with a single column
        prov <- matrix(da[,v])
      }else{prov <- da[,v]}
      ##################################
      alls <- unique(unlist(strsplit(prov,"")))
      alls <- alls[which(!is.na(alls))]
      ninds <- dim(prov)[1]
      fff <- apply(data.frame(alls),1,function(h){
        temp <- numeric(length = ninds)
        temp[grep(h,prov)]<-1
        #make sure is full rank
        
        return(temp)
      })#1 # assigning 1's
      #if(FULL){ # if user want to make sure only get the columns that will ensure full rank
      #  fff <- t(unique(t(fff)))
      #}
      colnames(fff) <- paste(mar.nam[i],alls, sep="/")
      
      mat.list[[i]] <- fff
      if(!silent){
        setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
      }
      
    }
    
    fin.mat <- do.call(cbind,mat.list)
    rownames(fin.mat) <- rownames(da)
    #############
    return(fin.mat)
  }
  
  ## remove all markers or columns that are all missing data
  all.na <- apply(data,2,function(x){length(which(is.na(x)))/length(x)})
  bad.na <- which(all.na==1)
  if(length(bad.na) > 0){
    data <- data[,-bad.na]
  }
  
  if(is.null(ref.alleles)){
    #############################
    if(by.allele){ ####&&&&&&&&&&&&&&&&&&&&&& use zero.one function
      ncolsData <- dim(data)[2]
      ncolsData <- max(ncolsData,round(ncolsData/20))
      user.code <- apply(data[,c(1:ncolsData),drop=FALSE], 2, function(x){q <- which(!is.na(x))[1];ss1 <- substr(x[q], start=1,stop=1);ss2 <- substr(x[q], start=2,stop=2);vv1 <-which(c(ss1,ss2)=="");if(length(vv1)>0){y <-1}else{y <- 0}; return(y)})
      
      AA <- sum(user.code, na.rm = TRUE)/length(user.code)
      if(AA > .9){ # means user is using single letter
        rnd <- rownames(data)
        data <- apply(data,2,gbs.to.bisnp);#W2[1:5,1:5]
        rownames(data) <- rnd
      }
      M <- zero.one(data)
      
    }else{ ###&&&&&&&&&&&&&&&&&&&&&&&&
      n.g <- apply(data,2,function(x){length(table(x))})
      bad <- which(n.g > 3)
      if(length(bad) == dim(data)[2]){
        stop("Error. All your markers are multiallelic. This function requires at least one bi-allelic marker\n")
      }
      
      # tells you which markers have double letter code, i.e. TT instead of T
      # 1: has only one letter
      # 0: has two letters
      ncolsData <- dim(data)[2]
      ncolsData <- max(ncolsData,round(ncolsData/20))
      user.code <- apply(data[,c(1:ncolsData), drop=FALSE], 2, function(x){q <- which(!is.na(x))[1];ss1 <- substr(x[q], start=1,stop=1);ss2 <- substr(x[q], start=2,stop=2);vv1 <-which(c(ss1,ss2)=="");if(length(vv1)>0){y <-1}else{y <- 0}; return(y)})
      AA <- sum(user.code, na.rm = TRUE)/length(user.code)
      if(AA > .9){
        rrn <- rownames(data)
        
        message("Converting GBS or single-letter code to biallelic code\n")
        if(silent){
          data <- apply(data, 2,gbs.to.bisnp)
        }else{
          data <- apply_pb(data, 2,gbs.to.bisnp) 
        }
        rownames(data) <- rrn
        data <- as.data.frame(data)
      }
      #### apply with progress bar ######
      s1 <- rownames(data)
      s2 <- colnames(data)
      data <- as.data.frame(t(data))
      rownames(data) <- s2
      colnames(data) <- s1
      bases <- c("A", "C", "G", "T","l","m","n","p","h","k","-","+","e","f","g","a","b","c","d")
      ## get reference allele function
      get.ref <- function(x, format) {
        if (format == "numeric") {
          ref.alt <- c(0, 1)
        }
        if (format == "AB") {
          ref.alt <- c("A", "B")
        }
        if (format == "ATCG") {
          y <- paste(na.omit(x), collapse = "")
          ans <- apply(array(bases), 1, function(z, y) {
            length(grep(z, y, fixed = T))
          }, y)
          if (sum(ans) > 2) {
            ref.alt <- (bases[which(ans == 1)])[1:2]
            #stop("Error in genotype matrix: More than 2 alleles")
          }
          if (sum(ans) == 2) {
            ref.alt <- bases[which(ans == 1)]
          }
          if (sum(ans) == 1) {
            ref.alt <- c(bases[which(ans == 1)], NA)
          }
        }
        return(ref.alt)
      }
      
      get.multi <- function(x, format) {
        if (format == "numeric") {
          ref.alt <- c(0, 1)
        }
        if (format == "AB") {
          ref.alt <- c("A", "B")
        }
        if (format == "ATCG") {
          y <- paste(na.omit(x), collapse = "")
          ans <- apply(array(bases), 1, function(z, y) {
            length(grep(z, y, fixed = T))
          }, y)
          if (sum(ans) > 2) {
            ref.alt <- TRUE
          }
          if (sum(ans) == 2) {
            ref.alt <- FALSE
          }
          if (sum(ans) == 1) {
            ref.alt <- FALSE
          }
        }
        return(ref.alt)
      }
      
      ####################################
      ## convert to matrix format
      ####################################
      markers <- as.matrix(data)
      ####################################
      # get reference alleles
      ####################################
      message("Obtaining reference alleles\n")
      if(silent){
        tmp <- apply(markers, 1, get.ref, format=format)
      }else{
        tmp <- apply_pb(markers, 1, get.ref, format=format) 
      }
      
      if(multi){ # if markers with multiple alleles should be removed
        message("Checking for markers with more than 2 alleles. If found will be removed.\n")
        if(silent){
          tmpo <- apply(markers, 1, get.multi, format = format)
        }else{
          tmpo <- apply_pb(markers, 1, get.multi, format = format) 
        }
        ###&&&&&&&&&&&& HERE WE MUST INSERT THE NEW FUNCTIONALITY, WHERE WE DETECTED MULTIPLE ALLELES
        multi.allelic <- which(!tmpo) # good markers
        markers <- markers[multi.allelic,,drop=FALSE]
        tmp <- tmp[, multi.allelic,drop=FALSE]
      }
      
      Ref <- tmp[1, ]
      Alt <- tmp[2, ]
      ####################################
      ## bind reference allele and markers and convert to numeric format based on the 
      # reference/alternate allele found
      ####################################
      message("Converting to numeric format\n")
      if(silent){
        M <- apply(cbind(Ref, markers), 1, function(x) {
          y <- gregexpr(pattern = x[1], text = x[-1], fixed = T)
          ans <- as.integer(lapply(y, function(z) {
            ifelse(z[1] < 0, ploidy, ploidy - length(z))
          }))
          return(ans)
        })
      }else{
        M <- apply_pb(cbind(Ref, markers), 1, function(x) {
          y <- gregexpr(pattern = x[1], text = x[-1], fixed = T)
          ans <- as.integer(lapply(y, function(z) {
            ifelse(z[1] < 0, ploidy, ploidy - length(z))
          }))
          return(ans)
        })
      }
      gid.geno <- s1 #colnames(geno)
      rownames(M) <- gid.geno
      ####################################
      # identify bad markers
      ####################################
      bad <- length(which(!is.element(na.omit(M), 0:ploidy)))
      if (bad > 0) {
        stop("Invalid marker calls.")
      }
      
    }
    #rownames(M) <- rownames(data)
    ####################################
    rownames(tmp) <- c("Alt","Ref")
  }else{# user provides reference alleles and just want a conversion
    
    common.mark <- intersect(colnames(data), colnames(ref.alleles))
    data <- data[,common.mark, drop=FALSE]
    tmp <- ref.alleles[,common.mark, drop=FALSE]; #rownames(refa) <- c("Alt","Ref")
    message("Converting to numeric format\n")
    M <- apply_pb(data.frame(1:ncol(data)),1,function(k){
      x <- as.character(data[,k])
      x2 <- strsplit(x,"")
      x3 <- unlist(lapply(x2,function(y){length(which(y == tmp[2,k]))}))
      return(x3)
    })
    #M <- M-1
    colnames(M) <- colnames(data)
    
  }
  
  ####################################
  # by column or markers calculate MAF
  ####################################
  message("Calculating minor allele frequency (MAF)\n")
  if(silent){
    MAF <- apply(M, 2, function(x) {
      AF <- mean(x, na.rm = T)/ploidy
      MAF <- ifelse(AF > 0.5, 1 - AF, AF)
    })
  }else{
    MAF <- apply_pb(M, 2, function(x) {
      AF <- mean(x, na.rm = T)/ploidy
      MAF <- ifelse(AF > 0.5, 1 - AF, AF)
    })
  }
  ####################################
  # which markers have MAF > 0, JUST GET THOSE
  ####################################
  polymorphic <- which(MAF > maf)
  M <- M[, polymorphic, drop=FALSE]
  ####################################
  # function to impute markers with the mode
  ####################################
  
  # time to impute
  if(imp){
    missing <- which(is.na(M))
    if (length(missing) > 0) {
      message("Imputing missing data with mode \n")
      if(silent){
        M <- apply(M, 2, impute.mode)
      }else{
        M <- apply_pb(M, 2, impute.mode)
      }
    }
  }else{
    message("Imputation not required. Be careful using non-imputed matrices in mixed model solvers\n")
  }
  ## ploidy 2 needs to be adjusted to -1,0,1
  # if(ploidy == 2){
  #   M <- M - 1
  # }
  
  return(list(M=M,ref.alleles=tmp))
}

build.HMM <- function(M1,M2, custom.hyb=NULL, return.combos.only=FALSE, separator=":"){
  # build hybrid marker matrix
  
  if(!is.null(custom.hyb)){
    pheno <- custom.hyb
    found <- length(which(colnames(pheno) %in% c("Var1","Var2","hybrid")))
    if(found != 3){
      stop("Column names Var1, Var2, hybrid need to be present when you provide \n       a data table to customize the hybrid genotypes to be build.\n", call. = FALSE)
    }
    return.combos.only=FALSE
  }else{
    a <- rownames(M1)
    b <- rownames(M2)
    pheno <- expand.grid(a,b)
    pheno <- pheno[!duplicated(t(apply(pheno, 1, sort))),]
    pheno$hybrid <- paste(pheno$Var1, pheno$Var2, sep=separator)
  }
  
  if(!return.combos.only){
    # check that marker matrices are in -1,0,1 format
    checkM1 <- c(length(which(M1 == -1)),length(which(M1 == 1)),length(which(M1 == 2)))
    checkM2 <- c(length(which(M2 == -1)),length(which(M2 == 1)),length(which(M2 == 2)))
    
    checkM1[which(checkM1 > 0)] <- 1
    checkM2[which(checkM2 > 0)] <- 1
    
    if(all(checkM1 == c(1,1,0))){ # homo markers were coded correctly as -1,1
    }else if(all(checkM1 == c(0,1,0)) | all(checkM1 == c(1,0,0))){ # homo markers were coded as 0 1
      warning("Either -1 or 1 alleles not detected in M1, we assume you have coded homozygotes \n       as 0 and 1 instead of -1 and 1. We'll fix it.\n")
    }else if(all(checkM1 == c(0,0,1))){ # homo markers were coded as 0 2
      warning("Either -1 or 1 alleles not detected in M1, we assume you have coded homozygotes \n       as 0 and 2 instead of -1 and 1. We'll fix it.\n")
    }
    
    if(all(checkM2 == c(1,1,0))){ # homo markers were coded correctly as -1,1
    }else if(all(checkM2 == c(0,1,0)) | all(checkM2 == c(1,0,0))){ # homo markers were coded as 0 1
      warning("Either -1 or 1 alleles not detected in M2, we assume you have coded homozygotes \n       as 0 and 1 instead of -1 and 1. We'll fix it.\n")
    }else if(all(checkM2 == c(0,0,1))){ # homo markers were coded as 0 2
      warning("Either -1 or 1 alleles not detected in M2, we assume you have coded homozygotes \n       as 0 and 2 instead of -1 and 1. We'll fix it.\n")
    }
    
    ## add markers coming from parents M1
    Z1 <- model.matrix(~Var1-1,pheno);dim(Z1); 
    colnames(Z1) <- gsub("Var1","",colnames(Z1))
    M1 <- M1[colnames(Z1),]
    #M1[1:4,1:4]; Z1[1:4,1:4]; 
    ## add markers coming from parents M2
    Z2 <- model.matrix(~Var2-1,pheno);dim(Z2); 
    colnames(Z2) <- gsub("Var2","",colnames(Z2))
    M2 <- M2[colnames(Z2),]
    #M2[1:4,1:4]; Z2[1:4,1:4];  
    
    ## create the 
    # Z3 <- model.matrix(~hybrid-1,pheno);dim(Z3);
    # colnames(Z3) <- gsub("hybrid","",colnames(Z3))
    # hyb.names <- colnames(Z3)[as.vector(apply(Z3,1,function(x){which(x==1)}))] # names of hybrids
    hyb.names <- pheno$hybrid
    ## marker matrix for hybrids one for each parent
    message(paste("Building hybrid marker matrix for",nrow(Z1),"hybrids\n"))
    
    # M1 <- as(M1, Class="dgCMatrix")
    # M2 <- as(M2, Class="dgCMatrix")
    # Z1 <- as(Z1, Class="dgCMatrix")
    # Z2 <- as(Z2, Class="dgCMatrix")
    
    message("Extracting M1 contribution\n")
    if(all(checkM1 == c(1,1,0))){ # homo markers were coded correctly as -1,1
      Md <- Z1 %*% M1;  # was already converted to -1,1
    }else if(all(checkM1 == c(0,1,0)) | all(checkM1 == c(1,0,0))){ # homo markers were coded as 0 1
      Md <- 2*Z1 %*% M1 - 1;  # 2*Z.dent %*% M.dent - 1   # convert to -1,1
    }else if(all(checkM1 == c(0,0,1))){ # homo markers were coded as 0 2
      Md <- Z1 %*% M1 - 1;  # Z.dent %*% M.dent - 1   # convert to -1,1
    }
    
    message("Extracting M2 contribution\n")
    if(all(checkM2 == c(1,1,0))){ # homo markers were coded correctly as -1,1
      Mf <- Z2 %*% M2;  # was already converted to -1,1
    }else if(all(checkM2 == c(0,1,0)) | all(checkM2 == c(1,0,0))){ # homo markers were coded as 0 1
      Mf <- 2*Z2 %*% M2 - 1;  # 2*Z.dent %*% M.dent - 1   # convert to -1,1
    }else if(all(checkM2 == c(0,0,1))){ # homo markers were coded as 0 2
      Mf <- Z2 %*% M2 - 1;  # Z.dent %*% M.dent - 1   # convert to -1,1
    }
    
    ## marker matrix coded as additive -1,0,1
    Mdf <- (Md + Mf)*(1/2) # normal marker matrix for the hybrids
    rownames(Mdf) <- hyb.names
    #hist(Mdf)
    
    ## dominance matrix for hybrids (0,1 coded)
    Delta <- 1/2*(1 - Md * Mf) #performs element wise multiplication = Hadamard product
    rownames(Delta) <- hyb.names
    #hist(Delta)
    message("Done!!\n")
    return(list(HMM.add=Mdf, HMM.dom=Delta, data.used=pheno))
    
  }else{
    return(list(HMM.add=NA, HMM.dom=NA, data.used=pheno))
  }
}

A.mat <- function (X, min.MAF = NULL) 
{
  
  X <- as.matrix(X)
  n <- nrow(X)
  frac.missing <- apply(X, 2, function(x) {
    length(which(is.na(x)))/n
  })
  missing <- max(frac.missing) > 0
  freq <- apply(X + 1, 2, function(x) {
    mean(x, na.rm = missing)
  })/2
  MAF <- apply(rbind(freq, 1 - freq), 2, min)
  if (is.null(min.MAF)) {
    min.MAF <- 1/(2 * n)
  }
  max.missing <- 1 - 1/(2 * n)
  markers <- which((MAF >= min.MAF) & (frac.missing <= max.missing))
  m <- length(markers)
  var.A <- 2 * mean(freq[markers] * (1 - freq[markers]))
  one <- matrix(1, n, 1)
  mono <- which(freq * (1 - freq) == 0)
  X[, mono] <- 2 * tcrossprod(one, matrix(freq[mono], length(mono), 
                                          1)) - 1
  freq.mat <- tcrossprod(one, matrix(freq[markers], m, 1))
  W <- X[, markers] + 1 - 2 * freq.mat
  A <- tcrossprod(W)/var.A/m
  return(A)
}

bbasis <- function (x, xl, xr, ndx, deg) 
{
  tpower <- function(x, t, p) {
    (x - t)^p * (x > t)
  }
  dx <- (xr - xl)/ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1)/(gamma(deg + 1) * dx^deg)
  B <- (-1)^(deg + 1) * P %*% t(D)
  B
}

tps <- function (columncoordinates, rowcoordinates, nsegments=NULL, 
                 minbound=NULL, maxbound=NULL, degree = c(3, 3), penaltyord = c(2, 2), 
                 nestorder = c(1, 1), asreml = "grp", eigenvalues = "include", 
                 method = "Lee", stub = NULL) 
{
  if (missing(columncoordinates)) 
    stop("columncoordinates argument must be set")
  if (missing(rowcoordinates)) 
    stop("rowcoordinates argument must be set")
  col <- columncoordinates
  nuc <- length(col)
  col.match <- match(columncoordinates, col)
  row <- sort(unique(rowcoordinates))
  nur <- length(row)
  row.match <- match(rowcoordinates, row)
  nv <- length(columncoordinates)
  if (is.null(minbound)) {
    cminval <- min(col)
    rminval <- min(row)
  } else {
    cminval <- min(c(minbound[1], min(col)))
    if (length(minbound) < 2) {
      rminval <- min(c(minbound[1], min(row)))
    }
    else {
      rminval <- min(c(minbound[2], min(row)))
    }
  }
  if (is.null(maxbound)) {
    cmaxval <- max(col)
    rmaxval <- max(row)
  }
  else {
    cmaxval <- max(c(maxbound[1], max(col)))
    if (length(maxbound) < 2) {
      rmaxval <- max(c(maxbound[1], max(row)))
    }
    else {
      rmaxval <- max(c(maxbound[2], max(row)))
    }
  }
  if (is.null(nsegments)) {
    nsegcol <- nuc - 1
    nsegrow <- nur - 1
  }
  else {
    nsegcol <- max(c(nsegments[1], 2))
  }
  if (length(nsegments) < 2) {
    nsegrow <- max(c(nsegments[1], 2))
  }
  else {
    nsegrow <- max(c(nsegments[2], 2))
  }
  nestcol <- floor(nestorder[1])
  if (length(nestorder) < 2) 
    nestrow <- floor(nestorder[1])
  else nestrow <- floor(nestorder[2])
  nsncol <- 0
  if (nestcol > 1) {
    if (nsegcol%%nestcol != 0) 
      warning("Column nesting ignored: number of column segments must be a multiple of nesting order")
    else nsncol <- nsegcol/nestcol
  }
  nsnrow <- 0
  if (nestrow > 1) {
    if (nsegrow%%nestrow != 0) 
      warning("Row nesting ignored: number of row segments must be a multiple of nesting order")
    else nsnrow <- nsegrow/nestrow
  }
  Bc <- bbasis(col, cminval, cmaxval, nsegcol, degree[1])
  nc <- ncol(Bc)
  if (length(degree) < 2) 
    degr <- degree[1]
  else degr <- degree[2]
  Br <- bbasis(row, rminval, rmaxval, nsegrow, degr)
  nr <- ncol(Br)
  if (nsncol > 0) {
    Bcn <- bbasis(col, cminval, cmaxval, nsncol, degree[1])
    ncn <- ncol(Bcn)
  }
  else ncn <- nc
  if (nsnrow > 1) {
    Brn <- bbasis(row, rminval, rmaxval, nsnrow, degr)
    nrn <- ncol(Brn)
  }
  else nrn <- nr
  diff.c <- penaltyord[[1]]
  Dc <- diff(diag(nc), diff = diff.c)
  svd.c <- svd(crossprod(Dc))
  nbc <- nc - diff.c
  U.Zc <- svd.c$u[, c(1:nbc)]
  U.Xc <- svd.c$u[, -c(1:nbc)]
  L.c <- sqrt(svd.c$d[c(1:nbc)])
  diagc <- L.c^2
  BcU <- Bc %*% U.Zc
  BcX <- Bc %*% U.Xc
  BcULi <- BcU %*% diag(1/L.c)
  if ("include" %in% eigenvalues) {
    BcZmat.df <- as.data.frame(BcULi)
    BcZmat <- BcULi
  }
  else {
    BcZmat.df <- as.data.frame(BcU)
    BcZmat <- BcU
  }
  BcZmat.df$TP.col <- col
  mat1c <- matrix(rep(1, nuc), nrow = nuc)
  BcXadj <- BcX - mat1c %*% t(mat1c) %*% BcX/nuc
  Xfc <- (svd(crossprod(BcXadj)))$u[, c(ncol(BcXadj):1)]
  BcX <- BcX %*% Xfc
  if (BcX[1, 1] < 0) 
    BcX[, 1] <- -1 * BcX[, 1]
  if (BcX[1, 2] > 0) 
    BcX[, 2] <- -1 * BcX[, 2]
  if (nsncol > 0) {
    Dcn <- diff(diag(ncn), diff = diff.c)
    svd.cn <- svd(crossprod(Dcn))
    nbcn <- ncn - diff.c
    U.Zcn <- svd.cn$u[, c(1:nbcn)]
    U.Xcn <- svd.cn$u[, -c(1:nbcn)]
    L.cn <- sqrt(svd.cn$d[c(1:nbcn)])
    BcnU <- Bcn %*% U.Zcn
    BcnX <- Bcn %*% U.Xcn
  }
  else {
    nbcn <- nbc
    BcnU <- BcU
    L.cn <- L.c
  }
  if (length(penaltyord) < 2) {
    diff.r <- penaltyord[1]
  }
  else {
    diff.r <- penaltyord[2]
  }
  Dr <- diff(diag(nr), diff = diff.r)
  svd.r <- svd(crossprod(Dr))
  nbr <- nr - diff.r
  U.Zr <- svd.r$u[, c(1:nbr)]
  U.Xr <- svd.r$u[, -c(1:nbr)]
  L.r <- sqrt(svd.r$d[c(1:nbr)])
  diagr <- L.r^2
  BrU <- Br %*% U.Zr
  BrX <- Br %*% U.Xr
  BrULi <- BrU %*% diag(1/L.r)
  if ("include" %in% eigenvalues) {
    BrZmat.df <- as.data.frame(BrULi)
    BrZmat <- BrULi
  }
  else {
    BrZmat.df <- as.data.frame(BrU)
    BrZmat <- BrU
  }
  BrZmat.df$TP.row <- row
  mat1r <- matrix(rep(1, nur), nrow = nur)
  BrXadj <- BrX - mat1r %*% t(mat1r) %*% BrX/nur
  Xfr <- (svd(crossprod(BrXadj)))$u[, c(ncol(BrXadj):1)]
  BrX <- BrX %*% Xfr
  if (BrX[1, 1] < 0) 
    BrX[, 1] <- -1 * BrX[, 1]
  if (BrX[1, 2] > 0) 
    BrX[, 2] <- -1 * BrX[, 2]
  if (nsnrow > 0) {
    Drn <- diff(diag(nrn), diff = diff.r)
    svd.rn <- svd(crossprod(Drn))
    nbrn <- nrn - diff.r
    U.Zrn <- svd.rn$u[, c(1:nbrn)]
    U.Xrn <- svd.rn$u[, -c(1:nbrn)]
    L.rn <- sqrt(svd.rn$d[c(1:nbrn)])
    BrnU <- Brn %*% U.Zrn
    BrnX <- Brn %*% U.Xrn
  }
  else {
    nbrn <- nbr
    BrnU <- BrU
    L.rn <- L.r
  }
  A <- 10^(floor(log10(max(row))) + 1)
  row.index <- rep(row, times = nuc)
  col.index <- rep(col, each = nur)
  index <- A * col.index + row.index
  C.R <- A * columncoordinates + rowcoordinates
  BcrZ1 <- BcnU[col.match, ] %x% matrix(rep(1, nbrn), nrow = 1, 
                                        ncol = nbrn)
  BcrZ2 <- matrix(rep(1, nbcn), nrow = 1, ncol = nbcn) %x% 
    BrnU[row.match, ]
  BcrZ <- BcrZ1 * BcrZ2
  diagrx <- rep(L.cn^2, each = nbrn)
  diagcx <- rep(L.rn^2, times = nbcn)
  if ("Lee" %in% method) {
    diagcr <- diagrx + diagcx
  }
  if ("Wood" %in% method) {
    diagcr <- diagrx * diagcx
  }
  if (!("Lee" %in% method) & !("Wood" %in% method)) {
    stop("Invalid setting of method argument")
  }
  BcrZLi <- BcrZ %*% diag(1/sqrt(diagcr))
  if ("include" %in% eigenvalues) {
    BcrZmat.df <- as.data.frame(BcrZLi)
    BcrZmat <- BcrZLi
  }
  else {
    BcrZmat.df <- as.data.frame(BcrZ)
    BcrZmat <- BcrZ
  }
  BcrZmat.df$TP.CxR <- C.R
  tracelist <- list()
  for (i in 1:diff.c) {
    nm <- paste0("Xc", i, ":Zr")
    tempmat <- (BcX[col.match, i] %x% matrix(rep(1, nbr), 
                                             nrow = 1)) * BrZmat[row.match, ]
    if ("include" %in% eigenvalues) 
      tempmatsc <- tempmat
    else tempmatsc <- tempmat * (rep(1, nv) %*% matrix((1/diagr), 
                                                       nrow = 1))
    tracelist[nm] <- sum(tempmatsc * tempmat)
  }
  for (i in 1:diff.r) {
    nm <- paste0("Zc:Xr", i)
    tempmat <- BcZmat[col.match, ] * (matrix(rep(1, nbc), 
                                             nrow = 1) %x% BrX[row.match, i])
    if ("include" %in% eigenvalues) 
      tempmatsc <- tempmat
    else tempmatsc <- tempmat * (rep(1, nv) %*% matrix((1/diagc), 
                                                       nrow = 1))
    tracelist[nm] <- sum(tempmatsc * tempmat)
  }
  if ("include" %in% eigenvalues) 
    tracelist["Zc:Zr"] <- sum(BcrZmat * BcrZmat)
  else {
    tempmatsc <- BcrZmat * (rep(1, nv) %*% matrix((1/diagcr), 
                                                  nrow = 1))
    tracelist["Zc:Zr"] <- sum(tempmatsc * BcrZmat)
  }
  # outdata <- as.data.frame(data)
  outdata <- data.frame(TP.col=columncoordinates)
  outdata$TP.row <- rowcoordinates
  outdata$TP.CxR <- C.R
  BcrX1 <- BcX[col.match, ] %x% matrix(rep(1, diff.r), nrow = 1)
  BcrX2 <- matrix(rep(1, diff.c), nrow = 1) %x% BrX[row.match, 
  ]
  BcrX <- BcrX1 * BcrX2
  fixed <- list()
  fixed$col <- data.frame(row.names = C.R)
  for (i in 1:diff.c) {
    c.fixed <- paste("TP.C", ".", i, sep = "")
    outdata[c.fixed] <- BcX[col.match, i]
    fixed$col[c.fixed] <- BcX[col.match, i]
  }
  fixed$row <- data.frame(row.names = C.R)
  for (i in 1:diff.r) {
    r.fixed <- paste("TP.R", ".", i, sep = "")
    outdata[r.fixed] <- BrX[row.match, i]
    fixed$row[r.fixed] <- BrX[row.match, i]
  }
  ncolX <- diff.c * diff.r
  fixed$int <- data.frame(row.names = C.R)
  for (i in 1:ncolX) {
    cr.fixed <- paste("TP.CR", ".", i, sep = "")
    outdata[cr.fixed] <- BcrX[, i]
    fixed$int[cr.fixed] <- BcrX[, i]
  }
  if (!missing(stub)) {
    cname <- paste0("BcZ", stub, ".df")
    rname <- paste0("BrZ", stub, ".df")
    crname <- paste0("BcrZ", stub, ".df")
  }
  else {
    cname <- "BcZ.df"
    rname <- "BrZ.df"
    crname <- "BcrZ.df"
  }
  mbftext <- paste0("list(TP.col=list(key=c(\"TP.col\",\"TP.col\"),cov=\"", 
                    cname, "\"),")
  mbftext <- paste0(mbftext, "TP.row=list(key=c(\"TP.row\",\"TP.row\"),cov=\"", 
                    rname, "\"),")
  mbftext <- paste0(mbftext, "TP.CxR=list(key=c(\"TP.CxR\",\"TP.CxR\"),cov=\"", 
                    crname, "\"))")
  mbflist <- eval(parse(text = mbftext))
  if ("grp" %in% asreml) {
    grp <- list()
    listnames <- list()
    start <- length(outdata)
    start0 <- start
    scale <- 1
    j <- 1
    for (i in 1:diff.c) {
      nm0 <- paste0(names(fixed$col[i]), "_frow")
      listnames[j] <- nm0
      for (k in 1:nbr) {
        nm <- paste0(nm0, "_", k)
        outdata[nm] <- scale * fixed$col[[i]] * BrZmat[row.match, 
                                                       k]
      }
      grp[[j]] <- seq(from = start + 1, to = start + nbr, 
                      by = 1)
      start <- start + nbr
      j <- j + 1
    }
    for (i in 1:diff.r) {
      nm0 <- paste0(names(fixed$row[i]), "_fcol")
      listnames[j] <- nm0
      for (k in 1:nbc) {
        nm <- paste0(nm0, "_", k)
        outdata[nm] <- scale * fixed$row[[i]] * BcZmat[col.match, 
                                                       k]
      }
      grp[[j]] <- seq(from = start + 1, to = start + nbc, 
                      by = 1)
      start <- start + nbc
      j <- j + 1
    }
    m <- 0
    nm0 <- "TP_fcol_frow"
    listnames[j] <- nm0
    for (k in 1:(nbrn * nbcn)) {
      nm <- paste0(nm0, "_", k)
      outdata[nm] <- scale * BcrZmat[, k]
    }
    grp[[j]] <- seq(from = start + 1, to = start + (nbcn * 
                                                      nbrn), by = 1)
    end <- start + (nbcn * nbrn)
    j <- j + 1
    listnames[j] <- "All"
    grp[[j]] <- seq(from = start0 + 1, to = end, by = 1)
    grp <- structure(grp, names = listnames)
  }
  if ("sepgrp" %in% asreml) {
    grp <- list()
    listnames <- list()
    start <- length(outdata)
    nm0 <- "TP_C"
    listnames[1] <- nm0
    for (i in 1:diff.c) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- fixed$col[[i]]
    }
    grp[[1]] <- seq(from = start + 1, to = start + diff.c, 
                    by = 1)
    start <- start + diff.c
    nm0 <- "TP_R"
    listnames[2] <- nm0
    for (i in 1:diff.r) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- fixed$row[[i]]
    }
    grp[[2]] <- seq(from = start + 1, to = start + diff.r, 
                    by = 1)
    start <- start + diff.r
    nm0 <- "TP_fcol"
    listnames[3] <- nm0
    for (k in 1:nbc) {
      nm <- paste0(nm0, "_", k)
      outdata[nm] <- BcZmat[col.match, k]
    }
    grp[[3]] <- seq(from = start + 1, to = start + nbc, by = 1)
    start <- start + nbc
    nm0 <- "TP_frow"
    listnames[4] <- nm0
    for (k in 1:nbr) {
      nm <- paste0(nm0, "_", k)
      outdata[nm] <- BrZmat[row.match, k]
    }
    grp[[4]] <- seq(from = start + 1, to = start + nbr, by = 1)
    start <- start + nbr
    grp <- structure(grp, names = listnames)
    nm0 <- "TP_fcol_frow"
    listnames[5] <- nm0
    for (k in 1:(nbrn * nbcn)) {
      nm <- paste0(nm0, "_", k)
      outdata[nm] <- BcrZmat[, k]
    }
    grp[[5]] <- seq(from = start + 1, to = start + (nbcn * 
                                                      nbrn), by = 1)
    grp <- structure(grp, names = listnames)
  }
  if ("own" %in% asreml) {
    grp <- list()
    listnames <- list()
    listnames[1] <- "All"
    start <- length(outdata)
    nm0 <- "Xc_Zr"
    Xc_Zr <- (BcX[col.match, ] %x% matrix(rep(1, nbr), nrow = 1)) * 
      (matrix(rep(1, diff.c), nrow = 1) %x% BrZmat[row.match, 
      ])
    nXc_Zr <- ncol(Xc_Zr)
    for (i in 1:nXc_Zr) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- Xc_Zr[, i]
    }
    nm0 <- "Zc_Xr"
    Zc_Xr <- (BcZmat[col.match, ] %x% matrix(rep(1, diff.r), 
                                             nrow = 1)) * (matrix(rep(1, nbc), nrow = 1) %x% BrX[row.match, 
                                             ])
    nZc_Xr <- ncol(Zc_Xr)
    for (i in 1:nZc_Xr) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- Zc_Xr[, i]
    }
    nm0 <- "Zc_Zr"
    Zc_Zr <- BcrZmat
    nZc_Zr <- ncol(Zc_Zr)
    for (i in 1:nZc_Zr) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- Zc_Zr[, i]
    }
    grp[[1]] <- seq(from = start + 1, to = start + nXc_Zr + 
                      nZc_Xr + nZc_Zr, by = 1)
    grp <- structure(grp, names = listnames)
  }
  res <- list()
  res$data <- outdata
  res$mbflist <- mbflist
  res[["BcZ.df"]] <- BcZmat.df
  res[["BrZ.df"]] <- BrZmat.df
  res[["BcrZ.df"]] <- BcrZmat.df
  res[["All"]] <- as.matrix(outdata[,grp$All])
  res$dim <- c(diff.c = diff.c, nbc = nbc, nbcn = nbcn, diff.r = diff.r, 
               nbr = nbr, nbrn = nbrn)
  res$trace <- tracelist
  if ("grp" %in% asreml) 
    res$grp <- grp
  if ("sepgrp" %in% asreml) 
    res$grp <- grp
  if ("own" %in% asreml) 
    res$grp <- grp
  if ("mbf" %in% asreml) 
    res$grp <- NULL
  if (!("include" %in% eigenvalues)) 
    res$eigen <- list(diagc = diagc, diagr = diagr, diagcr = diagcr)
  res
}

