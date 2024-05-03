

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

rrc <- function(x=NULL, H=NULL, nPC=2, returnGamma=FALSE, cholD=TRUE){
  if(is.null(x) ){stop("Please provide the x argument.", call. = FALSE)}
  if(is.null(H) ){stop("Please provide the x argument.", call. = FALSE)}
  # these are called PC models by Meyer 2009, GSE. This is a reduced rank implementation
  # we produce loadings, the Z*L so we can use it to estimate factor scores in mmec()
  Y <- apply(H,2, sommer::imputev)
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
    return(list(Gamma=Gamma, H0=H0, Sigma=Sigma))
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
