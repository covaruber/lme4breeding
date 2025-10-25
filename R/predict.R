predict.lmeb <- function(object, hyperTable=NULL, classify=NULL, usePEV=FALSE, ...)  {
  
  if(is.null(classify)){
    stop("Please provide the classify argument to build the D matrix.", call. = FALSE )
  }
  '%!in%' <- function(x,y)!('%in%'(x,y))
  if(is.null(hyperTable)){ # default rules for the hypertable
    message(magenta("hyperTable argument not provided. Building a hyper table based on the classify argument. Please check the output slot 'hyperTable' to ensure that the different effects have been included and average as you expected."))
    hyperTable <- Dtable(object)
    # if the user wants a simple averaging  (no include) we add 1s to all rows and 'average'
    hyperTable$average[which(hyperTable$type == "fixed")]=1
    # if the group hypertable matches perfectly the classify we do 'include'
    perfect <- which(hyperTable$group %in% classify)
    hyperTable$include[perfect]=1
    hyperTable$average[perfect]=0
    # if the group hypertable includes an intercept we do 'include'
    hyperTable$include[which(hyperTable$group %in% "(Intercept)")]=1
    # if the group hypertable matches imperfectly the classify we do 'include' and 'average'
    imperfect <- which( unlist( lapply(as.list(hyperTable$group ), function(x){
      xx <- strsplit(x,split=":")[[1]]
      myMatch <- length(which( xx %in% classify ))
      return(myMatch)
    })
    ) > 0)
    imperfect <- setdiff(imperfect, perfect)
    if(length(imperfect) > 0){
      hyperTable$include[imperfect]=1
      hyperTable$average[imperfect]=1
    }
  }
  # get all information from the model
  BLUP <- ranef(object, condVar=TRUE, includeCVM=TRUE)
  # get D table information
  mapCondVar <- attr(BLUP, which="mapCondVar")
  # get inverse of coefficient matrix
  if(usePEV){
    condVarMat <- getMME(object=object, vc=VarCorr(object) )$Ci
  }else{
    condVarMat <- attr(BLUP, which="condVarMat")
  }
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
  for(iRow in 1:nrow(hyperTable)){ # iRow=2
    iVar <- hyperTable[iRow,"variable"]
    iGroup <- hyperTable[iRow,"group"]
    # if we want to 'include' (nested we decide if we average or not)
    if(hyperTable[iRow,"include"]>0){
      v <- which(mapCondVar[,"variable"]==iVar ) # match the intercept
      m <- which( unlist(lapply(as.list(mapCondVar[,"group"]), function(x){length(which( c( strsplit(x, split = ":")[[1]], x )== iGroup ))} )) > 0 ) # match the slope
      v <- intersect(v,m) # columns involving in one way or another the slope within this intecept: hyperTable[iRow,"group"]
      for (jRow in 1:nrow(D)) { # jRow=3
        ## direct match (same variable/intercept and group/slope)
        w <- which( unlist( lapply(as.list(mapCondVar[,"level"]), function(x){
          length( which(strsplit(x, split = ":")[[1]] %in%
                          c(
                            rownames(D)[jRow], 
                            levsOr[jRow],
                            strsplit(rownames(D)[jRow], split = ":")[[1]],
                            strsplit(levsOr[jRow], split = ":")[[1]],
                            "(Intercept)"
                          )
          ) )
        } ) ) > 0 )
        
        myMatch <- intersect(v,w)
        if (length(myMatch) > 0) {
          D[jRow, myMatch] = 1
        }
        
        ## indirect match
      }
      # in addition to include we ask to average
      if(hyperTable[iRow,"average"]>0){
        nReps <- max(apply(D[,v,drop=FALSE],1,function(x){length(which(x>0))}))
        D[, v] = D[, v]/nReps
      }
    }
    # if simple averaging we just add 1s to all rows first, 
    # then we divide over the number of replicates for that effect
    if(hyperTable[iRow,"average"]>0 & hyperTable[iRow,"include"]==0){
      v <- which(mapCondVar[,"variable"]==iVar & mapCondVar[,"group"]==iGroup )
      D[, v] = 1
      nReps <- max(apply(D[,v,drop=FALSE],1,function(x){length(which(x>0))}))
      if(hyperTable[iRow,"type"] == "fixed"){
        # if averaging effect we need to check if intercept exist and add it
        nReps <- nReps + length(which(mapCondVar[,"level"] %in% "(Intercept)"))
      }
      D[, v, drop=FALSE] =  D[, v, drop=FALSE]/nReps
    }
    ## end of rules
  }
  # compute the predicted values and std errors
  predicted.value <- D %*% b
  vcov <- D %*% condVarMat %*% t(D)
  std.error <- sqrt(diag(vcov))
  pvals <- data.frame(id = rownames(D),
                      predicted.value = predicted.value[,1], 
                      std.error = std.error)
  # compile results
  ans <- list(pvals=pvals, b=b, condVarMat=condVarMat, D=D, mapCondVar=mapCondVar, hyperTable=hyperTable, classify=classify )
  return(ans)
}
