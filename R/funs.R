#### "relmat" class methods


lmebreed <-
    function(formula, data, family = NULL, REML = TRUE, relmat = list(), 
             addmat=list(),
             control = list(), start = NULL, verbose = FALSE, 
             subset, weights, na.action, offset, contrasts = NULL,
             model = TRUE, x = TRUE, ...)
{
    gaus <- FALSE
    if (is.null(family)) {
        gaus <- TRUE
    } else {
        ## copied from glm()
        if (is.character(family)) 
            family <- get(family, mode = "function", envir = parent.frame())
        if (is.function(family)) 
            family <- family()
        if (!inherits(family, "family")) stop("unknown family type")
        gaus <- family$family == "gaussian" && family$link == "identity"
    }
    mc <- match.call()
    lmerc <- mc                         # create a call to lmer
    # lmerc$formula <- formula; lmerc$data <- data
    lmerc[[1]] <- if (gaus) as.name("lmer") else as.name("glmer")
    lmerc$relmat <- NULL
    lmerc$addmat <- NULL
    if (!gaus) lmerc$REML <- NULL

    if (!length(relmat) & !length(addmat))  {
      return(eval.parent(lmerc))
    }            # call [g]lmer instead
    
    stopifnot(is.list(relmat),        # check the relmat argument
              length(names(relmat)) == length(relmat),
              all( sapply(relmat, inherits, what = c("relmat","matrix","dtCMatrix"))  ))
    
    lmf <- eval(lmerc, parent.frame()) 

    
    relfac <- relmat          # copy the relmat list for relfactor
    pnms <- names(relmat)
    pnms2 <- names(addmat)
    pp <- lmf@pp
    resp <- lmf@resp
    fl <- lmf@flist
    stopifnot(all(pnms %in% names(fl)))
    asgn <- attr(fl, "assign")
    Zt <- pp$Zt
    ## replace additional matrices
    for (i in seq_along(addmat)) {
      tn0 <- which(match(pnms2[i], names(fl)) == asgn)
      for(j in 1:length(tn0)){ # diagonal and unstructured models require to multiple all matrices by the same relfactor
        ind <- (lmf@Gp)[tn0[j]:(tn0[j]+1L)]
        rowsi <- (ind[1]+1L):ind[2]
        provZt <- t(addmat[[i]])
        Zt[rowsi,] <- provZt[rownames(Zt[rowsi,]),]
      }
    }
    ## use the relfactors
    for (i in seq_along(relmat)) {
        tn <- which(match(pnms[i], names(fl)) == asgn)
        for(j in 1:length(tn)){ # diagonal and unstructured models require to multiple all matrices by the same relfactor
          ind <- (lmf@Gp)[tn[j]:(tn[j]+1L)]
          rowsi <- (ind[1]+1L):ind[2]
          pick <- intersect( rownames(Zt), rownames(relfac[[i]])  )
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
    cls <- if (gaus) "lmerlmebreed" else "glmerlmebreed"
    ans <- do.call(new, list(Class=cls, relfac=relfac,
                             frame=mm@frame, flist=mm@flist, cnms=mm@cnms, Gp=mm@Gp,
                             theta=mm@theta, beta=mm@beta,u=mm@u,lower=mm@lower,
                             devcomp=mm@devcomp, pp=mm@pp,resp=mm@resp,optinfo=mm@optinfo))
    ans@call <- evalq(mc)
    ans
}

setMethod("ranef", signature(object = "lmebreed"),
          function(object, postVar = FALSE, drop = FALSE, whichel = names(ans), relmat = TRUE, ...)
      {
          if ((postVar <- as.logical(postVar)) && (relmat <- as.logical(relmat)))
              stop("code for applying relmat and posterior variances not yet written")
          ans <- ranef(as(object, "merMod"), postVar, drop = FALSE)
          ans <- ans[whichel]
          if (relmat) {
              if (postVar)
                  stop("postVar and relmat cannot both be true")
              rf <- object@relfac
              for (nm in names(rf)) {
                  dm <- data.matrix(ans[[nm]])
                  cn <- colnames(dm)
                  rn <- rownames(dm)
                  dm <- as.matrix(rf[[nm]] %*% dm)
                  colnames(dm) <- cn
                  rownames(dm) <- rn
                  ans[[nm]] <- data.frame(dm, check.names = FALSE)
              }
          }
          if (drop)
              ans <- lapply(ans, function(el)
                        {
                            if (ncol(el) > 1) return(el)
                            pv <- drop(attr(el, "postVar"))
                            el <- drop(as.matrix(el))
                            if (!is.null(pv))
                                attr(el, "postVar") <- pv
                            el
                        })
          ans
      })


setMethod("fitted", signature(object = "lmebreed"),
          function(object, ...) {
              stop("fitted() not applicable to lmebreed objects")
          })


setMethod("residuals", signature(object = "lmebreed"),
          function(object, ...) {
              stop("residuals() not applicable to lmebreed objects")
          })

