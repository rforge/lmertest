### NOT USED
Dev2 <- function(rho, vec.matr, nll = FALSE) {
  ### Deviance of a LMM as a function of the variance-covariance
  ### parameters.
  ### nll: should the negative log-likelihood rather than the deviance
  ### (= 2*nll) be returned?  
  sigma2 <- vec.matr[[1]]
  Lambda <- makeLambda2(rho, c(sqrt(vec.matr[[1]]), sqrt(vec.matr[-1]/vec.matr[[1]])))
  #Lambdat@x[] <<- theta[Lind] 
  Ut  <-  crossprod(Lambda,rho$Zt)
  L <- Cholesky(tcrossprod(Ut), LDL = FALSE, Imult = 1)
  cu <- solve(L, solve(L, Ut %*% rho$y, sys = "P"), sys = "L")
  RZX <- solve(L, solve(L, Ut %*% rho$X, sys = "P"), sys = "L")
  RX <- chol(rho$XtX - crossprod(RZX))
  cb <- solve(t(RX),crossprod(rho$X,rho$y)- crossprod(RZX, cu))
  beta <- solve(RX, cb)
  u <- solve(L,solve(L,cu - RZX %*% beta, sys="Lt"), sys="Pt")
  fitted <- as.vector(crossprod(Ut, u) + rho$X %*% beta)
  ## evaluate using dnorm?
  prss <- sum(c(rho$y - fitted, as.vector(u))^2)
  ## rho$prss <- prss
  n <- length(fitted); p <- ncol(RX)
  ## ML deviance:
  dev <- as.vector(n * log(2 * pi * sigma2) + prss / sigma2 +
                     c(2 * determinant(L)$modulus))
  if(rho$REML) ## REML deviance:
    dev <- dev + as.vector(c(2 * determinant(RX)$modulus) -
                             p * log(2 * pi * sigma2))
  if(nll) ## return negative log-likelihood rather than deviance?
    dev <- dev/2
  return(as.vector(dev))
}


### NOT USED
###########################################################################
#calculates Lambda matrix (see lmer theory)
###########################################################################
makeLambda2 <- function(rho,vec.matr)
{  
  #if there is correlation between intercept and slope in random term
  if(rho$corr.intsl)
  {
    Lambda <- matrix(nrow=0,ncol=0)
    for(i in 1:length(rho$nlev))
    {
      
      lambda1 <- vec.matr[which(rho$param$vec.num==i)+1]
      #the correlation between intercept and slope is present
      if(length(lambda1)>1)
      {
        #one random coefficient
        #ST.full <- c(lambda1[1:2],0,lambda1[length(lambda1)])
        #ST <- matrix(ST.full,nrow=2,ncol=2)
        #multiple random coefficients
        ST <- matrix(0,nrow=rho$param$STdim[i],ncol=rho$param$STdim[i])
        ST[lower.tri(ST, diag=TRUE)] <- lambda1
        Lambda  <-  bdiag(Lambda,kronecker(ST, Diagonal(rho$nlev[i])))
        # a new one
        #Lambda <- bdiag(Lambda,kronecker(Diagonal(rho$nlev[i]), ST))
        
      }
      else
      {
        Lambda <- bdiag(Lambda,kronecker(lambda1, Diagonal(rho$nlev[i])))
      }
    }
  }
  if(!(rho$corr.intsl))
  { 
    Lambda <- Diagonal(x=rep.int(vec.matr[-1],rho$nlev))
  }
  
  return(Lambda)
}


## modified devfun2 from lme4: REML criterion was added
devfun3 <- function (fm, useSc, signames, reml = TRUE) 
{
  stopifnot(is(fm, "merMod"))
  #fm <- refitML(fm)
  #basedev <- deviance(fm)
  vlist <- sapply(fm@cnms, length)
  #sig <- sigma(fm)
  #stdErr <- unname(coef(summary(fm, ddf = "lme4"))[, 2])
  pp <- fm@pp$copy()
  isCor <- checkCorr(fm)
  #if (useSc) {
  #  opt <- Cv_to_Vv(pp$theta, n = vlist, s = sig)
  #  names(opt) <- if (signames) {
  #    c(sprintf(".sig%02d", seq(length(opt) - 1)), ".sigma")
  #  }
    #else {
    #  c(tnames(fm, old = FALSE, prefix = c("sd", "cov")), 
    #    "sigma2")
    #}
  #}
  #else {
  #  opt <- Cv_to_Sv(pp$theta, n = vlist)
  #  names(opt) <- if (signames) {
  #    sprintf(".sig%02d", seq_along(opt))
  #  }
  #  else {
  #    tnames(fm, old = FALSE, prefix = c("sd", "cor"))
  #  }
  #}
  #opt <- c(opt, fixef(fm))
  resp <- fm@resp$copy()
  np <- length(pp$theta)
  nf <- length(fixef(fm))
  if (!isGLMM(fm)) 
    np <- np + 1L
  n <- nrow(pp$V)
  if (isLMM(fm)) {
    ans <- function(pars) {
      stopifnot(is.numeric(pars), length(pars) == np)
      ## TO FIX: when the correlations are present then not execute the next line
      if(!isCor)
        pars[which(pars < 0)] <- 0
      sigma2 <- pars[np]
      thpars <- Vv_to_Cv(pars, n = vlist, s = sqrt(sigma2))
      .Call("lmer_Deviance", pp$ptr(), resp$ptr(), thpars, PACKAGE = "lme4")
      sigsq <- sigma2
      dev <- pp$ldL2() + (resp$wrss() + pp$sqrL(1))/sigsq + n * 
        log(2 * pi * sigsq)      
      if(reml){
        p <- ncol(pp$RX())
        dev <- dev + 2*determinant(pp$RX())$modulus - p * log(2 * pi * sigsq)              
      }
      return(dev)     
    }
  }
  #attr(ans, "optimum") <- opt
  #attr(ans, "basedev") <- basedev
  attr(ans, "thopt") <- pp$theta
  attr(ans, "isCor") <- isCor
  #attr(ans, "stderr") <- stdErr
  class(ans) <- "devfun3"
  ans
}


## modified devfun3: now depends on theta and not covariate parameters
devfun4 <- function (fm, useSc, signames, reml = TRUE) 
{
  stopifnot(is(fm, "merMod"))
 
  vlist <- sapply(fm@cnms, length)
 
  pp <- fm@pp$copy()
  isCor <- checkCorr(fm)
  
  resp <- fm@resp$copy()
  np <- length(pp$theta)
  nf <- length(fixef(fm))
  if (!isGLMM(fm)) 
    np <- np + 1L
  n <- nrow(pp$V)
  if (isLMM(fm)) {
    ans <- function(thpars) {
      stopifnot(is.numeric(thpars), length(thpars) == np)
      ## TO FIX: when the correlations are present then not execute the next line
      #if(!isCor)
      #  pars[which(pars < 0)] <- 0
      #sigma2 <- thpars[np]^2
      #thpars <- Vv_to_Cv(pars, n = vlist, s = sqrt(sigma2))
      .Call("lmer_Deviance", pp$ptr(), resp$ptr(), thpars[-np], PACKAGE = "lme4")
      sigsq <- thpars[np]^2
      dev <- pp$ldL2() + (resp$wrss() + pp$sqrL(1))/sigsq + n * 
        log(2 * pi * sigsq)      
      if(reml){
        p <- ncol(pp$RX())
        dev <- dev + 2*determinant(pp$RX())$modulus - p * log(2 * pi * sigsq)              
      }
      return(dev)     
    }
  }
  #attr(ans, "optimum") <- opt
  #attr(ans, "basedev") <- basedev
  attr(ans, "thopt") <- pp$theta
  attr(ans, "isCor") <- isCor
  #attr(ans, "stderr") <- stdErr
  class(ans) <- "devfun4"
  ans
}



devCritFun <- function(object, REML = NULL) 
{
  if (isTRUE(REML) && !isLMM(object)) 
    stop("can't compute REML deviance for a non-LMM")
  cmp <- object@devcomp$cmp
  if (is.null(REML) || is.na(REML[1])) 
    REML <- isREML(object)
  if (REML) {
    if (isREML(object)) {
      cmp[["REML"]]
    }
    else {
      lnum <- log(2 * pi * cmp[["pwrss"]])
      n <- object@devcomp$dims[["n"]]
      nmp <- n - length(object@beta)
      ldW <- sum(log(weights(object, method = "prior")))
      -ldW + cmp[["ldL2"]] + cmp[["ldRX2"]] + nmp * (1 + 
                                                       lnum - log(nmp))
    }
  }
  else {
    if (!isREML(object)) {
      cmp[["dev"]]
    }
    else {
      n <- object@devcomp$dims[["n"]]
      lnum <- log(2 * pi * cmp[["pwrss"]])
      ldW <- sum(log(weights(object, method = "prior")))
      -ldW + cmp[["ldL2"]] + n * (1 + lnum - log(n))
    }
  }
}



## calc vcov of fixed effects based on var cor parameters of random ## effects
vcovJSS <- function(fm)
{
  stopifnot(is(fm, "merMod"))
  #fm <- refitML(fm)
  #basedev <- deviance(fm)
  vlist <- sapply(fm@cnms, length)
  #sig <- sigma(fm)
  #stdErr <- unname(coef(summary(fm, ddf = "lme4"))[, 2])
  pp <- fm@pp$copy()
  isCor <- checkCorr(fm)
#   if (useSc) {
#     opt <- Cv_to_Vv(pp$theta, n = vlist, s = sig)
#     names(opt) <- if (signames) {
#       c(sprintf(".sig%02d", seq(length(opt) - 1)), ".sigma")
#     }
#     #else {
#     #  c(tnames(fm, old = FALSE, prefix = c("sd", "cov")), 
#     #    "sigma2")
#     #}
#   }
#   else {
#     opt <- Cv_to_Sv(pp$theta, n = vlist)
#     names(opt) <- if (signames) {
#       sprintf(".sig%02d", seq_along(opt))
#     }
#     else {
#       tnames(fm, old = FALSE, prefix = c("sd", "cor"))
#     }
#   }
  #opt <- c(opt, fixef(fm))
  resp <- fm@resp$copy()
  np <- length(pp$theta)
  nf <- length(fixef(fm))
  if (!isGLMM(fm)) 
    np <- np + 1L
  n <- nrow(pp$V)
  if (isLMM(fm)) {
    ans <- function(Lc, pars) {
      stopifnot(is.numeric(pars), length(pars) == np)
      ## TO FIX: when the correlations are present then not execute the next line
      if(!isCor)
        pars[which(pars < 0)] <- 0
      sigma2 <- pars[np]
      thpars <- Vv_to_Cv(pars, n = vlist, s = sqrt(sigma2))
      .Call("lmer_Deviance", pp$ptr(), resp$ptr(), thpars, PACKAGE = "lme4")      
      vcov_out <- sigma2 * tcrossprod(pp$RXi()) #chol2inv(pp$RX())
      #return(as.matrix(vcov_out))
      return(as.matrix(Lc %*% as.matrix(vcov_out) %*% t(Lc)))        
    }
  } 
  #attr(ans, "optimum") <- opt
  #attr(ans, "basedev") <- basedev
  #attr(ans, "thopt") <- pp$theta
  #attr(ans, "stderr") <- stdErr
  attr(ans, "isCor") <- isCor
  class(ans) <- "vcovJSS"
  ans
}

## calc vcov of fixed effects based on theta parameters of random ## effects
vcovJSStheta <- function(fm)
{
  stopifnot(is(fm, "merMod"))

  vlist <- sapply(fm@cnms, length)
  
  pp <- fm@pp$copy()
  isCor <- checkCorr(fm)
  
  resp <- fm@resp$copy()
  np <- length(pp$theta)
  nf <- length(fixef(fm))
  if (!isGLMM(fm)) 
    np <- np + 1L
  n <- nrow(pp$V)
  if (isLMM(fm)) {
    ans <- function(Lc, thpars) {
      stopifnot(is.numeric(thpars), length(thpars) == np)
      ## TO FIX: when the correlations are present then not execute the next line
      #if(!isCor)
      #  pars[which(pars < 0)] <- 0
      sigma2 <- thpars[np]^2
      #thpars <- Vv_to_Cv(pars, n = vlist, s = sqrt(sigma2))
      .Call("lmer_Deviance", pp$ptr(), resp$ptr(), thpars[-np], PACKAGE = "lme4")      
      vcov_out <- sigma2 * tcrossprod(pp$RXi()) #chol2inv(pp$RX())
      #return(as.matrix(vcov_out))
      return(as.matrix(Lc %*% as.matrix(vcov_out) %*% t(Lc)))        
    }
  } 
  #attr(ans, "optimum") <- opt
  #attr(ans, "basedev") <- basedev
  #attr(ans, "thopt") <- pp$theta
  #attr(ans, "stderr") <- stdErr
  attr(ans, "isCor") <- isCor
  class(ans) <- "vcovJSStheta"
  ans
}





# vv <- as.data.frame(VarCorr(fm2.ML))  ## need ML estimates!
# 
# 
# ## name for consistency with profile object below
# vnames <- c(sprintf(".sig%02d", 1), ".sigma")
# pars <- setNames(vv[c(1,  2), "sdcor"], vnames)
# 
# vv[c(1,  2), "sdcor"]
# Dev(rho, rho$param$vec.matr)
# 
# pars <- vv[c(2,  1), "vcov"]
# Dev2(rho, pars)
# hh1 <- hessian(function(x) Dev2(rho,x), pars)
# vv2 <- 2 * solve(hh1)
