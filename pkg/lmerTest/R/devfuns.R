
## NOT USED NOW 
## calculated deviance based on variance-covariance parameters
# ## modified devfun2 from lme4: REML criterion was added
# devfun3 <- function (fm, useSc, signames, reml = TRUE) 
# {
#   stopifnot(is(fm, "merMod"))
#   #fm <- refitML(fm)
#   #basedev <- deviance(fm)
#   vlist <- sapply(fm@cnms, length)
#   #sig <- sigma(fm)
#   #stdErr <- unname(coef(summary(fm, ddf = "lme4"))[, 2])
#   pp <- fm@pp$copy()
#   isCor <- checkCorr(fm)
#   #if (useSc) {
#   #  opt <- Cv_to_Vv(pp$theta, n = vlist, s = sig)
#   #  names(opt) <- if (signames) {
#   #    c(sprintf(".sig%02d", seq(length(opt) - 1)), ".sigma")
#   #  }
#     #else {
#     #  c(tnames(fm, old = FALSE, prefix = c("sd", "cov")), 
#     #    "sigma2")
#     #}
#   #}
#   #else {
#   #  opt <- Cv_to_Sv(pp$theta, n = vlist)
#   #  names(opt) <- if (signames) {
#   #    sprintf(".sig%02d", seq_along(opt))
#   #  }
#   #  else {
#   #    tnames(fm, old = FALSE, prefix = c("sd", "cor"))
#   #  }
#   #}
#   #opt <- c(opt, fixef(fm))
#   resp <- fm@resp$copy()
#   np <- length(pp$theta)
#   nf <- length(fixef(fm))
#   if (!isGLMM(fm)) 
#     np <- np + 1L
#   n <- nrow(pp$V)
#   if (isLMM(fm)) {
#     ans <- function(pars) {
#       stopifnot(is.numeric(pars), length(pars) == np)
#       ## TO FIX: when the correlations are present then not execute the next line
#       if(!isCor)
#         pars[which(pars < 0)] <- 0
#       sigma2 <- pars[np]
#       thpars <- Vv_to_Cv(pars, n = vlist, s = sqrt(sigma2))
#       .Call("lmer_Deviance", pp$ptr(), resp$ptr(), thpars, PACKAGE = "lme4")
#       sigsq <- sigma2
#       dev <- pp$ldL2() + (resp$wrss() + pp$sqrL(1))/sigsq + n * 
#         log(2 * pi * sigsq)      
#       if(reml){
#         p <- ncol(pp$RX())
#         dev <- dev + 2*determinant(pp$RX())$modulus - p * log(2 * pi * sigsq)              
#       }
#       return(dev)     
#     }
#   }
#   #attr(ans, "optimum") <- opt
#   #attr(ans, "basedev") <- basedev
#   attr(ans, "thopt") <- pp$theta
#   attr(ans, "isCor") <- isCor
#   #attr(ans, "stderr") <- stdErr
#   class(ans) <- "devfun3"
#   #ans
# }


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

  attr(ans, "thopt") <- pp$theta
  attr(ans, "isCor") <- isCor
  class(ans) <- "devfun4"
  ans
}




## NOT USED NOW
## calc vcov of fixed effects based on var cor parameters of random ## effects
vcovJSS <- function(fm)
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
    ans <- function(Lc, pars) {
      stopifnot(is.numeric(pars), length(pars) == np)
      ## TO FIX: when the correlations are present then not execute the next line
      if(!isCor)
        pars[which(pars < 0)] <- 0
      sigma2 <- pars[np]
      thpars <- Vv_to_Cv(pars, n = vlist, s = sqrt(sigma2))
      .Call("lmer_Deviance", pp$ptr(), resp$ptr(), thpars, PACKAGE = "lme4")      
      vcov_out <- sigma2 * tcrossprod(pp$RXi()) #chol2inv(pp$RX())
      return(as.matrix(Lc %*% as.matrix(vcov_out) %*% t(Lc)))        
    }
  } 

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

      sigma2 <- thpars[np]^2

      .Call("lmer_Deviance", pp$ptr(), resp$ptr(), thpars[-np], PACKAGE = "lme4")      
      vcov_out <- sigma2 * tcrossprod(pp$RXi()) 

      return(as.matrix(Lc %*% as.matrix(vcov_out) %*% t(Lc)))        
    }
  } 

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
