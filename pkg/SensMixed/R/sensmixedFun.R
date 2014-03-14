##############################################################################
# performs  analysis of sensory data
##############################################################################
sensmixedFun <- function(attributes, Prod_effects, replication, individual, data, 
                         product_structure = 3, error_structure="No_Rep",
                         MAM=FALSE, alpha.random = 0.1, alpha.fixed = 0.05 )
{
  ##product_structure=1  (default structure) : Analysis of main fixed effects
  ###product_structure=2 : Main effects AND all 2-factor interactions. 
  ###product_structure=3 : Full factorial model with ALL possible fixed effects
  ##error_structure 
  

  
  if(!is.null(replication) && error_structure!="No-Rep")
    random <- list(individual=individual, replication=replication)
  else
    random <- individual
  isRandReduce <- TRUE
  isFixReduce <- FALSE
  isLsmeans <- TRUE
  
  
  
  #check if there are correlations between intercepts and slopes
  checkCorr <- function(model)
  {
    corr.intsl <- FALSE
    lnST <- length(getME(model, "ST"))
    for(i in 1:lnST)
    {    
      if(nrow(getME(model, "ST")[[i]])>1)
        corr.intsl <- TRUE
    } 
    return(corr.intsl) 
  }
  
  
  #resultFULL <- vector("list",length(attributes))
  nbvar <- length(attributes)
  
  # create the initial model
  model <- createLMERmodel(structure=list(product_structure=product_structure, 
                           error_structure=error_structure), 
                           data, attributes[1],
                           fixed = list(Product=Prod_effects, Consumer=NULL),
                           random = random, FALSE)
  
  if(checkCorr(model))
    isRandReduce <- FALSE
  
  #number of effects in the model
  fixedrand <- .fixedrand(model)
  # number of fixed effects
  nfixe <- length(fixedrand$fixedeffs)

  # number of random effects
  nrand <- length(fixedrand$randeffs)
  
  n <- nfixe + nrand
  
  # Matrix that will contain the p values
  #pvalue <- matrix(NA, n, nbvar)
  #colnames(pvalue) <- attributes
  pvalueF <- matrix(NA, nfixe, nbvar)
  colnames(pvalueF) <- attributes
  pvalueChi <- matrix(NA, nrand, nbvar)
  colnames(pvalueChi) <- attributes
  
  ## TODO: change rownames for pvaues according to elimRand.R
  ##       should work for slopes as well
  #rownames(pvalue) <- c(fixedrand$fixedeffs,  substr(fixedrand$randeffs, 5, nchar(fixedrand$randeffs)))  
  rownames(pvalueF) <- fixedrand$fixedeffs
  rownames(pvalueChi) <- fixedrand$randeffs

  Fval <- matrix(0, nfixe, nbvar)
  colnames(Fval) <- attributes
  rownames(Fval) <- rownames(pvalueF)
  Chi <- matrix(0, nrand, nbvar)
  colnames(Chi) <- attributes
  rownames(Chi) <- rownames(pvalueChi)
  
  ### using parallel
  func <- local({
    #data
    refit 
    model
    data
    
    function(x)
    {
      #attach(data)
      #print(paste("Calculating for", x,"...", sep=" "))
      #m <- eval(substitute(refit(object=model, newresp=resp, 
      #                              rename.response = TRUE), 
      #                               list(resp=as.name(x))))
      assign("x", x, envir=environment(formula(model)))
      m <- refit(object=model, newresp=x, 
                  rename.response = TRUE)
      s <- step(m, reduce.fixed = FALSE, reduce.random = TRUE, 
                alpha.random = 0.1, alpha.fixed = 0.05, lsmeans.calc=FALSE,
                difflsmeans.calc = FALSE) 
      #detach(data)
      return(s)
    }})
  
  
  ncpus <- detectCores()
  cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))             
  parallel::clusterExport(cl, varlist=c( "refit", "model", "data" ),
                          envir=environment())
  if(RNGkind()[1L] == "L'Ecuyer-CMRG")
    parallel::clusterSetRNGStream(cl)
  res <- parallel::parLapply(cl, data[,attributes], func)
  parallel::stopCluster(cl) 
  

  
  #### fill the results
  for(i in 1:length(attributes))
  {    
    #fill pvalues for fixed effects
    calc <- res[[i]]$anova.table
    pvalueF[rownames(calc),i] <- calc[,6]
    
    #fill pvalues for random effects
    calcrand <- res[[i]]$rand.table    
    pvalueChi[rownames(calcrand),i] <- calcrand[,4]
    
    # fill F and Chi values
    Fval[rownames(calc),i] <- calc[,5]
    Chi[rownames(calcrand),i] <- calcrand[,1]    
  }
  
  pvalueF[is.na(pvalueF)] <- 1
  pvalueChi[is.na(pvalueChi)] <- 1
  

  return(list(fixed = list(Fval=Fval, pvalueF=pvalueF), random = list(Chi=Chi, pvalueChi=pvalueChi)))
  
}
