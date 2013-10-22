##############################################################################
# performs  analysis of sensory data
##############################################################################
sensOldmixedFun <- function(attributes, Prod_effects, replication, individual, data, product_structure = 3, error_structure="No_Rep", alpha.random = 0.1, alpha.fixed = 0.05 )
{
  #product_structure=1  (default structure) : Analysis of main fixed effects
  #product_structure=2 : Main effects AND all 2-factor interactions. 
  #product_structure=3 : Full factorial model with ALL possible fixed effects
  #error_structure=
  
  attach(data)
  
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
  
  
  resultFULL <- vector("list",length(attributes))
  nbvar <- length(attributes)
  
  # create the initial model
  model <- createLMERmodel(structure=list(product_structure=product_structure, error_structure=error_structure), data, attributes[1], fixed = list(Product=Prod_effects, Consumer=NULL), random = random, FALSE)
  
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
  pvalue <- matrix(NA, n, nbvar)
  colnames(pvalue) <- attributes
  rownames(pvalue) <- c(fixedrand$fixedeffs, fixedrand$randeffs)  
  
  # Matrix that will contains the statistics
  #(F value for fixed effect and Chi2 for random effect)
  FChi <- matrix(0, n, nbvar)
  colnames(FChi) <- attributes
  rownames(FChi) <- rownames(pvalue)
  

  #### fill the results
  for(i in 1:length(attributes))
  {
        if(i > 1)
        {
          #m <- createLMERmodel(structure=list(product_structure=product_structure, error_structure=error_structure), data, attributes[i], fixed = list(Product=Prod_effects, Consumer=NULL), random = random, FALSE)
          #j <- attributes[i]
          #m <- refit(object=model, newresp = data[, j], rename.response = TRUE)
          #m <- eval(substitute(refit(object=model, newresp=resp, rename.response = TRUE), list(resp = data[,attributes[i]])))
          #m <- eval(substitute(refit(object=model, newresp=resp, rename.response = TRUE), list(resp=as.name(attributes[i]))))
          m <- createLMERmodel(structure=list(product_structure=product_structure, error_structure=error_structure), data, attributes[i], fixed = list(Product=Prod_effects, Consumer=NULL), random = random, FALSE)
          
        } 
        else 
           m <- model
    
    
        print(paste("Calculating for", attributes[i],"...", sep=" "))
        
        t <- step(m, reduce.fixed = isFixReduce, reduce.random = isRandReduce, alpha.random = alpha.random, alpha.fixed = alpha.fixed, lsmeans.calc=FALSE, difflsmeans.calc = FALSE)
    
    
    #fill pvalues for fixed effects
    calc <- t$anova.table
    pvalue[rownames(calc),i] <- calc[,6]
    
    #fill pvalues for random effects
    calcrand <- t$rand.table
    rownames(calcrand) <- substr(rownames(calcrand), 2, nchar(rownames(calcrand))-1)
    pvalue[rownames(calcrand),i] <- calcrand[,4]
    
    # fill F and Chi values
    FChi[rownames(calc),i] <- calc[,5]
    FChi[rownames(calcrand),i] <- calcrand[,1]
    
  }
  
  pvalue[is.na(pvalue)] <- 1
  
  
  detach(data)
  return(list(FChi = FChi, pvalue = pvalue))
  
}
