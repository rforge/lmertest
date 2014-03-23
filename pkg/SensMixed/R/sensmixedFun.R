##############################################################################
# performs  analysis of sensory data
##############################################################################
sensmixedFun <- function(attributes, Prod_effects, replication, individual, data, 
                         product_structure = 3, error_structure="No_Rep",
                         MAM=FALSE, parallel=TRUE, alpha.random = 0.1, alpha.fixed = 0.05)
{
  ## product_structure=1  (default structure) : Analysis of main fixed effects
  ## product_structure=2 : Main effects AND all 2-factor interactions. 
  ## product_structure=3 : Full factorial model with ALL possible fixed effects
  ## error_structure 
  

  
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
  
  ## create the initial model
  model.init <- createLMERmodel(structure=list(product_structure=product_structure, 
                           error_structure=error_structure), 
                           data, attributes[1],
                           fixed = list(Product=Prod_effects, Consumer=NULL),
                           random = random, corr=FALSE, MAM)
  model <- if(MAM) model.init$model.anova else model.init
  model.lsm <- if(MAM) model.init$model.lsmeans else model.init
  ## use Per's function
#   if(MAM){
#     if(length(Prod_effects)>1)
#       ##model.frame()    
#     else
#       data[,c(random$individual, Prod_effects)]  
#   }
#   if(MAM){
#     an.model <- anova(model$model.anova, type=1)
#     if(length(Prod_effects) > 1)
#       lsm.model <- lsmeans::lsmeans(model$model.lsmeans, pairwise ~ prod)
#     else 
#       lsm.model <-  eval(substitute(lsmeans::lsmeans(object=model$model.lsmeans, 
#                                                      pairwise ~ prod), 
#                                     list(prod=as.name(Prod_effects))))   
#     return(list(an.model = an.model, lsm.model=lsm.model))
#   }
#   
 
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
  rownames(pvalueF)[unlist(lapply(rownames(pvalueF), function(x) grepl(":x",x)))] <- "Scaling"
  rownames(pvalueChi) <- fixedrand$randeffs

  Fval <- matrix(0, nfixe, nbvar)
  colnames(Fval) <- attributes
  rownames(Fval) <- rownames(pvalueF)
  Chi <- matrix(0, nrand, nbvar)
  colnames(Chi) <- attributes
  rownames(Chi) <- rownames(pvalueChi)
  #print(formula(model))
  #print(anova(model, type=1))
  ### using parallel
  if(parallel){
    respar <- tryCatch({
    func <- local({
      #data
      refit
      model
      model.lsm
      data
      
      function(new.resp)
      {
        #attach(data)
        #print(paste("Calculating for", x,"...", sep=" "))
        #m <- eval(substitute(refit(object=model, newresp=resp, 
        #                              rename.response = TRUE), 
        #                               list(resp=as.name(x))))
        #return(anova(as(model, "merModLmerTest"), type=1))
        assign("new.resp", new.resp, envir=environment(formula(model)))        
        if(MAM){
          data.an <- model.frame(model)
          fo <- paste(formula(model))
          #fo[2] <- "new.resp"
          #fo[3] <- paste(Prod_effects, collapse="*")
          data.an$new.resp <- new.resp
          lm.pred <- lm(as.formula(paste("new.resp", "~", 
                                         paste(Prod_effects, collapse="*"), sep="")),
                        data=data.an)
          data.an$x <- scale(predict(lm.pred), scale=FALSE)
          m <- lmer(as.formula(paste("new.resp", fo[1], fo[3], sep="")), 
                    data=data.an)
          #st <- step(m, lsmeans.calc=FALSE, difflsmeans.calc=FALSE, 
          #           reduce.fixed=FALSE)
          #anova.table <- anova(st$model, type=1)
          rand.table <- rand(m)$rand.table
          anova.table <- anova(m,  type=1)          
          if(any(is.nan(anova.table[, "Pr(>F)"])))
            anova.table <- anova(m, ddf="Kenward-Roger", type=1)     
          rownames(anova.table)[unlist(lapply(rownames(anova.table), function(y) grepl(":x", y)))] <- "Scaling"
          
          ## update model for lsmeans
#           m.lsm <- refit(object=model.lsm, newresp=new.resp, 
#                          rename.response = TRUE)        
#           fo.lsm <- paste(formula(model.lsm))
#           fo.lsm[2] <- "new.resp"
#           data.lsm <- model.frame(model.lsm)
#           data.lsm$new.resp <- new.resp 
#           m.lsm <- lmer(as.formula(paste(fo.lsm[2], fo.lsm[1], fo.lsm[3], sep="")), 
#                         data=data.lsm)
#           if(length(Prod_effects) > 1)
#             lsmeans.table <- lsmeans::lsmeans(m.lsm, pairwise ~ prod)
#           else 
#             lsmeans.table <-  eval(substitute(lsmeans::lsmeans(object=m.lsm, 
#                                                                pairwise ~ prod), 
#                                               list(prod=as.name(Prod_effects))))  
          return(list(anova.table=anova.table, rand.table=rand.table,
                      lsmeans.table=anova.table))        
        }
        m <- refit(object=model, newresp=new.resp, 
                   rename.response = TRUE)
        #anova(as(m,"merModLmerTest"), ddf="Kenward-Roger", type=1)    
        s <- step(m, reduce.fixed = FALSE, reduce.random = TRUE, 
                  alpha.random = 0.1, alpha.fixed = 0.05, lsmeans.calc=FALSE,
                  difflsmeans.calc = FALSE) 
        #detach(data)
        return(s)
      }})
    ncpus <- detectCores()
    cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))             
    parallel::clusterExport(cl, varlist=c( "refit", "model", "model.lsm", "data" ),
                            envir=environment())
    if(RNGkind()[1L] == "L'Ecuyer-CMRG")
      parallel::clusterSetRNGStream(cl)
    res <- parallel::parLapply(cl, data[,attributes], func)
    parallel::stopCluster(cl) 
    res
    },  error = function(e) { NULL })
    if(is.null(respar)){
      message(" \n ERROR in parallel has occurred, an unparallized version is used instead \n")     
      return(sensmixedFun(attributes, Prod_effects, replication, individual, data, 
                                      product_structure, error_structure,
                                      MAM, parallel=FALSE, alpha.random, alpha.fixed ))
    }
  }
  else{
    stepAllattr <- function(new.resp){    
      assign("new.resp", new.resp, envir=environment(formula(model)))
      #assign("new.resp", new.resp, envir=environment(formula(model.init$model.lsmeans)))      
      if(MAM){
        data.an <- model.frame(model)
        fo <- paste(formula(model))
        #fo[2] <- "new.resp"
        #fo[3] <- paste(Prod_effects, collapse="*")
        data.an$new.resp <- new.resp
        lm.pred <- lm(as.formula(paste("new.resp", "~", 
                                       paste(Prod_effects, collapse="*"), sep="")),
                      data=data.an)
        data.an$x <- scale(predict(lm.pred), scale=FALSE)
        m <- lmer(as.formula(paste("new.resp", fo[1], fo[3], sep="")), 
                      data=data.an)
        
       # st <- step(m, lsmeans.calc=FALSE, difflsmeans.calc=FALSE, 
       #                  reduce.fixed=FALSE)
        rand.table <- rand(m)$rand.table
        #fo <- paste(formula(m))
        #data.an <- model.frame(m)
        #data.an$Colourbalance <- TVbo$Colourbalance 
        #m <- lmer(as.formula(paste("Colourbalance", fo[1], fo[3], sep="")), 
        #                                       data=data.an)
        anova.table <- anova(m, type=1)
        if(any(is.nan(anova.table[, "Pr(>F)"])))
          anova.table <- anova(m,  ddf="Kenward-Roger", type=1)
        rownames(anova.table)[unlist(lapply(rownames(anova.table), function(y) grepl(":x", y)))] <- "Scaling"
        
        ## update model for lsmeans
        m.lsm <- refit(object=model.lsm, newresp=new.resp, 
                       rename.response = TRUE)        
        fo.lsm <- paste(formula(model.lsm))
        fo.lsm[2] <- "new.resp"
        data.lsm <- model.frame(model.lsm)
        data.lsm$new.resp <- new.resp 
        data.lsm$x <- data.an$x 
        m.lsm <- lmer(as.formula(paste(fo.lsm[2], fo.lsm[1], fo.lsm[3], sep="")), 
                      data=data.lsm)
         if(length(Prod_effects) > 1)
            lsmeans.table <- lsmeans::lsmeans( m.lsm, pairwise ~ prod)
         else 
            lsmeans.table <-  eval(substitute(lsmeans::lsmeans(object=m.lsm, 
                                                           pairwise ~ prod), 
                                          list(prod=as.name(Prod_effects))))  
        return(list(anova.table=anova.table, rand.table=rand.table,lsmeans.table=lsmeans.table))        
      }
      m <- refit(object=model, newresp=new.resp, 
                 rename.response = TRUE) 
      s <- step(m, reduce.fixed = FALSE, reduce.random = TRUE, 
                alpha.random = 0.1, alpha.fixed = 0.05, lsmeans.calc=FALSE,
                difflsmeans.calc = FALSE) 
      s
    }
    res <- lapply(data[,attributes], stepAllattr)
    #res <- lapply(attributes, stepAllattr)
  }
  
  

  
  #### fill the results
  for(i in 1:length(attributes))
  {    
    #fill pvalues for fixed effects
    calc <- res[[i]]$anova.table
    ## if the reduced is lm model
    if("Residuals" %in% rownames(res[[i]]$anova.table))
      pvalueF[rownames(calc)[-nrow(calc)], i] <- calc[-nrow(calc), 5]
    else
      pvalueF[rownames(calc), i] <- calc[, 6]
    
    #fill pvalues for random effects
    calcrand <- res[[i]]$rand.table    
    pvalueChi[rownames(calcrand),i] <- calcrand[, "p.value"]
    
    # fill F and Chi values
    if("Residuals" %in% rownames(res[[i]]$anova.table))
      Fval[rownames(calc)[-nrow(calc)],i] <- calc[-nrow(calc), 4]
    else
      Fval[rownames(calc),i] <- calc[,5]
    Chi[rownames(calcrand),i] <- calcrand[,"Chi.sq"]    
  }
  
  pvalueF[is.na(pvalueF)] <- 1
  pvalueChi[is.na(pvalueChi)] <- 1
  

  return(list(fixed = list(Fval=Fval, pvalueF=pvalueF), random = list(Chi=Chi, pvalueChi=pvalueChi)))
  
}
