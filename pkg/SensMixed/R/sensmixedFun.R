
##############################################################################
# performs  analysis of sensory data
##############################################################################
sensmixedFun <- function(attributes = NULL, Prod_effects, replication = NULL, 
                         individual, data, product_structure = 3, 
                         error_structure = "No_Rep", 
                         MAM = FALSE, MAM_PER = FALSE, adjustedMAM = FALSE, 
                         alpha_conditionalMAM = 1, calc_post_hoc = TRUE, 
                         parallel=TRUE, 
                         reduce.random=TRUE, alpha.random = 0.1, 
                         alpha.fixed = 0.05, interact.symbol = interact.symbol)
  
{
  ## product_structure=1  (default structure) : Analysis of main fixed effects
  ## product_structure=2 : Main effects AND all 2-factor interactions. 
  ## product_structure=3 : Full factorial model with ALL possible fixed effects
  ## error_structure 
  if(is.null(attributes))
    attributes <- colnames(data)[!sapply(data, 
                                         function(x) "factor" %in% class(x))]
  
  if(class(attributes)=="integer")
    attributes <- colnames(data)[attributes]
  
  if(length(attributes) < 7)
    parallel <- FALSE
  
  
  
  if(MAM_PER)
  {
      return(runMAM(data, Prod_effects, individual, attributes, adjustedMAM, 
                    alpha_conditionalMAM))   
  }

  
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
  
  

  nbvar <- length(attributes)
  
  
  
  ## create the initial model
  suppressMessages(model.init <- createLMERmodel(structure = 
                                  list(product_structure=product_structure, 
                           error_structure = error_structure), 
                           data, attributes[1],
                           fixed = list(Product = Prod_effects, Consumer = NULL),
                           random = random, corr=FALSE, MAM, 
                           mult.scaling = FALSE, calc_post_hoc = calc_post_hoc))
  model <- if(MAM) model.init$model.anova else model.init
  if(calc_post_hoc)
    model.lsm <- if(MAM) model.init$model.lsmeans else model.init
  
 
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
  pvalueF <- matrix(NA, nfixe, nbvar)
  colnames(pvalueF) <- attributes
  pvalueChi <- matrix(NA, nrand, nbvar)
  colnames(pvalueChi) <- attributes
  
  ## TODO: change rownames for pvaues according to elimRand.R
  ##       should work for slopes as well
  rownames(pvalueF) <- fixedrand$fixedeffs
  pvalueF <- .renameScalingTerm(pvalueF, Prod_effects) 
  
  rownames(pvalueChi) <- fixedrand$randeffs

  Fval <- matrix(0, nfixe, nbvar)
  colnames(Fval) <- attributes
  rownames(Fval) <- rownames(pvalueF)
  Chi <- matrix(0, nrand, nbvar)
  colnames(Chi) <- attributes
  rownames(Chi) <- rownames(pvalueChi)
 
  ### using parallel
  if(parallel){
    respar <- tryCatch({
    funcNOMAM <- local({
      #data
      refit
      step
      model
      data
      reduce.random
      alpha.random
      alpha.fixed
      calc_post_hoc
      summary
     
      function(new.resp.private.sensmixed)
      {
        assign("new.resp.private.sensmixed", new.resp.private.sensmixed, 
               envir=environment(formula(model)))
        suppressMessages(m <- refit(object=model, 
                                    newresp=new.resp.private.sensmixed, 
                                    rename.response = TRUE))
        suppressMessages(s <- step(m, reduce.fixed = FALSE, 
                                   reduce.random = reduce.random, 
                                   alpha.random = alpha.random, 
                                   alpha.fixed = alpha.fixed, 
                                   lsmeans.calc=TRUE,
                                   difflsmeans.calc = TRUE))        
        
        if(calc_post_hoc){
#           ## calculate averaged squared d-prime  
#           sigma <- summary(s$model, "lme4")$sigma
#           rows <- sapply(rownames(s$diffs.lsmeans.table), 
#                          function(x) strsplit(x, " ")[[1]][1]) 
#           s$anova.table$dprimeav <- rep(1, nrow(s$anova.table))
#           for(eff in rownames(s$anova.table)){
#             dp <- s$diffs.lsmeans.table[which(rows %in% eff),1]/sigma
#             av.dp <- sqrt(sum(dp^2)/length(dp))
#             s$anova.table[eff, "dprimeav"] <- av.dp
          
#           }
          #s$anova.table <- .calcAvDprime(s$model, s$anova.table, 
         #                                 s$diffs.lsmeans.table, s$lsmeans.table)  
        }       
        
        s        
      }})
    funcMAM <- local({
      #step     
      #data
      #createLMERmodel
      
      function(attr)
      {
        model.init <- suppressMessages(createLMERmodel(
          structure = list(product_structure=product_structure,
                           error_structure=error_structure), 
          data = data, response = attr, fixed = list(Product = Prod_effects, 
                                                     Consumer = NULL),
          random = random, corr = FALSE, MAM = TRUE, 
          mult.scaling = FALSE, 
          calc_post_hoc = calc_post_hoc))
        
         model.an <- model.init$model.anova
         model.lsm <- model.init$model.lsmeans
         #return(model.an)
         
         st <- suppressMessages(step(as(model.an,"merModLmerTest"), 
                                     fixed.calc=FALSE))
        rand.table <- st$rand.table        
        
        if(reduce.random){
          anova.table <- suppressMessages(anova(as(st$model, "merModLmerTest"), 
                                                type = 1))
          if(length(which(anova.table[, "Pr(>F)"] == "NaN") > 0))
            anova.table <- suppressMessages(anova(as(st$model, "merModLmerTest"), 
                                                  type = 1, ddf="Kenward-Roger")) 
        }
        else{
          anova.table <- suppressMessages(anova(model.an, type = 1))
          if(length(which(anova.table[, "Pr(>F)"] == "NaN") > 0))
            anova.table <- suppressMessages(anova(as(model.an, "merModLmerTest"), 
                                                  type = 1, ddf="Kenward-Roger")) 
        }
        
        anova.table <- .renameScalingTerm(anova.table, Prod_effects) 
        
        
        if(length(Prod_effects) > 1)
          lsmeans.table <- lsmeans::.old.lsmeans(model.lsm, pairwise ~ prod)
        else 
          lsmeans.table <-  eval(substitute(lsmeans::.old.lsmeans(object = model.lsm, 
                                                             pairwise ~ prod), 
                                            list(prod=as.name(Prod_effects))))  
        return(list(anova.table = anova.table, rand.table = rand.table,
                    lsmeans.table = lsmeans.table)) 
        #attach(data)
        #print(paste("Calculating for", x,"...", sep=" "))
        #m <- eval(substitute(refit(object=model, newresp=resp, 
        #                              rename.response = TRUE), 
        #                               list(resp=as.name(x))))
        #return(anova(as(model, "merModLmerTest"), type=1))
        #         assign("new.resp", new.resp, envir=environment(formula(model)))        
        #         if(MAM){
        #           model.an <- refit(object=model, newresp=new.resp, 
        #                             rename.response = TRUE)
        #           data.an <- model.frame(model.an)
        #           fo <- paste(formula(model.an))
        #           #fo[2] <- "new.resp"
        #           #fo[3] <- paste(Prod_effects, collapse="*")
        #           #data.an$new.resp <- new.resp
        #           lm.pred <- lm(as.formula(paste("new.resp", "~", 
        #                                          paste(Prod_effects, collapse="*"), sep="")),
        #                         data=data.an)
        #           data.an$x <- scale(predict(lm.pred), scale=FALSE)
        #           m <- lmer(as.formula(paste("new.resp", fo[1], fo[3], sep="")), 
        #                     data=data.an)
        #           st <- step(m, fixed.calc=FALSE)
        #           #anova.table <- anova(st$model, type=1)
        #           #rand.table <- rand(m)$rand.table
        #           
        #           rand.table <- st$rand.table
        #           if(reduce.random)
        #             anova.table <- anova(as(st$model, "merModLmerTest"), type=1)
        #           else
        #             anova.table <- anova(m, type=1)
        #           #if(any(is.nan(anova.table[, "Pr(>F)"])))
        #           #  anova.table <- anova(m, ddf="Kenward-Roger", type=1)     
        #           rownames(anova.table)[unlist(lapply(rownames(anova.table), 
        #                                               function(y) grepl(":x", y)))] <- "Scaling"
        #           
        #           ## update model for lsmeans
        # #           m.lsm <- refit(object=model.lsm, newresp=new.resp, 
        # #                          rename.response = TRUE)        
        # #           fo.lsm <- paste(formula(model.lsm))
        # #           fo.lsm[2] <- "new.resp"
        # #           data.lsm <- model.frame(model.lsm)
        # #           data.lsm$new.resp <- new.resp 
        # #           m.lsm <- lmer(as.formula(paste(fo.lsm[2], fo.lsm[1], fo.lsm[3], sep="")), 
        # #                         data=data.lsm)
        # #           if(length(Prod_effects) > 1)
        # #             lsmeans.table <- lsmeans::lsmeans(m.lsm, pairwise ~ prod)
        # #           else 
        # #             lsmeans.table <-  eval(substitute(lsmeans::lsmeans(object=m.lsm, 
        # #                                                                pairwise ~ prod), 
        # #                                               list(prod=as.name(Prod_effects))))  
        #           return(list(anova.table=anova.table, rand.table=rand.table,
        #                       lsmeans.table=anova.table))        
        #         }
        #         m <- refit(object=model, newresp=new.resp, 
        #                    rename.response = TRUE)
        #         #anova(as(m,"merModLmerTest"), ddf="Kenward-Roger", type=1)    
        #         s <- step(m, reduce.fixed = FALSE, reduce.random = reduce.random, 
        #                   alpha.random = 0.1, alpha.fixed = 0.05, lsmeans.calc=FALSE,
        #                   difflsmeans.calc = FALSE) 
        #         #detach(data)
        #         return(s)
      }})
    ncpus <- detectCores()
    cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))             
    parallel::clusterExport(cl, varlist=c("refit", "data", "createLMERmodel",
                                          "createFormulaAllFixRand",
                                          "checkComb",
                                          "checkNumberInteract",
                                          "checkZeroCell",
                                          "nrandTerms",
                                          "fixedFormula",
                                          "substring.location", "step",
                                          ".renameScalingTerm",
                                          "isLMM",
                                          "fixef", "anova", 
                                          ".calcAvDprime",
                                          "getPureInter",
                                          ".calcPureDiffs",
                                          ".getIndTermsContained"
                                          ),
                                           #"model", "model.lsm", "data" ),
                            envir=environment())
    if(RNGkind()[1L] == "L'Ecuyer-CMRG")
      parallel::clusterSetRNGStream(cl)
    if(!MAM)
      res <- parallel::parLapply(cl, data[,attributes], funcNOMAM)
    else
      res <- parallel::parLapply(cl, attributes, funcMAM)
    parallel::stopCluster(cl) 
    res
    },  error = function(e) { NULL })
    if(is.null(respar)){
      message(" \n WARNING: error in parallel has occurred, an unparallized version is used instead \n")    
      return(sensmixedFun(attributes = attributes, Prod_effects = Prod_effects, 
                          replication = replication, individual = individual, 
                          data = data, product_structure = product_structure, 
                          error_structure = error_structure, MAM = MAM, 
                          MAM_PER = MAM_PER, adjustedMAM = adjustedMAM, 
                          alpha_conditionalMAM = alpha_conditionalMAM, 
                          calc_post_hoc = calc_post_hoc, parallel = FALSE, 
                          reduce.random = reduce.random, 
                          alpha.random = alpha.random, alpha.fixed = alpha.fixed, 
                          interact.symbol = interact.symbol))
    }
  }
  else{    
    if(!MAM){      
      res <- llply(data[,attributes], .stepAllAttrNoMAM, model, reduce.random,
                   alpha.random, alpha.fixed, calc_post_hoc, .progress="text")
    }
    else{
      res <- llply(attributes, .stepAllAttrMAM, 
                    product_structure = product_structure, 
                    error_structure = error_structure,
                    data = data, Prod_effects = Prod_effects, random = random,
                    reduce.random = reduce.random, alpha.random = alpha.random, 
                    alpha.fixed = alpha.fixed, 
                    mult.scaling = FALSE, 
                    calc_post_hoc = calc_post_hoc, .progress="text")    
    } 
  }
  
  
 
  
  #### fill the results
  finalres <- tryCatch({
    if(calc_post_hoc){
      post_hoc <- vector(mode = "list", length = length(res))      
      if(MAM)
        names(post_hoc) <- attributes
      else{
        names(post_hoc) <- names(res)
        ## create d prime output
        dprimeav <- Fval
      }
    }
      
    for(i in 1:length(attributes))
    {    
      res[[i]]$response <- attributes[i]
      
      #fill pvalues for fixed effects
      calc <- res[[i]]$anova.table
      ## if the reduced is lm model
      if("Residuals" %in% rownames(res[[i]]$anova.table))
        pvalueF[rownames(calc)[-nrow(calc)], i] <- calc[-nrow(calc), 5]
      else
        pvalueF[rownames(calc), i] <- calc[, "Pr(>F)"]
      
      #fill pvalues for random effects
      calcrand <- res[[i]]$rand.table    
      pvalueChi[rownames(calcrand),i] <- calcrand[, "p.value"]
      
      # fill F and Chi values
      if("Residuals" %in% rownames(res[[i]]$anova.table))
        Fval[rownames(calc)[-nrow(calc)],i] <- calc[-nrow(calc), 4]
      else
        Fval[rownames(calc),i] <- calc[,"F.value"]
      Chi[rownames(calcrand),i] <- calcrand[,"Chi.sq"] 
      
      ## fill differences of lsmeans
      if(calc_post_hoc){
        if(MAM)
          post_hoc[[i]] <- res[[i]]$lsmeans.table[[2]]
        else{
          post_hoc[[i]] <- res[[i]]$diffs.lsmeans.table
          dprimeav[rownames(calc), i] <- calc[, "dprimeav"]
        }
      }
    }
    
    pvalueF[is.na(pvalueF)] <- 1
    pvalueChi[is.na(pvalueChi)] <- 1  
      

    }, error = function(e) { NULL })
  if(is.null(finalres) && parallel){
    message(" \n WARNING: error in parallel has occurred: cannot call Kenward-Roger 
            in parallel. instead use unparallelized version \n")     
    return(sensmixedFun(attributes = attributes, Prod_effects = Prod_effects, 
                        replication = replication, individual = individual, 
                        data = data, product_structure = product_structure, 
                        error_structure = error_structure, MAM = MAM, 
                        MAM_PER = MAM_PER, adjustedMAM = adjustedMAM, 
                        alpha_conditionalMAM = alpha_conditionalMAM, 
                        calc_post_hoc = calc_post_hoc, parallel = FALSE, 
                        reduce.random = reduce.random, 
                        alpha.random = alpha.random, alpha.fixed = alpha.fixed, 
                        interact.symbol = interact.symbol))
  }
  ## change the output
  ## output for the random effects
  #tr_rand <- .changeOutput(Chi, pvalueChi, TRUE)
  
   
  
  
  if(MAM){   
    ind.scaling <- grepl("Scaling", rownames(pvalueF))
    pvalueScaling <- pvalueF[ind.scaling, , drop=FALSE]
    pvalueF <- pvalueF[!ind.scaling, , drop=FALSE]
    FScaling <- Fval[ind.scaling, , drop=FALSE]
    Fval <- Fval[!ind.scaling, , drop=FALSE]   
    #tr_fixed <- .changeOutput(Fval, pvalueF, FALSE)
    #tr_scaling <- .changeOutput(FScaling, pvalueScaling, FALSE)    
    if(calc_post_hoc)
      return(list(fixed = list(Fval = Fval, pvalueF = pvalueF), 
                random = list(Chi = Chi, pvalueChi = pvalueChi), 
                scaling = list(FScaling = FScaling, 
                              pvalueScaling = pvalueScaling), 
                post_hoc = post_hoc))
    else
      return(list(fixed = list(Fval = Fval, pvalueF = pvalueF), 
                  random = list(Chi = Chi, pvalueChi = pvalueChi), 
                  scaling = list(FScaling = FScaling, 
                                 pvalueScaling = pvalueScaling)))
  }
  #tr_fixed <- .changeOutput(Fval, pvalueF, FALSE)
  
  if(calc_post_hoc)    
      return(list(fixed = list(Fval = Fval, dprimeav = dprimeav, pvalueF = pvalueF), 
                  random = list(Chi = Chi, pvalueChi = pvalueChi), 
                  post_hoc = post_hoc, step_res = res))  
  return(list(fixed = list(Fval = Fval, pvalueF = pvalueF), 
              random = list(Chi = Chi, pvalueChi = pvalueChi), step_res = res))  
}
