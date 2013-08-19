totalAnovaRandLsmeans <- function(model, ddf="Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.05, reduce.fixed = TRUE, reduce.random = TRUE, lsmeans.calc = TRUE, difflsmeans.calc=TRUE,  isTotal=FALSE, isAnova=FALSE, isRand=FALSE, isLSMEANS=FALSE, isDiffLSMEANS=FALSE, isTtest=FALSE, test.effs=NULL, method.grad="simple")
{
  #errors in specifying the parameters
  if(!isRand && !(ddf %in% c("Satterthwaite","Kenward-Roger")))
  {
    print ('Error: parameter ddf is wrongly specified')
    stop()
  }
  if(!isRand && !(type %in% c(1,3)))
  {
    print ('Error: parameter type is wrongly specified')
    stop()
  } 
  
  data <- model.frame(model)#summary(model,"lme4")@frame 
  
  
  #update contrasts
  mm <- model.matrix(model)
  l <- attr(mm,"contrasts")
  contr <- l
  ### change contrasts for F tests calculations
  #list of contrasts for factors
  if( isAnova || isTotal )
  {    
    if( length(which(unlist(contr)!="contr.SAS")) > 0 )
    {
      names.facs <- names(contr)
      l <- as.list(rep("contr.SAS",length(names.facs)))
      names(l) <- names(contr)
      #warning(" \nmodel has been refitted with contrasts=contr.SAS \n")
      #model<-update(model,.~., data=data, contrasts=l)
      model <- updateModel(model, .~., getME(model, "is_REML"), l)     
    }    
  }
  else
  {
    #update model to mer class
    model <- updateModel(model, .~., getME(model, "is_REML"), l)
  }
  
  
  
  
  #not to show the warnings  
  #options(warn=-1) 
  
  result <- NULL
  anova.table <- NULL
  
  result$response <- rownames(attr(terms(model),"factors"))[1]#names(attr(terms(model),"dataClasses")[1])
  
  #model<-lmer(formula=formula, data=data)
  
  #update model
  # change unordered contrasts to contr.SAS
  # change REML to TRUE
  #options(contrasts=c(unordered="contr.SAS", ordered="contr.poly"))
  
  #model<-update(model, REML=TRUE)
  if( isRand || isTotal || (ddf=="Kenward-Roger" && (isTotal || isAnova)) )
    model<-
    if (getME(model, "is_REML") == 1)
    {
      model
    }
  else
  {
    warning("\n model has been refitted with REML=TRUE \n")
    updateModel(model, .~., reml=TRUE, l)
  }
  
  mf.final <- update.formula(formula(model),formula(model)) 
  
  
  #update data 
  #eliminate missing values
  ### old code ############
  #if(ncol(get_all_vars(mf.final, data))>ncol(data))
  #{
  #   #mframe<-model.frame(mf.final, data=data, na.action=na.pass)
  #   #data$response<-mframe[,1]
  #   data$response<-data[,1]
  #   data<-data[,-1]
  #   fm<-paste(mf.final)
  #   fm[2]<-"response"
  #   mf.final<- as.formula(paste(fm[2],fm[1],fm[3], sep=""))
  #   mf.final<-update.formula(mf.final,mf.final)
  #}
  
  #data<-get_all_vars(mf.final, data)
  ###############################################################
  
  data <- data[complete.cases(data),]
  
  
  #if(!is.null(l))
  #   model<-eval(substitute(lmer(mf.final, data=data, REML=model@dims[["REML"]], contrasts=l),list(mf.final=mf.final)))
  #else
  #   model<-eval(substitute(lmer(mf.final, data=data, REML=model@dims[["REML"]]),list(mf.final=mf.final)))
  #if(!is.null(l))
  #   model<-update(model, mf.final, data=data, REML=model@dims[["REML"]], contrasts=l)
  #else
  #   model<-update(model, mf.final, data=data, REML=model@dims[["REML"]])
  
  model <- updateModel(model, mf.final, getME(model, "is_REML"), l)
  
  # check if there are no fixed effects
  #if(length(attr(delete.response(terms(model)),"term.labels"))==0)
  # {
  #  fm<-getFormula(model, withRand=TRUE)
  #  fm<-as.formula(paste(fm[2],fm[1],paste("1+", fm[3], sep=""), sep=""))
  #  model<-eval(substitute(lmer(mf.final, data=data),list(mf.final=mf.final)))
  #}
  
  
  
  # save the call of the model              
  result$call <- model@call
  
  #check if there are correlations between intercept and slope
  result$corr.intsl <- checkCorr(model)  
  
  
  #Reduction of random effects
  #if(isRand || isTotal )
  #{     
  # perform reduction of random effects
  #if(!isRandReduce)
  #  result.rand<-elimRandLmer(model, data, 1)
  #else
  #  result.rand<-elimRandLmer(model, data, alpha.rand)
  
  if(isRand || isTotal)
  {
    result.rand <- elimZeroVarOrCorr(model, data, l)
    model <- result.rand$model      
  }    
  
  
  #save results for fixed effects for model with only fixed effects
  if(class(model) == "lm" | class(model) == "gls")
  {
    result <- saveResultsFixModel(result, model)
    result$rand.table=NULL
    return(result)
  }
  
  
  
  #if(!isRandReduce)
  #  result.rand <- elimRandEffs(model, data, 1)
  
  #analysis of the random part  
  if(isRand || isTotal)
  {
    if(isRand)
      reduce.random <- FALSE
    #if(reduce.random)
    result.rand <- elimRandEffs(model, data, alpha.random, reduce.random, l)  
    #if(!reduce.random && isRand)
    #  result.rand <- elimRandEffs(model, data, 1)
    
    model <- result.rand$model
    #convert rand table to data frame
    rt <- as.data.frame(result.rand$TAB.rand)
    rt$Chi.DF <- as.integer(rt$Chi.DF)
    if(!is.null(rt$elim.num))
      rt$elim.num <- as.integer(rt$elim.num)
    
    result$rand.table <- rt
    if(isRand)
      return(result)
  }
  
  
  
  
  #}
  
  
  
  #save results for fixed effects for model with only fixed effects
  if(class(model) == "lm" | class(model) == "gls")
    return(saveResultsFixModel(result, model))
  
  
  #perform reduction of fixed effects for model with mixed effects
  stop = FALSE
  is.first.anova <- TRUE
  is.first.sign <- TRUE  
  
  
  
  while(!stop)
  {      
    
    # if there are no fixed terms
    if(nrow(anova(model, ddf="lme4"))==0)
    {
      if(is.null(anova.table))
      {
        
        if(isLSMEANS || isDiffLSMEANS)
        {
          lsmeans.summ <-  matrix(ncol=7,nrow=0)
          colnames(lsmeans.summ) <- c("Estimate","Standard Error", "DF", "t-value", "Lower CI", "Upper CI", "p-value")
          lsmeans.summ <- as.data.frame(lsmeans.summ)
          if(isLSMEANS)
            result$lsmeans.table <- lsmeans.summ
          if(isDiffLSMEANS)
            result$diffs.lsmeans.table <- lsmeans.summ
          return(result)
        }
        if(isTtest)
        {
          # save lmer outcome in rho environmental variable
          rho <- rhoInit(model)     
          
          # calculate asymptotic covariance matrix A
          h  <-  hessian(function(x) Dev(rho,x), rho$param$vec.matr)
          rho$A <- 2*solve(h)
          #rho$A <- 2*ginv(h)
          
          tsummary <- calculateTtest(rho, diag(rep(1,length(rho$fixEffs))),length(rho$fixEffs), method.grad)
          result$ttest <- list(df=tsummary[,"df"], tvalue=tsummary[,"t value"], tpvalue=tsummary[,"p-value"])
        }
        result$model <- model
        result$anova.table <- anova(model, ddf="lme4")
        return(result)        
        
      }          
      break
    }        
    
    
    # save lmer outcome in rho environmental variable
    rho <- rhoInit(model)     
    
    # calculate asymptotic covariance matrix A
    h  <-  hessian(function(x) Dev(rho,x), rho$param$vec.matr)
    rho$A  <-  2*solve(h)
    #rho$A  <-  2*ginv(h)
    
    #Check if A is positive-definite
    isposA <- all(eigen(rho$A)$values>0)      
    if(!isposA)
    {
      print("Asymptotic covariance matrix A is not positive!")
      #previous code: return ERROR
      #print("ERROR: asymptotic covariance matrix A is not positive!")
      #result$model <- model
      #TABs <- emptyAnovaLsmeansTAB()
      #result$anova.table <- TABs$TAB.fixed
      #result$lsmeans.table <- TABs$TAB.lsmeans
      #result$diffs.lsmeans.table <- TABs$TAB.lsmeans
      #return(result)
    }
    
    
    #calculate ttest and p-values for summary
    if(isTtest)
    {
      tsummary <- calculateTtest(rho, diag(rep(1,length(rho$fixEffs))), length(rho$fixEffs), method.grad)
      result$ttest <- list(df=tsummary[,"df"], tvalue=tsummary[,"t value"], tpvalue=tsummary[,"p-value"])
      return(result)
    }
    
    
    #calculate lsmeans of differences of LSMEANS of the final model
    if(isLSMEANS || isDiffLSMEANS)
    {
      if(isLSMEANS)
      {
        lsmeans.tab <- calcLSMEANS(model, data, rho, alpha.fixed, test.effs=test.effs, method.grad=method.grad, lsmeansORdiff=TRUE, l)
        result$lsmeans.table <- lsmeans.tab$summ.data
        result$diffs.lsmeans.table <- NULL
      }
      if(isDiffLSMEANS)
      {
        lsmeans.tab <- calcLSMEANS(model, data, rho, alpha.fixed, test.effs=test.effs, method.grad=method.grad, lsmeansORdiff=FALSE, l)
        result$diffs.lsmeans.table <- lsmeans.tab$summ.data
        result$lsmeans.table <- NULL
      }
      return(result)
    }
    
    
    # Calculate  F-test with Satterthwaite's approximation
    # create X design matrix for fixed effects
    X.design.list <- createDesignMat(model,data)
    X.design <- X.design.list$X.design
    names.design.withLevels <- X.design.list$names.design.withLevels
    
    # define full coefficients for model
    #coefs.real<-fullCoefs(model, data)
    #coefs.real<-getDummyCoefs(model, data)
    
    
    #save full coefficients in rho
    #rho$s.test<-coefs.real
    nums.dummy.coefs <- getNumsDummyCoefs(model, data, l)
    rho$nums.zeroCoefs <- nums.dummy.coefs$nums.zeroCoefs
    rho$nums.Coefs <- nums.dummy.coefs$nums.Coefs
    fullCoefs <- rep(0, ncol(X.design))
    fullCoefs[rho$nums.Coefs] <- rho$fixEffs
    
    
    #define the terms that are to be tested
    test.terms <- attr(terms(model),"term.labels")
    
    #initialize anova table
    if(is.first.anova)
    {
      anova.table <- initAnovaTable(model, reduce.fixed)
      is.first.anova <- FALSE
      elim.num <- 1
    }
    
    # calculate general set matrix for type 3 hypothesis
    if( type==3 )
      L <- calcGeneralSetForHypothesis(X.design, rho)  
    
    
    # calculate type 1 hypothesis matrices for each term
    if( type==1 )
    {
      #X <- model.matrix(model)
      X <- X.design
      p <- ncol(X)
      XtX <- crossprod(X)
      U <- doolittle(XtX)$U
      d <- diag(U)
      for(i in 1:nrow(U))
        if(d[i] > 0) U[i, ] <- U[i, ] / d[i]

    }
 
    #calculate ss F value ddf and p value for each term 
    resultFpvalueSS <- lapply(test.terms, calcFpvalueMAIN, L=L, X.design=X.design, fullCoefs=fullCoefs, model=model, rho=rho, ddf=ddf, type=type)  
    #fill anova table
    anova.table <- fillAnovaTable(resultFpvalueSS,  anova.table)
    
   
    if(!reduce.fixed)             
      break     
    else
    {
      resNSelim <- elimNSFixedTerm(model, anova.table, data, alpha.fixed, elim.num, l)
      if(is.null(resNSelim))
        break
      else
      {
        model <- resNSelim$model
        mf.final <- update.formula(formula(model),formula(model))
        #model <- eval(substitute(lmer(mf.final, data=data),list(mf.final=mf.final)))
        
        model <- updateModel(model, mf.final, getME(model, "is_REML"), l)
        
        anova.table <- resNSelim$anova.table
        elim.num <- elim.num+1
      }        
    }      
    
  }
  
  
  #convert anova table to data frame
  anova.table <- as.data.frame(anova.table)
  anova.table$NumDF <- as.integer(anova.table$NumDF)
  if(!is.null(anova.table$elim.num))
    anova.table$elim.num <- as.integer(anova.table$elim.num)
  
  if(isTotal || isAnova)
  {
    result$anova.table <- anova.table
    if(isAnova)
      return(result)
  }
  
  #if in step function least squares means of diffs of LSMEANS are required
  if(lsmeans.calc)
  {
    lsmeans.tab <- calcLSMEANS(model, data, rho, alpha.fixed, test.effs=test.effs, method.grad=method.grad, lsmeansORdiff=TRUE, l)
    result$lsmeans.table <- lsmeans.tab$summ.data
  }
  else
  {
    result$lsmeans.table <- NULL
  }
  if(difflsmeans.calc)
  {
    lsmeans.tab <- calcLSMEANS(model, data, rho, alpha.fixed, test.effs=test.effs, method.grad=method.grad, lsmeansORdiff=FALSE, l)
    result$diffs.lsmeans.table <- lsmeans.tab$summ.data
  }
  else
  {
    result$diffs.lsmeans.table <- NULL
  }
  
  tsummary <- calculateTtest(rho, diag(rep(1,length(rho$fixEffs))), length(rho$fixEffs), method.grad)
  result$ttest <- list(df=tsummary[,"df"], tvalue=tsummary[,"t value"], tpvalue=tsummary[,"p-value"])
  
  
  #update final model
  mf.final <- update.formula(formula(model),formula(model))
  #if(!is.null(l))
  #  model<-eval(substitute(lmer(mf.final, data=data, REML=model@dims[["REML"]], contrasts=l),list(mf.final=mf.final)))
  #else
  #  model<-eval(substitute(lmer(mf.final, data=data, REML=model@dims[["REML"]]),list(mf.final=mf.final)))
  model <- updateModel(model, mf.final, getME(model, "is_REML"), contr)
  
  #save model
  model <- as(model,"merLmerTest")
  model@t.pval <- result$ttest$tpvalue
  result$model <- model
  return(result)
}


# generic functions for total analysis on mixed models
#totalAnalysis <- function(model, data,  alpha.rand = 0.05, alpha.fix = 0.05, isFixReduce = FALSE, isRandReduce = FALSE, test.effs=NULL, plot=FALSE, ...) UseMethod("totalAnalysis")


step <- function(model, ddf="Satterthwaite", type=3, alpha.random = 0.1, alpha.fixed = 0.05, reduce.fixed = TRUE, reduce.random = TRUE, lsmeans.calc=TRUE, difflsmeans.calc=TRUE, test.effs=NULL, method.grad="simple",...)
{  
  result <- totalAnovaRandLsmeans(model=model, ddf=ddf , type=type,  alpha.random=alpha.random, alpha.fixed=alpha.fixed, reduce.fixed=reduce.fixed, reduce.random=reduce.random, lsmeans.calc=lsmeans.calc, difflsmeans.calc=difflsmeans.calc, isTotal=TRUE, isTtest=FALSE, test.effs=test.effs, method.grad=method.grad)
  class(result) <- "step"
  result
}

### UNUSED function
#totalAnalysis.formula <- function(formula, data, ...)
#{
#  model <- lmer(formula=formula, data=data)
#  resAnalysis<-totalAnalysis.default(model, data, ...)
#  resAnalysis$call<-match.call()
#  resAnalysis
#}

print.step <- function(x, ...)
{
  #cat("Call:\n")
  #print(x$call)
  
  if(!is.null(x$rand.table))
  {
    cat("\nRandom effects:\n")  
    printCoefmat(x$rand.table, digits=3 , dig.tst=1  ,tst.ind=c(which(colnames(x$rand.table)=="Chi.DF"),which(colnames(x$rand.table)=="elim.num")), P.values=TRUE, has.Pvalue=TRUE)
  }
  
  if(nrow(x$anova.table) != 0)
  {
    if(class(x$model) == "lm" | class(x$model) == "gls")
    {
      cat("\nFixed effects:\n")
      print(x$anova.table)
      cat("\nLeast squares means:\n")
      print(x$lsmeans.table) 
      cat("\nFinal model:\n")
      print(x$model)
      return()
    }
    else
    {
      cat("\nFixed effects:\n")
      printCoefmat(x$anova.table, dig.tst=1, tst.ind=c(1,2), cs.ind=3, digits=3 ,P.values = TRUE, has.Pvalue=TRUE)
      if(!is.null(x$lsmeans.table))
      {
        cat("\nLeast squares means:\n")
        printCoefmat(x$lsmeans.table, dig.tst=1  ,tst.ind=c(1:(which(colnames(x$lsmeans.table)=="Estimate")-1),which(colnames(x$lsmeans.table)=="DF")), digits=3 ,P.values = TRUE, has.Pvalue=TRUE)
      }
      if(!is.null(x$diffs.lsmeans.table))
      {
        cat("\n Differences of LSMEANS:\n")
        printCoefmat(x$diffs.lsmeans.table, dig.tst=1  ,tst.ind=c(1:(which(colnames(x$diffs.lsmeans.table)=="Estimate")-1),which(colnames(x$diffs.lsmeans.table)=="DF")), digits=3 ,P.values = TRUE, has.Pvalue=TRUE)
      }
      
    }    
  }
  else
    print(x$anova.table)
  cat("\nFinal model:\n")
  print(x$model@call)  
}


plot.step <- function(x, ...)
{
  if(!is.null(x$lsmeans.table) && nrow(x$lsmeans.table)>0)
    plotLSMEANS(x$lsmeans.table, x$response, "LSMEANS")     
  if(!is.null(x$diffs.lsmeans.table) && nrow(x$diffs.lsmeans.table)>0)
    plotLSMEANS(x$diffs.lsmeans.table, x$response, "DIFF of LSMEANS")
}

lmer <-
  function(formula, data, family = NULL, REML = TRUE,
           control = list(), start = NULL, verbose = FALSE, doFit = TRUE,
           subset, weights, na.action, offset, contrasts = NULL,
           model = TRUE, x = TRUE, ...)
  {
    mc <- match.call()
    mc[[1]] <- quote(lme4::lmer)
    model <- eval.parent(mc)
    model <- as(model,"merLmerTest")
    #tryCatch(  { result = glm( y~x , family = binomial( link = "logit" ) ) } , error = function(e) { print("test") } )
    #t.pval <- tryCatch( {totalAnovaRandLsmeans(model=model, ddf="Satterthwaite", isTtest=TRUE)$ttest$tpvalue}, error = function(e) { NULL })
    #if(!is.null(t.pval))
    #{
    #  model@t.pval  <- t.pval
    #}
    #else
    #{
    #  model <- as(model,"lmerMod")
    #}
    return(model)
  }


setMethod("anova", signature(object="merLmerTest"),
          function(object,...)  
          {
            mCall <- match.call(expand.dots = TRUE)
            dots <- list(...)
            modp <- if (length(dots))
              sapply(dots, is, "merLmerTest") | sapply(dots, is, "mer") | sapply(dots, is, "lm") else logical(0)
            if (any(modp)) {
              return(callNextMethod())
            }
            else
            {
              cnm <- callNextMethod()
              #if(!is.null(ddf)&& ddf=="lme4") 
              #  return(cnm)
              if(("ddf" %in% names(dots)) && dots$"ddf"=="lme4") 
                return(cnm)
              if("ddf" %in% names(dots))
                ddf <- dots$"ddf"
              else
                ddf <- "Satterthwaite"  
              if("method.grad" %in% names(dots))
                method.grad <- dots$"method.grad"
              else
                method.grad <- "simple"
              if("type" %in% names(dots))
                type <- dots$"type"
              else
                type <- 3
{
  table <- cnm
  an.table <- tryCatch({totalAnovaRandLsmeans(model=object, ddf=ddf, type=type, isAnova=TRUE, reduce.random=FALSE, reduce.fixed=FALSE, method.grad=method.grad)$anova.table}, error = function(e) { NULL })
  if(!is.null(an.table))
  {
    rnames <- rownames(table)
    if(nrow(an.table)>0)
    {
      table <- as.data.frame(cbind(table$Df, an.table$"Sum Sq", an.table$"Mean Sq", an.table[,"F.value"], an.table[,"DenDF"], an.table[,"Pr(>F)"]))
      colnames(table) <- c("Df", "Sum Sq", "Mean Sq", "F value", "Denom", "Pr(>F)")
      dimnames(table) <- list(rnames,
                              c("Df", "Sum Sq", "Mean Sq", "F value", "Denom", "Pr(>F)"))
    }
    else
      table <- an.table
    attr(table, "heading") <- paste("Analysis of Variance Table with ", ddf, " approximation for degrees of freedom")
  }
  
  class(table) <- c("anova", "data.frame")
  return(table)
}  
            }
            
          })
setMethod("summary", signature(object = "merLmerTest"),
    function(object, ddf="Satterthwaite", ...)
    {
      cl <- callNextMethod()
      if(!is.null(ddf) && ddf=="lme4") return(cl)
      else
      {
         #coefs.satt <- cbind(cl@coefs,totalAnovaRandLsmeans(model=object, ddf=ddf, isTtest=TRUE)$ttest$tpvalue) 
         t.pval <- tryCatch( {totalAnovaRandLsmeans(model=object, ddf="Satterthwaite", isTtest=TRUE)$ttest$tpvalue}, error = function(e) { NULL })
         coefs.satt <- cbind(cl$coefficients, t.pval) 
         cl$coefficients <- coefs.satt
         
#            coefs.satt <- cbind("df"= result$ttest$df, "p value"= result$ttest$tpvalue)
         colnames(cl$coefficients)[4] <- "Pr(>|t|)"
      }
      return(cl)
      #return(as(cl,"summary.merLmerTest"))
    }
) 



#randTAB.default<-function(model, data, ...)
rand <- function(model, ...)
{
  result <- totalAnovaRandLsmeans(model=model, isRand=TRUE, reduce.random=FALSE)  
  res <- list(rand.table=result$rand.table, isCorr = result$corr.intsl)
  class(res) <- "rand"
  res
}

print.rand <- function(x, ...)
{
  
  cat("Analysis of Random effects Table:\n")
  #if(x$isCorr)
  #  print(x$rand.table)
  #else
  if(!is.null(x$rand.table))
    printCoefmat(x$rand.table, digits=3 , dig.tst=1  ,tst.ind=c(which(colnames(x$rand.table)=="Chi.DF"),which(colnames(x$rand.table)=="elim.num")), P.values=TRUE, has.Pvalue=TRUE)        
}



# generic functions for LSMEANS
#lsmeansTAB <- function(model, data, test.effs=NULL, plot=FALSE , ...) UseMethod("lsmeansTAB")


#lsmeansTAB.default<-function(model, data, test.effs=NULL, plot=FALSE , ...)
lsmeans <- function(model, test.effs=NULL, method.grad="simple", ...)
{
  result <- totalAnovaRandLsmeans(model=model, ddf="Satterthwaite", isLSMEANS=TRUE, test.effs=test.effs, reduce.random=FALSE, reduce.fixed=FALSE, method.grad=method.grad)  
  res <- list(lsmeans.table=result$lsmeans.table, response=result$response)
  class(res) <- "lsmeans"
  res 
}

print.lsmeans <- function(x, ...)
{
  
  cat("Least Squares Means table:\n")
  printCoefmat(data.matrix(x$lsmeans.table), dig.tst=1, tst.ind=c(1:(which(colnames(x$lsmeans.table)=="Estimate")-1),which(colnames(x$lsmeans.table)=="DF")), digits=3 , P.values=TRUE, has.Pvalue=TRUE)       
}

plot.lsmeans <- function(x, ...)
{
  
  #plots for LSMEANS
  if(!is.null(x$lsmeans.table) && nrow(x$lsmeans.table)>0)
    plotLSMEANS(x$lsmeans.table, x$response, "LSMEANS")     
}

difflsmeans <- function(model, test.effs=NULL, method.grad="simple", ...)
{
  result <- totalAnovaRandLsmeans(model=model, ddf="Satterthwaite", isDiffLSMEANS=TRUE, test.effs=test.effs, reduce.random=FALSE, reduce.fixed=FALSE, method.grad=method.grad)  
  res <- list(diffs.lsmeans.table=result$diffs.lsmeans.table, response=result$response)
  class(res) <- "difflsmeans"
  res 
}

print.difflsmeans <- function(x, ...)
{
  
  cat("Differences of LSMEANS:\n")
  printCoefmat(data.matrix(x$diffs.lsmeans.table), dig.tst=1, tst.ind=c(1:(which(colnames(x$diffs.lsmeans.table)=="Estimate")-1),which(colnames(x$diffs.lsmeans.table)=="DF")), digits=3 , P.values=TRUE, has.Pvalue=TRUE)
  
}

plot.difflsmeans <- function(x, ...)
{
  
  #plots for DIFF of LSMEANS
  if(!is.null(x$diffs.lsmeans.table) && nrow(x$diffs.lsmeans.table)>0)
    plotLSMEANS(x$diffs.lsmeans.table, x$response, "DIFF of LSMEANS")   
}








#setGeneric("sigma", function(object, ...) standardGeneric("sigma"))

#setMethod("sigma", signature(object = "mer"),
#          function (object, ...) {
#              dd <- object@dims
#        if(!dd[["useSc"]]) return(1)
#	      object@deviance[[if(dd[["REML"]]) "sigmaREML" else "sigmaML"]]
#          })

#setMethod("summary", signature(object = "summary.mer"), function(object) object)

#setMethod("coef", signature(object = "summary.mer"),
#          function(object, ...) object@coefs)

#setMethod("summary", signature(object = "mer"),
#    function(object, ddf="Satterthwaite", ...)
#      {
#          
#          
#          #calculate df, t-value and p-values
#          coefs.satt<-NULL
#          if(ddf=="Satterthwaite")
#          {
#            #added code: refit model with contr.SAS
#            #options(contrasts=c(unordered="contr.SAS", ordered="contr.poly"))
#            #warning("\n model has been refitted with contrasts=contr.SAS \n")
#            #object<-update(object, .~., object@dims[["REML"]])
# 
#            result <- totalAnovaRandLsmeans(model=object, ddf=ddf, isTtest=TRUE)  
#            coefs.satt <- cbind("df"= result$ttest$df, "p value"= result$ttest$tpvalue)
#          }
#          
#                    
#          REML <- object@dims[["REML"]]
#          fcoef <- fixef(object)
#          vcov <- vcov(object)
#          corF <- vcov@factors$correlation
#          dims <- object@dims
#          coefs <- cbind("Estimate" = fcoef, "Std. Error" = corF@sd) #, DF = DF)
#          llik <- logLik(object, REML)
#          dev <- object@deviance
#          mType <- if((non <- as.logical(length(object@V)))) "NMM" else "LMM"
#          if (gen <- as.logical(length(object@muEta)))
#              mType <- paste("G", mType, sep = '')
#          mName <- switch(mType, LMM = "Linear", NMM = "Nonlinear",
#                          GLMM = "Generalized linear",
#                          GNMM = "Generalized nonlinear")
#	  method <- {
#	      if (mType == "LMM")
#		  if(REML) "REML" else "maximum likelihood"
#	      else
#		  paste("the", if(dims[["nAGQ"]] == 1) "Laplace" else
#			"adaptive Gaussian Hermite",
#			"approximation")
#	  }
#          AICframe <- data.frame(AIC = AIC(llik), BIC = BIC(llik),
#                                 logLik = as.vector(llik),
#                                 deviance = dev[["ML"]],
#                                 REMLdev = dev[["REML"]],
#                                 row.names = "")
#          if (is.na(AICframe$REMLdev)) AICframe$REMLdev <- NULL
#          varcor <- VarCorr(object)
#          REmat <- formatVC(varcor)
#          if (is.na(attr(varcor, "sc")))
#              REmat <- REmat[-nrow(REmat), , drop = FALSE]
#
#          if (nrow(coefs) > 0) {
#              if (!dims[["useSc"]]) {
#                  coefs <- coefs[, 1:2, drop = FALSE]
#                  stat <- coefs[,1]/coefs[,2]
#                  pval <- 2*pnorm(abs(stat), lower = FALSE)
#                  coefs <- cbind(coefs, "z value" = stat, "Pr(>|z|)" = pval)
#              } else {
#                  stat <- coefs[,1]/coefs[,2]
#                  ##pval <- 2*pt(abs(stat), coefs[,3], lower = FALSE)
#                  coefs <- cbind(coefs, "t value" = stat) #, "Pr(>|t|)" = pval)
#                  #added code: inclusion of p-values, df and t-values SAS
#                  coefs <- cbind(coefs, coefs.satt)
#
#              }
#          } ## else : append columns to 0-row matrix ...
#          new("summary.mer",
#              object,
#              methTitle = paste(mName, "mixed model fit by", method),
#              logLik = llik,
#              ngrps = sapply(object@flist, function(x) length(levels(x))),
#              sigma = sigma(object),
#              coefs = coefs,
#              vcov = vcov,
#              REmat = REmat,
#              AICtab= AICframe
#              )
#      })## summary()


### code from lme4 package
## This is modeled a bit after  print.summary.lm :
#printMer <- function(x, digits = max(3, getOption("digits") - 3),
#                     correlation = TRUE, symbolic.cor = FALSE,
#                     signif.stars = getOption("show.signif.stars"), ...)
#{
#    so <- summary(x)
#    REML <- so@dims[["REML"]]
#    llik <- so@logLik
#    dev <- so@deviance
#    dims <- x@dims
#
#    cat(so@methTitle, "\n")
#    if (!is.null(x@call$formula))
#        cat("Formula:", deparse(x@call$formula),"\n")
#    if (!is.null(x@call$data))
#        cat("   Data:", deparse(x@call$data), "\n")
#    if (!is.null(x@call$subset))
#        cat(" Subset:", deparse(x@call$subset),"\n")
#    print(so@AICtab, digits = digits)
#
#    cat("Random effects:\n")
#    print(so@REmat, quote = FALSE, digits = digits, ...)
#
#    ngrps <- so@ngrps
#    cat(sprintf("Number of obs: %d, groups: ", dims[["n"]]))
#    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
#    cat("\n")
#    if (is.na(so@sigma))
#  cat("\nEstimated scale (compare to 1):",
#            sqrt(exp(so@deviance[["lr2"]])/so@dims[["n"]]), "\n")
#    if (nrow(so@coefs) > 0) {
#  ### here the code is changed in order to put p-values
#  cat("\nFixed effects:\n")
#  if(ncol(so@coefs)==5)
#  printCoefmat(so@coefs, zap.ind = 3, #, tst.ind = 4
#		     digits = digits, signif.stars = signif.stars, has.Pvalue=TRUE ,P.values = TRUE)
#  else
#  printCoefmat(so@coefs, zap.ind = 3, #, tst.ind = 4
#  	     digits = digits, signif.stars = signif.stars)

#  if(correlation) {
#	    corF <- so@vcov@factors$correlation
#	    if (!is.null(corF)) {
#		p <- ncol(corF)
#		if (p > 1) {
#		    rn <- rownames(so@coefs)
#		    rns <- abbreviate(rn, minlength=11)
#		    cat("\nCorrelation of Fixed Effects:\n")
#		    if (is.logical(symbolic.cor) && symbolic.cor) {
#			corf <- as(corF, "matrix")
#			dimnames(corf) <- list(rns,
#					       abbreviate(rn, minlength=1, strict=TRUE))
#			print(symnum(corf))
#		    }
#		    else {
#			corf <- matrix(format(round(corF@x, 3), nsmall = 3),
#				       ncol = p, dimnames = list(rns,
#					       abbreviate(rn, minlength=6)))
#			corf[!lower.tri(corf)] <- ""
#			print(corf[-1, -p, drop=FALSE], quote = FALSE)
#		    }
#		}
#	    }
#	}
#   }
#    invisible(x)
#}

#setMethod("print", "mer", printMer)
#setMethod("show", "mer", function(object) printMer(object))
#####setMethod("print", "summary.mer", printMer)

#setMethod("deviance", signature(object = "summary.mer"), function(object) object@deviance)
#setMethod("logLik", signature(object = "summary.mer"), function(object) object@logLik)
#setMethod("vcov", signature(object = "summary.mer"), function(object) object@vcov)
#setMethod("summary", signature(object = "summary.mer"), function(object) object)

#setClass("summary.mer",                 # Additional slots in a summary object
#         representation(           
#  		methTitle = "character",
#			logLik= "logLik",
#			ngrps = "integer",
#			sigma = "numeric", # scale, non-negative number
#			coefs = "matrix",
#			vcov = "dpoMatrix",
#			REmat = "matrix",
#			AICtab= "data.frame"),
#        contains = "mer")

#################################################################################################
##### summary for the newest version of lme4
#################################################################################################


### "format()" the 'VarCorr' matrix of the random effects -- for show()ing
# formatVC <- function(varc, digits = max(3, getOption("digits") - 2),
#                      useScale) {
#   sc <- unname(attr(varc, "sc"))
#   recorr <- lapply(varc, attr, "correlation")
#   reStdDev <- c(lapply(varc, attr, "stddev"), if(useScale) list(Residual = sc))
#   reLens <- unlist(c(lapply(reStdDev, length)))
#   nr <- sum(reLens)
#   reMat <- array('', c(nr, 4),
#                  list(rep.int('', nr),
#                       c("Groups", "Name", "Variance", "Std.Dev.")))
#   reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
#   reMat[,2] <- c(unlist(lapply(varc, colnames)), if(useScale) "")
#   reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
#   reMat[,4] <- format(unlist(reStdDev), digits = digits)
#   if (any(reLens > 1)) {
#     maxlen <- max(reLens)
#     corr <-
#       do.call("rBind",
#               lapply(recorr,
#                      function(x) {
#                        x <- as(x, "matrix")
#                        cc <- format(round(x, 3), nsmall = 3)
#                        cc[!lower.tri(cc)] <- ""
#                        nr <- dim(cc)[1]
#                        if (nr >= maxlen) return(cc)
#                        cbind(cc, matrix("", nr, maxlen-nr))
#                      }))[, -maxlen, drop = FALSE]
#     if (nrow(corr) < nrow(reMat))
#       corr <- rbind(corr, matrix("", nrow = nrow(reMat) - nrow(corr), ncol = ncol(corr)))
#     colnames(corr) <- rep.int("", ncol(corr))
#     colnames(corr)[1] <- "Corr"
#     cbind(reMat, corr)
#   } else reMat
# }
# 
# ## This is modeled a bit after  print.summary.lm :
# ## Prints *both*  'mer' and 'merenv' - as it uses summary(x) mainly
printMerenv <- function(x, digits = max(3, getOption("digits") - 3),
                        correlation = NULL, symbolic.cor = FALSE,
                        signif.stars = getOption("show.signif.stars"), ...)
{
  print("asdfg")
  return()
  so <- summary(x)
  cat(sprintf("%s ['%s']\n",so$methTitle, class(x)))
  if (!is.null(f <- so$family)) {
    cat(" Family:", f)
    if (!(is.null(ll <- so$link))) cat(" (", ll, ")")
    cat("\n")
  }
  ## FIXME: commenting out for now, restore after release?
  ## cat("Scaled residuals:\n")
  ## print(summary(residuals(x,type="pearson",scaled=TRUE)),digits=digits)
  if (!is.null(cc <- so$call$formula))
    cat("Formula:", deparse(cc),"\n")
  ## if (!is.null(so$family)) {
  ##     cat("Family: ",so$family,
  ##         " (link=",so$link,")\n",
  ##         sep="")
  ## }
  if (!is.null(cc <- so$call$data))
    cat("   Data:", deparse(cc), "\n")
  if (!is.null(cc <- so$call$subset))
    cat(" Subset:", deparse(asOneSidedFormula(cc)[[2]]),"\n")
  cat("\n")
  tab <- so$AICtab
  if (length(tab) == 1 && names(tab) == "REML")
    cat("REML criterion at convergence:", round(tab, 4), "\n")
  else print(round(so$AICtab, 4))
  cat("\nRandom effects:\n")
  print(formatVC(so$varcor, digits = digits, useScale = so$useScale),
        quote = FALSE, digits = digits, ...)
  
  ngrps <- so$ngrps
  cat(sprintf("Number of obs: %d, groups: ", so$devcomp$dims[["n"]]))
  cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
  cat("\n")
  p <- nrow(so$coefficients)
  if (p > 0) {
    cat("\nFixed effects:\n")
    printCoefmat(so$coefficients, zap.ind = 3, #, tst.ind = 4
                 digits = digits, signif.stars = signif.stars)
    if(!is.logical(correlation)) { # default
      correlation <- p <= 20
      if(!correlation) {
        nam <- deparse(substitute(x)) # << TODO: improve if this is called from show()
        cat(sprintf(paste("\nCorrelation matrix not shown by default, as p = %d > 20.",
                          "Use print(%s, correlation=TRUE)  or",
                          "    vcov(%s)	 if you need it\n", sep="\n"),
                    p, nam, nam))
      }
    }
    if(correlation) {
      if(is.null(VC <- so$vcov)) VC <- vcov(x)
      corF <- VC@factors$correlation
      if (is.null(corF)) {
        cat("\nCorrelation of Fixed Effets is not available\n")
      }
      else {
        p <- ncol(corF)
        if (p > 1) {
          rn <- rownames(so$coefficients)
          rns <- abbreviate(rn, minlength=11)
          cat("\nCorrelation of Fixed Effects:\n")
          if (is.logical(symbolic.cor) && symbolic.cor) {
            corf <- as(corF, "matrix")
            dimnames(corf) <- list(rns,
                                   abbreviate(rn, minlength=1, strict=TRUE))
            print(symnum(corf))
          }
          else {
            corf <- matrix(format(round(corF@x, 3), nsmall = 3),
                           ncol = p,
                           dimnames = list(rns, abbreviate(rn, minlength=6)))
            corf[!lower.tri(corf)] <- ""
            print(corf[-1, -p, drop=FALSE], quote = FALSE)
          }
        }
      }
    }
  }
  invisible(x)
}## printMerenv()
# 
# ##' @importFrom stats vcov
# ##' @S3method vcov summary.merMod
# vcov.summary.merLmerTest <- function(object, correlation = TRUE, ...) {
#   if(is.null(object$vcov)) stop("logic error in summary of merlmerTest object")
#   object$vcov
# }
# 
# ##' @S3method print merMod
 print.merLmerTest <- printMerenv
# 
# ##' @exportMethod show
# #setMethod("show",  "merLmerTest", function(object) printMerenv(object))
# 
# 
# ##' @S3method print summary.merMod
# print.summary.merLmerTest <- printMerenv
# 
# ##' @S3method summary merMod
# summary.merLmerTest <- function(object, ...)
# {
#   resp <- object@resp
#   devC <- object@devcomp
#   dd <- devC$dims
#   cmp <- devC$cmp
#   useSc <- as.logical(dd["useSc"])
#   sig <- sigma(object)
#   REML <- isREML(object)
#   
#   link <- fam <- NULL
#   if(is(resp, "glmResp")) {
#     fam <- resp$family$family
#     link <- resp$family$link
#   }
#   coefs <- cbind("Estimate" = fixef(object),
#                  "Std. Error" = sig * sqrt(diag(object@pp$unsc())))
#   if (nrow(coefs) > 0) {
#     coefs <- cbind(coefs, coefs[,1]/coefs[,2], deparse.level=0)
#     colnames(coefs)[3] <- paste(if(useSc) "t" else "z", "value")
#     ## additional code for p-values
#     if(isLMM(object)) {
#       coefs.satt <- cbind(coefs,object@t.pval) 
#       coefs <- coefs.satt    
#       colnames(coefs)[4] <- "Pr(>|t|)"
#     }
#     
#     if (isGLMM(object)) {
#       coefs <- cbind(coefs,2*pnorm(abs(coefs[,3]),lower.tail=FALSE))
#       colnames(coefs)[4] <- c("Pr(>|z|)")
#     }
#   }
#   mName <- paste(switch(1L + dd[["GLMM"]] * 2L + dd[["NLMM"]],
#                         "Linear", "Nonlinear",
#                         "Generalized linear", "Generalized nonlinear"),
#                  "mixed model fit by",
#                  if(REML) "REML" else "maximum likelihood")
#   llik <- logLik(object)   # returns NA for a REML fit - maybe change?
#   AICstats <- {
#     if (REML) cmp["REML"] # do *not* show likelihood stats here
#     else {
#       c(AIC = AIC(llik), BIC = BIC(llik), logLik = c(llik),
#         deviance = deviance(object))
#     }
#   }
#   ## FIXME: You can't count on object@re@flist,
#   ##        nor compute VarCorr() unless is(re, "reTrms"):
#   varcor <- VarCorr(object)
#   # use S3 class for now
#   structure(list(methTitle=mName, devcomp=devC,
#                  isLmer=is(resp, "lmerResp"), useScale=useSc,
#                  logLik=llik, family=fam, link=link,
#                  ngrps=sapply(object@flist, function(x) length(levels(x))),
#                  coefficients=coefs, sigma=sig,
#                  vcov=vcov(object, correlation=TRUE, sigm=sig),
#                  varcor=varcor, # and use formatVC(.) for printing.
#                  AICtab=AICstats, call=object@call
#   ), class = "summary.merLmerTest")
# }
# 
# ##' @S3method summary summary.merMod
# summary.summary.merLmerTest <- function(object, varcov = FALSE, ...) {
#   if(varcov && is.null(object$vcov))
#     object$vcov <- vcov(object, correlation=TRUE, sigm = object$sigma)
#   object
# }
