totalAnovaRandLsmeans<-function(model, ddf, type=3, alpha.random = 0.1, alpha.fixed = 0.05, reduce.fixed = TRUE, reduce.random = TRUE, lsmeans.calc = TRUE, difflsmeans.calc=TRUE,  isTotal=FALSE, isAnova=FALSE, isRand=FALSE, isLSMEANS=FALSE, isDiffLSMEANS=FALSE, test.effs=NULL, method.grad="Richardson")
{
  #errors in specifying the parameters
  if(!isRand && !(ddf %in% c("Satterthwaite","Kenward-Roger")))
  {
    print ('Error: parameter ddf is wrongly specified')
    stop()
  }
  if(!isRand && type!=3)
  {
    print ('Error: parameter type is wrongly specified')
    stop()
  }  
  
  data<-summary(model)@frame
  #not to show the warnings  
  options(warn=-1)    
    
  result<-NULL
  anova.table<-NULL
  
  result$response<-names(attr(terms(model),"dataClasses")[1])
  
  #model<-lmer(formula=formula, data=data)
  
  #update model
  # change unordered contrasts to contr.SAS
  # change REML to TRUE
  options(contrasts=c(unordered="contr.SAS", ordered="contr.poly"))
  
  #model<-update(model, REML=TRUE)
  model<-
        if (model@dims['REML'] == 1)
        {
            model
        }
        else
        {
  	    warning("\n model has been refitted with REML=TRUE \n")
            update(model,.~.,REML=TRUE)
        }
  mf.final<-update.formula(formula(model),formula(model))
  
   
  
  #update data 
  #eliminate missing values
  if(ncol(get_all_vars(mf.final, data))>ncol(data))
  {
     #mframe<-model.frame(mf.final, data=data, na.action=na.pass)
     #data$response<-mframe[,1]
     data$response<-data[,1]
     data<-data[,-1]
     fm<-paste(mf.final)
     fm[2]<-"response"
     mf.final<- as.formula(paste(fm[2],fm[1],fm[3], sep=""))
     mf.final<-update.formula(mf.final,mf.final)
  }
  
  data<-get_all_vars(mf.final, data)
  data<-data[complete.cases(data),]
  
  
  
  model<-eval(substitute(lmer(mf.final, data=data),list(mf.final=mf.final)))
  
  
  
    # check if there are no fixed effects
  #if(length(attr(delete.response(terms(model)),"term.labels"))==0)
 # {
  #  fm<-getFormula(model, withRand=TRUE)
  #  fm<-as.formula(paste(fm[2],fm[1],paste("1+", fm[3], sep=""), sep=""))
  #  model<-eval(substitute(lmer(mf.final, data=data),list(mf.final=mf.final)))
  #}
   
  
   
  # save the call of the model              
  result$call<-model@call
  
  #check if there are correlations between intercept and slope
  result$corr.intsl<-checkCorr(model)
  
  
  
  
  #Reduction of random effects
  #if(isRand || isTotal )
  #{     
    # perform reduction of random effects
    #if(!isRandReduce)
    #  result.rand<-elimRandLmer(model, data, 1)
    #else
    #  result.rand<-elimRandLmer(model, data, alpha.rand)
    
    result.rand<-elimZeroVarOrCorr(model, data)
    model<-result.rand$model
    #save results for fixed effects for model with only fixed effects
    if(class(model) == "lm" | class(model) == "gls")
    {
      result<-saveResultsFixModel(result, model)
      result$rand.table=NULL
      return(result)
    }
    
    #if(!isRandReduce)
    #  result.rand<-elimRandEffs(model, data, 1)
  
    #analysis of the random part  
    if(isRand)
      reduce.random<-FALSE
    #if(reduce.random)
      result.rand<-elimRandEffs(model, data, alpha.random, reduce.random)  
    #if(!reduce.random && isRand)
    #  result.rand<-elimRandEffs(model, data, 1)
  
    model<-result.rand$model
       
    #convert rand table to data frame
    rt<-as.data.frame(result.rand$TAB.rand)
    rt$Chi.DF<-as.integer(rt$Chi.DF)
    if(!is.null(rt$elim.num))
      rt$elim.num<-as.integer(rt$elim.num)
  
    result$rand.table<-rt
    if(isRand)
      return(result)
    
  #}
  
  
  
  #save results for fixed effects for model with only fixed effects
  if(class(model) == "lm" | class(model) == "gls")
    return(saveResultsFixModel(result, model))
  
    
  #perform reduction of fixed effects for model with mixed effects
  stop = FALSE
  is.first.anova<-TRUE
  is.first.sign<-TRUE
  
  
  while(!stop)
  {
  
      
    
      # if there are no fixed terms
      if(nrow(anova(model))==0)
      {
        if(is.null(anova.table))
        {
          
          if(isLSMEANS || isDiffLSMEANS)
          {
             lsmeans.summ<- matrix(ncol=7,nrow=0)
             colnames(lsmeans.summ)<-c("Estimate","Standard Error", "DF", "t-value", "Lower CI", "Upper CI", "p-value")
             lsmeans.summ<-as.data.frame(lsmeans.summ)
             if(isLSMEANS)
               result$lsmeans.table<-lsmeans.summ
             if(isDiffLSMEANS)
               result$diffs.lsmeans.table<-lsmeans.summ
             return(result)
          }
          result$model<-model
          result$anova.table<-anova(model)
          return(result)        
          
        }          
        break
      }
        
          
      
      # save lmer outcome in rho environmental variable
      rho<-rhoInit(model)     
      
      # calculate asymptotic covariance matrix A
      h <- hessian(function(x) Dev(rho,x), rho$param$vec.matr)
      rho$A <- 2*solve(h)
      #rho$A <- 2*ginv(h)
      
      #Check if A is positive-definite
      isposA<-all(eigen(rho$A)$values>0)
      
      if(!isposA)
      {
        print("ERROR: asymptotic covariance matrix A is not positive!")
        result$model<-model
        TABs<-emptyAnovaLsmeansTAB()
        result$anova.table<-TABs$TAB.fixed
        result$lsmeans.table<-TABs$TAB.lsmeans
        result$diffs.lsmeans.table<-TABs$TAB.lsmeans
        return(result)
      }
      
      #calculate lsmeans of differences of LSMEANS of the final model
      if(isLSMEANS || isDiffLSMEANS)
      {
        if(isLSMEANS)
        {
          lsmeans.tab<-calcLSMEANS(model, data, rho, alpha.fixed, test.effs=test.effs, method.grad=method.grad, lsmeansORdiff=TRUE)
          result$lsmeans.table<-lsmeans.tab$summ.data
          result$diffs.lsmeans.table<-NULL
        }
        if(isDiffLSMEANS)
        {
          lsmeans.tab<-calcLSMEANS(model, data, rho, alpha.fixed, test.effs=test.effs, method.grad=method.grad, lsmeansORdiff=FALSE)
          result$diffs.lsmeans.table<-lsmeans.tab$summ.data
          result$lsmeans.table<-NULL
        }
        return(result)
      }
      
          
      # Calculate  F-test with Satterthwaite's approximation
      # create X design matrix for fixed effects
      X.design <- createDesignMat(model,data)
       
      
      # define full coefficients for model
      #coefs.real<-fullCoefs(model, data)
      #coefs.real<-getDummyCoefs(model, data)
      
      
      #save full coefficients in rho
      #rho$s.test<-coefs.real
      nums.dummy.coefs <- getNumsDummyCoefs(model, data)
      rho$nums.zeroCoefs <- nums.dummy.coefs$nums.zeroCoefs
      rho$nums.Coefs <- nums.dummy.coefs$nums.Coefs
      
     
      
      #define the terms that are to be tested
      test.terms <- attr(terms(model),"term.labels")
      
      #initialize anova table
      if(is.first.anova)
      {
        anova.table<-initAnovaTable(model, reduce.fixed)
        is.first.anova<-FALSE
        elim.num<-1
      }
      
      # calculate general set matrix for type 3 hypothesis
      L <- calcGeneralSetForHypothesis(X.design, rho)
      
      for(i in 1:length(test.terms))
      {
                 
         result.fstat <- calcFpvalue(test.terms[i], L, model, rho, ddf, method.grad=method.grad)   
         if(!is.null(result.fstat))
         {
            anova.table[which(rownames(anova.table) == test.terms[i]),2]<-result.fstat$denom#round(result.fstat$denom,2)
            anova.table[which(rownames(anova.table) == test.terms[i]),3]<-result.fstat$Fstat#round(result.fstat$Fstat,2)
            anova.table[which(rownames(anova.table) == test.terms[i]),which(colnames(anova.table)=="p.value")]<-result.fstat$pvalue#round(result.fstat$pvalue,4)
            #if(!is.first.anova)
            #{
            #  anova.model<-anova(model)
            #  anova.table[which(rownames(anova.table) == test.terms[i]),1:4]<-anova.model[which(rownames(anova.model) == test.terms[i]),]
            #}
              
         }
       
      }
     
      if(!reduce.fixed)
        break
      else
      {
        resNSelim <- elimNSFixedTerm(model, anova.table, data, alpha.fixed, elim.num)
        if(is.null(resNSelim))
          break
        else
        {
          model <- resNSelim$model
          mf.final <- update.formula(formula(model),formula(model))
          model <- eval(substitute(lmer(mf.final, data=data),list(mf.final=mf.final)))
          anova.table <- resNSelim$anova.table
          elim.num <- elim.num+1
        }        
      }      
    
  }
  
  #convert anova table to data frame
  anova.table<-as.data.frame(anova.table)
  anova.table$NumDF<-as.integer(anova.table$NumDF)
  if(!is.null(anova.table$elim.num))
    anova.table$elim.num<-as.integer(anova.table$elim.num)
  
  if(isTotal || isAnova)
  {
    result$anova.table <- anova.table
    if(isAnova)
      return(result)
  }
  
  #if in step function least squares means of diffs of LSMEANS are required
  if(lsmeans.calc)
  {
    lsmeans.tab <- calcLSMEANS(model, data, rho, alpha.fixed, test.effs=test.effs, method.grad=method.grad, lsmeansORdiff=TRUE)
    result$lsmeans.table <- lsmeans.tab$summ.data
  }
  else
  {
    result$lsmeans.table <- NULL
  }
  if(difflsmeans.calc)
  {
    lsmeans.tab <- calcLSMEANS(model, data, rho, alpha.fixed, test.effs=test.effs, method.grad=method.grad, lsmeansORdiff=FALSE)
    result$diffs.lsmeans.table <- lsmeans.tab$summ.data
  }
  else
  {
    result$diffs.lsmeans.table <- NULL
  }
  
  #update model
  mf.final <- update.formula(formula(model),formula(model))
  model <- eval(substitute(lmer(mf.final, data=data),list(mf.final=mf.final)))
  
  #save model
  result$model <- model
  return(result)
}


# generic functions for total analysis on mixed models
#totalAnalysis <- function(model, data,  alpha.rand = 0.05, alpha.fix = 0.05, isFixReduce = FALSE, isRandReduce = FALSE, test.effs=NULL, plot=FALSE, ...) UseMethod("totalAnalysis")


step <- function(model, ddf="Satterthwaite", type=3, alpha.random = 0.1, alpha.fixed = 0.05, reduce.fixed = TRUE, reduce.random = TRUE, lsmeans.calc=TRUE, difflsmeans.calc=TRUE, test.effs=NULL, method.grad="Richardson",...)
{  
  result <- totalAnovaRandLsmeans(model=model, ddf=ddf , type=type,  alpha.random=alpha.random, alpha.fixed=alpha.fixed, reduce.fixed=reduce.fixed, reduce.random=reduce.random, lsmeans.calc=lsmeans.calc, difflsmeans.calc=difflsmeans.calc, isTotal=TRUE, test.effs=test.effs, method.grad=method.grad)
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
  cat("Call:\n")
  print(x$call)
  
  if(!is.null(x$rand.table))
  {
    cat("\nRandom effects:\n")  
    printCoefmat(x$rand.table, digits=3 , dig.tst=1  ,tst.ind=c(which(colnames(x$rand.table)=="Chi.DF"),which(colnames(x$rand.table)=="elim.num")), P.values=TRUE, has.Pvalue=TRUE)
  }
  
  if(nrow(x$anova.table)!=0)
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


plot.step<-function(x, ...)
{
  if(!is.null(x$lsmeans.table))
    plotLSMEANS(x$lsmeans.table, x$response, "LSMEANS")     
  if(!is.null(x$diffs.lsmeans.table))
    plotLSMEANS(x$diffs.lsmeans.table, x$response, "DIFF of LSMEANS")
}

# generic functions for ANOVA
#anovaTAB <- function(model, data, ...) UseMethod("anovaTAB")


#anovaTAB.default<-function(model, data, ...)
Anova <- function(model, ddf="Satterthwaite", type=3, method.grad="Richardson",...)
{
  result <- totalAnovaRandLsmeans(model=model, ddf=ddf, type=type, isAnova=TRUE, reduce.random=FALSE, reduce.fixed=FALSE, method.grad=method.grad)  
  res <- list(anova.table=result$anova.table)
  class(res) <- "Anova"
  res
}

print.Anova <- function(x, ...)
{

  if(nrow(x$anova.table)!=0)
  {
    cat("Analysis of Variance Table:\n")    
    printCoefmat(x$anova.table, dig.tst=1, tst.ind=c(1,2), cs.ind=3, digits=3 , P.values=TRUE, has.Pvalue=TRUE)
  }
  else
    print(x$anova.table)
}


# generic functions for random effects
#randTAB <- function(model, data, ...) UseMethod("randTAB")


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
lsmeans <- function(model, test.effs=NULL, method.grad="Richardson", ...)
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
    plotLSMEANS(x$lsmeans.table, x$response, "LSMEANS")     
}

difflsmeans <- function(model, test.effs=NULL, method.grad="Richardson", ...)
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
    plotLSMEANS(x$diffs.lsmeans.table, x$response, "DIFF of LSMEANS")   
}
