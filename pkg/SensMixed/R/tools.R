## check if the data is balanced
isbalanced <- function(data)
{
  suppressWarnings(!is.list(replications(~ . , data)))
}

### function checks  if there are zero cells in a factor term
checkZeroCell <- function(data, factors)
{
  t <- table(data[, match(factors, names(data))])
  if(length(which(t==0))>0)
  {
    message(paste("Some of the combinations of ", paste(factors,collapse=":"), " has no data, therefore this combination will not be part of the initial model"))
    cat("\n")
    return(TRUE)
  }
  
  return(FALSE)
}

### checks if the number of levels for an interaction term is equal to number of observations
checkNumberInteract <- function(data, factors)
{
  ## returns TRUE if number of levels is equal to nrow of data
  
  nlev <- 1
  for(i in 1:length(factors))
  {    
    if(!is.factor(data[, match(factors[i], names(data))]))
      next()
    nlev <- nlev * nlevels(data[, match(factors[i], names(data))])
  }
  if(nlev >= nrow(data))
  {
    warning.str <- "Number of levels for "
    if(length(factors) > 1)
      warning.str <- c(warning.str," interaction ", sep=" ")
    #for(i in length(factors))
    #    warning.str <- paste(warning.str, factors[i],sep=" ")  
    warning.str <- c(warning.str, paste(factors,collapse=":"), " is more or equal to the number of observations in data", sep=" ")    
    message(warning.str)
    cat("\n")
    return(TRUE)
  }
  return(FALSE)
}


### Function converts variables to factors
convertToFactors <- function(data, facs)
{
  #convert effects to factors
  for(fac in facs)
    data[,fac] <- as.factor(data[,fac])
  data
}

## create formula with only fixed terms
fixedFormula <- function(fmodel)
{
  terms.fm <- attr(terms.formula(fmodel),"term.labels")
  ind.rand.terms <- which(unlist(lapply(terms.fm,function(x) substring.location(x, "|")$first))!=0)
  terms.fm[ind.rand.terms] <- unlist(lapply(terms.fm[ind.rand.terms],function(x) paste("(",x,")",sep="")))
  fm <- paste(fmodel)
  fm[3] <- paste(terms.fm[-ind.rand.terms],collapse=" + ")
  if(fm[3]=="")
    fo <- as.formula(paste(fm[2],fm[1],1, sep=""))
  else
    fo <- as.formula(paste(fm[2],fm[1],fm[3], sep=""))
  return(fo)
}

### Create an lmer model
createLMERmodel <- function(structure, data, response, fixed, random, corr, MAM=FALSE)
{ 
  
  #construct formula for lmer model    
  mf.final <- createFormulaAllFixRand(structure, data, response, fixed, random, corr)    
  ## if MAM needs to be contructed
  if(MAM){
    if(length(fixed$Product)>1){
      data$prod <- interaction(data[, fixed$Product[1]], data[, fixed$Product[2]])
      mf.final.lsm <- createFormulaAllFixRand(structure, data, response, 
                                              list(Product="prod", Consumer=fixed$Consumer), random, corr)   
      #ff <- as.formula(paste(mf.final[2], mf.final[1], "prod", sep=""))
    }else{
      mf.final.lsm <- mf.final
    }
    
    ## create formulas for anova and lsmeans   
    ############################################################################
    ff <- fixedFormula(mf.final)   
    data$x <- scale(predict(lm(ff, data=data)), scale=FALSE)    
    ## for anova
    fm <- paste(mf.final)
    fm[3] <- paste(fm[3], paste(random$individual, "x", sep=":"), sep=" + ")
    fo.anova <- as.formula(paste(fm[2], fm[1], fm[3], sep=""))    
    ## for lsmeans
    fm.lsm <- paste(mf.final.lsm)
    fm.lsm[3] <- paste(fm.lsm[3],  paste(random$individual, "x", sep=":"),
                       "x", sep=" + ")
    fo.lsm <- as.formula(paste(fm.lsm[2], fm.lsm[1], fm.lsm[3], sep=""))
    
    
    ## create models for anova and lsmeans
    ############################################################################
    
    ## for anova
    model.anova <- lmerTest::lmer(fo.anova, data)
    #anova(model.anova, type=1) ## TODO: compare with SAS
    #st.anova <- step(model.anova, lsmeans.calc=FALSE, difflsmeans.calc=FALSE, 
    #reduce.fixed=FALSE)
    
    ## change contrasts for lsmeans to be contr.sum
    if(length(fixed$Product)==1){
      mm <- model.matrix(model.anova)
      l <- attr(mm, "contrasts")
      contr <- l
      names.facs <- names(contr)
      l <- as.list(rep("contr.sum", length(names.facs)))
      names(l) <- names(contr)
    }else{
      l <- as.list(rep("contr.sum", 2))
      names(l) <- c("prod", random$individual)       
    }   
    
    ## model for lsmeans
    model.lsmeans <- lmerTest::lmer(fo.lsm, data, contrasts=l)
    return(list(model.anova=model.anova, model.lsmeans=model.lsmeans))
    #summaryBy(Coloursaturation ~ prod , data)
    #st.lsmeans <- step(model.lsmeans, lsmeans.calc=FALSE, difflsmeans.calc=FALSE, reduce.fixed=FALSE)
    #newm <- lmerTest::lmer(formula(st.lsmeans$model), data=data, contrasts=l)
    #lsmeans::lsmeans(object=model.lsmeans,  
    #                 pairwise ~ prod)
    #if(length(fixed$Product)==1)
    #  eval(substitute(lsmeans::lsmeans(object=model.lsmeans,  
    #                                 pairwise ~ prod), 
    #                                 list(prod=as.name(fixed$Product))))
    #else
      #eval(substitute(lsmeans::lsmeans(object=model.lsmeans,  
      #                               pairwise ~ prod), 
      #              list(prod=as.name(paste(fixed$Product, collapse=":")))))
    
   #lsmeans::lsmeans(model.lsmeans, pairwise ~ TVset:Picture)
  }else{
    model <- lmerTest::lmer(mf.final, data) 
    return(model)
  }
  
    
  #model <- as(model,"mer")
  #model <- update(model)
  
  #mf.final <- update.formula(formula(model),formula(model))
  #model <- eval(substitute(lmer(mf.final, data=data),list(mf.final=mf.final)))
  #model <- update(model, data=data ,REML=TRUE)
  
  return(model)
}

# check an interaction term for validity
checkComb <- function(data, factors)
{
  return(checkNumberInteract(data,factors) || checkZeroCell(data, factors))
}

.fixedrand <- function(model)
{
  effs <- attr(terms(formula(model)), "term.labels")
  neffs <- length(effs)
  randeffs <- effs[grep(" | ", effs)]
  randeffs <- sapply(randeffs, function(x) substring(x, 5, nchar(x)))
  fixedeffs <- effs[!(effs %in% names(randeffs))]
  return(list(randeffs=randeffs, fixedeffs=fixedeffs))
}

.fillpvalues <- function(x, pvalue)
{
  pvalue[rownames(x$anova.table),x$response] <- x$anova.table[,6]
  pvalue
}