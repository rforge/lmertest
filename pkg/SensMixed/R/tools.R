## check if the data is balanced
isbalanced <- function(data)
{
  suppressWarnings(!is.list(replications(~ . , data)))
}

runMAM <- function(data, Prod_effects, individual, attributes, adjustedMAM=FALSE, 
                   alpha_conditionalMAM=1){
  if(length(attributes) < 2)
    stop("number of attributes for MAM should be more than 1")
  if(length(Prod_effects) > 1)
    stop("should be one-way product structure")  
  dataMAM <- data[, c(individual, Prod_effects)]
  dataMAM$replication <- rep(0, nrow(data))
  dataMAM[, 1] <- as.factor(dataMAM[, 1])
  dataMAM[, 2] <- as.factor(dataMAM[, 2])
  if(nlevels(dataMAM[, 2]) < 3)
    stop("There MUST be at least 3 products")
  ## create a rep factor
  assprod <- interaction(dataMAM[, 1], dataMAM[, 2])
  t <- table(assprod)
  if(length(unique(t))!=1)
    stop("data is unbalanced")
  for(i in 1:length(names(t)))
    dataMAM$replication[assprod==names(t)[i]] <- 1:unique(t)
  dataMAM <- cbind(dataMAM, data[, attributes])
  return(MAManalysis(dataMAM, adjustedMAM, alpha_conditionalMAM))
}


### function checks  if there are zero cells in a factor term
checkZeroCell <- function(data, factors)
{
  t <- table(data[, match(factors, names(data))])
  if(length(which(t==0))>0)
  {
    message(paste("Some of the combinations of ", paste(factors,collapse=":"), 
                  " has no data, therefore this combination will not be part of the initial model"))
    cat("\n")
    return(TRUE)
  }
  
  return(FALSE)
}

## checks if the number of levels for an interaction term 
## is equal to number of observations
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
    warning.str <- c(warning.str, paste(factors,collapse=":"), 
                     " is more or equal to the number of observations in data", 
                     sep=" ")    
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
  ind.rand.terms <- which(unlist(lapply(terms.fm,
                                        function(x) substring.location(x, "|")$first))!=0)
  terms.fm[ind.rand.terms] <- unlist(lapply(terms.fm[ind.rand.terms],
                                            function(x) paste("(",x,")",sep="")))
  fm <- paste(fmodel)
  fm[3] <- paste(terms.fm[-ind.rand.terms],collapse=" + ")
  if(fm[3]=="")
    fo <- as.formula(paste(fm[2],fm[1],1, sep=""))
  else
    fo <- as.formula(paste(fm[2],fm[1],fm[3], sep=""))
  return(fo)
}

### Create an lmer model
createLMERmodel <- function(structure, data, response, fixed, random, corr, 
                            MAM=FALSE, mult.scaling = FALSE, 
                            calc_post_hoc = FALSE)
{ 
  
  #construct formula for lmer model    
  mf.final <- createFormulaAllFixRand(structure, data, response, fixed, random, 
                                      corr)    
  ## if MAM needs to be contructed
  if(MAM){
    if(length(fixed$Product)>1){
      data$prod <- interaction(data[, fixed$Product[1]], 
                               data[, fixed$Product[2]])
      mf.final.lsm <- createFormulaAllFixRand(structure, data, response, 
                                              list(Product="prod", 
                                                   Consumer=fixed$Consumer), 
                                              random, corr)   
    }else{
      mf.final.lsm <- mf.final
    }
    
    ## create formulas for anova and lsmeans   
    ############################################################################
   
    
    ff <- fixedFormula(mf.final)
    if(length(fixed$Product)>1 && mult.scaling){
      prods <- paste("x", fixed$Product, sep="")
      data[, prods] <- lapply(paste(ff[2], ff[1], fixed$Product), 
                                function(formulas)  
                                scale(predict(lm(as.formula(formulas), 
                                                 data=data)), scale=FALSE) )
    }
    
    # create x out of predicted values from lm
    data$x <- rep(NA, nrow(data))
    x.prd <- scale(predict(lm(ff, data=data)), scale=FALSE)
    notNA <- rownames(x.prd)
    data[notNA, "x"] <- x.prd
    
    ## for anova
    fm <- paste(mf.final)
    if(is.list(random)){
      if(length(fixed$Product)>1 && mult.scaling)
        fm[3] <- paste(fm[3], paste(random$individual, prods, sep=":", 
                                    collapse=" + "), sep =" + ")
      else
        fm[3] <- paste(fm[3], paste(random$individual, "x", sep=":"), sep=" + ")
    }
    else{
      if(length(fixed$Product)>1 && mult.scaling)
        fm[3] <- paste(fm[3], paste(random, prods, sep=":", 
                                    collapse=" + "), sep =" + ")
      else
        fm[3] <- paste(fm[3], paste(random, "x", sep=":"), sep=" + ")      
    }
    fo.anova <- as.formula(paste(fm[2], fm[1], fm[3], sep=""))    
    
    ## for lsmeans
    fm.lsm <- paste(mf.final.lsm)
    if(is.list(random))
      fm.lsm[3] <- paste(fm.lsm[3],  paste(random$individual, "x", sep=":"),
                       "x", sep=" + ")
    else
      fm.lsm[3] <- paste(fm.lsm[3],  paste(random, "x", sep=":"),
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
      if(is.list(random))
        names(l) <- c("prod", random$individual) 
      else
        names(l) <- c("prod", random) 
    }   
    
    ## model for lsmeans
    if(calc_post_hoc){
      model.lsmeans <- lmerTest::lmer(fo.lsm, data, contrasts=l)
      return(list(model.anova=model.anova, model.lsmeans=model.lsmeans))
    }
    else return(list(model.anova=model.anova))
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

.renameScalingTerm <- function(tableWithScaling, Prod_effects){
#   if(length(Prod_effects)>1){    
#     xprods <- paste("x", Prod_effects, sep="") 
#     for(indProd in 1:length(xprods)){
#       rownames(tableWithScaling)[unlist(lapply(rownames(tableWithScaling), 
#                                           function(x) 
#                                             grepl(xprods[indProd], x)))] <-
#         paste("Scaling", substring(xprods[indProd], 2, nchar(xprods[indProd])),
#               sep=" ")
#     }
#   }
#   else
    rownames(tableWithScaling)[unlist(lapply(rownames(tableWithScaling), 
                                        function(x) grepl(":x", x)))] <- 
    "Scaling"
  tableWithScaling
}


###############################################################################
# get terms contained  - from lmerTest package
###############################################################################
getIndTermsContained <- function(allterms, ind.hoi)
{
  
  terms.hoi.split <- strsplit(allterms[ind.hoi],":")
  ind.terms.contain <- NULL
  #check which of the terms are contained in the highest order terms
  for(i in (1:length(allterms))[-ind.hoi]) 
  {
    isContained<-FALSE
    for(j in 1:length(terms.hoi.split))
    {
      #if the term is contained in some of the highest order interactions then 
      #we cannot test it for significance
      if(length(which(unlist(strsplit(allterms[i],":")) %in% terms.hoi.split[[j]] == FALSE))==0)
      {
        isContained <- TRUE
        break
      }                
    }
    if(isContained)
      ind.terms.contain <- c(ind.terms.contain,i)
    
  }
  # if there are no terms that are contained in the maximum order effects
  # then compare all the terms between each other for the maximum p value
  if( is.null(ind.terms.contain) )
    return(NULL)
  return(ind.terms.contain)
}

###############################################################################
# get terms contained 
###############################################################################
.getIndTermsContained <- function(allterms, ind.hoi)
{
  
  terms.hoi.split <- strsplit(allterms[ind.hoi],":")
  ind.terms.contain <- NULL
  #check which of the terms are contained in the highest order terms
  for(i in (1:length(allterms))[-ind.hoi]) 
  {
    isContained<-FALSE
    for(j in 1:length(terms.hoi.split))
    {
      #if the term is contained in some of the highest order interactions then 
      #we cannot test it for significance
      if(length(which(unlist(strsplit(allterms[i],":")) %in% terms.hoi.split[[j]] == FALSE))==0)
      {
        isContained <- TRUE
        break
      }                
    }
    if(isContained)
      ind.terms.contain <- c(ind.terms.contain,i)
    
  }
  # if there are no terms that are contained in the maximum order effects
  # then compare all the terms between each other for the maximum p value
  if( is.null(ind.terms.contain) )
    return(NULL)
  return(ind.terms.contain)
}


## get the pure lsmeans for an interaction term
getPureInter <- function(lsm.table, anova.table, eff){

  rows.lsm <- sapply(rownames(lsm.table), 
                     function(x) strsplit(x, " ")[[1]][1]) 
  pure.inter.lsm <- lsm.table[which(rows.lsm %in% eff), ]
  
  contained.effs <- 
    rownames(anova.table)[.getIndTermsContained(rownames(anova.table), 
                                               which(rownames(anova.table) 
                                                     == eff))]
  ## deltas for 3 way interactions
  if( length(unlist(strsplit(eff,":"))) == 3 ){
    ind.inteffs <- grep(":", contained.effs)
    ##plust main effs
    main.effs <- contained.effs[-ind.inteffs]
    ##minus the interactions
    contained.effs <- contained.effs[ind.inteffs]
  }

  p1 <- pure.inter.lsm[ , 1:which(colnames(pure.inter.lsm)=="Estimate") ]
  for(ceff in contained.effs){
    p1  <- merge(p1, 
                 lsm.table[rows.lsm == ceff, 
                           c(unlist(strsplit(ceff,":")), "Estimate")], 
                 by = unlist(strsplit(ceff,":")))
    p1[, "Estimate.x"] <- p1[, "Estimate.x"] - p1[, "Estimate.y"]
    colnames(p1)[which(colnames(p1) == "Estimate.x")] <- "Estimate"
    p1 <- p1[ ,- which(colnames(p1) == "Estimate.y")]        
  }
  ## plus the main effs
  if( length(unlist(strsplit(eff,":"))) == 3 ){
    for(ceff in main.effs){
      p1  <- merge(p1, 
                   lsm.table[rows.lsm == ceff, 
                             c(unlist(strsplit(ceff,":")), "Estimate")], 
                   by = unlist(strsplit(ceff,":")))
      p1[, "Estimate.x"] <- p1[, "Estimate.x"] + p1[, "Estimate.y"]
      colnames(p1)[which(colnames(p1) == "Estimate.x")] <- "Estimate"
      p1 <- p1[ ,- which(colnames(p1) == "Estimate.y")]        
    }
  }
  p1
}

## calculate pure diffs
.calcPureDiffs <- function(pureinter){
  puredifs <- matrix(0, ncol=nrow(pureinter), nrow = nrow(pureinter))
  for (i in 1:nrow(pureinter)) for (j in 1:i) 
    puredifs[i,j] <- pureinter[i, "Estimate"] -  pureinter[j, "Estimate"]  
  puredifs
}

## calculate average d prime from the step function
.calcAvDprime <- function(model, anova.table, dlsm.table, lsm.table){
  sigma <- summary(model, "lme4")$sigma
  rows <- sapply(rownames(dlsm.table), 
                 function(x) strsplit(x, " ")[[1]][1]) 
  anova.table$dprimeav <- rep(1, nrow(anova.table))
  
  for(eff in rownames(anova.table)){
    #lsm <- lsmeans::.old.lsmeans( model , pairwise ~ Track:SPL:Car)
    #dp <- lsm[[2]][, 1]/sigma
    
    
    ## for interaction  - 
    #eff <- attr(terms(model),"term.labels")[3]
    
    #lsm <- lsmeans::.old.lsmeans(model,  as.formula(paste("pairwise ~ ", eff)), 
    #                             lf = TRUE)
    #fixef()
    #split.eff  <-  unlist(strsplit(eff,":"))
    #if( length(split.eff) > 1 ){
      
    pureinter <- getPureInter(lsm.table, anova.table, eff)
    puredifs <- .calcPureDiffs(pureinter) 
      
      
    #all.equal(sum(pure.inter.pairs^2)/(12), sum(puredifs^2)/(12), tol = 1e-4)
    dp <- puredifs / sigma      
    av.dp <- sqrt(sum(dp^2)/(nrow(dp)*(nrow(dp)-1)/2))
    #}
    #else{
     # lsm.eff <- getPureInter(lsm.table, eff)
      
      #dp <- dlsm.table[which(rows %in% eff), 1] / sigma
      #av.dp <- sqrt(sum(dp^2)/length(dp))
      
    #}
    anova.table[eff, "dprimeav"] <- av.dp 
  }
  anova.table
  
  #m.bo <- lme4::lmer(att1 ~ Track*SPL*Car + (1|Assessor) + (1|SPL:Assessor) + 
  #                     (1|Track:SPL:Assessor) + (1|Car:SPL:Assessor), 
  #                   data = sound_data_balanced)
  #lsm <- lsmeans::.old.lsmeans(m.bo,  pairwise ~ Track:SPL:Car, lf = TRUE)
  #lsm[[1]] %*% fixef(model)
  #with(sound_data_balanced, tapply(att1, factor(Track:SPL:Car), mean))
}


## step function for NO MAM
.stepAllAttrNoMAM <- function(new.resp.private.sensmixed, 
                              model = model,
                              reduce.random = reduce.random, 
                              alpha.random = alpha.random, 
                              alpha.fixed = alpha.fixed, 
                              calc_post_hoc = calc_post_hoc){
  assign("new.resp.private.sensmixed", new.resp.private.sensmixed, 
         envir=environment(formula(model)))
  suppressMessages(m <- refit(object=model, newresp = new.resp.private.sensmixed, 
             rename.response = TRUE))
  suppressMessages(st <- step(m, reduce.fixed = FALSE, 
                             reduce.random = reduce.random, 
                             alpha.random = alpha.random, 
                             alpha.fixed = alpha.fixed, 
                             lsmeans.calc = TRUE,
                             difflsmeans.calc = calc_post_hoc))
  
  if(calc_post_hoc)
    st$anova.table <- .calcAvDprime(st$model, st$anova.table, 
                                    st$diffs.lsmeans.table, st$lsmeans.table)    
  st
}

## step function for MAM
.stepAllAttrMAM <- function(attr, product_structure, error_structure,
                            data, Prod_effects, random,
                            reduce.random = reduce.random, 
                            alpha.random = alpha.random, 
                            alpha.fixed = alpha.fixed, 
                            mult.scaling = mult.scaling, 
                            calc_post_hoc = calc_post_hoc){
  model.init <- suppressMessages(createLMERmodel(structure = 
                                  list(product_structure = product_structure, 
                                       error_structure = error_structure), 
                                  data = data, response = attr,
                                  fixed = list(Product = Prod_effects, 
                                             Consumer=NULL),
                                random = random, corr = FALSE, MAM = TRUE,
                                mult.scaling = mult.scaling, 
                                calc_post_hoc = calc_post_hoc))
  model.an <- model.init$model.anova
  model.lsm <- model.init$model.lsmeans
  

  st <- suppressMessages(step(model.an, fixed.calc = FALSE))
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
 
  if(calc_post_hoc){
    if(length(Prod_effects) > 1)
      lsmeans.table <- lsmeans::.old.lsmeans( model.lsm, pairwise ~ prod)
    else 
      lsmeans.table <-  eval(substitute(lsmeans::.old.lsmeans(object=model.lsm, 
                                                         pairwise ~ prod), 
                                        list(prod=as.name(Prod_effects)))) 
    return(list(anova.table=anova.table, rand.table=rand.table,
                lsmeans.table=lsmeans.table)) 
  }
 
  return(list(anova.table=anova.table, rand.table=rand.table)) 
}

calc.cols <- function(x)
{
  if(x<0.001) 
    return("red") 
  if(x<0.01) 
    return("orange") 
  if(x<0.05) 
    return("yellow") 
  return("grey")
}



## UNUSED function
.plotFixedPartsSensmixed <- function(Fval, pvalueF, cex=2, interact.symbol){
    #x11()
    #plot.new()
   # layout(matrix(c(rep(2,2),3,rep(2,2),3,rep(2,2),3, rep(1,2),3), 4, 3, 
   #               byrow = TRUE))
    layout(matrix(c(rep(2,2),3,rep(2,2),3,rep(2,2),3, rep(1,2),3), 4, 3, 
                  byrow = TRUE), 
           heights=c(0.4, 1 , 1.4), widths = c(2,2,4.3))
    
    #Fval <- resSensMixed$fixed$Fval
    #pvalueF <- resSensMixed$fixed$pvalueF
    
    #### plots for F value
    cex.gr <- cex
    names.fixed <- rownames(Fval)# [inds.fixed]
    if(!interact.symbol==":"){        
      names.fixed <- sapply(names.fixed, change.inter.symbol, interact.symbol)
    }
    
    ylim <- c(0, max(sqrt(Fval)) + 0.5)
    
    names.fixed.effs <- LETTERS[1:nrow(Fval)]
    names.fixed.effs.legend <- paste(names.fixed.effs, collapse="")
    #plot(x=bp1[1,], y=rep(1,15), type="n", axes=F, xlab="", ylab="")
    plot.new()
    #else if (is.matrix(FChi.fvalue))
    #  ylim <- c(0, max(apply(FChi.fvalue,2,sum)))
    for(i in 1:ncol(Fval))
    {
      
      fvals <- matrix(0, nrow(Fval), ncol(Fval))
      fvals[,i] <- sqrt(Fval[,i])
      
      
      if(i == ncol(Fval))
      {
        if(length(colnames(Fval)) > 10)
          cex.names <- cex - 0.2
        else
          cex.names <- cex
        bp <- barplot(fvals, col= unlist(lapply(pvalueF[,i], calc.cols)), 
                      ylim=ylim, las=2, main=expression(paste("Barplot for ",
                                                              sqrt(F), 
                                                              " values"
                      )), names.arg =
                        colnames(Fval), las=2, cex.names=cex.names, beside=TRUE, 
                      add=TRUE, xpd=FALSE, cex.main=cex, cex.axis=cex)       
        #text(x=bp, y=rep(0.5, ncol(Fval)), names.fixed.effs, font=1, cex=cex)
        if(sqrt(max(Fval)) > 8)
          text(x=bp, y=sqrt(Fval)+0.3, names.fixed.effs, font=1, cex=cex)
        else
          text(x=bp, y=sqrt(Fval)+0.15, names.fixed.effs, font=1, cex=cex)
        
        plot.new()        
        legend("right", names.fixed, pch=names.fixed.effs.legend,  
               bty="n", pt.lwd=cex, pt.cex=cex, text.font=1, cex=cex)
        legend("topright", c("ns","p < 0.05", "p < 0.01", "p < 0.001"), pch=15, 
               col=c("grey","yellow","orange","red"), title="Significance", 
               bty="n", cex=cex, text.font=1)
      }
      else{
        if(i==1)
          barplot(fvals, col=unlist(lapply(pvalueF[,i], calc.cols)), 
                  ylim=ylim, axes=FALSE, las=2, cex.names=cex, 
                  beside=TRUE)    
        else
          barplot(fvals, col=unlist(lapply(pvalueF[,i], calc.cols)), 
                  ylim=ylim, axes=FALSE, las=2, cex.names=cex, 
                  beside=TRUE, add=TRUE)          
      }       
      if(i < ncol(Fval))
        par(new=TRUE)
    }

} 

change.inter.symbol <- function(x, interact.symbol){
  if(grepl(":", x)){
    symb.loc <- substring.location(x, ":")
    spl.effs <- strsplit(x,":")[[1]]
    x <- paste(spl.effs, collapse=interact.symbol)
    return(x)
  }
  x
}

.changeOutput <- function(vals, pvals, isRand){
  colnames.out <- rownames(vals)
  names <- colnames(vals)
  tr <- vector("list", length(colnames.out))
  
  for(i in 1:length(colnames.out)){       
    tr[[i]] <- createTexreg(
      coef.names = names, se=vals[i,],
      coef = vals[i,],
      pvalues = pvals[i,], isRand=isRand)
  }
    
  names(tr) <- colnames.out
  return(tr)
}
# 
# .convertOutputToMatrix <- function(result){
#   resSensMixed$random
#   Chi <- matrix(0, nrow = length(result), ncol = nrow(result[[1]]))
# }

.plotBars <- function(val, pval, title, plotLegend = TRUE, plotLetters = TRUE, 
                      reduceNames = TRUE, cex = 2, cex.main = 2, 
                      ylim = ylim, 
                      names.effs = NULL, 
                      names.effs.legend = NULL){
  if(reduceNames)
    names.arg <- sapply(colnames(val), 
                        function(x) paste(substring(x,1,7),"..", sep=""))
  else
    names.arg <- colnames(val)
  if(plotLegend)
    cex.names <- cex
  else
    cex.names <- cex - 0.7
  for(i in 1:ncol(val))
  {
    vals <- matrix(0, nrow(val), ncol(pval))
    vals[,i] <- val[,i]      
    if(i == ncol(val)){
      if(length(colnames(val)) > 10)
        cex.names <- cex.names - 0.2
      else
        cex.names <- cex.names
      bp <- barplot(vals, col = unlist(lapply(pval[,i], calc.cols)), 
                    ylim=ylim, las = 2, 
                    main = title,
                    names.arg = names.arg, las = 2,
                    cex.names = cex.names, beside = TRUE, add = TRUE, 
                    xpd = FALSE, cex.main = cex.main, cex.axis = cex - 0.3) 
      
      if(plotLetters){
        if(max(val) > 8)
          text(x = bp, y = val + 0.3, names.effs, font = 1, cex = cex - 0.3)
        else
          text(x = bp, y = val + 0.15, names.effs, font = 1, cex = cex - 0.1)
      }
      
      #text(x=bp, y=rep(0.5, ncol(Chi)), names.rand.effs, font=1, cex=cex)
      
      if(plotLegend){
        plot.new()        
        legend("right", rownames(val), pch = names.effs.legend,  
               bty="n", pt.lwd = cex, pt.cex = cex, text.font = 1, 
               cex = cex - 0.4)
        legend("topright", c("ns","p < 0.05", "p < 0.01", "p < 0.001"), pch = 15, 
               col = c("grey","yellow","orange","red"), title = "Significance", 
               bty = "n", cex = cex, text.font = 1)
      }      
      
    }
    else{
      if(i==1)
        barplot(vals,col = unlist(lapply(pval[,i], calc.cols)),
                ylim = ylim, axes = FALSE, las = 2, cex.names = cex, 
                cex.axis = cex, beside = TRUE)    
      else
        barplot(vals,col = unlist(lapply(pval[,i], calc.cols)), 
                ylim = ylim, axes = FALSE, las = 2, cex.names = cex, 
                cex.axis = cex, beside=TRUE, add=TRUE)         
      
    }
    if(i < ncol(val))
      par(new = TRUE)
  }
}

.plotSensMixed <- function(val, pval, title, mult = FALSE, sep = FALSE,
                           cex = 2,                           
                           interact.symbol = ":"){
  ylim <- c(0, max(val) + 0.5)
  
  ## change the interaction symbol
  if(!interact.symbol == ":")      
    rownames(pval) <- rownames(val) <-  sapply(rownames(val), change.inter.symbol, 
                           interact.symbol)  
 
  
  ## multiple plots
  if(mult){
    reduceNames <- TRUE
    neff <- nrow(val)
    if(sep){
      layout(matrix(c(rep(2,2),3,rep(2,2),3,rep(2,2),3, rep(1,2),3), 4, 3, 
                    byrow = TRUE), 
             heights=c(0.4, 1 , 1.4), widths = c(2,2,4.3))
      reduceNames <- FALSE
    }else{
      if(neff < 2)
        layout( matrix(1:2, 1, 2, byrow=TRUE)) 
      else if(neff < 4)
        layout(matrix(1:4, 2, 2, byrow=TRUE))            
      else if(neff < 5)
        layout(cbind(matrix(1:4, 2, 2, byrow=TRUE), 5:6),
               heights=c(1, 1)) 
      else if(neff < 7)
        layout(cbind(matrix(1:6, 3, 2, byrow=TRUE), 7:9))   
      else if(neff < 10)
        layout(cbind(matrix(1:9, 3, 3, byrow=TRUE), 10:12))  
    }
    
    for(eff in rownames(pval)){
      if(sep){
        plot.new()
        .plotBars(val[eff, , drop=FALSE], pval[eff, , drop=FALSE], 
                  title = eff, plotLegend = FALSE, 
                  plotLetters = FALSE, reduceNames = reduceNames, cex = cex, 
                  cex.main = cex - 0.7, ylim = ylim)
        plot.new()
        legend("right", c("ns","p < 0.05", "p < 0.01", "p < 0.001"), pch=15, 
               col=c("grey","yellow","orange","red"), title="Significance", 
               bty="n", cex = cex - 0.5, text.font=1)
      }
      else{
        .plotBars(val[eff, , drop=FALSE], pval[eff, , drop=FALSE], 
                  title = eff, plotLegend = FALSE, 
                  plotLetters = FALSE, reduceNames = reduceNames, cex = cex, 
                  cex.main = cex - 0.7, ylim = ylim)
      }      
     }
     if(!sep){
       plot.new()
       legend("right", c("ns","p < 0.05", "p < 0.01", "p < 0.001"), pch=15, 
              col=c("grey","yellow","orange","red"), title="Significance", 
              bty="n", cex = cex - 0.5, text.font=1)
     }        
    }else{
        layout(matrix(c(rep(2,2),3,rep(2,2),3,rep(2,2),3, rep(1,2),3), 4, 3, 
                      byrow = TRUE), 
               heights=c(0.4, 1 , 1.4), widths = c(2,2,4.3))
        names.effs <- LETTERS[1:nrow(val)]
        names.effs.legend <- paste(names.effs, collapse="")
        
        plot.new()
        .plotBars(val, pval, title = title, 
                  plotLegend = TRUE, plotLetters = TRUE, reduceNames = FALSE, 
                  cex = cex, cex.main = cex - 0.2, ylim = ylim,
                  names.effs = names.effs, 
                  names.effs.legend = names.effs.legend)
      }
}

.changeConsmixedOutputForDoc <- function(table, name.pval){  
  table[, name.pval] <- gsub("<", "&lt ", table[, name.pval])
  table  
}

## output for the sensmixed
.createDocOutputSensmixed <- function(x, file = NA, bold = FALSE, append = TRUE){
  colnames.out.rand <- rownames(x$rand$Chi)
  names <- colnames(x$rand$Chi)
  tr_rand <- vector("list", length(colnames.out.rand))
  
  for(i in 1:length(colnames.out.rand)){       
    tr_rand[[i]] <- createTexreg(
      coef.names = names, se=x$rand$Chi[i,],
      coef = x$rand$Chi[i,],
      pvalues = x$rand$pvalueChi[i,], isRand=TRUE    
    )     
  } 
  
  
  ## output for the fixed effects
  colnames.out.fixed <- rownames(x$fixed$Fval)
  names <- colnames(x$fixed$Fval)
  tr <- vector("list", length(colnames.out.fixed))
  
  for(i in 1:length(colnames.out.fixed)){       
    tr[[i]] <- createTexreg(
      coef.names = names, se=x$fixed$Fval[i,],
      coef = x$fixed$Fval[i,],
      pvalues = x$fixed$pvalueF[i,],
      isRand=FALSE
    )     
  }
  
  if("scaling" %in% names(x)){
    ## output for the scaling  effects if presented
    colnames.out.scaling <- rownames(x$scaling$FScaling)
    names <- colnames(x$scaling$FScaling)
    tr_scal <- vector("list", length(colnames.out.scaling))
    
    for(i in 1:length(colnames.out.scaling)){       
      tr_scal[[i]] <- createTexreg(
        coef.names = names, se=x$scaling$FScaling[i,],
        coef = x$scaling$FScaling[i,],
        pvalues = x$scaling$pvalueScaling[i,],
        isRand=FALSE
      )     
    }
    regres <- list(lrand = tr_rand, lfixed = tr, lscale = tr_scal)
  }
  else 
    regres <- list(lrand = tr_rand, lfixed = tr)
  
  #  wdGet()
  #  funny<-function(){
  #    c <- plot(x, mult = TRUE)
  #    print(c)
  #  }
  #  wdPlot(plotfun=funny,method="bitmap", height = 10 , width = 10)
  #  
  
  if(bold)
    stars <- numeric(0)
  else
    stars <- c(0.001, 
               0.01, 0.05)

  htmlreg(regres, 
          file = file, inline.css = FALSE, 
          doctype = FALSE, html.tag = FALSE, head.tag = FALSE, 
          body.tag = FALSE,
          custom.model.names =list(
            custom.model.names.rand = colnames.out.rand, 
            custom.model.names.fixed = colnames.out.fixed), 
          caption = list(
            caption.rand="Likelihood ration test for the random effects",
            caption.fixed="F-test for the fixed effects"), bold=bold,
          stars=stars, append = append)
  
  
  if(!is.null(x$post_hoc)){
    if("scaling" %in% names(x))
      name.pval <- "p.value"
    else
      name.pval <- "p-value"
    #names(x$post_hoc)
    sink(file = file, append = append)
    for(i in 1:length(x$post_hoc)){
      x$post_hoc[[i]][,  name.pval] <- format.pval(x$post_hoc[[i]][, name.pval],
                                                  digits=3, eps=1e-3)
      x$post_hoc[[i]] <- .changeConsmixedOutputForDoc(x$post_hoc[[i]],  name.pval)
      if("scaling" %in% names(x))
        xt.posthoc <- xtable(x$post_hoc[[i]], align="lccccc",
                           display=c("s","f", "f", "d", "f", "s"))
      else
        xt.posthoc <- xtable(x$post_hoc[[i]], align="lccccccc",
                             display=c("s","f", "f", "d", "f", "f", "f", "s"))
      caption(xt.posthoc) <- 
        paste("Post-hoc for the attribute ", names(x$post_hoc)[i])
      print(xt.posthoc, caption.placement="top", table.placement="H",
            sanitize.text.function=function(x){x}, size="\\small", type = "html")
    }    
    sink()
  }
  
}

## output for the consmixed
.createDocOutputConsmixed <- function(x, file = NA, bold = FALSE, append = TRUE){
  sink(file = file, append = append)
  
  ## tests for the random effects
  x$rand.table[, "p.value"] <- format.pval(x$rand.table[,"p.value"],
                                           digits=3, eps=1e-3)
  x$rand.table <- .changeConsmixedOutputForDoc(x$rand.table, "p.value")
  if("elim.num" %in% colnames(x$rand.table))
    xt.rand <- xtable(x$rand.table, align="lcccc", 
                      display=c("s","f","d","s","s"))
  else
    xt.rand <- xtable(x$rand.table, align="lccc", 
                      display=c("s","f","d","s"))
  caption(xt.rand) <- "Likelihood ratio tests for the random-effects
  and their order of elimination"
  print(xt.rand, caption.placement="top", table.placement="H",
        sanitize.text.function=function(x){x}, size="\\small", type = "html")
  
  ## tests for the fixed effects
  x$anova.table[, "Pr(>F)"] <- format.pval(x$anova.table[,"Pr(>F)"],
                                           digits=3, eps=1e-3)
  x$anova.table <- .changeConsmixedOutputForDoc(x$anova.table, "Pr(>F)")
  if("elim.num" %in% colnames(x$anova.table)) 
    xt.anova <- xtable(x$anova.table, align="lccccccc",
                       display=c("s","f", "f", "d", "f", "f", "s", "s"))     
  else
    xt.anova <- xtable(x$anova.table, align="lcccccc",
                       display=c("s","f", "f", "d", "f", "f","s"))
  caption(xt.anova) <- 
    "F-tests for the fixed-effects and their order of elimination"
  
  
  print(xt.anova, caption.placement="top", table.placement="H",
        sanitize.text.function=function(x){x}, size="\\small", type = "html")
  
  ## post hoc output
  x$diffs.lsmeans.table[, "p-value"] <- 
    format.pval(x$diffs.lsmeans.table[,"p-value"], digits=3, eps=1e-3)
  x$diffs.lsmeans.table <- 
    .changeConsmixedOutputForDoc(x$diffs.lsmeans.table, "p-value")    
  xt.lsmeans <- xtable(x$diffs.lsmeans.table, align="lccccccc",
                       display=c("s","f", "f", "f", "f", "f","f", "s"))
  caption(xt.lsmeans) <- 
    "Differences of Least Squares Means"
  print(xt.lsmeans, caption.placement="top", table.placement="H",
        sanitize.text.function=function(x){x}, size="\\small", type = "html")
  sink()
}