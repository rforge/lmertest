###################################################################################################################
## Different plots for SensMixed package 
###################################################################################################################
### plot results
plotSensMixed <- function(resSensMixed, fixed=FALSE, rand=FALSE)
{
  
  
  dens <- function(x)
  {
    if(x < 0.01) 
      return(dens=500)    
    if(x < 0.05) 
      return(dens=100) 
    if(x < 0.1)
      return(dens=50)
    return(dens=10)
  }
  
  col.bars.F <- function(x)
  {
    gr <- gray.colors(3)
    if(x < 0.05)
      return(gr[1])
    if(x < 0.1)
      return(gr[2])
    return(gr[3])
  }
  
  col.bars.Chi <- function(x)
  {
    gr <- gray.colors(2)
    if(x<0.05) 
      return(gr[1]) 
    return(gr[2])
  }
  
  #######################################################################################################
  #######  Plot of random effects. Run and look at results before and then run for the fixed effects plot
  #######################################################################################################
  if(rand || (!rand && !fixed)){
    
    #### plots as in lmerTest  paper
    #par(mar=c(15,4,4,2))
    plot.new()
    #inds.rand <- rownames(resSensMixed$random$Chi)#grep(" | ", rownames(FChi))
    #FChi.chi <- FChi[inds.rand,]
    
    Chi <- sqrt(resSensMixed$random$Chi)
    pvalueChi <- resSensMixed$random$pvalueChi
     
    names.random <- rownames(Chi) #rev(unlist(lapply(rownames(FChi.chi), function(x) substring2(x, 4, nchar(x)-1))))
    if(length(names.random)==1)
      col.bars <- "pink"
    else if(length(names.random)==2)
      col.bars <- c("pink", "blue")
    else
      col.bars <- brewer.pal(nrow(Chi), "Set3")#gray.colors(5)#gray.colors(5, start=0, end=1)
    #if (is.vector(FChi.chi))
      ylim <- c(0,max(Chi))
    #else if (is.matrix(FChi.chi))
    #  ylim <- c(0, max(apply(FChi.chi,2,sum)))
    #ang <- c(0, 45, 90, 120, 160)
    cex.gr <- 0.7
    for(i in 1:ncol(Chi))
    {
      Chisqs <- matrix(0, nrow(Chi), ncol(Chi))
      Chisqs[,i] <- Chi[,i]    
      if(i == ncol(Chi))
        barplot(Chisqs,col=col.bars, ylim=ylim, density=unlist(lapply(pvalueChi[,i],dens)), names.arg = colnames(Chi), las=2,  main=expression(paste("Barplot for ", sqrt(chi^2),  " values of the LRT")), cex.names=cex.gr, beside=TRUE) 
      else
        barplot(Chisqs,col=col.bars, ylim=ylim, density=unlist(lapply(pvalueChi[,i],dens)), axes=FALSE, las=2, cex.names=cex.gr, beside=TRUE)
      if(i < ncol(Chi))
        par(new=TRUE)
    }
    #############################################################################################################################################################################
    legend(locator(1), names.random,  col=col.bars, pch=15,  bty="n", cex=cex.gr)
    legend(locator(1), rev(c("p-value < 0.01", "p-value < 0.05", "p-value < 0.1", "p-value >= 0.1")),  col="black", density=rev(c(500, 100, 50, 10)), bty="n", cex=cex.gr)
  }
  
 
  
  #######################################################################################################
  #######  Plot of fixed effects
  #######################################################################################################
  if(fixed || (!rand && !fixed)){
    plot.new()
    #library(RColorBrewer)
    
    
    Fval <- resSensMixed$fixed$Fval
    pvalueF <- resSensMixed$fixed$pvalueF
      
    #### plots for F value
    #inds.fixed <- (1:nrow(Fval))[-inds.rand]
    #FChi.fvalue <- FChi[inds.fixed,]
    cex.gr <- 0.7
    names.fixed <- rownames(Fval)# [inds.fixed]
    if(length(names.fixed)==1)
      col.bars <- "blue"
    else if(length(names.fixed)==2)
      col.bars <- c("blue", "red")
    else 
      col.bars <- brewer.pal(nrow(Fval), "Set1")
    #par(mar=c(15,4,4,2))
    #col.bars<-gray.colors(4, start=0.7, end=1)
    #if (is.vector(FChi.fvalue))
      ylim <- c(0, max(sqrt(Fval)))
    #else if (is.matrix(FChi.fvalue))
    #  ylim <- c(0, max(apply(FChi.fvalue,2,sum)))
    for(i in 1:ncol(Fval))
    {
        
      fvals <- matrix(0, nrow(Fval), ncol(Fval))
      fvals[,i] <- sqrt(Fval[,i])
  
      if(i == ncol(Fval))
        barplot(fvals, col=col.bars, ylim=ylim, density=unlist(lapply(pvalueF[,i],dens)), names.arg = colnames(Fval), las=2, main=expression(paste("Barplot for ", sqrt(F), " values")), cex.names=cex.gr, beside=TRUE)
      else
        barplot(fvals, col=col.bars, ylim=ylim, density=unlist(lapply(pvalueF[,i],dens)), axes=FALSE, las=2, cex.names=cex.gr, beside=TRUE)
      if(i < ncol(Fval))
        par(new=TRUE)
    }
    legend(locator(1), names.fixed,  col=rev(col.bars), pch=15,  bty="n", cex=cex.gr)
    legend(locator(1), rev(c("p-value < 0.01", "p-value < 0.05", "p-value < 0.1", "p-value >= 0.1")),  col="black", density=rev(c(500, 100, 50, 10)), bty="n", cex=cex.gr)
  }
}
