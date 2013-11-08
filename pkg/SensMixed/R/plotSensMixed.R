###################################################################################################################
## Different plots for SensMixed package 
###################################################################################################################
### plot results
plotSensMixed <- function(resSensMixed)
{
  
  library(RColorBrewer)
  FChi <- resSensMixed$FChi
  pvalue <- resSensMixed$pvalue
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
  #### plots as in lmerTest  paper
  #par(mar=c(15,4,4,2))
 plot.new()
  inds.rand <- grep(" | ", rownames(FChi))
  FChi.chi <- FChi[inds.rand,]
  names.random <- rownames(FChi)[inds.rand]#rev(unlist(lapply(rownames(FChi.chi), function(x) substring2(x, 4, nchar(x)-1))))
  if(length(inds.rand)==1)
    col.bars <- "pink"
  else if(length(inds.rand)==2)
    col.bars <- c("pink", "blue")
  else
    col.bars <- brewer.pal(nrow(FChi.chi), "Set3")#gray.colors(5)#gray.colors(5, start=0, end=1)
  #if (is.vector(FChi.chi))
    ylim <- c(0,max(FChi.chi))
  #else if (is.matrix(FChi.chi))
  #  ylim <- c(0, max(apply(FChi.chi,2,sum)))
  #ang <- c(0, 45, 90, 120, 160)
  cex.gr <- 0.7
  for(i in 1:ncol(FChi))
  {
    if(is.vector(FChi.chi))
    {
      Chisqs <- rep(0,ncol(FChi))
      Chisqs[i] <- FChi.chi[i]
    }      
    else if (is.matrix(FChi.chi))
    {
      Chisqs <- matrix(0,nrow(FChi.chi),ncol(FChi.chi))
      Chisqs[,i] <- FChi.chi[,i]
    }
    if(i == ncol(FChi))
      barplot(Chisqs,col=col.bars, ylim=ylim, density=unlist(lapply(pvalue[inds.rand,i],dens)), names.arg = colnames(FChi.chi), las=2,  main="Barplot for Chi square values of the LRT", cex.names=cex.gr, beside=TRUE) 
    else
      barplot(Chisqs,col=col.bars, ylim=ylim, density=unlist(lapply(pvalue[inds.rand,i],dens)), axes=FALSE, las=2, cex.names=cex.gr, beside=TRUE)
    if(i < ncol(FChi))
      par(new=TRUE)
  }
  #############################################################################################################################################################################
  legend(locator(1), names.random,  col=rev(col.bars), pch=15,  bty="n", cex=cex.gr)
  legend(locator(1), rev(c("p-value < 0.01", "p-value < 0.05", "p-value < 0.1", "p-value >= 0.1")),  col="black", density=rev(c(500, 100, 50, 10)), bty="n", cex=cex.gr)
  
  
 plot.new()
  
  #######################################################################################################
  #######  Plot of fixed effects
  #######################################################################################################
  
  #### plots for F value
  inds.fixed <- (1:nrow(FChi))[-inds.rand]
  FChi.fvalue <- FChi[inds.fixed,]
  cex.gr <- 0.7
  names.fixed <- rownames(FChi)[inds.fixed]
  if(length(inds.fixed)==1)
    col.bars <- "blue"
  else if(length(inds.fixed)==2)
    col.bars <- c("blue", "red")
  else 
    col.bars <- brewer.pal(nrow(FChi.fvalue), "Set1")
  #par(mar=c(15,4,4,2))
  #col.bars<-gray.colors(4, start=0.7, end=1)
  #if (is.vector(FChi.fvalue))
    ylim <- c(0, max(FChi.fvalue))
  #else if (is.matrix(FChi.fvalue))
  #  ylim <- c(0, max(apply(FChi.fvalue,2,sum)))
  for(i in 1:ncol(FChi))
  {
    if(is.vector(FChi.fvalue))
    {
      Chisqs <- rep(0,ncol(FChi))
      Chisqs[i] <- FChi.fvalue[i]
    }      
    else if (is.matrix(FChi.fvalue))
    {
      Chisqs <- matrix(0,nrow(FChi.fvalue),ncol(FChi.fvalue))
      Chisqs[,i] <- FChi.fvalue[,i]
    }
    if(i == ncol(FChi))
      barplot(Chisqs, col=col.bars, ylim=ylim, density=unlist(lapply(pvalue[inds.fixed,i],dens)), names.arg = colnames(FChi), las=2, main="Barplot for F values", cex.names=cex.gr, beside=TRUE)
    else
      barplot(Chisqs,col=col.bars, ylim=ylim, density=unlist(lapply(pvalue[inds.fixed,i],dens)), axes=FALSE, las=2, cex.names=cex.gr, beside=TRUE)
    if(i < ncol(FChi))
      par(new=TRUE)
  }
  legend(locator(1), names.fixed,  col=rev(col.bars), pch=15,  bty="n", cex=cex.gr)
  legend(locator(1), rev(c("p-value < 0.01", "p-value < 0.05", "p-value < 0.1", "p-value >= 0.1")),  col="black", density=rev(c(500, 100, 50, 10)), bty="n", cex=cex.gr)
  
}
