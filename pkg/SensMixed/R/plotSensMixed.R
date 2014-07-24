###################################################################################################################
## Different plots for SensMixed package 
###################################################################################################################
### plot results
plotSensMixed <- function(resSensMixed, mult = FALSE, sep = FALSE, cex = 2, 
                          interact.symbol = ":", isFixed = TRUE, isRand = TRUE)
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
  

  

 Chi <- sqrt(resSensMixed$random$Chi)
 pvalueChi <- resSensMixed$random$pvalueChi
 Fval <- sqrt(resSensMixed$fixed$Fval)
 pvalueF <- resSensMixed$fixed$pvalueF  
 if(mult == FALSE){
   if(isRand)
     .plotSensMixed(Chi, pvalueChi, title =  expression(paste("Barplot for ", 
                                                            sqrt(chi^2))), 
                    mult=FALSE, sep = FALSE, interact.symbol = interact.symbol, 
                    cex = cex)
   if(isFixed)
     .plotSensMixed(Fval, pvalueF, title = expression(paste("Barplot for ",
                                                            sqrt(F), " values")), 
                    mult=FALSE, sep = FALSE, interact.symbol = interact.symbol, 
                    cex = cex)
 }
 else{
   if(isRand)
     .plotSensMixed(Chi, pvalueChi, mult = TRUE, sep = sep,
                    interact.symbol = interact.symbol, cex = cex)
   if(isFixed)
     .plotSensMixed(Fval, pvalueF, mult = TRUE, sep = sep,
                    interact.symbol = interact.symbol, cex = cex)
 }
   
 
}
