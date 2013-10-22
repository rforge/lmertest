sensOldmixed <- function(attributes, Prod_effects, replication, individual, data, product_structure = 3, error_structure="No_Rep", alpha.random = 0.1, alpha.fixed = 0.05, ...)
{  
  result <- sensOldmixedFun(attributes, Prod_effects, replication, individual, data, product_structure = product_structure, error_structure=error_structure, alpha.random = alpha.random, alpha.fixed = alpha.fixed)
  class(result) <- "sensOldmixed"
  result
}

print.sensOldmixed <- function(x, ...)
{
  cat("\n matrix of F and Chi square values:\n")
  print(round(x$FChi,2))
  cat("\n matrix of p-values:\n")
  res <- apply(x$pvalue,2, format.pval, digits=2)
  rownames(res) <- rownames(x$pvalue)
  print(res)
} 

sensmixed <- function(attributes, Prod_effects, replication, individual, data, product_structure = 3, error_structure="No_Rep", alpha.random = 0.1, alpha.fixed = 0.05, ...)
{  
  result <- sensmixedFun(attributes, Prod_effects, replication, individual, data, product_structure = product_structure, error_structure=error_structure, alpha.random = alpha.random, alpha.fixed = alpha.fixed)
  class(result) <- "sensmixed"
  result
}

print.sensmixed <- function(x, ...)
{
  cat("\n matrix of F and Chi square values:\n")
  print(round(x$FChi,2))
  cat("\n matrix of p-values:\n")
  res <- apply(x$pvalue,2, format.pval, digits=2)
  rownames(res) <- rownames(x$pvalue)
  print(res)
}  

plot.sensmixed <- function(x, ...)
{
  plotSensMixed(x)
}

consmixed <- function(response, Prod_effects, Cons_effects=NULL, Cons, data, structure = 3, alpha.random = 0.1, alpha.fixed = 0.05, ...)
{  
  result <- consmixedFun(response=response, Prod_effects, Cons_effects=Cons_effects, Cons, data, structure = structure, alpha.random = alpha.random, alpha.fixed = alpha.fixed)
  class(result)<-"consmixed"
  result
}

print.consmixed <- function(x, ...)
{
  if(!is.null(x$rand.table))
  {
    cat("\nRandom effects:\n") 
    x$rand.table[,"p.value"] <- format.pval(x$rand.table[,"p.value"], digits=4, eps=1e-7)
    x$rand.table[,"Chi.sq"] <- round(x$rand.table[,"Chi.sq"],2)
    print(x$rand.table)
    #printCoefmat(x$rand.table, digits=3 , dig.tst=1  ,tst.ind=which(colnames(x$rand.table)=="Chi.DF"), P.values=TRUE, has.Pvalue=TRUE, na.print = "KEEP")
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
      x$anova.table[,"Pr(>F)"] <- format.pval(x$anova.table[,"Pr(>F)"], digits=4, eps=1e-7)
      x$anova.table[,c("Sum Sq","Mean Sq", "F.value")] <- round(x$anova.table[,c("Sum Sq","Mean Sq", "F.value")],4)
      x$anova.table[,"DenDF"] <- round(x$anova.table[,"DenDF"],2)
      print(x$anova.table) 
      #printCoefmat(x$anova.table, dig.tst=3, tst.ind=3, cs.ind=3, digits=3 ,P.values = TRUE, has.Pvalue=TRUE)
      if(!is.null(x$lsmeans.table))
      {
        cat("\nLeast squares means:\n")
        printCoefmat(x$lsmeans.table, dig.tst=3 ,tst.ind=c(1:(which(colnames(x$lsmeans.table)=="Estimate")-1),which(colnames(x$lsmeans.table)=="DF")), digits=3 ,P.values = TRUE, has.Pvalue=TRUE)
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