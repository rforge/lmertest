\name{sensmixed}
\alias{sensmixed}

\title{Automated model selection process for each attribute of sensory data}

\description{
Constructs the biggest possible models for the selected attributes and reduces them to the best by principle of parcimony models. First elimination of random effects is performed following by elimination of fixed effects.
The LRT test is used for testing random terms, F-type hypothesis test is used for testing fixed terms
}
\usage{
sensmixed(attributes, Prod_effects, replication, 
individual, data, product_structure = 3, error_structure="No_Rep",
alpha.random = 0.1, alpha.fixed = 0.05, ...)
}
\arguments{
  \item{attributes}{
    
     vector with names of sensory attributes
}

\item{Prod_effects}{
    
     names of the variables related to the product
}

\item{replication}{
    
     names of the replication column in the data, if present
}

\item{individual}{
    
    name of the column in the data that represent assessors
}

\item{data}{
    
     data frame (data from sensory studies)
}

\item{product_structure}{
    
    one of the values in c(1,2,3). 1: only main effects will enter the initial biggest model. 2: main effects and 2-way interaction. 3: all main effects and all possible interaction 
}

\item{error_structure}{
    
    one of the values in c("No_Rep", "2-WAY", "3-WAY"). "No_Rep" and "2-WAY" - assessor effect and all possible interactions between assessor and Product_effects. "3-WAY" - assessor and replicate effect and interaction between them and interaction between them and Product_effects
}

\item{alpha.random}{
    
     significance level for elimination of the random part (for LRT test)
}

\item{alpha.fixed}{
    
     significance level for elimination of the fixed part (for F test)
}

  \item{\dots}{
     other potential arguments.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{FCHi}{matrix with Chi square values from LRT test and F values form F-type test for the selected attributes}
  \item{pvalue}{matrix withp-values for all effects for the selected attributes}
  }
  
\author{
Alexandra Kuznetsova, Per Bruun Brockhoff, Rune Haubo Bojesen Christensen
}

\examples{


#import SensMixed package
library(SensMixed)

#import TVbo data from lmerTest package
data(TVbo)

#run automated selection process
sensmixed(c("Coloursaturation","Colourbalance"),
Prod_effects=c("TVset", "Picture"), replication="Repeat", 
individual="Assessor", data=TVbo)

\dontrun{
#plot F and Chi square values
#plot(result)
}

\dontrun{
result <- sensmixed(names(TVbo)[5:ncol(TVbo)],
Prod_effects=c("TVset", "Picture"),
replication="Repeat", individual="Assessor", data=TVbo)

plot(result)
}
}