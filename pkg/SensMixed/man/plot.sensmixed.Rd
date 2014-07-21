\name{plot}
\alias{plot.sensmixed}

\title{
function creates plots for the sensmixed object 
}
\description{
function creates barplots for the  square roots of F statistics and 
square roots of chi square values for all attributes
}
\usage{
\method{plot}{sensmixed} (x, mult = FALSE, cex = 2,  interact.symbol = ":", ...)
}

\arguments{
  \item{x}{
object of class sensmixed
}
 \item{mult}{
logical. Should multiple plots be plotted, that is barplots 
for each effect in a separate plot
}
\item{cex}{
The magnification to be used
}
\item{interact.symbol}{
The symbol to be used for the interaction effects
}
 \item{\dots}{
     other potential arguments.
}
}

\value{
 NULL is returned
}

\author{
Alexandra Kuznetsova
}




\examples{
res <- sensmixed(c("Coloursaturation", "Colourbalance"),
                  Prod_effects=c("TVset"), 
                  individual="Assessor", data=TVbo, MAM=TRUE, 
                  reduce.random=FALSE)
plot(res)
plot(res, mult = TRUE)
plot(res, interact.symbol = " x ")
}
