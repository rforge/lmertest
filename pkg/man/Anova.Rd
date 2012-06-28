\name{Anova}
\alias{Anova}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Calculates Analysis of Variance Table for mixed linear model of object of class lmer.
}
\description{
 Performs type 3 hypothesis (definition from SAS proc mixed) for fixed effects part of a linear mixed effects model.
 The F statistics is calculated based on Satterthwaite's approximation for denominator degrees of freedom. 
}
\usage{
Anova(model, ddf="Satterthwaite", type=3, method.grad="Richardson", ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
    
     linear mixed effects model (lmer object from lme4 package).
}

 \item{ddf}{
    
     approximation for denominator degrees of freedom. By default Satterthwaite's approximation. ddf=Kenward-Roger calculates Kenward-Roger approximation
}

\item{type}{
   type of hypothesis to be tested. Currently there is only one choice type=3, which means the hypothesis of type 3.
}


\item{method.grad}{
   approximation method for the grad function, which is used in calculation of denominator degrees of freedom. Could be "simple" or "Richardson". the default one is "Richardson".
}

  \item{\dots}{
     other potential arguments.
}
}
\details{
  The definition of type 3 hypothesis is taken from SAS software. The options for the contrasts should be "contr.SAS". Calculation of Kenward-Roger approximation is taken from pbkrtest package.
}
\value{
 Produces a data frame, which contains F tests  for model's fixed terms.
}
\references{
Tests of Hypotheses in Fixed-Effects Linear Models Copyright 1978 by SAS Institute Inc., Cary, NC, USA
}
\author{
 Alexandra Kuznetsova, Per Bruun Brockhoff, Rune Haubo Bojesen Christensen
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{step}}, \code{\link{rand}}, \code{\link{lsmeans}}
}
\examples{

#import lme4 package and lmerTest package
library(lme4)
library(lmerTest)

#import ham data from lmerTest package
data(ham)

#specify lmer model
m1 <- lmer(Informed.liking ~ Gender*Information*Product + (1|Product:Consumer) +(1|Information:Consumer) + (1|Consumer), data=ham)

# anova table produced by lmerTest package
Anova(m1)

# anova table with "simple" computational method for calculating gradient
# for denominator degrees of freedom
Anova(m1, method.grad="simple")

# using Kenward-Roger approximation for denominator degrees of freedom
Anova(m1, ddf="Kenward-Roger")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
