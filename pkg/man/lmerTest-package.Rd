\name{lmerTest-package}
\alias{lmerTest-package}
\alias{lmerTest}
\docType{package}
\title{
 The package performs different kinds of tests on lmer objects(of lme4 package), such as F tests of type 3 hypothesis for the fixed part, 
log-likelihood ratio tests for the random part, least squares means and differences of least squares means for the factors of the fixed part with corresponding plots. The package also provides with a function step, that preforms backward elimination of non-significant effects, starting from the random effects, and then fixed ones. 
}
\description{
 The package provides Anova function, that gives data frame similar to what gives lmer function(of lme4 package) but with p-values calculated from F statistics of type 3 hypothesis. There are two options for denominator degrees of freedom of F statistics: "Satterthwaite" and "Kenward-Roger". The calculation of Kenward-Roger's approximation is taken from pbkrtest package, the calculation of Satterthwaite's approximation is taken from SAS proc mixed theory (see reference).The type 3 hypothesis (marginal) is calculated according to SAS theory(SAS Institute Inc., 1978). The tests on random effects are performed using log-likelihood ratio tests with one degree of freedom, which means, testing one effect in a time.
}

\details{
\tabular{ll}{
Package: \tab lmerTest\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-01-10\cr
License: \tab GPL\cr
}

The calculation of statistics for the fixed part was developed according to SAS Proc Mixed Theory(see reference). Therefore the contrasts matrix for the fixed effects should be contr.SAS. The option for the contrasts is automatically changed, and the model is updated. 
In the elimination of the random part the models with REML=TRUE are compared, therefor if the model with REML=FALSE is specified, then it is updated automatically to REML=TRUE

}
\author{
Alexandra Kuznetsova <alku@imm.dtu.dk>, Per Bruun Brockhoff, Rune Haubo Bojesen Christensen 
}
\references{
Tests of Hypotheses in Fixed-Effects Linear Models Copyright 1978 by SAS Institute Inc., Cary, NC, USA.
J.H. Goodnight. General Linear Model Procedure S.A.S. Institute, Inc.
Approximations to Distributions of Test Statistics in Complex Mixed Linear Models Using SAS Proc MIXED
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation ~~
%%~~ directory ~~
%%\keyword{ package }
%%\seealso{
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%%}
\examples{
options(contrasts=c(unordered="contr.SAS", ordered="contr.poly"))

#import lme4 package and lmerTest package
library(lme4)
library(lmerTest)

#import ham data from lmerTest package
data(ham)


m <- lmer(Informed.liking ~ Gender*Information*Product +(1|Consumer), data=ham)
step(m)

data(TVbo)
m <- lmer(Cutting ~ TVset*Picture+(1|Assessor)
+ (1|Assessor:Picture) + (1|Assessor:TVset), data=TVbo)

step(m, test.effs="TVset")
Anova(m)
rand(m)
lsmeans(m)


data(carrots)
m2 <- lmer(Preference~sens2+sens1+BITTER+Homesize+
sens2:Homesize+(1+sens2|Consumer), data=carrots)

rand(m2)
an<-Anova(m2)
step(m2)

}
