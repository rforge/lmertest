\name{lmer}
\alias{lmer}
\title{Fit Linear Mixed-Effects Models}
\description{
  Fit a linear mixed model
}

\details{

  This \code{lmer} function is an overloaded function of \code{lmer} of \code{mer} class (or \code{merMod} class for \code{lme4} package version > "0.999999-3") with an additional slot t.pval, which gives p-values for the t test of fixed terms with Satterthwaite approximation for degees of freedom.
  
}
\value{
  An object of class \code{"\linkS4class{merLmerTest}"} (or \code{"\linkS4class{merModLmerTest}"})
}
\seealso{
  \code{\linkS4class{merLmerTest}} class, \code{\linkS4class{merMod}} class, \code{\linkS4class{merModLmerTest}} class
}
\examples{

library(lmerTest)

## linear mixed models
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)

# anova table the same as of class mer but with additional F statistics and  
# p-values calculated based on Satterthwaite's approximations
anova(fm1)

# anova table the same as of class mer (merMod) 
anova(fm1, ddf="lme4")

# gives summary of merLmerTest class. The same as of class summary.mer but with 
# additional p-values calculated based on Satterthwate's approximations
summary(fm1)

}
\keyword{models}
\keyword{methods}
