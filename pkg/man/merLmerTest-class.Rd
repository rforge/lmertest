\name{merLmerTest-class}
\title{Mixed Model Representations}
\docType{class}
% Classes
\alias{merLmerTest-class}
%\alias{coef.mer-class}
\alias{summary.merLmerTest-class}
%\alias{ranef.mer-class}
% Methods
\alias{anova,merLmerTest-method}

\alias{summary,merLmerTest-method}
\alias{summary,summary.merLmerTest-method}

\description{The \code{merLmerTest} \emph{contains} mer class of lme4 package and overloads \code{anova} and \code{summary} functions.  The \code{summary.merLmerTest} class represents the summary of objects of class \code{merLmerTest}.


  The \code{"summary.merLmerTest"} class \emph{contains} the \code{"summary.mer"}
}

\section{Objects from the Class}{Objects can be created via the
  \code{\link{lmer}}  functions.
} 
 
\section{Methods}{
anova \code{signature(object = "merLmerTest")}: returns the data frame with F test of type 3 hypothesis for each term in formula and corresponding p-values. For multiple arguments, model comparison statistics. The definition of type 3 hypothesis comes from SAS theory proc mixed.
      \code{\link{anova}}
}
 
\seealso{
  \code{\link{lmer}()}
}
\examples{
(m <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject),
             data = sleepstudy))
print(m, digits = 10, corr = FALSE) # more precision; no 
anova(m)
anova(m, ddf="Kenward-Roger")
anova(m, ddf="lme4")
}
\keyword{classes}
