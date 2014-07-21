\name{plot.consmixed}
\alias{plot.consmixed}

\title{
plots the post-hoc for the consmixed object
}
\description{
plots the least squares means and differences of least squares means 
together with the confidence intervals for the fixed effects
}
\usage{
\method{plot}{consmixed} (x, ...)
}

\arguments{
  \item{x}{
 object of class consmixed
}
 \item{\dots}{
     other potential arguments.
}
}

\value{
returns NULL
}

\author{
Alexandra Kuznetsova 
}


\examples{
res <- consmixed(response="Informed.liking", 
Prod_effects= c("Product","Information"), 
Cons_effects=c("Gender","Age"), Cons = "Consumer", data =ham, structure=1)

plot(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
