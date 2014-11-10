require(lmerTest)

m <- lmer(Coloursaturation ~ TVset*Picture+
            (1|Assessor)+(1|Assessor:TVset), data=TVbo)

f1 <- function(x){
  m1 <- refit(object=m, newresp=x, 
                 rename.response = TRUE)
  m1
}

f2 <- function(x){
  m2 <- as(refit(object=m, newresp=x, 
           rename.response = TRUE), "merModLmerTest")
  m2
}

f4 <- function(x){
 
  assign("x", x, env=environment(formula(m)))
  rf <- refit(object=m, newresp=x, 
              rename.response = TRUE)
  step(rf)
}

tools::assertError(update(f1(TVbo$Colourbalance)))
tools::assertError(step(f2(TVbo$Colourbalance)))
lapply(TVbo[, 7, drop=FALSE], f4)

## after the assignment the error disappears
## Why?! seems like x becomes attached... 
## is it OK to do like that within a package?
update(f1(TVbo$Colourbalance))
step(f2(TVbo$Colourbalance))
