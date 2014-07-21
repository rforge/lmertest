require(SensMixed)

res <- sensmixed(c("Coloursaturation", "Colourbalance"),
                                  Prod_effects=c("TVset"), 
                                  individual="Assessor", 
                 data=TVbo, MAM=TRUE, reduce.random=FALSE)
res_MAM <- sensmixed(c("Coloursaturation", "Colourbalance"),
                                    Prod_effects=c("TVset"), 
                                    individual="Assessor", data=TVbo, 
                     MAM_PER=TRUE)
TOL <- 1e-3
stopifnot(all.equal(res_MAM[[3]][, , 1][2:3,"F"], c(res$fixed$Fval[, 1],
                                                    res$scaling$FScaling[, 1]), 
                    tol=TOL, check.attributes = FALSE))
