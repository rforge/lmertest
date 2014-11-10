require(lmerTest)

modelCarrots <- lmer(Preference ~  sens2*sens1*Homesize*Age
                                    + (1 | product) + (1 + sens1 + sens2 | Consumer),
                                     data=carrots)

## the results for the rand function differ from the step
## because of the update function - changing the contrasts to contr.SAS
rnd <- rand(modelCarrots)


stopifnot(all.equal(rnd$rand.table[,"p.value"], c(0.0000289113, 0.4852719557,
                                                  0.0186101913)))

stp <- step(modelCarrots)


stopifnot(all.equal(stp$rand.table[,"p.value"], c(8.828516e-01, 3.130212e-05,
                                                  1.632085e-02), tol=1e-4))