require(SensMixed)

tools::assertError(res2 <- sensmixed(c("Noise", "Elasticeffect"),
                  Prod_effects = c("TVset"), replication="Repeat", 
                  individual="Assessor", data=TVbo, parallel = FALSE))