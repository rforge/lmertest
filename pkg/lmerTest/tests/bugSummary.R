require(lmerTest)

load(system.file("testdata","bugSummaryData.RData", package="lmerTest"))

lmer3 <- lmer(cog_Mid ~ (1|item) + (1 + vowel|speaker) + 
                Language*vowel*sex, 
              data=data1.frame, REML = FALSE, na.action = na.omit)

## no more errors if the C code from lme4 is used
#tools::assertCondition(s <- summary(lmer3))

lmer4 <- lmer(cog_Mid ~ (1|item) + (vowel - 1|speaker) + Language*vowel*sex, 
              data = data1.frame)

## anova does not work in this case A is not positiv definite
## SAS seems to work
anova(lmer4)

