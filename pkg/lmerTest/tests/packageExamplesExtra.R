require(lmerTest)

## from merModLmerTest
m <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject),
          data = sleepstudy)

anova(m, type=1)

anova(m, ddf="Kenward-Roger")

## from lmerTest
m <- lmer(Informed.liking ~ Gender+Information+Product +(1|Consumer), data=ham)

anova(m, ddf = "Kenward-Roger")

## from anova methods
m.ham <- lmer(Informed.liking ~ Product*Information*Gender 
              + (1|Consumer), data = ham)

anova(m.ham, type = 1)

anova(m.ham, ddf = "lme4")

fm2 <- lmer(Preference ~ sens2 + I(sens1^2)  +
              (1+sens2|Consumer), data=carrots)

## from lmer
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)

anova(fm1, ddf="Kenward-Roger")

# anova table the same as of class merMod
anova(fm1, ddf="lme4")

anova(fm1, fm2)
