if(packageVersion("lme4") <= "0.999999-3")
{
  merLmerTest <- setClass("merLmerTest", representation(t.pval="numeric") , contains = "mer")
  summary.merLmerTest <- setClass("summary.merLmerTest", contains = "summary.mer")
}
if(packageVersion("lme4") > "0.999999-3")
{
  merModLmerTest <- setClass("merModLmerTest",  contains = c("merMod", "lmerMod"))
}
  


