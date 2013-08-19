merLmerTest <- setClass("merLmerTest", representation(t.pval="numeric") , contains = c("merMod", "lmerMod"))
#summary.merLmerTest <- setClass("summary.merLmerTest", contains = "summary.merMod")