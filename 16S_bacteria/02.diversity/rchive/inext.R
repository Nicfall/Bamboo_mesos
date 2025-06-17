library("iNEXT")

setwd("~/nicolagk@hawaii.edu - Google Drive/My Drive/Bamboo_mesos/Bamboo_mesos/02.diversity")

sum.lar <- readRDS("otutab.larvae.sum.rds")
sum.wat <- readRDS("otutab.water.sum.rds")

inext.lar <- iNEXT(sum.lar,datatype="abundance")
inext.wat <- iNEXT(sum.wat,datatype="abundance")

saveRDS(inext.lar,file="inext.larvae.rds")
saveRDS(inext.wat,file="inext.water.rds")