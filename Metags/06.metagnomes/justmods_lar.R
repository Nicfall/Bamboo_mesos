library("glmmTMB")
library("bbmle")

setwd("~/nicolagk@hawaii.edu - Google Drive/My Drive/Bamboo_mesos/Bamboo_mesos/06.metagnomes")

df.larv.long <- readRDS("df.larv.long.rds")

str(df.larv.long)
df.larv.long$rowid <- as.factor(df.larv.long$rowid)
df.larv.long$Microbe_treatment <- as.factor(df.larv.long$Microbe_treatment)
df.larv.long$Food_type <- as.factor(df.larv.long$Food_type)
#df.larv.long$food.microbes <- as.factor(df.larv.long$food.microbes)
df.larv.long$Mesocosm_id <- as.factor(df.larv.long$Mesocosm_id)
str(df.larv.long)

larv.mod.nb1 <- glmmTMB(value~offset(log(sum))+Microbe_treatment+Food_type+(1|variable)+(1|variable:Microbe_treatment)+(1|variable:Food_type),family=nbinom1,data=df.larv.long)
summary(larv.mod.nb1)

larv.mod.nb2 <- glmmTMB(value~offset(log(sum))+Microbe_treatment+Food_type+(1|variable)+(1|variable:Microbe_treatment)+(1|variable:Food_type),family=nbinom2,data=df.larv.long)
summary(larv.mod.nb2)

AICtab(larv.mod.nb1,larv.mod.nb2)

saveRDS(larv.mod.nb1,file="larv.mod.nb1.full.rds")
saveRDS(larv.mod.nb2,file="larv.mod.nb2.full.rds")

#larv.mod.nb1 <- readRDS("larv.mod.nb1.full.rds")

larv.mod.nb1.nomicr <- update(larv.mod.nb1,.~.-(1|variable:Microbe_treatment))

#anova(larv.mod.nb1,larv.mod.nb1.nomicr)

larv.mod.nb1.nofood <- update(larv.mod.nb1,.~.-(1|variable:Food_type))

#anova(larv.mod.nb1,larv.mod.nb1.nomic)

saveRDS(larv.mod.nb1.nomicr,file="larv.mod.nb1.nomicrobes.rds")
saveRDS(larv.mod.nb1.nofood,file="larv.mod.nb1.nofood.rds")

