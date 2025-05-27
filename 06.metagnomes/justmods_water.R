library("glmmTMB")
library("bbmle")

setwd("~/nicolagk@hawaii.edu - Google Drive/My Drive/Bamboo_mesos/Bamboo_mesos/06.metagnomes")

df.water.long <- readRDS("df.water.long.rds")
#df.larv.long <- readRDS("df.larv.long.rds")

str(df.water.long)
df.water.long$rowid <- as.factor(df.water.long$rowid)
df.water.long$Microbe_treatment <- as.factor(df.water.long$Microbe_treatment)
df.water.long$Food_type <- as.factor(df.water.long$Food_type)
#df.water.long$food.microbes <- as.factor(df.water.long$food.microbes)
df.water.long$Mesocosm_id <- as.factor(df.water.long$Mesocosm_id)
str(df.water.long)

water.mod.nb1 <- glmmTMB(value~offset(log(sum))+Microbe_treatment+Food_type+(1|variable)+(1|variable:Microbe_treatment)+(1|variable:Food_type),family=nbinom1,data=df.water.long)
summary(water.mod.nb1)

water.mod.nb2 <- glmmTMB(value~offset(log(sum))+Microbe_treatment+Food_type+(1|variable)+(1|variable:Microbe_treatment)+(1|variable:Food_type),family=nbinom2,data=df.water.long)
summary(water.mod.nb2)

AICtab(water.mod.nb1,water.mod.nb2)

saveRDS(water.mod.nb1,file="water.mod.nb1.full.rds")
saveRDS(water.mod.nb2,file="water.mod.nb2.full.rds")

#water.mod.nb1 <- readRDS("water.mod.nb1.full.rds")

water.mod.nb1.nomicr <- update(water.mod.nb1,.~.-(1|variable:Microbe_treatment))

#anova(water.mod.nb1,water.mod.nb1.nomic)

water.mod.nb1.nofood <- update(water.mod.nb1,.~.-(1|variable:Food_type))

#anova(water.mod.nb1,water.mod.nb1.nomic)

saveRDS(water.mod.nb1.nomicr,file="water.mod.nb1.nomicrobes.rds")
saveRDS(water.mod.nb1.nofood,file="water.mod.nb1.nofood.rds")


