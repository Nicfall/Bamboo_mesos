## Larvae

str(df.larv.long)
df.larv.long$rowid <- as.factor(df.larv.long$rowid)
df.larv.long$Microbe_treatment <- as.factor(df.larv.long$Microbe_treatment)
df.larv.long$Food_type <- as.factor(df.larv.long$Food_type)
#df.larv.long$food.microbes <- as.factor(df.larv.long$food.microbes)
df.larv.long$Mesocosm_id <- as.factor(df.larv.long$Mesocosm_id)
str(df.larv.long)

##full model
larv.mod.nb1 <- glmmTMB(value~offset(log(sum))+Microbe_treatment*Food_type+(1|variable)+(1|variable:Microbe_treatment)+(1|variable:Food_type)+(1|variable:Microbe_treatment:Food_type),family=nbinom1,dispformula=~Microbe_treatment*Food_type,data=df.larv.long)
summary(larv.mod.nb1)

larv.mod.nb2 <- glmmTMB(value~offset(log(sum))+Microbe_treatment*Food_type+(1|variable)+(1|variable:Microbe_treatment)+(1|variable:Food_type)+(1|variable:Microbe_treatment:Food_type),family=nbinom2,dispformula=~Microbe_treatment*Food_type,data=df.larv.long)
summary(larv.mod.nb2)

AICtab(larv.mod.nb1,larv.mod.nb2)

larv.mod.nb1.resid <- simulateResiduals(fittedModel = larv.mod.nb1, plot = T)
larv.mod.nb2.resid <- simulateResiduals(fittedModel = larv.mod.nb2, plot = T)

##dispersal formula checks
larv.mod.nb1.disnoint <- glmmTMB(value~offset(log(sum))+Microbe_treatment*Food_type+(1|variable)+(1|variable:Microbe_treatment)+(1|variable:Food_type)+(1|variable:Microbe_treatment:Food_type),family=nbinom1,dispformula=~Microbe_treatment+Food_type,data=df.larv.long)

anova(larv.mod.nb1.disnoint,larv.mod.nb1) #sig*

##removing 3 way interaction
larv.mod.nb1.nofm <- update(larv.mod.nb1,.~.-(1|variable:Microbe_treatment:Food_type))
larv.mod.nb1.nofm
anova(larv.mod.nb1,larv.mod.nb1.nofm) #ns

##removing just food
larv.mod.nb1.nofood <- update(larv.mod.nb1,.~.-(1|variable:Food_type))
anova(larv.mod.nb1,larv.mod.nb1.nofood) #super sig***

##removing just microbes
larv.mod.nb1.nomicr <- update(larv.mod.nb1,.~.-(1|variable:Microbe_treatment))
anova(larv.mod.nb1,larv.mod.nb1.nomicr) #super sig***

##checking out nbinom2
larv.mod.nb2.nofm <- update(larv.mod.nb2,.~.-(1|variable:Microbe_treatment:Food_type))
anova(larv.mod.nb2,larv.mod.nb2.nofm)
##super sig

###save conditional modes
larv.mod.ranef <- as.data.frame(ranef(larv.mod.nb1))
#write.csv(larv.mod.ranef, "larv.mod.ranef.csv")
