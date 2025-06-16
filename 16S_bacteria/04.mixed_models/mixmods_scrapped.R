ps.exp <- subset_samples(ps.clean,Mesocosm_type=="Experiment"&Microbe_treatment!="AHK"&Exp_day!=20&Microbe_treatment!="ECO"&Microbe_treatment!="ALL")
ps.exp #60 samples

##water samples
ps.exp.wtr <- subset_samples(ps.exp,Sample_type=="Water")
ps.exp.wtr #30 samples

##larvae samples
ps.exp.lar <- subset_samples(ps.exp,Sample_type=="Larvae")
ps.exp.lar #30 samples

##water
ps.exp.wtr.trim <- filter_taxa(ps.exp.wtr, function(x) sum(x > 1) > (0.20*length(x)), TRUE)
ps.exp.wtr.trim ##62 taxa
tail(sort(sample_sums(ps.exp.wtr.trim),decreasing=T))

##larvae
#ps.exp.lar.trim <- filter_taxa(ps.exp.lar, function(x) sum(x > 1) > (0.20*length(x)), TRUE)
#ps.exp.lar.trim #44 taxa 
##with all and eco

#ps.exp.all.lar.trim <- filter_taxa(ps.exp.all.lar, function(x) sum(x > 1) > (0.15*length(x)), TRUE)

##water
otu.exp.wtr <- data.frame(ps.exp.wtr@otu_table)
sum(otu.exp.wtr == 0) / (ncol(otu.exp.wtr) * nrow(otu.exp.wtr)) * 100
#93.19066%

otu.exp.wtr.trim <- data.frame(ps.exp.wtr.trim@otu_table)
sum(otu.exp.wtr.trim == 0) / (ncol(otu.exp.wtr.trim) * nrow(otu.exp.wtr.trim)) * 100
#54.94624%

##larvae
otu.exp.lar <- data.frame(ps.exp.lar@otu_table)
sum(otu.exp.lar == 0) / (ncol(otu.exp.lar) * nrow(otu.exp.lar)) * 100
#95.32425%

otu.exp.lar.trim <- data.frame(ps.exp.lar.trim@otu_table)
sum(otu.exp.lar.trim == 0) / (ncol(otu.exp.lar.trim) * nrow(otu.exp.lar.trim)) * 100
#61.81818%
### Which ASVs to subset from fasta

otus.lar.trim <- colnames(otu.exp.lar.trim)
otus.wtr.trim <- colnames(otu.exp.wtr.trim)

#write.table(otus.lar.trim, file="otus.larvae.trim.txt", append = FALSE, sep = "/n", row.names = FALSE, col.names = FALSE,quote=FALSE)
#write.table(otus.wtr.trim, file="otus.water.trim.txt", append = FALSE, sep = "/n", row.names = FALSE, col.names = FALSE,quote=FALSE)

### Just EMO, MMO, and MEM 

```{r}
otu.wtr <- data.frame(ps.exp.wtr.trim@otu_table)
otu.lar <- data.frame(ps.exp.lar.trim@otu_table)

#total read count colums
otu.wtr$sum <- rowSums(otu.wtr)
otu.lar$sum <- rowSums(otu.lar)

otu.wtr$Short_label <- row.names(otu.wtr)
otu.lar$Short_label <- row.names(otu.lar)

#bring in the metadata again
samdf.wtr <- data.frame(ps.exp.wtr.trim@sam_data)
samdf.lar <- data.frame(ps.exp.lar.trim@sam_data)

##add extra metadata
#meso.data <- read.csv("../bamboomesos_metadata.csv")

#samdf.wtr.more <- merge(samdf.wtr,meso.data)
#samdf.lar.more <- merge(samdf.lar,meso.data)

#samdf.wtr$food.microbes <- paste0(samdf.wtr$Food_type,"_",samdf.wtr$Microbe_treatment)
#samdf.lar$food.microbes <- paste0(samdf.lar$Food_type,"_",samdf.lar$Microbe_treatment)
# 
# #join the data sets by sample name
df.water <- plyr::join(otu.wtr, samdf.wtr, by="Short_label", type="left")
df.larv <- plyr::join(otu.lar, samdf.lar, by="Short_label", type="left")
colnames(df.water)

##getting rid of some irrelevant metadat
colnames(df.water)
df.water1 <- df.water %>%
  dplyr::select(-Short_label,-Longer_name,-Sample_type,-Num_larvae,-Date_collected,-Exp_day,-Mesocosm_treatment,-Mesocosm_type,-Treatment_notes,-Notes,-Raw_reads,-Num_larvae_d00,-Num_larvae_d08,-Num_larvae_d34,-Num_larvae_d40,-Total_pupae,-Proportion_pupated,-Biomass_day30,-lib_size_clean,-is.neg)
colnames(df.water1)

df.water.long <- melt(df.water1, id.vars=c("sum", "Mesocosm_id", "Microbe_treatment","Food_type", "Biomass_day40"))

df.water.long$rowid <- 1:nrow(df.water.long)

df.larv1 <- df.larv %>%
  dplyr::select(-Short_label,-Longer_name,-Sample_type,-Num_larvae,-Date_collected,-Exp_day,-Mesocosm_treatment,-Mesocosm_type,-Treatment_notes,-Notes,-Raw_reads,-Num_larvae_d00,-Num_larvae_d08,-Num_larvae_d34,-Num_larvae_d40,-Total_pupae,-Proportion_pupated,-Biomass_day30,-lib_size_clean,-is.neg)
colnames(df.larv1)

#df.larv2 <- df.larv1[df.larv1$sum!=0,]
#df.larv2$sum==0

df.larv.long <- melt(df.larv1, id.vars=c("sum", "Mesocosm_id", "Microbe_treatment","Food_type", "Biomass_day40"))

df.larv.long$rowid <- 1:nrow(df.larv.long)

##save
#saveRDS(df.water.long,file="df.water.long.rds")
#saveRDS(df.larv.long,file="df.larv.long.rds")
```

df.water.long <- readRDS("df.water.long.rds")
df.larv.long <- readRDS("df.larv.long.rds")

## Water

```{r}
str(df.water.long)
df.water.long$rowid <- as.factor(df.water.long$rowid)
df.water.long$Microbe_treatment <- as.factor(df.water.long$Microbe_treatment)
df.water.long$Food_type <- as.factor(df.water.long$Food_type)
#df.water.long$food.microbes <- as.factor(df.water.long$food.microbes)
df.water.long$Mesocosm_id <- as.factor(df.water.long$Mesocosm_id)
str(df.water.long)

##full model
water.mod.nb1 <- glmmTMB(value~offset(log(sum))+Microbe_treatment*Food_type+(1|variable)+(1|variable:Microbe_treatment)+(1|variable:Food_type)+(1|variable:Microbe_treatment:Food_type),family=nbinom1,dispformula=~variable+Microbe_treatment+Food_type,data=df.water.long)
summary(water.mod.nb1)

water.mod.nb2 <- glmmTMB(value~offset(log(sum))+Microbe_treatment*Food_type+(1|variable)+(1|variable:Microbe_treatment)+(1|variable:Food_type)+(1|variable:Microbe_treatment:Food_type),family=nbinom2,dispformula=~variable+Microbe_treatment+Food_type,data=df.water.long)
summary(water.mod.nb2)

AICtab(water.mod.nb1,water.mod.nb2)

water.mod.nb1.resid <- simulateResiduals(fittedModel = water.mod.nb1, plot = T)
water.mod.nb2.resid <- simulateResiduals(fittedModel = water.mod.nb2, plot = T)

##dispersal formula checks
water.mod.nb1.disnovar <- glmmTMB(value~offset(log(sum))+Microbe_treatment*Food_type+(1|variable)+(1|variable:Microbe_treatment)+(1|variable:Food_type)+(1|variable:Microbe_treatment:Food_type),family=nbinom1,dispformula=~Microbe_treatment+Food_type,data=df.water.long)

anova(water.mod.nb1.disnovar,water.mod.nb1) #sig***

##removing 3 way interaction
water.mod.nb1.nofm <- update(water.mod.nb1,.~.-(1|variable:Microbe_treatment:Food_type))
anova(water.mod.nb1,water.mod.nb1.nofm)
##super sig 0.000242

water.mod.nb2.nofm <- update(water.mod.nb2,.~.-(1|variable:Microbe_treatment:Food_type))
anova(water.mod.nb2,water.mod.nb2.nofm)
##also super sig

###save conditional modes
water.mod.ranef <- as.data.frame(ranef(water.mod.nb1))
#write.csv(water.mod.ranef, "water.mod.ranef.csv")
#saveRDS(water.mod.nb1,file="water.mod.nb1.rds")
```

### Nbinom1 vs. nbinom2 things

Following this wonderful resource, on page 16

Ben Bolker, Mollie Brooks, Beth Gardner, Cleridy Lennert, Mihoko Minami, October 23, 2012, Owls example: a zero-inflated, generalized linear mixed model for count data. https://groups.nceas.ucsb.edu/non-linear-modeling/projects/owls/WRITEUP/owls.pdf

```{r}
mvtab <- ddply(df.water.long,
               .(variable:Microbe_treatment:Food_type),
               summarise,
               callmean=mean(value),
               callvar=var(value))

q1 <- qplot(callmean,callvar,data=mvtab)
print(q1+
        ## linear (quasi-Poisson/NB1) fit
        geom_smooth(method="lm",formula=y~x-1,colour="pink")+
        ## smooth (loess)
        geom_smooth(colour="red")+
        ## semi-quadratic (NB2/LNP)
        geom_smooth(method="lm",formula=y~I(x^2)+offset(x)-1,colour="purple")+
        ## Poisson (v=m)
        geom_abline(a=0,b=1,lty=2))
```

## Larvae

```{r}
str(df.larv.long)
df.larv.long$rowid <- as.factor(df.larv.long$rowid)
df.larv.long$Microbe_treatment <- as.factor(df.larv.long$Microbe_treatment)
df.larv.long$Food_type <- as.factor(df.larv.long$Food_type)
#df.larv.long$food.microbes <- as.factor(df.larv.long$food.microbes)
df.larv.long$Mesocosm_id <- as.factor(df.larv.long$Mesocosm_id)
str(df.larv.long)

##binomial
larv.mod.bin <- glmmTMB(cbind(value, sum-value)~Microbe_treatment*Food_type+(1|variable)+(1|variable:Microbe_treatment)+(1|variable:Food_type)+(1|variable:Microbe_treatment:Food_type)+(1|rowid),family=binomial,data=df.larv.long)
summary(larv.mod.bin)

##full model
larv.mod.nb1 <- glmmTMB(value~offset(log(sum))+Microbe_treatment*Food_type+(1|variable)+(1|variable:Microbe_treatment)+(1|variable:Food_type)+(1|variable:Microbe_treatment:Food_type),disp=~variable+Food_type+Microbe_treatment,family=nbinom1,data=df.larv.long)
summary(larv.mod.nb1)
#saveRDS(larv.mod.nb1,file="larv.mod.nb1.rds")

larv.mod.nb2 <- glmmTMB(value~offset(log(sum))+Microbe_treatment*Food_type+(1|variable)+(1|variable:Microbe_treatment)+(1|variable:Food_type)+(1|variable:Microbe_treatment:Food_type),disp=~variable+Food_type+Microbe_treatment,family=nbinom2,data=df.larv.long)
summary(larv.mod.nb2)

AICtab(larv.mod.nb1,larv.mod.nb2,larv.mod.bin)

larv.mod.bin.resid <- simulateResiduals(fittedModel = larv.mod.bin, plot = T)
larv.mod.nb1.resid <- simulateResiduals(fittedModel = larv.mod.nb1, plot = T)
larv.mod.nb2.resid <- simulateResiduals(fittedModel = larv.mod.nb2, plot = T)

##dispersal formula checks
larv.mod.nb1.disnovar <- glmmTMB(value~offset(log(sum))+Microbe_treatment*Food_type+(1|variable)+(1|variable:Microbe_treatment)+(1|variable:Food_type)+(1|variable:Microbe_treatment:Food_type),family=nbinom1,dispformula=~Microbe_treatment+Food_type,data=df.larv.long)

anova(larv.mod.nb1.disnovar,larv.mod.nb1) #sig***

##removing 3 way interaction
larv.mod.nb1.nofm <- update(larv.mod.nb1,.~.-(1|variable:Microbe_treatment:Food_type))
summary(larv.mod.nb1.nofm)
anova(larv.mod.nb1,larv.mod.nb1.nofm) #sig **

# ##removing just food
# larv.mod.nb1.nofood <- update(larv.mod.nb1,.~.-(1|variable:Food_type))
# anova(larv.mod.nb1,larv.mod.nb1.nofood) #super sig***
# 
# ##removing just microbes
# larv.mod.nb1.nomicr <- update(larv.mod.nb1,.~.-(1|variable:Microbe_treatment))
# anova(larv.mod.nb1,larv.mod.nb1.nomicr) #super sig***

##checking out nbinom2
larv.mod.nb2.nofm <- update(larv.mod.nb2,.~.-(1|variable:Microbe_treatment:Food_type))
anova(larv.mod.nb2,larv.mod.nb2.nofm)
##sig**

##checking out binomial
larv.mod.bin.nofm <- update(larv.mod.bin,.~.-(1|variable:Microbe_treatment:Food_type))
summary(larv.mod.bin.nofm)
anova(larv.mod.bin,larv.mod.bin.nofm) #sig***

###save conditional modes
larv.mod.ranef <- as.data.frame(ranef(larv.mod.nb1))
#write.csv(larv.mod.ranef, "larv.mod.ranef.csv")

# library("ggstats")
# ggcoef_multicomponents(larv.mod.nb1)
```

### Nbinom1 vs. nbinom2 things

```{r}
##note: some errors/warnings if the treatments aren't factors for some reason
mvtab <- ddply(df.larv.long,
               .(variable:Microbe_treatment:Food_type),
               summarise,
               callmean=mean(value),
               callvar=var(value))

q1 <- qplot(callmean,callvar,data=mvtab)
print(q1+
        ## linear (quasi-Poisson/NB1) fit
        geom_smooth(method="lm",formula=y~x-1,colour="pink")+
        ## smooth (loess)
        geom_smooth(colour="red")+
        ## semi-quadratic (NB2/LNP)
        geom_smooth(method="lm",formula=y~I(x^2)+offset(x)-1,colour="purple")+
        ## Poisson (v=m)
        geom_abline(a=0,b=1,lty=2))#+
#ylim(0,1.75e8)+
#xlim(0,20000)
```


##dispersal formula checks
#larv.mod.nb1.disnovar <- glmmTMB(value~offset(log(sum))+Microbe_treatment*Food_type+(1|variable)+(1|variable:Microbe_treatment)+(1|variable:Food_type)+(1|variable:Microbe_treatment:Food_type),family=nbinom1,dispformula=~Microbe_treatment+Food_type,data=df.all.larv.long)

#anova(larv.mod.nb1.disnovar,larv.mod.nb1) #sig***

# ##removing just food
# water.mod.nb1.nofood <- update(water.mod.nb1,.~.-(1|variable:Food_type))
# anova(water.mod.nb1,water.mod.nb1.nofood) #super sig***
# 
# ##removing just microbes
# water.mod.nb1.nomicr <- update(water.mod.nb1,.~.-(1|variable:Microbe_treatment))
# anova(water.mod.nb1,water.mod.nb1.nomicr) #super sig***

##checking out nbinom2
water.mod.nb2.nofm <- update(water.mod.nb2,.~.-(1|variable:Microbe_treatment:Food_type))
anova(water.mod.nb2,water.mod.nb2.nofm)
##sig**

##checking out binomial
water.mod.bin.nofm <- update(water.mod.bin,.~.-(1|variable:Microbe_treatment:Food_type))
summary(water.mod.bin.nofm)
anova(water.mod.bin,water.mod.bin.nofm) #sig***

###save conditional modes
water.mod.ranef <- as.data.frame(ranef(water.mod.nb1))
#water.mod.ranef <- as.data.frame(ranef(water.mod.nb1))
#write.csv(water.mod.ranef, "water.mod.ranef.all.csv")
#saveRDS(water.mod.nb1,file="water.mod.nb1.all.rds")

# water.mod.ranef <- as.data.frame(ranef(water.mod.nb2))
# write.csv(water.mod.ranef, file="water.mod.ranef.all.NB2.csv")
# 
# water.mod.ranef <- as.data.frame(ranef(water.mod.nb1.nofm))
# write.csv(water.mod.ranef, file="water.mod.ranef.all.no3way.csv")
# library("ggstats")
# ggcoef_multicomponents(water.mod.nb1)