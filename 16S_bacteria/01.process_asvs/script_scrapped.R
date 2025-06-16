## ASVs left after quality filtering to subset from fasta - didn't do yet

```{r}
ps.exp1 <- subset_samples(ps.clean,Mesocosm_type=="Experiment")
ps.exp1
ps.exp2 <- subset_samples(ps.exp1,Exp_day!=20)
ps.exp2
ps.exp3 <- subset_samples(ps.exp2,Microbe_treatment!="AHK")
ps.exp3
ps.exp <- prune_taxa(taxa_sums(ps.exp3)!=0,ps.exp3)
ps.exp
otu.exp <- data.frame(ps.exp@otu_table)

otu.exp.list <- colnames(otu.exp)
#write.table(otu.exp.list, file="otus.experiment.clean.txt", append = FALSE, sep = "/n", row.names = FALSE, col.names = FALSE,quote=FALSE)

tax.exp <- data.frame(ps.exp@tax_table)
#write.csv(tax.exp,file="otus.experiment.taxa.csv")
```

# Rarefy

Doesn't make sense to rarefy a bunch of these samples

```{r, eval=F}
setwd("~/Library/CloudStorage/GoogleDrive-nicolagk@hawaii.edu/My Drive/Bamboo_mesos/Bamboo_mesos/01.process_asvs")
ps.clean.decontam <- readRDS("bamboomesos_ps.lulu.clean.decontam.rds")
ps.clean.decontam #514 taxa, 169 samples

sample_data(ps.clean.decontam)$lib_size_clean <- sample_sums(ps.clean.decontam)

samdf.clean <- data.frame(ps.clean.decontam@sam_data)
#write.csv(samdf.clean,file="samdf.clean.csv")

samdf.low <- subset(samdf.clean,lib_size_clean<20000)
samdf.low

ggplot(data=samdf.low, aes(x=Longer_name, y=lib_size_clean, color=Mesocosm_type))+
  geom_point()+
  theme(axis.text.x=element_text(angle=45,hjust=1))
##all the lowest ones that aren't controls are in the E. coli only treatment which makes sense

df <- data.frame(ps.clean.decontam@otu_table) # Put sample_data into a ggplot-friendly data.frame
df$libsize <- rowSums(df)

df2 <- df[order(df$libsize),]

df2$Index <- seq(nrow(df2))
##all of them:
ggplot(data=df2, aes(x=Index, y=libsize))+
  geom_point()

head(df2$libsize,n=25) 
tail(df2$libsize,n=50)

##higher ones
ggplot(data=df2, aes(x=Index, y=libsize))+
  geom_point()+
  ylim(120000,270000)+
  xlim(130,173)
#there's definitely a natural break between 130950 and 144273
#10% of that would be 13,095
#or between 166362 and 195031
#so 10% of that would be 16k

##lower ones
ggplot(data=df2, aes(x=Index, y=libsize))+
  geom_point()+
  ylim(0,50000)+
  xlim(0,50)

set.seed(2939)

#ps.rare <- rarefy_even_depth(ps.clean.decontam,sample.size=10999,replace=FALSE)
#saveRDS(ps.rare,"bamboomesos_ps.lulu.clean.decontam.rare11k.rds")
##7 samples removed, 6 OTUs removed

ps.rare <- rarefy_even_depth(ps.clean.decontam,sample.size=6950,replace=FALSE)
#saveRDS(ps.rare,"bamboomesos_ps.lulu.clean.decontam.rare6.9k.rds")
##12 samples removed, 14 OTUs removed

##removing the same samples from the other ps object for continuity with comparisons
#ps.rare <- readRDS("bamboomesos_ps.lulu.clean.decontam.rare11k.rds")
#ps.rare #156 samples, 500 taxa

ps.clean.decontam

samdf.clean <- data.frame(ps.clean.decontam@sam_data)

#samdf.low <- subset(samdf.clean,lib_size_clean<10999)
samdf.low <- subset(samdf.clean,lib_size_clean<6950)

samdf.low
rownames(samdf.low)

otu.clean <- data.frame(ps.clean.decontam@otu_table)

otu.less <- subset(otu.clean,!rownames(otu.clean)%in% c(rownames(samdf.low)))
ps.clean.decontam.copy <- ps.clean.decontam
ps.clean.decontam.copy@otu_table <- otu_table(otu.less,taxa_are_rows = F)
ps.clean.decontam.copy

ps.less1 <- phyloseq(otu_table(otu.less, taxa_are_rows=FALSE), 
                     sample_data(samdf.clean), 
                     tax_table(as.matrix(taxa)))
ps.less1

ps.less <- prune_taxa(taxa_sums(ps.less1)>0,ps.less1)
ps.less #500 taxa, 157 samples

#saveRDS(ps.less,"bamboomesos_ps.lulu.clean.decontam.less.rds")
```



## Distance from center for metagenomics

```{r}
library("dplyr")

df.segs <- segs %>%
  group_by(all.trts) %>%
  mutate(distance_from_center = sqrt((NMDS1 - oNMDS1)^2 + (NMDS2 - oNMDS2)^2))

df.segs

df.outliers <- df.segs %>%
  group_by(all.trts) %>%
  slice(which.max(distance_from_center))

##checking the outliers are true visually
ggplot(df.gg.exp.nmds.horn, aes(x=NMDS1, y=NMDS2,shape=sample.food,color=Microbe_treatment,label=Short_label))+
  scale_shape_manual(values=c(17,15,24,22),name="Sample type - diet",labels=c("Larvae - bamboo","Larvae - standard","Water - bamboo","Water - standard"))+
  theme_cowplot()+
  geom_segment(data=segs, mapping = aes(xend = oNMDS1, yend = oNMDS2),alpha=0.5)+ # spiders
  geom_point(size=2,fill="white")+
  scale_color_manual(values=new.colors,name="Microbes",labels=c("Lab","Mos.","Env.","Mos.+Env.","Mos.+Env.+Lab"))+
  annotate("text",x=-1.8,y=-1.9,label="2D stress = 0.15")+
  geom_text(vjust=0,hjust=0)

write.csv(df.outliers,file="comp.outliers.csv")
```