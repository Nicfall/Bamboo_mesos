
# Script ends here

# Reading SQM results into R

## Setup

```{r}
#BiocManager::install("SQMtools")
library("SQMtools")

##3 samples test runs:
#setwd("~/nicolagk@hawaii.edu - Google Drive/My Drive/Bamboo_mesos/metagenomes/SQM_zips_3samps")

##real deal:
setwd("~/nicolagk@hawaii.edu - Google Drive/My Drive/Bamboo_mesos/metagenomes/SQM_zips")
```

## Doing it

```{r}
# Define the directory containing the .zip files
##tests:
#zip_dir <- "~/nicolagk@hawaii.edu - Google Drive/My Drive/Bamboo_mesos/metagenomes/SQM_zips_3samps"

##real deal:
zip_dir <- "~/nicolagk@hawaii.edu - Google Drive/My Drive/Bamboo_mesos/metagenomes/SQM_zips"

# Get a list of all .zip files
zip_files <- list.files(zip_dir, pattern = "\\.zip$", full.names = TRUE)
head(zip_files)

# Loop through each .zip file
for (zip_file in zip_files){
  # Extract the base name (without extension) to use as the object name
  obj_name <- tools::file_path_sans_ext(basename(zip_file))
  # Assign the result of loadSQM() to a variable with the extracted name
  assign(obj_name, loadSQM(zip_file), envir = .GlobalEnv)
}

##selecting my R objects that start with these 2 capital letters
object_names <- ls(pattern = "^[LW]")

for (obj in object_names) {
  # Construct the new variable name
  new_name <- paste0(obj, ".bac")
  # Print the name for debugging
  print(paste("Creating:", new_name))
  # Use get() to retrieve the actual object
  obj_data <- get(obj)
  # Apply the function to the object
  result <- subsetTax(obj_data, "superkingdom", "Bacteria")
  # Assign the result to a new variable with the dynamically created name
  assign(new_name, result, envir = .GlobalEnv)
}

# Get all objects that end with ".bac" (matching our dynamically created ones)
bac_objects <- ls(pattern = "\\.bac$")

# Retrieve the actual objects using mget()
bac_list <- mget(bac_objects)

# Combine them into one using combineSQMlite()
bactab.3samps <- do.call(combineSQMlite, bac_list)



```



# Full tables

```{r}
library("SQMtools")

setwd("~/nicolagk@hawaii.edu - Google Drive/My Drive/Bamboo_mesos/metagenomes")
```

## Bacteria

```{r}
bactab.all <- readRDS("bactab.all.rds")

bactab.ks <- data.frame(bactab.all[["functions"]][["KEGG"]][["abund"]])
bactab.knames <- data.frame(bactab.all[["misc"]][["KEGG_names"]])
##rename the column
colnames(bactab.knames) <- c("KEGG_names")
##remove kegg names that aren't in our data frame
bactab.knames.inc <- bactab.knames[row.names(bactab.knames) %in% row.names(bactab.ks),]

setdiff(row.names(bactab.ks), row.names(bactab.knames))  # Rows in bactab.ks but not in bactab.knames
setdiff(row.names(bactab.knames), row.names(bactab.ks))  # Rows in bactab.knames but not in bactab.ks

write.csv(bactab.ks)

##checking out just one example
lb5.bac.ks <- data.frame(lb5.bac[["functions"]][["KEGG"]][["abund"]])
lb5.bac.knames <- data.frame(lb5.bac[["misc"]][["KEGG_names"]])

setdiff(row.names(lb5.bac.ks), row.names(lb5.bac.knames))  # Rows in bactab.ks but not in bactab.knames
setdiff(row.names(lb5.bac.knames), row.names(lb5.bac.ks))  # Rows in bactab.knames but not in bactab.ks

setdiff(row.names(lb5.bac.ks), row.names(bactab.ks[,"LB5",drop=F]))
setdiff(row.names(bactab.ks[,"LB5",drop=F]), row.names(lb5.bac.ks))

bactab.ks.lb5 <- bactab.ks[!is.na(bactab.ks[,"LB5"]) & bactab.ks[,"LB5"] != 0, ]
lb5.bac.ks.no0 <- lb5.bac.ks[!is.na(lb5.bac.ks[,"LB5", drop=F]) & lb5.bac.ks[,"LB5", drop=F] != 0,,drop=F]

#bactab.ks.lb5 <- bactab.ks[bactab.ks[,"LB5"]!=0,]
setdiff(row.names(bactab.ks.lb5), row.names(lb5.bac.ks.no0))
setdiff(row.names(lb5.bac.ks.no0), row.names(bactab.ks.lb5))

#bactab.ks.lb5[,"LB5"] == lb5.bac.ks.no0
##supposed to check if any don't match
any(bactab.ks.lb5[,"LB5"] != lb5.bac.ks.no0) 
```

I don't know why, but there are lots of gaps in the full bacteria table of the full kegg names, so going to do it this way instead:

```{r}
download.file("https://rest.kegg.jp/list/ko", "kegg_ko_list.txt")

lines <- readLines("kegg_ko_list.txt")
length(lines)
kegg_data <- read.delim("kegg_ko_list.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)

colnames(kegg_data) <- c("KEGG_ID", "Description")

bactab.kegg <- data.frame(KEGG_ID=row.names(bactab.ks))

bactab.kegg.info <- merge(bactab.kegg,kegg_data,by="KEGG_ID")

setdiff(bactab.kegg.info$KEGG_ID, row.names(bactab.ks))
setdiff(row.names(bactab.ks), bactab.kegg.info$KEGG_ID)


setdiff(bactab.kegg$KEGG_ID, row.names(bactab.ks))


# Example: Find the name for "K18069"
kegg_data[kegg_data$KEGG_ID == "ko:K18069", ]

row.names(bactab.ks)

#write.csv(bactab.kegg.info,"bactab.allsamps.kegg.info.csv")
#write.csv(bactab.ks,"bactab.allsamps.counts.csv")

keggpaths.all <- data.frame(paths=bactab.all[["misc"]][["KEGG_paths"]])
library(stringr)

paths.split <- data.frame(str_split_fixed(keggpaths.all$paths, pattern="; ",n=17))
paths.split[56,]
keggpaths.all[56,]

keggpaths.all.split1 <- cbind(keggpaths.all,paths.split)

keggpaths.all.split <- cbind(KEGG_ID = rownames(keggpaths.all.split1), keggpaths.all.split1)

#write.csv(keggpaths.all.split,"bactab.allsamps.keggpaths.csv")

```

```{r}
library("SQMtools")

setwd("~/nicolagk@hawaii.edu - Google Drive/My Drive/Bamboo_mesos/metagenomes/SQM_zips_3samps")

#LB5 <- readRDS("LB5.rds")
#WB9 <- readRDS("WB9.rds")
#LL11 <- readRDS("LL11.rds")

lb5 <- loadSQM("LB5.zip")
wb9 <- loadSQM("WB9.zip")
ll11 <- loadSQM("LL11.zip")

lb5.bac <- subsetTax(lb5, "superkingdom", "Bacteria")
wb9.bac <- subsetTax(wb9, "superkingdom", "Bacteria")
ll11.bac <- subsetTax(ll11, "superkingdom", "Bacteria")

bactab.3samps <- combineSQMlite(lb5.bac,wb9.bac,ll11.bac)
bactab.3samps


lb5.kpaths <- data.frame(paths=lb5.bac[["misc"]][["KEGG_paths"]])
lb5.kinfo <- data.frame(paths=lb5.bac[["misc"]][["KEGG_names"]])

head(lb5.kinfo[!row.names(lb5.kinfo) %in% row.names(lb5.kpaths),,drop=F])
any(row.names(lb5.kinfo) != row.names(lb5.kpaths))

lb5.bac.counts <- lb5.bac[["functions"]][["KEGG"]][["abund"]]

head(lb5.bac.counts[!row.names(lb5.bac.counts) %in% row.names(lb5.kpaths),,drop=F])

pathlist.3samps <- bactab.3samps[["misc"]][["KEGG_paths"]]
pathlist.lb5 <- lb5.bac[["misc"]][["KEGG_paths"]]

head(pathlist.3samps)
head(pathlist.3samps)

head(lb5.kinfo[!pathlist.3samps %in% row.names(lb5.kpaths),,drop=F])

pathlist.lb5.unique <- unique(pathlist.lb5)
pathlist.lb5.names.unique <- unique(names(pathlist.lb5))

```

Doing this on the cluster because it's toooo many entries

```{r}

kegg.paths <- readRDS("~/Downloads/KEGG_paths_list.rds")

# Initialize an empty data frame
unique_paths_df <- data.frame(KO = character(), Pathway = character(), stringsAsFactors = FALSE)

# Loop through each sample in kegg.paths
for (sample in names(kegg.paths)) {
  # Extract the named vector of pathways (KO IDs are names)
  pathways_list <- kegg.paths[[sample]]
  
  # Loop through each KO entry
  for (ko in names(pathways_list)) {
    # Get the pathway string and split into individual pathways
    split_paths <- unlist(strsplit(pathways_list[[ko]], "; ", fixed = TRUE))
    
    # Create a temporary dataframe for this KO
    temp_df <- data.frame(KO = ko, Pathway = unique(split_paths), stringsAsFactors = FALSE)
    
    # Combine with the main dataframe
    unique_paths_df <- rbind(unique_paths_df, temp_df)
  }
}

# Remove duplicates in case the same KO-pathway pair appears multiple times
unique_paths_df <- unique(unique_paths_df)

# View the table
print(unique_paths_df)

# Save to CSV if needed
write.csv(unique_paths_df, "unique_kegg_paths_with_KO.csv", row.names = FALSE)

```

# Back in R

```{r}
library("SQMtools")

setwd("~/nicolagk@hawaii.edu - Google Drive/My Drive/Bamboo_mesos/metagenomes/SQM_zips")

bac.counts <- read.csv("bactab.allsamps.counts.csv",row.names=1)
kegg.paths <- read.csv("unique_kegg_paths_with_KO.csv")
kegg.names <- read.csv("unique_kegg_names_with_KO.csv")
colnames(kegg.names) <- c("KO","names")

##they all match
any(kegg.paths$KO != kegg.names$KO)

bac.kos <- row.names(bac.counts)

kegg.paths.kos <- kegg.paths[kegg.paths$KO %in% bac.kos,]
kegg.names.kos <- kegg.names[kegg.names$KO %in% bac.kos,]

kegg.infos <- merge(kegg.names.kos,kegg.paths.kos,by="KO")

#write.csv(kegg.infos,file="bactab.kegginfos.csv")
```

# Fungi? 

```{r}
library("SQMtools")

setwd("~/nicolagk@hawaii.edu - Google Drive/My Drive/Bamboo_mesos/metagenomes/SQM_zips_3samps")

#LB5 <- readRDS("LB5.rds")
#WB9 <- readRDS("WB9.rds")
#LL11 <- readRDS("LL11.rds")

lb5 <- loadSQM("LB5.zip")
wb9 <- loadSQM("WB9.zip")
ll11 <- loadSQM("LL11.zip")

wb9.fungi <- subsetTax(wb9, "phylum", c("Ascomycota", "Basidiomycota", "Chytridiomycota","Microsporidia", "Mucoromycota", "Zoopagomycota", "Olpidiomycota"))

fungal_phyla <- c("Ascomycota", "Basidiomycota", "Chytridiomycota","Microsporidia", "Mucoromycota", "Zoopagomycota", "Olpidiomycota")

wb9.fungi <- do.call(rbind, lapply(fungal_phyla, function(phy) {
  subsetTax(wb9, "phylum", phy)
}))
lb5.fun1 <- subsetTax(wb9, "phylum", "Ascomycota")
lb5.fun1 <- subsetTax(wb9, "phylum", "Basidiomycota")

fungal_phyla <- c("Ascomycota", "Basidiomycota", "Chytridiomycota","Microsporidia", "Mucoromycota", "Zoopagomycota", "Olpidiomycota")

subsetTaxMulti = function (SQM, rank, taxa,
                           trusted_functions_only = FALSE,
                           ignore_unclassified_functions = FALSE, 
                           rescale_tpm = TRUE,
                           rescale_copy_number = TRUE)
{
  subs = lapply(taxa, FUN=function(tax) subsetTax(SQM, rank, tax, trusted_functions_only, ignore_unclassified_functions, rescale_tpm, rescale_copy_number))
  return(combineSQM(subs, tax_source = "contigs", trusted_functions_only = trusted_functions_only, ignore_unclassified_functions = ignore_unclassified_functions, rescale_tpm = rescale_tpm, rescale_copy_number = rescale_copy_number))
}

wb9.fun <- subsetTaxMulti(wb9,"phylum",fungal_phyla)

#### HERE? ####
LB5 <- loadSQM("LB5.zip")
WB9 <- loadSQM("WB9.zip")
LL11 <- loadSQM("LL11.zip")

# List objects that start with "L" or "W"
object_names <- ls(pattern = "^[LW]")

# Define the list of fungal phyla
fungal_phyla <- c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Microsporidia", "Mucoromycota", "Zoopagomycota", "Olpidiomycota")

# List objects that start with "L" or "W"
object_names <- ls(pattern = "^[LW]")

# Define the list of fungal phyla
fungal_phyla <- c("Ascomycota", "Basidiomycota", "Chytridiomycota", 
                  "Microsporidia", "Mucoromycota", "Zoopagomycota", "Olpidiomycota")

for (obj in object_names) {
  # Construct the new variable name
  new_name <- paste0(obj, ".fun")
  # Retrieve the object
  obj_data <- get(obj)
  # Extract phylum names from the SQM structure
  existing_phyla <- rownames(obj_data[["taxa"]][["phylum"]]$abun)  # Extract row names (phylum names)
  # Find which fungal phyla are present
  valid_fungal_phyla <- intersect(fungal_phyla, existing_phyla)
  # Print debugging information
  print(paste("Processing:", obj, "→", new_name))
  print(paste("Valid phyla found:", paste(valid_fungal_phyla, collapse = ", ")))
  # Only proceed if there are valid fungal phyla
  if (length(valid_fungal_phyla) > 0) {
    result <- subsetTaxMulti(obj_data, "phylum", valid_fungal_phyla)
  } else {
    result <- NULL  # Or handle empty cases differently
    message(paste("No fungal phyla found in", obj))
  }
  
  # Assign result to a new variable
  assign(new_name, result, envir = .GlobalEnv)
}



#lb5.euk <- subsetTax(lb5, "superkingdom", "Eukaryota")
#wb9.euk <- subsetTax(wb9, "superkingdom", "Eukaryota")
#ll11.euk <- subsetTax(ll11, "superkingdom", "Eukaryota")

euktab.3samps <- combineSQMlite(lb5.euk,wb9.euk,ll11.euk)


setwd("/mnt/lustre/koa/koastore/cmaiki_group/nicolagk/bamboomesos_metagnomes/mostsamps/rready/rproc")

euktab <- readRDS("euktab.allsamps.rds")
euk.phyla <- euktab[["taxa"]][["phylum"]][["abun"]]
names(euk.phyla)

tax.phyla <- euktab[["taxa"]][["phylum"]][["abund"]]
row.names(tax.phyla)
write.table(row.names(tax.phyla),file="tax.phyla.txt")

```

On cluster:
  
  ```{r}
library("SQMtools")

setwd("/mnt/lustre/koa/koastore/cmaiki_group/nicolagk/bamboomesos_metagnomes/mostsamps/rready/fun_proc")
zip_dir <- "/mnt/lustre/koa/koastore/cmaiki_group/nicolagk/bamboomesos_metagnomes/mostsamps/rready/fun_proc"


subsetTaxMulti = function (SQM, rank, taxa,
                           trusted_functions_only = FALSE,
                           ignore_unclassified_functions = FALSE, 
                           rescale_tpm = TRUE,
                           rescale_copy_number = TRUE)
{
  subs = lapply(taxa, FUN=function(tax) subsetTax(SQM, rank, tax, trusted_functions_only, ignore_unclassified_functions, rescale_tpm, rescale_copy_number))
  return(combineSQM(subs, tax_source = "contigs", trusted_functions_only = trusted_functions_only, ignore_unclassified_functions = ignore_unclassified_functions, rescale_tpm = rescale_tpm, rescale_copy_number = rescale_copy_number))
}

zip_dir <- "~/nicolagk@hawaii.edu - Google Drive/My Drive/Bamboo_mesos/metagenomes/SQM_zips_3samps"

##real deal:
#zip_dir <- "~/nicolagk@hawaii.edu - Google Drive/My Drive/Bamboo_mesos/metagenomes/SQM_zips"

# Get a list of all .zip files
zip_files <- list.files(zip_dir, pattern = "\\.zip$", full.names = TRUE)
head(zip_files)

# Loop through each .zip file
for (zip_file in zip_files){
  # Extract the base name (without extension) to use as the object name
  obj_name <- tools::file_path_sans_ext(basename(zip_file))
  # Assign the result of loadSQM() to a variable with the extracted name
  assign(obj_name, loadSQM(zip_file), envir = .GlobalEnv)
}

##selecting my R objects that start with these 2 capital letters
object_names <- ls(pattern = "^[LW]")
object_names

# Define the list of fungal phyla
fungal_phyla <- c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Microsporidia", "Mucoromycota", "Zoopagomycota", "Olpidiomycota")

for (obj in object_names) {
  # Construct the new variable name
  new_name <- paste0(obj, ".fun")
  # Retrieve the object
  obj_data <- get(obj)
  # Extract phylum names from the SQM structure
  existing_phyla <- rownames(obj_data[["taxa"]][["phylum"]]$abun)  # Extract row names (phylum names)
  # Find which fungal phyla are present
  valid_fungal_phyla <- intersect(fungal_phyla, existing_phyla)
  # Print debugging information
  print(paste("Processing:", obj, "→", new_name))
  print(paste("Valid phyla found:", paste(valid_fungal_phyla, collapse = ", ")))
  # Only proceed if there are valid fungal phyla
  if (length(valid_fungal_phyla) > 0) {
    result <- subsetTaxMulti(obj_data, "phylum", valid_fungal_phyla)
  } else {
    result <- NULL  # Or handle empty cases differently
    message(paste("No fungal phyla found in", obj))
  }
  
  # Assign result to a new variable
  assign(new_name, result, envir = .GlobalEnv)
}

# Get all objects that end with ".fun" (matching our dynamically created ones)
fun_objects <- ls(pattern = "\\.fun$")

# Retrieve the actual objects using mget()
fun_list <- mget(fun_objects)

# Combine them into one using combineSQMlite()
funtab.allsamps <- do.call(combineSQMlite, fun_list)

saveRDS(funtab.allsamps, file="funtab.allsamps.rds")
```

Worked locally but didn't work on the cluster, so I assume one or more of the samples were troublesome

No fungal phyla found in WL7
No fungal phyla found in WL8

```{r}
funtab1 <- readRDS("../SQM_zips/funtab.allsamps.rds")

funtab <- funtab1[["functions"]][["KEGG"]]$abund

#write.csv(funtab,file="funtab.counts.csv")
```


## Virus functions

```{r}
library("SQMtools")

setwd("/mnt/lustre/koa/koastore/cmaiki_group/nicolagk/bamboomesos_metagnomes/mostsamps/rready/vir_proc")
zip_dir <- "/mnt/lustre/koa/koastore/cmaiki_group/nicolagk/bamboomesos_metagnomes/mostsamps/rready"

# Get a list of all .zip files
zip_files <- list.files(zip_dir, pattern = "\\.zip$", full.names = TRUE)
head(zip_files)

# Loop through each .zip file
for (zip_file in zip_files){
  # Extract the base name (without extension) to use as the object name
  obj_name <- tools::file_path_sans_ext(basename(zip_file))
  # Assign the result of loadSQM() to a variable with the extracted name
  assign(obj_name, loadSQM(zip_file), envir = .GlobalEnv)
}

##selecting my R objects that start with these 2 capital letters
object_names <- ls(pattern = "^[LW]")
object_names

for (obj in object_names) {
  # Construct the new variable name
  new_name <- paste0(obj, ".vir")
  # Print the name for debugging
  print(paste("Creating:", new_name))
  # Use get() to retrieve the actual object
  obj_data <- get(obj)
  # Apply the function to the object
  result <- subsetTax(obj_data, "superkingdom", "Viruses")
  # Assign the result to a new variable with the dynamically created name
  assign(new_name, result, envir = .GlobalEnv)
}

# Get all objects that end with ".vir" (matching our dynamically created ones)
vir_objects <- ls(pattern = "\\.vir$")

# Retrieve the actual objects using mget()
vir_list <- mget(vir_objects)

# Combine them into one using combineSQMlite()
virtab.allsamps <- do.call(combineSQMlite, vir_list)

saveRDS(virtab.allsamps, file="virtab.allsamps.rds")
```

Make job

```{bash, eval=F}
nano rproc_virus.slurm
```

Text within:

```
#!/bin/bash
#SBATCH --job-name=rproc_virus
#SBATCH --partition=shared
##3 day max run time for public partitions, except 4 hour max runtime for the sandbox partition
#SBATCH --time=0-24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --error=%A.err
#SBATCH --output=%A.out ##%A = filled with job id
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --mail-user=nicolagk@hawaii.edu

module load lang/R/4.2.1-foss-2022a
Rscript rproc_virus.R
```

When I ran fungi before, I got an error when there weren't any fungi counts which makes sense, so might have to fix that again

Potential fixing:
  
  ```{r}
library("SQMtools")

setwd("/mnt/lustre/koa/koastore/cmaiki_group/nicolagk/bamboomesos_metagnomes/mostsamps/rready/vir_proc")
zip_dir <- "/mnt/lustre/koa/koastore/cmaiki_group/nicolagk/bamboomesos_metagnomes/mostsamps/rready"

# Get a list of all .zip files
zip_files <- list.files(zip_dir, pattern = "\\.zip$", full.names = TRUE)
head(zip_files)

# Loop through each .zip file
for (zip_file in zip_files){
  # Extract the base name (without extension) to use as the object name
  obj_name <- tools::file_path_sans_ext(basename(zip_file))
  # Assign the result of loadSQM() to a variable with the extracted name
  assign(obj_name, loadSQM(zip_file), envir = .GlobalEnv)
}

## Selecting R objects that start with these 2 capital letters
object_names <- ls(pattern = "^[LW]")
object_names

for (obj in object_names) {
  # Construct the new variable name
  new_name <- paste0(obj, ".vir")
  # Retrieve the object
  obj_data <- get(obj)
  # Check if "Viruses" is present in the superkingdom column
  if (any(rownames(obj_data[["taxa"]][["superkingdom"]]$abun) == "Viruses")) {
    print(paste("Processing:", obj, "→", new_name))
    # Subset for viruses
    result <- subsetTax(obj_data, "superkingdom", "Viruses")
    # Assign result to a new variable
    assign(new_name, result, envir = .GlobalEnv)
  } else {
    message(paste("No viruses found in", obj))
  }
}

# Get all objects that end with ".vir" (matching our dynamically created ones)
vir_objects <- ls(pattern = "\\.vir$")
# Retrieve the actual objects using mget()
vir_list <- mget(vir_objects)

# Combine them into one using combineSQMlite()
virtab.allsamps <- do.call(combineSQMlite, vir_list)

saveRDS(virtab.allsamps, file="virtab.allsamps.rds")

```


```{r}
table(top.wat.bam$Pathway)
head(top.wat.bam$Pathway)

library(dplyr)
library(ggplot2)

# Step 1: Count the number of times each pathway appears
top.wat.bam20 <- top.wat.bam %>%
  count(Pathway) %>%
  arrange(desc(n)) %>%
  slice_head(n = 20)

# Step 2: Reorder for plotting
top.wat.bam20.form <- top.wat.bam20 %>%
  mutate(Pathway = factor(Pathway, levels = rev(Pathway)))

#top.wat.bam20 <- top.wat.bam21.form[top.wat.bam21.form$Pathway!="Not Included in Pathway or Brite; Unclassified: metabolism; Enzymes with EC numbers",]

# Step 3: Plot
ggplot(top.wat.bam20.form, aes(x = Pathway, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 - water, bamboo diet",
       x = "Pathway",
       y = "Count (Number of Times Detected)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 7))

top.wat.sta20 <- top.wat.lar %>%
  count(Pathway) %>%
  arrange(desc(n)) %>%
  slice_head(n = 20)

# Step 2: Reorder for plotting
top.wat.sta20.form <- top.wat.sta20 %>%
  mutate(Pathway = factor(Pathway, levels = rev(Pathway)))

#top.wat.sta20 <- top.wat.sta21.form[top.wat.sta21.form$Pathway!="Not Included in Pathway or Brite; Unclassified: metabolism; Enzymes with EC numbers",]

# Step 3: Plot
ggplot(top.wat.sta20.form, aes(x = Pathway, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 - water, standard diet",
       x = "Pathway",
       y = "Count (Number of Times Detected)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 7))



```



## Top X% - larvae

```{r}
# top.mem.lar <- subset(larv.ranef.all, scale(MEM) > 2.326) #top 1%, one-tailed Z score test 
# top.mmo.lar <- subset(larv.ranef.all, scale(MMO) > 2.326) #top 1%, one-tailed Z score test 
# top.emo.lar <- subset(larv.ranef.all, scale(EMO) > 2.326) #top 1%, one-tailed Z score test 
# 
# top.bam.lar <- subset(larv.ranef.all, scale(BvL) > 1.645) #top 1%, one-tailed Z score test 
# top.lar.lar <- subset(larv.ranef.all, scale(BvL) < -1.645) #top 1%, one-tailed Z score test 
larv.ranef.all$sBvL <- c(scale(larv.ranef.all$BvL))

Lbamboo <- subset(larv.ranef.all, sBvL > 1.28 & MEM > 0) #top 10%, one-tailed Z score test
Lstandard <- subset(larv.ranef.all, sBvL < -1.28 & MMO > 0) #bottom 10%, one-tailed Z score test

#larv.ranef.all$sMMO <- scale(larv.ranef.all$MMO)
#larv.ranef.all$sMEM <- scale(larv.ranef.all$MEM)

#Lstandard.1 <- subset(larv.ranef.all, sBvL < -1.282 & MMO > 0) #bottom 10%, one-tailed Z score test
#Lstandard2.1 <- subset(larv.ranef.all, sBvL < -1.282 & MMO > 0  & MEM < 0) #bottom 10%, one-tailed Z score test
#Lstandard$kegg.name == Lstandard.1$kegg.name

#Lstandard <- subset(larv.ranef.all,  sBvL < -1.282 & MMO > 0) #bottom 10%, one-tailed Z score test
Lstandard <- subset(larv.ranef.all,  sBvL < -1.282 & MMO > 0  & MEM < 0) #bottom 10%, one-tailed Z score test
#Lbamboo <- subset(larv.ranef.all,  sBvL > 1.282 ) #bottom 10%, one-tailed Z score test

Lbamboo <- subset(larv.ranef.all,  sBvL > 1.282 & MMO > 0  & MEM < 0) #bottom 10%, one-tailed Z score test


# 
# top.lar <- subset(larv.ranef.all, scale(BvL) > 1.282 | scale(BvL) < -1.282) #top 10%, one-tailed Z score test 
# top.wat <- subset(larv.ranef.all, scale(BvL) > 1.282 | scale(BvL) < -1.282) #top 10%, one-tailed Z score test 
# Lbamboo <- subset(modes,  sBvL > 1.28 & sMEM>0) #top 10%, one-tailed Z score test
# Lstandard <- subset(modes,  sBvL < -1.28 & sMMO>0) #bottom 10%, one-tailed Z score test

top.lar.mem <- subset(larv.ranef.all, scale(MEM) > 2.326) #top 1%, one-tailed Z score test 
top.lar.mmo <- subset(larv.ranef.all, scale(MMO) > 2.326) #top 1%, one-tailed Z score test 
top.lar.emo <- subset(larv.ranef.all, scale(EMO) > 2.326) #top 1%, one-tailed Z score test 

top.lar.bam <- subset(larv.ranef.all, scale(Bamboo) > 2.326) #top 1%, one-tailed Z score test 
top.lar.sta <- subset(larv.ranef.all, scale(`Larval food`) > 2.326) #top 1%, one-tailed Z score test 

top.lar.uniq <- unique(c(top.lar.mem$kegg.name,top.lar.mmo$kegg.name,top.lar.emo$kegg.name,top.lar.bam$kegg.name))

larv.top1p <- subset(larv.ranef.all, larv.ranef.all$kegg.name %in% top.lar.uniq)

larv.top1p.sub <- larv.top1p[,c("kegg.name","MMO","EMO","MEM","Bamboo","Larval food")]

#cond.all.df

#subset(cond.df)

row.names(larv.top1p.sub) <- larv.top1p.sub$kegg.name
larv.top1p.sub2 <- larv.top1p.sub[,colnames(larv.top1p.sub)!="kegg.name"]
cond.top1p.mat <- as.matrix(larv.top1p.sub2)

library("pheatmap")
##full thing
pheatmap(cond.top1p.mat,
         #color = my_palette,
         #color = colorRampPalette(c("blue","white","red"))(100),
         #breaks = breaks,
         #color = colorRampPalette(brewer.pal(n=9,name="Purple-Blue"))(100)
         color = colorRampPalette(rev(hcl.colors(n=10,"Vik")))(100),
         scale="row",
         show_rownames=F,
         cluster_cols=F)
```

```{r}
table(top.bam.lar$Pathway)
head(top.bam.lar$Pathway)

library(dplyr)
library(ggplot2)

# Step 1: Count the number of times each pathway appears
top.lar.bam20 <- top.bam.lar %>%
  count(Pathway) %>%
  arrange(desc(n)) %>%
  slice_head(n = 20)

# Step 2: Reorder for plotting
top.lar.bam20.form <- top.lar.bam20 %>%
  mutate(Pathway = factor(Pathway, levels = rev(Pathway)))

#top.lar.bam20 <- top.lar.bam21.form[top.lar.bam21.form$Pathway!="Not Included in Pathway or Brite; Unclassified: metabolism; Enzymes with EC numbers",]

# Step 3: Plot
ggplot(top.lar.bam20.form, aes(x = Pathway, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 - larvae, bamboo diet",
       x = "Pathway",
       y = "Count (Number of Times Detected)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 7))

top.lar.sta20 <- top.lar.lar %>%
  count(Pathway) %>%
  arrange(desc(n)) %>%
  slice_head(n = 20)

# Step 2: Reorder for plotting
top.lar.sta20.form <- top.lar.sta20 %>%
  mutate(Pathway = factor(Pathway, levels = rev(Pathway)))

#top.lar.sta20 <- top.lar.sta21.form[top.lar.sta21.form$Pathway!="Not Included in Pathway or Brite; Unclassified: metabolism; Enzymes with EC numbers",]

# Step 3: Plot
ggplot(top.lar.sta20.form, aes(x = Pathway, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 - larvae, standard diet",
       x = "Pathway",
       y = "Count (Number of Times Detected)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 7))

```

# Combining genes lists

```{r}
library(dplyr)

# Get the unique kegg.names for each df
Lbamboo_only <- Lbamboo %>%
  filter(!(kegg.name %in% c(Lstandard$kegg.name)))

Lstandard_only <- Lstandard %>%
  filter(!(kegg.name %in% c(Lbamboo$kegg.name)))

Lbamboo_only$food_indicator <- "Bamboo"
Lstandard_only$food_indicator <- "Standard"

Lboth.top <- rbind(Lbamboo_only,Lstandard_only)
write.csv(Lboth.top,file="fxns.larvae.top.csv")

Wbamboo_only <- Wbamboo %>%
  filter(!(kegg.name %in% c(Lstandard$kegg.name,Lbamboo$kegg.name)))

Wstandard_only <- Wstandard %>%
  filter(!(kegg.name %in% c(Lstandard$kegg.name,Lbamboo$kegg.name)))

Wbamboo_only$food_indicator <- "Bamboo"
Wstandard_only$food_indicator <- "Standard"

Wboth.top <- rbind(Wbamboo_only,Wstandard_only)

Wbamboo$food_indicator <- "Bamboo"
Wstandard$food_indicator <- "Standard"

Wboth.top_notuniq <- rbind(Wstandard,Wbamboo)

write.csv(Wboth.top_notuniq,file="fxns.water.top.csv")
write.csv(Wboth.top,file="fxns.water.top_nosharing.csv")

```

# Checking out pathways

## Larvae

```{r}
library(dplyr)

top.lar2 <- top.lar %>%
  separate_rows(Pathway, sep = "\\|")

top.lar2$Pathway <- trimws(top.lar2$Pathway)

top.lar.bam <- subset(top.lar2,BvL > 0)
top.lar.sta <- subset(top.lar2,BvL < 0)

top20.lar.bam <- top.lar.bam %>%
  count(Pathway, sort = TRUE) %>%
  slice_max(n, n = 20)

ggplot(top20.lar.bam,aes(x=n,y=reorder(Pathway,n)))+
  geom_bar(stat="identity")+
  ggtitle("Larvae, bamboo diet")+
  theme_cowplot()
#theme(axis.text.x=element_text(angle=90,hjust=1),axis.title.x=element_blank())

top.lar.sta <- subset(top.lar2,BvL < 0)

top20.lar.sta <- top.lar.sta %>%
  count(Pathway, sort = TRUE) %>%
  slice_max(n, n = 20)

ggplot(top20.lar.sta,aes(x=n,y=reorder(Pathway,n)))+
  geom_bar(stat="identity")+
  ggtitle("Larvae, standard diet")+
  theme_cowplot()
#theme(axis.text.x=element_text(angle=90,hjust=1),axis.title.x=element_blank())

##plus looking at microbes
top.lar.bam.mem <- subset(top.lar2,BvL > 0 & MEM > 0)
top.lar.sta.mmo <- subset(top.lar2,BvL < 0 & MMO > 0)

top20.lar.bam.mem <- top.lar.bam.mem %>%
  count(Pathway, sort = TRUE) %>%
  slice_max(n, n = 20)

ggplot(top20.lar.bam.mem,aes(x=n,y=reorder(Pathway,n)))+
  geom_bar(stat="identity")+
  ggtitle("Larvae, bamboo diet, MEM > 0")+
  theme_cowplot()
#theme(axis.text.x=element_text(angle=90,hjust=1),axis.title.x=element_blank())

top20.lar.sta.mmo <- top.lar.sta.mmo %>%
  count(Pathway, sort = TRUE) %>%
  slice_max(n, n = 20)

ggplot(top20.lar.sta.mmo,aes(x=n,y=reorder(Pathway,n)))+
  geom_bar(stat="identity")+
  ggtitle("Larvae, standard diet, MMO > 0")+
  theme_cowplot()

top.lar.sta <- subset(top.lar2,BvL < 0)

top20.lar.sta <- top.lar.sta %>%
  count(Pathway, sort = TRUE) %>%
  slice_max(n, n = 20)

ggplot(top20.lar.sta,aes(x=n,y=reorder(Pathway,n)))+
  geom_bar(stat="identity")+
  ggtitle("Larvae, standard diet")+
  theme_cowplot()
#theme(axis.text.x=element_text(angle=90,hjust=1),axis.title.x=element_blank())



mean_BvL_by_pathway_filtered <- mean_BvL_by_pathway %>%
  filter(n_genes > 1)

top.lar.bam <- subset(mean_BvL_by_pathway_filtered, scale(mean_BvL) > 1.645) #top 1%, one-tailed Z score test 
top.lar.sta <- subset(mean_BvL_by_pathway_filtered, scale(mean_BvL) < -1.645) #top 1%, one-tailed Z score test 

write.csv(mean_BvL_by_pathway_filtered,file="pathways.all.csv")


write.csv(top.lar.bam,file="pathways.meanBvL.larvae.bamboo.csv")
write.csv(top.lar.sta,file="pathways.meanBvL.larvae.standard.csv")


top20.lar.bam <- top.lar.bam[order(-top.lar.bam$mean_BvL), ][1:20, ]
top20.lar.bam

ggplot(top20.lar.bam,aes(x=reorder(Pathway,mean_BvL),y=mean_BvL))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=90,hjust=1))

top20.lar.sta <- top.lar.sta[order(-top.lar.sta$mean_BvL), ][1:20, ]
top20.lar.sta
```

```{r}
library(dplyr)

larv.ranef.all2 <- larv.ranef.all %>%
  separate_rows(Pathway, sep = "\\|")

larv.ranef.all2$Pathway <- trimws(larv.ranef.all2$Pathway)

mean_BvL_by_pathway <- larv.ranef.all2 %>%
  filter(!is.na(Pathway) & Pathway != "") %>%
  select(Pathway, BvL) %>%
  group_by(Pathway) %>%
  summarise(
    mean_BvL = mean(BvL, na.rm = TRUE),
    n_genes = n()
  ) %>%
  arrange(desc(mean_BvL))

mean_BvL_by_pathway_filtered <- mean_BvL_by_pathway %>%
  filter(n_genes > 1)

top.lar.bam <- subset(mean_BvL_by_pathway_filtered, scale(mean_BvL) > 1.645) #top 1%, one-tailed Z score test 
top.lar.sta <- subset(mean_BvL_by_pathway_filtered, scale(mean_BvL) < -1.645) #top 1%, one-tailed Z score test 

write.csv(mean_BvL_by_pathway_filtered,file="pathways.all.csv")


write.csv(top.lar.bam,file="pathways.meanBvL.larvae.bamboo.csv")
write.csv(top.lar.sta,file="pathways.meanBvL.larvae.standard.csv")


top20.lar.bam <- top.lar.bam[order(-top.lar.bam$mean_BvL), ][1:20, ]
top20.lar.bam

ggplot(top20.lar.bam,aes(x=reorder(Pathway,mean_BvL),y=mean_BvL))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=90,hjust=1))

top20.lar.sta <- top.lar.sta[order(-top.lar.sta$mean_BvL), ][1:20, ]
top20.lar.sta
```

# Heat map in the style of the other one

```{r}
##unique IDs with the name
Lmaster$kegg.path <- paste0(Lmaster$KEGG_ID,"_",Lmaster$Name)

cond.all.df <- Lmaster[,c("kegg.path","MMO","EMO","MEM","Bamboo","Larval.food")]

#cond.all.df

#subset(cond.df)

row.names(cond.all.df) <- cond.all.df$kegg.path
cond.all.df2 <- cond.all.df[,colnames(cond.all.df)!="kegg.path"]
cond.all.mat <- as.matrix(cond.all.df2)

library("pheatmap")
##full thing
pheatmap(cond.all.mat,
         #color = my_palette,
         #color = colorRampPalette(c("blue","white","red"))(100),
         #breaks = breaks,
         #color = colorRampPalette(brewer.pal(n=9,name="Purple-Blue"))(100)
         color = colorRampPalette(rev(hcl.colors(n=10,"Vik")))(100),
         #scale="row",
         show_rownames=F,
         cluster_cols=F)

##need to rescale so 0 is in the middle

# Actual data limits
min_val <- min(cond.all.mat, na.rm = TRUE)
max_val <- max(cond.all.mat, na.rm = TRUE)

# Define number of colors on each side of 0 (you can adjust this)
n_neg <- round(100 * abs(min_val) / (abs(min_val) + max_val))  # e.g., ~33
n_pos <- 100 - n_neg                                           # e.g., ~67
##Vik
#neg_cols <- colorRampPalette(c(hcl.colors(n=10,"Vik")[1:5],"white"))(n_neg)
#pos_cols <- colorRampPalette(c("white",hcl.colors(n=10,"Vik")[6:10]))(n_pos)

neg_cols <- colorRampPalette(c("wheat4", "white"))(n_neg)
pos_cols <- colorRampPalette(c("white","mediumorchid4","black"))(n_pos)
#neg_cols <- colorRampPalette(rev(hcl.colors(10, "Temps"))[1:5])(n_neg)  # pink to white
#pos_cols <- colorRampPalette(c("white",rev(hcl.colors(9, "BrwnYl"))))(n_pos)  # white to green
my_palette <- c(neg_cols, pos_cols)
##Geyser

#neg_colors <- colorRampPalette(c("black", "white"))(n_neg)
#pos_colors <- colorRampPalette(c("white","#967969","black"))(n_pos)
#my_palette <- c(neg_colors, pos_colors)

# Define breaks so that 0 falls exactly at the transition
breaks <- seq(min_val, max_val, length.out = length(my_palette) + 1)

# Draw the heatmap
pheatmap(cond.all.mat,
         color = my_palette,
         #color = colorRampPalette(c("blue","white","red"))(100),
         breaks = breaks,
         #annotation_row = row.source,
         #annotation_colors = ann_colors,
         border_color = F,
         show_rownames=F,
         cluster_cols=F)

```
