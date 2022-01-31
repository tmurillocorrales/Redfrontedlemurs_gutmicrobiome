#### Eukaryotic 18S amplicon data ####

## Dominik Schneider, Dirk Berkelmann and Avril von Hoyningen-Huene contributed to these scripts

# Load libraries ####
library(ampvis2)
library(ape)
library(phangorn)
library(stringr)
library(reshape2)
library(rayshader)
library(ggdark)
require(GUniFrac)
library(tibble)
library(viridis)
library(scales)
library(dplyr)
library(tidyverse)

# Load data ####
myotutable = read.table("18S_BLCA_ASV_table.txt", sep = "\t", header=T, check.names=F, comment.char="", skip = 1, quote = "")
# Note for metadata do not use capital letters
mymetadata = read.table("mymetadata_18s.txt", sep="\t", header=T, comment.char = "", quote = "\"")

# Modify data for ampvis2####
# Get taxonomy from OTU table
row.names(myotutable)=myotutable$`#OTU ID`
tax_info = data.matrix(myotutable$taxonomy)
OTUID = row.names(myotutable)

# Construct a data frame from taxonomy info and separate by semicolon
tax_info = str_split_fixed(tax_info, '\\;',14)
tax_info = data.frame(tax_info, stringsAsFactors = FALSE)

# Taxomony was assigned with Bayesian LCA-based Taxonomic Classification Method (BLCA) -> https://github.com/qunfengdong/BLCA
# convert columns to numeric
tax_info$X2 = as.numeric(as.character(tax_info$X2))
tax_info$X4 = as.numeric(as.character(tax_info$X4))
tax_info$X6 = as.numeric(as.character(tax_info$X6))
tax_info$X8 = as.numeric(as.character(tax_info$X8))
tax_info$X10 = as.numeric(as.character(tax_info$X10))
tax_info$X12 = as.numeric(as.character(tax_info$X12))
tax_info$X14 = as.numeric(as.character(tax_info$X14))

# A confidence score threshold of 80% was used as recommended by the authors
# Apply threshold to set to unclassified
tax_info$X1[ tax_info$X2 < 80 ] = "Unclassified"
tax_info$X3[ tax_info$X4 < 80 ] = "Unclassified"
tax_info$X5[ tax_info$X6 < 80 ] = "Unclassified"
tax_info$X7[ tax_info$X8 < 80 ] = "Unclassified"
tax_info$X9[ tax_info$X10 < 80 ] = "Unclassified"
tax_info$X11[ tax_info$X12 < 80 ] = "Unclassified"
tax_info$X13[ tax_info$X14 < 80 ] = "Unclassified"
tax_info = data.frame(cbind(tax_info$X1,tax_info$X3,tax_info$X5,tax_info$X7,tax_info$X9,tax_info$X11,tax_info$X13))

#Renaming
tax_info[] = lapply(tax_info, gsub, pattern = "superkingdom:", replacement = "")
tax_info[] = lapply(tax_info, gsub, pattern = "phylum:", replacement = "")
tax_info[] = lapply(tax_info, gsub, pattern = "class:", replacement = "")
tax_info[] = lapply(tax_info, gsub, pattern = "order:", replacement = "")
tax_info[] = lapply(tax_info, gsub, pattern = "family:", replacement = "")
tax_info[] = lapply(tax_info, gsub, pattern = "genus:", replacement = "")
tax_info[] = lapply(tax_info, gsub, pattern = "species:", replacement = "")
tax_info[] = lapply(tax_info, gsub, pattern = "unclassified", replacement = "Unclassified")
tax_info$X1 = gsub(x = tax_info$X1, pattern = "Unclassified", replacement = "k__Unclassified")
tax_info$X2 = gsub(x = tax_info$X2, pattern = "Unclassified", replacement = "p__Unclassified")
tax_info$X3 = gsub(x = tax_info$X3, pattern = "Unclassified", replacement = "c__Unclassified")
tax_info$X4 = gsub(x = tax_info$X4, pattern = "Unclassified", replacement = "o__Unclassified")
tax_info$X5 = gsub(x = tax_info$X5, pattern = "Unclassified", replacement = "f__Unclassified")
tax_info$X6 = gsub(x = tax_info$X6, pattern = "Unclassified", replacement = "g__Unclassified")
tax_info$X7 = gsub(x = tax_info$X7, pattern = "Unclassified", replacement = "s__Unclassified")
tax_info[] = lapply(tax_info, gsub, pattern = "__Unclassified .*", replacement = "__Unclassified")

# Add column names
colnames(tax_info) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(tax_info) = OTUID

# Remove old taxonomy from OTU table
myotutable$taxonomy = NULL

# Add new taxonomy (merge both tables)
myotutable = cbind(myotutable, tax_info)
myotutable$tax_info = NULL

# Rename first column
colnames(myotutable)[1] = "OTU"
i = sapply(myotutable, is.factor)
myotutable[i] = lapply(myotutable[i], as.character)

# Combine the data to ampvis2 class
dataset_BLCA = amp_load(otutable = myotutable,
                        metadata = mymetadata)


dataset_BLCA
#400 samples, 6276 OTUs

#Remove data from April 2018 and samples with less than 6000 reads
dataset_BLCA <- amp_subset_samples(dataset_BLCA, !month_year %in% c("04_Apr18"), minreads = 6000)
#380 samples 6245 otus

# Pick out parasite phyla previously reported in animals and humans ##
#Load txt file with parasite phyla
parasites_phy = scan("para_phylum.txt", what="", sep="\n")
#Extract phyla from ampvis object
parasites_phy = amp_subset_taxa(dataset_BLCA, tax_vector = parasites_phy, normalise = F, remove = F)
#Remove common environmental contaminant "Enoplea".
parasites_phy = amp_subset_taxa(parasites_phy, tax_vector = c("Enoplea"), normalise = F, remove = T)
parasites_phy
#447 OTUS

### Merge data ####
#Merge all triplicates of each individual per month.
#IN ORDER FOR THE MERGE FUNCTION TO WORK THE SAMPLE ID COLUMN SHOULD BE TITLED SampleID.

#Use merge functions from ampvis2
dataset_merge <- amp_mergereplicates(dataset_BLCA, merge_var = "individual_month")

dataset_merge
#298 samples 6245

#Only parasite phyla of interest
parasites_phy
#380 samples otus 447
dataset_phy_merge <- amp_mergereplicates(parasites_phy, merge_var = "individual_month")
dataset_phy_merge
#298 samples OTUs 447

##PCOA ####
#PCOA showed in supplementary figure S1c
PCOA_phy = amp_ordinate(dataset_phy_merge,
                        num_threads = 4,
                        type = "PCOA",
                        distmeasure = "jaccard",
                        sample_color_by = "season",   
                        sample_shape_by = "group",
                        transform = "none",
                        sample_colorframe = FALSE,
                        sample_point_size = 4)

PCOA_phy


### Extract jaccard distance matrix from PCOA ##
jaccard = as.matrix(PCOA_phy$plot_env$distmatrix)
## Extract metadata from PCOA ###
metadata_pcoa = PCOA_phy$data

### Permanova #######
perma = adonis(jaccard ~ group + season + sex + age_months, 
               metadata_pcoa, permutations = 10000, strata = metadata_pcoa$individual)
perma
#Results:
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#group        3     2.062 0.68717  2.4434 0.02404 0.3855    
#season       1     0.977 0.97658  3.4725 0.01139 0.0002 ***
#sex          1     0.233 0.23340  0.8299 0.00272 0.5176    
#age_months   1     0.632 0.63215  2.2478 0.00737 0.5458    
#Residuals  291    81.838 0.28123         0.95447           
#Total      297    85.742                 1.00000

#Create dataframe with results
permanova_results = as.data.frame(perma$aov.tab)

#### Multiple testing correction ####
#Order pvalues from min to max
pvals = perma$aov.tab$`Pr(>F)`
pvals
#Set pvalues to numeric
pvals = as.numeric(pvals)

#Calculate BH adjusted pvalues
BH  = p.adjust(pvals, method = "BH")
#Add corrected pvalues to data frame
permanova_results$BH = BH
#Export results which are shown in supplementary table S7.
write.table(permanova_results, file = "permanova_jaccard.csv", sep = ",", row.names = TRUE, col.names = TRUE)

## Alpha diversity ####
# Calculate alpha diversity indices in ampvis2
alphadiv = amp_alphadiv(parasites_phy, measure= c("observed","shannon","simpson","invsimpson"), richness = TRUE)

# Fix rownames
rownames(alphadiv)=alphadiv$samplesid
alphadiv$samplesid= NULL

#Calculate mean and standard deviation for each group per month
stats = alphadiv %>% group_by(month_year, group) %>% summarise(mean.richness = mean(ObservedOTUs), sd.richness = sd(ObservedOTUs),
                                                               min.richness = min(ObservedOTUs), max.richness = max(ObservedOTUs))

#Calculate mean and standard deviation for the group 
stats_group = alphadiv %>% group_by(group) %>% summarise(mean.richness = mean(ObservedOTUs), sd.richness = sd(ObservedOTUs))

### Organize month data according to sampling timeline
alphadiv$month_year <- factor(alphadiv$month_year, 
                              levels = c("05_May18", "06_Jun18","07_July18",
                                         "08_Aug18","09_Sep18", "10_Oct18", "11_Nov18", "12_Dec18", 
                                         "01_Jan19", "02_Feb19", "03_Mar19","04_Apr19"))

view(alphadiv)

#Box plot to prepare figure 1d ##
d3=ggplot(alphadiv, aes(x= month_year, y=ObservedOTUs, fill=month_year)) +
  geom_boxplot(aes()) +
  theme_gray()+
  labs(y = "Parasite richness") + 
  scale_fill_manual(values = c("05_May18" = "#bc7369", "06_Jun18" = "#bc7369", "07_July18" = "#bc7369",
                               "08_Aug18" = "#bc7369", "09_Sep18" = "#bc7369", "10_Oct18" = "#bc7369",
                               "11_Nov18" = "#5ea07c", "12_Dec18" = "#5ea07c", "01_Jan19" = "#5ea07c",
                               "02_Feb19" = "#5ea07c", "03_Mar19" = "#5ea07c", "04_Apr19" = "#bc7369")) +
  theme_gray(base_size = 14)+ 
  facet_wrap( ~ group, ncol=4) +
  theme(axis.title.x = element_blank(),
        legend.position = "none", axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 11), 
        axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5), strip.text = element_text(size = 12))
d3

# Normalized data with GMPR####
library(GMPR)

# Add GMPR
gmpr.size.factor = GMPR(t(dataset_BLCA$abund), min_ct = 2, intersect_no = 4) # ignore error message in Windows, 
#it seems to lead to the same results as Linux (without error there)
otu.tab.norm = data.frame(t(t(dataset_BLCA$abund) / gmpr.size.factor), check.names = FALSE)
dataset_BLCA$GMPR = otu.tab.norm
dataset_BLCA$READ_COUNTS = dataset_BLCA$abund

# Replace read counts with GMPR
dataset_BLCA$abund = dataset_BLCA$GMPR
dataset_BLCA$abund = dataset_BLCA$READ_COUNTS

# Pick out parasite phyla previously reported in animals and humans ##
#Load txt file with parasite phyla
parasites_phy = scan("para_phylum.txt", what="", sep="\n")
#Extract phyla from ampvis object
parasites_phy = amp_subset_taxa(dataset_BLCA, tax_vector = parasites_phy, normalise = F, remove = F)
#Remove common environmental contaminant "Enoplea".
parasites_phy = amp_subset_taxa(parasites_phy, tax_vector = c("Enoplea"), normalise = F, remove = T)
parasites_phy
#447 OTUS

# Line charts ####
#Chunk of script to produce supplemental figure 1b. Panel figures were put together in inkscape.

#Prepare group A separately as is missing two sampling months (January and February)
## Groups B,F,J ####
#Subset data in ampvis2
group.bfj <- amp_subset_samples(parasites_phy, !group %in% c("A"))
#After: 330 samples and 441 OTUs

# Export data
amp_export_otutable(group.bfj, "table_linecharts_bfj", sep = "\t")
write.table(group.bfj$metadata, file = "metadata_linecharts_bfj.csv", sep = "\t")

# Load data
otu = read.table("table_linecharts_bfj.csv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata_linecharts_bfj.csv", row.names = 1, header = T)

# Check if rows and collumns match (very important, otherwise everything will be mixed in the aggregation step, must be TRUE)
all(rownames(metadata) %in% colnames(otu))

# Collapse/Aggregate samples by metadata column (also used for sorting samples!!!)
taxonomy = paste0(otu$Phylum,"; ",otu$Order,"; ",otu$Class)
unique(taxonomy)

otu_transposed = t(otu[,1:(dim(otu)[2]-7)])

# You don't want to aggregate your data? Use this line instead:
aggdata = otu_transposed

# "Restore" OTU table
otu_tax = as.data.frame(t(aggdata))

# Combine tables
otu_tax = cbind(otu_tax, taxonomy)

# Fix names...
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = "No blast hit; ", replacement = "No blast hit")

# Collapse phylogeny according selected level
taxa_level = aggregate(otu_tax[,1:dim(otu_tax)[2]-1], by = list(otu_tax$taxonomy), FUN = sum)

# Sanity check
colSums(taxa_level[,2:dim(taxa_level)[2]])

# Add row names
rownames(taxa_level) = taxa_level[,1]
taxa_level = taxa_level[,2:dim(taxa_level)[2]]

# Transform to relative abundance
taxa_level = scale(taxa_level, center=F, scale=colSums(taxa_level))
taxa_level = taxa_level*100
colSums(taxa_level)

# Remove, Group low abundant phyla
# Pick only phyla >= 0.25%
percentage = 0.25
#            ^ number to filter phylum/class/order abundance

abundant.phyla = apply(taxa_level, 1, function(x){perc = sum(x>percentage)/length(x); perc>0})
rare.phyla = apply(taxa_level, 1, function(x){perc = sum(x<percentage)/length(x); perc>((100-percentage)/100)})
filtered_table = taxa_level[abundant.phyla,]
filtered_rare_table = taxa_level[rare.phyla,]

# Sustract >= phyla from 100% to get percentage of rare_taxa
rare_taxa = 100-colSums(filtered_table)

# Add rare_taxa as column to table
filtered_table = rbind(filtered_table, rare_taxa)

# Show only top taxa
TOP = 3 # change to top X
sub = filtered_table[ order(rowMeans(filtered_table), decreasing = T), ]
sub = sub[1:TOP,]

# Add metadata to combine samples by
table_final = rbind(year = as.character(metadata$year), month = as.character(metadata$month), 
                    individual = as.character(metadata$individual), group = as.character(metadata$group),sub)

table_final = data.frame(t(table_final))

filtered_table_melt = melt(table_final, id = c("year","month","individual", "group"))

#Check for presence of NAs
which(is.na.data.frame(filtered_table_melt))

#Rename columns
names(filtered_table_melt)[names(filtered_table_melt) == "variable"] = "Taxon"
names(filtered_table_melt)[names(filtered_table_melt) == "value"] = "Abundance"

#Set abundance to numeric
filtered_table_melt$Abundance = as.numeric(filtered_table_melt$Abundance)
str(filtered_table_melt)

#Rename month column to month_year
names(filtered_table_melt)[names(filtered_table_melt)=="month"] = "month_year"

filtered_table_melt[c("month", "year")] = str_split_fixed(filtered_table_melt$month, "_", 2)

#Calculate means and standard deviations
stats_bfj = filtered_table_melt %>% group_by(month_year, group, Taxon) %>% summarise(mean = mean(Abundance), sd = sd(Abundance))
stats_bfj_year = filtered_table_melt %>% group_by(group, Taxon) %>% summarise(mean = mean(Abundance), sd = sd(Abundance))

#Order the month variable according to the sampling timeline
filtered_table_melt$month_year <- as.factor(filtered_table_melt$month_year, 
                                         levels = c("05_May18", "06_Jun18","07_July18",
                                                    "08_Aug18","09_Sep18", "10_Oct18", 
                                                    "11_Nov18", "12_Dec18", 
                                                    "01_Jan19", "02_Feb19", "03_Mar19","04_Apr19"))

filtered_table_melt$month <- factor(filtered_table_melt$month, 
                                    levels = c("05", "06","07",
                                               "08","09", "10", 
                                               "11", "12", 
                                               "01", "02", "03","04"))

#Convert month to numeric
filtered_table_melt$month= as.numeric(filtered_table_melt$month)

# Create lineplot
lineplot = ggplot(filtered_table_melt, aes(x = month, y = Abundance, fill = Taxon, colour = Taxon, group = month))

# all together
lineplot +
  geom_smooth(method="loess", span = 0.3, fullrange = F, se = T, aes(group = Taxon, fill = Taxon), 
              level = 0.95, size = 1) +
  theme_gray(base_size = 14)  +  
  scale_x_continuous(breaks = 1:length(unique(filtered_table_melt$month)), labels = levels(unique(filtered_table_melt$month))) +
  scale_fill_viridis_d(option = "C", begin= 0, end = 0.8) +
  scale_color_viridis_d(option = "C", begin= 0, end = 0.8) +
  labs(y = "Relative abundance")

# Facet_wrap by group
lineplot +
  geom_smooth(method="loess", span = 0.3, fullrange = F, se = T, aes(group = Taxon, fill = Taxon),
              level = 0.95, size = 1) +
  scale_x_continuous(breaks = 1:length(unique(filtered_table_melt$month)), 
                     labels = levels(unique(filtered_table_melt$month))) +
  scale_fill_viridis_d(option = "C", begin= 0, end = 0.8) +
  scale_color_viridis_d(option = "C", begin= 0, end = 0.8) +
  labs(y = "Relative abundance")+
  theme(axis.text.x = element_text(size = 10 , angle = 90, vjust = 0.4, hjust = 1,face = "bold"),
        axis.text.y = element_text(size = 10), panel.grid.major = element_line(size = 0.1)) +
  facet_wrap( ~ group, scales = "free_y", ncol = 4) +
  theme(axis.title.x=element_blank(),  axis.title.y = element_blank(), strip.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12), legend.position = "none")

## Group A ####
#Group A was plotted separately as it was not sampled in January and February
#Subset data in ampvis2
group.a <- amp_subset_samples(parasites_phy, group %in% c("A"))
#50 samples and 206 OTUs

# Export data
amp_export_otutable(group.a, "table_linecharts_A", sep = "\t")
write.table(group.a$metadata, file = "metadata_linecharts_A.csv", sep = "\t")

# Load data
otu = read.table("table_linecharts_A.csv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata_linecharts_A.csv", row.names = 1, header = T)

# Check if rows and collumns match (very important, otherwise everything will be mixed in the aggregation step, must be TRUE)
all(rownames(metadata) %in% colnames(otu))

# Collapse/Aggregate samples by metadata column (also used for sorting samples)
taxonomy = paste0(otu$Phylum,"; ",otu$Order,"; ",otu$Class)
unique(taxonomy)
otu_transposed = t(otu[,1:(dim(otu)[2]-7)])
aggdata = otu_transposed

# "Restore" OTU table
otu_tax = as.data.frame(t(aggdata))

# Combine tables
otu_tax = cbind(otu_tax, taxonomy)

# Fix names...
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = "No blast hit; ", replacement = "No blast hit")
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = "Bacteria; ", replacement = "")

# Collapse phylogeny according selected level
# Defined above (Order in this case)
taxa_level = aggregate(otu_tax[,1:dim(otu_tax)[2]-1], by = list(otu_tax$taxonomy), FUN = sum)

# Sanity check
colSums(taxa_level[,2:dim(taxa_level)[2]])

# Add row names
rownames(taxa_level) = taxa_level[,1]
taxa_level = taxa_level[,2:dim(taxa_level)[2]]

# Transform to relative abundance
taxa_level = scale(taxa_level, center=F, scale=colSums(taxa_level))
taxa_level = taxa_level*100
colSums(taxa_level)

## Remove, Group low abundant phyla
# Pick only phyla >= 0.25%
percentage = 0.25
#            ^ number to filter phylum/class/order abundance

abundant.phyla = apply(taxa_level, 1, function(x){perc = sum(x>percentage)/length(x); perc>0})
rare.phyla = apply(taxa_level, 1, function(x){perc = sum(x<percentage)/length(x); perc>((100-percentage)/100)})
filtered_table = taxa_level[abundant.phyla,]
filtered_rare_table = taxa_level[rare.phyla,]

# Sustract >= phyla from 100% to get percentage of rare_taxa
rare_taxa = 100-colSums(filtered_table)

# Add rare_taxa as column to table
filtered_table = rbind(filtered_table, rare_taxa)

# Show only top 3
TOP = 3 # change to top X
sub = filtered_table[ order(rowMeans(filtered_table), decreasing = T), ]
sub = sub[1:TOP,]

# Add metadata to combine samples by
table_final = rbind(year = as.character(metadata$year), month = as.character(metadata$month), 
                    individual = as.character(metadata$individual), group = as.character(metadata$group),sub)

table_final = data.frame(t(table_final))

filtered_table_melt = melt(table_final, id = c("year","month","individual", "group"))

#Check for presence of NAs
which(is.na.data.frame(filtered_table_melt))

#Rename columns
names(filtered_table_melt)[names(filtered_table_melt) == "variable"] = "Taxon"
names(filtered_table_melt)[names(filtered_table_melt) == "value"] = "Abundance"

#Set abundance to numeric
filtered_table_melt$Abundance = as.numeric(filtered_table_melt$Abundance)
str(filtered_table_melt)

#Rename month column to month_year
names(filtered_table_melt)[names(filtered_table_melt)=="month"] = "month_year"

filtered_table_melt[c("month", "year")] = str_split_fixed(filtered_table_melt$month, "_", 2)

#Calculate means and standard deviations
stats_a = filtered_table_melt %>% group_by(month_year, group, Taxon) %>% summarise(mean = mean(Abundance), sd = sd(Abundance))
stats_a_year = filtered_table_melt %>% group_by(group, Taxon) %>% summarise(mean = mean(Abundance), sd = sd(Abundance))

#Order the month variable according to the sampling timeline
filtered_table_melt$month_year <- as.factor(filtered_table_melt$month_year, 
                                            levels = c("05_May18", "06_Jun18","07_July18",
                                                       "08_Aug18","09_Sep18", "10_Oct18", 
                                                       "11_Nov18", "12_Dec18", 
                                                       "01_Jan19", "02_Feb19", "03_Mar19","04_Apr19"))

filtered_table_melt$month <- factor(filtered_table_melt$month, 
                                    levels = c("05", "06","07",
                                               "08","09", "10", 
                                               "11", "12", 
                                               "01", "02", "03","04"))

#Convert month to numeric
filtered_table_melt$month= as.numeric(filtered_table_melt$month)

# Create lineplot
lineplot = ggplot(filtered_table_melt, aes(x = month, y = Abundance, fill = Taxon, colour = Taxon, group = month))

# all together
lineplot +
  geom_smooth(method="loess", span = 0.3, fullrange = F, se = T, aes(group = Taxon, fill = Taxon), 
              level = 0.95, size = 1) +
  theme_gray(base_size = 14)  +  
  scale_x_continuous(breaks = 1:length(unique(filtered_table_melt$month)), labels = levels(unique(filtered_table_melt$month))) +
  scale_fill_viridis_d(option = "C", begin= 0, end = 0.8) +
  scale_color_viridis_d(option = "C", begin= 0, end = 0.8) +
  labs(y = "Relative abundance")

# Bar charts ####
# Sort your metadata table and OTU table by sample name

###### Full dataset ####
#Prepare data from normalized counts
# Export data
amp_export_otutable(dataset_BLCA, "table_barcharts", sep = "\t")
write.table(dataset_BLCA$metadata, file = "metadata_barcharts.csv", sep = "\t")

# Load data
otu = read.table("table_barcharts.csv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata_barcharts.csv", row.names = 1, header = T)

# Check if rows and collumns match (very important, otherwise everything will be mixed in the aggregation step, must be TRUE)
all(rownames(metadata) %in% colnames(otu))

# Collapse/Aggregate samples by metadata column (also used for sorting samples)
taxonomy = paste0(otu$Kingdom,"; ",otu$Phylum,"; ", otu$Class)
unique(taxonomy)

otu_transposed = t(otu[,1:(dim(otu)[2]-7)])

aggdata = aggregate(otu_transposed, by=list(metadata$month_year), FUN=sum)
#                                                     ^ Argument for aggregating and sorting

#Renaming
rownames(aggdata) = aggdata$Group.1
aggdata$Group.1 = NULL

# "Restore" OTU table
otu_tax = as.data.frame(t(aggdata))

# Combine tables
otu_tax = cbind(otu_tax, taxonomy)

# Fix names...
# No Blast Hit
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ; ", replacement = "")
# Separator
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ; ", replacement = "; ")
# Remove accession nr. identity and evalue
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = "~.*", replacement = "")

# Collapse phylogeny according selected level
taxa_level = aggregate(otu_tax[,1:dim(otu_tax)[2]-1], by = list(otu_tax$taxonomy), FUN = sum)

# Add row names
rownames(taxa_level) = taxa_level[,1]
taxa_level = taxa_level[,2:dim(taxa_level)[2]]

# Transform to relative abundance
taxa_level = scale(taxa_level, center=F, scale=colSums(taxa_level))
taxa_level = taxa_level*100
colSums(taxa_level)

# Remove group low abundant phyla as rare taxa
# Pick only phyla >= 2%
percentage = 2
#            ^ number to filter phylum/class/order abundance (1 = everything below 1% will be summarized as "rare taxa" in this case) ####

abundant.phyla = apply(taxa_level, 1, function(x){perc = sum(x>percentage)/length(x); perc>0})
rare.phyla = apply(taxa_level, 1, function(x){perc = sum(x<percentage)/length(x); perc>((100-percentage)/100)})
filtered_table = taxa_level[abundant.phyla,]
filtered_rare_table = taxa_level[rare.phyla,]

# Sustract >= phyla from 100% to get percentage of rare_taxa
rare_taxa = 100-colSums(filtered_table)

# Add rare_taxa as column to table
filtered_table = rbind(filtered_table, rare_taxa)

#Prepare table with percentages
write.table(filtered_table, "percentagesbarplots_alldataset.csv", sep = ",", row.names = TRUE, col.names = TRUE)

#Bar charts####
#Chunk of script to produce figure S1a and supplementary table S6. 
#Bar charts were created separately for each group and panel figures were prepared in inkscape.

#### Group A ####
#Subset samples in ampvis2
group.a <- amp_subset_samples(dataset_BLCA, group %in% c("A"))
#After: 50 samples and 3607 OTUs

# Export data
amp_export_otutable(group.a, "table_barcharts_a", sep = "\t")
write.table(group.a$metadata, file = "metadata_barcharts_a.csv", sep = "\t")

# Load data
otu = read.table("table_barcharts_a.csv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata_barcharts_a.csv", row.names = 1, header = T)

# Check if rows and collumns match (very important, otherwise everything will be mixed in the aggregation step, must be TRUE)
all(rownames(metadata) %in% colnames(otu))

# Collapse/Aggregate samples by metadata column (also used for sorting samples)
taxonomy = paste(otu$Kingdom,';',otu$Phylum,';',otu$Class)
unique(taxonomy)
otu_transposed = t(otu[,1:(dim(otu)[2]-7)])

aggdata = aggregate(otu_transposed, by=list(metadata$month_year), FUN=sum)
#                                                     ^ Argument to aggregate and sort

#Renaming
rownames(aggdata) = aggdata$Group.1
aggdata$Group.1 = NULL

# "Restore" OTU table
otu_tax = as.data.frame(t(aggdata))

# Combine tables
otu_tax = cbind(otu_tax, taxonomy)

# Fix names...
# No Blast Hit
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ; ", replacement = "")
# Separator
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ; ", replacement = "; ")
# Remove accession nr. identity and evalue
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = "~.*", replacement = "")

# Collapse phylogeny according selected level
# Defined above (Order in this case)
taxa_level = aggregate(otu_tax[,1:dim(otu_tax)[2]-1], by = list(otu_tax$taxonomy), FUN = sum)

# Add row names
rownames(taxa_level) = taxa_level[,1]
taxa_level = taxa_level[,2:dim(taxa_level)[2]]

# Transform to relative abundance
taxa_level = scale(taxa_level, center=F, scale=colSums(taxa_level))
taxa_level = taxa_level*100
colSums(taxa_level)

# Group low abundant phyla as rare taxa
# Pick only phyla >= 2%
percentage = 2  
#            ^ number to filter phylum/class/order abundance

abundant.phyla = apply(taxa_level, 1, function(x){perc = sum(x>percentage)/length(x); perc>0})
rare.phyla = apply(taxa_level, 1, function(x){perc = sum(x<percentage)/length(x); perc>((100-percentage)/100)})
filtered_table = taxa_level[abundant.phyla,]
filtered_rare_table = taxa_level[rare.phyla,]

# Sustract >= phyla from 100% to get percentage of rare_taxa
rare_taxa = 100-colSums(filtered_table)

# Add rare_taxa as column to table
filtered_table = rbind(filtered_table, rare_taxa)

# Melt for ggplot
filtered_table_melt = melt(filtered_table)

colnames(filtered_table_melt)[1] = "Taxon"
colnames(filtered_table_melt)[2] = "Month_year"
colnames(filtered_table_melt)[3] = "Abundance"

# Stacked bar charts (ggplot)
barplot = ggplot(filtered_table_melt, aes(x=Month_year, y = Abundance, fill = Taxon))
barplot = barplot + geom_bar(stat = "identity", position = "stack")

# Define colors for each Taxon 
colors = rownames(filtered_table)

#Color palette with 2% rare taxa

Taxa_color=rep('black',nrow(filtered_table))
Taxa_color[rownames(filtered_table)=="Alveolata; Ciliophora; Litostomatea"]='#d15c85'
Taxa_color[rownames(filtered_table)=="Archaeplastida; Streptophyta; Embryophyceae"]='#936aca'
Taxa_color[rownames(filtered_table)=="Excavata; Metamonada; Trichomonadea"]='#6faa41'
Taxa_color[rownames(filtered_table)=="k__Unclassified"]='#d44258'
Taxa_color[rownames(filtered_table)=="Opisthokonta; Fungi; Ascomycota"]='#c49c32'
Taxa_color[rownames(filtered_table)=="Opisthokonta; Fungi; Basidiomycota"]='#7c8bc0'
Taxa_color[rownames(filtered_table)=="Opisthokonta; Metazoa; Arthropoda"]='#5ea164'
Taxa_color[rownames(filtered_table)=="Opisthokonta; Metazoa; Craniata"]='#c8643b'
#Taxa_color[rownames(filtered_table)=="Opisthokonta; Metazoa; Mollusca"]='#9c854a' 
Taxa_color[rownames(filtered_table)=="Opisthokonta; Metazoa; Nematoda"]='#47ae9b'
Taxa_color[rownames(filtered_table)=="rare_taxa"]='#c27271'

# Custom colors
barplot + 
  theme_light() +
  scale_x_discrete(limits=c(order_x_by = c("05_May18", "06_Jun18", "07_July18", "08_Aug18", "09_Sep18", 
                                           "10_Oct18", "11_Nov18", "12_Dec18",
                                           "03_Mar19", "04_Apr19")))+
  xlab("Month")+
  ggtitle("Group A")+
  scale_y_continuous(name="Relative abundance (%)", breaks=seq(0,100, 10)) +
  theme(axis.text.x = element_text(size = 14 , angle = 90, vjust = 0.4, hjust = 1), 
        axis.text.y = element_text(size = 14), panel.grid.major = element_line(size = 0.1),
        legend.title = element_text(size = 18), legend.text = element_text(size = 14),
        axis.title.x = element_text(size =18), axis.title.y = element_text(size =18),
        title = element_text(size = 18)) +
  scale_fill_manual(values=Taxa_color) + geom_bar(size=0.25, stat = "identity", position = "stack", colour="black")

#### Group B ####
#Subset samples in ampvis2
group.b <- amp_subset_samples(dataset_BLCA, group %in% c("B"))
#83 samples and 4785  OTUs

# Export data
amp_export_otutable(group.b, "table_barcharts_b", sep = "\t")
write.table(group.b$metadata, file = "metadata_barcharts_b.csv", sep = "\t")

# Load data
otu = read.table("table_barcharts_b.csv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata_barcharts_b.csv", row.names = 1, header = T)

# Check if rows and collumns match (very important, otherwise everything will be mixed in the aggregation step, must be TRUE)
all(rownames(metadata) %in% colnames(otu))

# Collapse/Aggregate samples by metadata column (also used for sorting samples!!!)
taxonomy = paste(otu$Kingdom,';',otu$Phylum,';',otu$Class)
unique(taxonomy)
otu_transposed = t(otu[,1:(dim(otu)[2]-7)])

aggdata = aggregate(otu_transposed, by=list(metadata$month_year), FUN=sum)
#                                                     ^ Argument to aggregate and sort data

#Renaming
rownames(aggdata) = aggdata$Group.1
aggdata$Group.1 = NULL

# "Restore" OTU table
otu_tax = as.data.frame(t(aggdata))

# Combine tables
otu_tax = cbind(otu_tax, taxonomy)

# Fix names...
# No Blast Hit
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ; ", replacement = "")
# Separator
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ; ", replacement = "; ")
# Remove accession nr. identity and evalue
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = "~.*", replacement = "")

# Collapse phylogeny according selected level
# Defined above (Order in this case)
taxa_level = aggregate(otu_tax[,1:dim(otu_tax)[2]-1], by = list(otu_tax$taxonomy), FUN = sum)

# Add row names
rownames(taxa_level) = taxa_level[,1]
taxa_level = taxa_level[,2:dim(taxa_level)[2]]

# Transform to relative abundance
taxa_level = scale(taxa_level, center=F, scale=colSums(taxa_level))
taxa_level = taxa_level*100
colSums(taxa_level)

# Group low abundant taxon as rare taxa
# Pick only phyla >= 2%
percentage = 2  
#            ^ number to filter phylum/class/order abundance

abundant.phyla = apply(taxa_level, 1, function(x){perc = sum(x>percentage)/length(x); perc>0})
rare.phyla = apply(taxa_level, 1, function(x){perc = sum(x<percentage)/length(x); perc>((100-percentage)/100)})
filtered_table = taxa_level[abundant.phyla,]
filtered_rare_table = taxa_level[rare.phyla,]

#Sustract >= phyla from 100% to get percentage of rare_taxa
rare_taxa = 100-colSums(filtered_table)

# Add rare_taxa as column to table
filtered_table = rbind(filtered_table, rare_taxa)

# Melt for ggplot
filtered_table_melt = melt(filtered_table)

colnames(filtered_table_melt)[1] = "Taxon"
colnames(filtered_table_melt)[2] = "Month_year"
colnames(filtered_table_melt)[3] = "Abundance"

# Stacked bar charts (ggplot)
barplot = ggplot(filtered_table_melt, aes(x=Month_year, y = Abundance, fill = Taxon))
barplot = barplot + geom_bar(stat = "identity", position = "stack")

# Define colors for each taxon ##
colors = rownames(filtered_table)

# ABUNDANT taxa color list #

#Color palette with 2% rare taxa
Taxa_color=rep('black',nrow(filtered_table))
Taxa_color[rownames(filtered_table)=="Alveolata; Ciliophora; Litostomatea"]='#d15c85'
Taxa_color[rownames(filtered_table)=="Archaeplastida; Streptophyta; Embryophyceae"]='#936aca'
Taxa_color[rownames(filtered_table)=="Excavata; Metamonada; Trichomonadea"]='#6faa41'
Taxa_color[rownames(filtered_table)=="k__Unclassified"]='#d44258'
Taxa_color[rownames(filtered_table)=="Opisthokonta; Fungi; Ascomycota"]='#c49c32'
Taxa_color[rownames(filtered_table)=="Opisthokonta; Fungi; Basidiomycota"]='#7c8bc0'
Taxa_color[rownames(filtered_table)=="Opisthokonta; Metazoa; Arthropoda"]='#5ea164'
Taxa_color[rownames(filtered_table)=="Opisthokonta; Metazoa; Craniata"]='#c8643b'
#Taxa_color[rownames(filtered_table)=="Opisthokonta; Metazoa; Mollusca"]='#9c854a' 
Taxa_color[rownames(filtered_table)=="Opisthokonta; Metazoa; Nematoda"]='#47ae9b'
Taxa_color[rownames(filtered_table)=="rare_taxa"]='#c27271'

# Custom colors
barplot + 
  theme_light() +
  scale_x_discrete(limits=c(order_x_by = c("05_May18", "06_Jun18", "07_July18", "08_Aug18", "09_Sep18", 
                                           "10_Oct18", "11_Nov18", "12_Dec18", "01_Jan19", "02_Feb19",
                                           "03_Mar19", "04_Apr19")))+
  ggtitle("Group B")+
  xlab("Month")+
  scale_y_continuous(name="Relative abundance (%)", breaks=seq(0,100, 10)) +
  theme(axis.text.x = element_text(size = 14 , angle = 90, vjust = 0.4, hjust = 1), 
        axis.text.y = element_text(size = 14), panel.grid.major = element_line(size = 0.1),
        legend.title = element_text(size = 18), legend.text = element_text(size = 14),
        axis.title.x = element_text(size =18), axis.title.y = element_text(size =18),
        title = element_text(size = 18)) +
  scale_fill_manual(values=Taxa_color) + geom_bar(size=0.25, stat = "identity", position = "stack", colour="black")

#### Group F ####
#Subset samples in ampvis2
group.f <- amp_subset_samples(dataset_BLCA, group %in% c("F"))
#After: 146 samples and 5090 OTUs

# Export data
amp_export_otutable(group.f, "table_barcharts_f", sep = "\t")
write.table(group.f$metadata, file = "metadata_barcharts_f.csv", sep = "\t")

# Load data
otu = read.table("table_barcharts_f.csv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata_barcharts_f.csv", row.names = 1, header = T)

# Check if rows and collumns match (very important, otherwise everything will be mixed in the aggregation step, must be TRUE)
all(rownames(metadata) %in% colnames(otu))

# Collapse/Aggregate samples by metadata column (also used for sorting samples)
taxonomy = paste(otu$Kingdom,';',otu$Phylum,';',otu$Class)
unique(taxonomy)

otu_transposed = t(otu[,1:(dim(otu)[2]-7)])

aggdata = aggregate(otu_transposed, by=list(metadata$month_year), FUN=sum)
#                                                     ^ Argument to aggregate and sort data

#Renaming
rownames(aggdata) = aggdata$Group.1
aggdata$Group.1 = NULL

# "Restore" OTU table
otu_tax = as.data.frame(t(aggdata))

# Combine tables
otu_tax = cbind(otu_tax, taxonomy)

# Fix names...
# No Blast Hit
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ; ", replacement = "")
# Separator
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ; ", replacement = "; ")
# Remove accession nr. identity and evalue
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = "~.*", replacement = "")

# Collapse phylogeny according selected level
# Defined above (Order in this case)
taxa_level = aggregate(otu_tax[,1:dim(otu_tax)[2]-1], by = list(otu_tax$taxonomy), FUN = sum)

# Add row names
rownames(taxa_level) = taxa_level[,1]
taxa_level = taxa_level[,2:dim(taxa_level)[2]]

# Transform to relative abundance
taxa_level = scale(taxa_level, center=F, scale=colSums(taxa_level))
taxa_level = taxa_level*100
colSums(taxa_level)

# Group low abundant taxon as rare taxa 
# Pick only phyla >= 2%
percentage = 2  
#            ^ number to filter phylum/class/order abundance

abundant.phyla = apply(taxa_level, 1, function(x){perc = sum(x>percentage)/length(x); perc>0})
rare.phyla = apply(taxa_level, 1, function(x){perc = sum(x<percentage)/length(x); perc>((100-percentage)/100)})
filtered_table = taxa_level[abundant.phyla,]
filtered_rare_table = taxa_level[rare.phyla,]

# Sustract >= phyla from 100% to get percentage of rare_taxa
rare_taxa = 100-colSums(filtered_table)

# Add rare_taxa as column to table
filtered_table = rbind(filtered_table, rare_taxa)

#Export table with percentages
write.table(filtered_table, file="Percentagesbarcharts_phyla_f.csv", sep = ",", row.names = TRUE,  col.names = NA)

# Melt for ggplot
filtered_table_melt = melt(filtered_table)

colnames(filtered_table_melt)[1] = "Taxon"
colnames(filtered_table_melt)[2] = "Month_year"
colnames(filtered_table_melt)[3] = "Abundance"

# Stacked bar charts (ggplot)
barplot = ggplot(filtered_table_melt, aes(x=Month_year, y = Abundance, fill = Taxon))
barplot = barplot + geom_bar(stat = "identity", position = "stack")

# Bar plot ##
# Define colors for each taxon
colors = rownames(filtered_table)

#Color palette with 2% rare taxa
Taxa_color=rep('black',nrow(filtered_table))
Taxa_color[rownames(filtered_table)=="Alveolata; Ciliophora; Litostomatea"]='#d15c85'
Taxa_color[rownames(filtered_table)=="Archaeplastida; Chlorophyta; Trebouxiophyceae"]='#6d80d8'
Taxa_color[rownames(filtered_table)=="Archaeplastida; Streptophyta; Embryophyceae"]='#936aca'
Taxa_color[rownames(filtered_table)=="Excavata; Metamonada; Trichomonadea"]='#6faa41'
Taxa_color[rownames(filtered_table)=="k__Unclassified"]='#d44258'
Taxa_color[rownames(filtered_table)=="Opisthokonta; Fungi; Ascomycota"]='#c49c32'
Taxa_color[rownames(filtered_table)=="Opisthokonta; Fungi; Basidiomycota"]='#7c8bc0'
Taxa_color[rownames(filtered_table)=="Opisthokonta; Metazoa; Arthropoda"]='#5ea164'
Taxa_color[rownames(filtered_table)=="Opisthokonta; Metazoa; Craniata"]='#c8643b'
Taxa_color[rownames(filtered_table)=="Opisthokonta; Metazoa; Mollusca"]='#9c854a' 
Taxa_color[rownames(filtered_table)=="Opisthokonta; Metazoa; Nematoda"]='#47ae9b'
Taxa_color[rownames(filtered_table)=="rare_taxa"]='#c27271'

barplot + 
  theme_light() +
  scale_x_discrete(limits=c(order_x_by = c("05_May18", "06_Jun18", "07_July18", "08_Aug18", "09_Sep18", 
                                           "10_Oct18", "11_Nov18", "12_Dec18", "01_Jan19", "02_Feb19",
                                           "03_Mar19", "04_Apr19")))+
  ggtitle("Group F")+
  xlab("Month")+
  scale_y_continuous(name="Relative abundance (%)", breaks=seq(0,100, 10)) +
  theme(axis.text.x = element_text(size = 14 , angle = 90, vjust = 0.4, hjust = 1), 
        axis.text.y = element_text(size = 14), panel.grid.major = element_line(size = 0.1),
        legend.title = element_text(size = 18), legend.text = element_text(size = 14),
        axis.title.x = element_text(size =18), axis.title.y = element_text(size =18),
        title = element_text(size = 18)) +
  scale_fill_manual(values=Taxa_color) + geom_bar(size=0.25, stat = "identity", position = "stack", colour="black")

#### Group J ####
#Subset samples in ampvis2
group.j <- amp_subset_samples(dataset_BLCA, group %in% c("J"))
#After: 101 samples and 4594 OTUs

# Export data
amp_export_otutable(group.j, "table_barcharts_j", sep = "\t")
write.table(group.j$metadata, file = "metadata_barcharts_j.csv", sep = "\t")

# Load data
otu = read.table("table_barcharts_j.csv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata_barcharts_j.csv", row.names = 1, header = T)

# Check if rows and collumns match (very important, otherwise everything will be mixed in the aggregation step, must be TRUE)
all(rownames(metadata) %in% colnames(otu))

# Collapse/Aggregate samples by metadata column (also used for sorting samples)
taxonomy = paste(otu$Kingdom,';',otu$Phylum,';',otu$Class)
unique(taxonomy)
otu_transposed = t(otu[,1:(dim(otu)[2]-7)])

aggdata = aggregate(otu_transposed, by=list(metadata$month_year), FUN=sum)
#                                                     ^ Argument to aggregate and sort samples
#renaming...
rownames(aggdata) = aggdata$Group.1
aggdata$Group.1 = NULL

# "Restore" OTU table
otu_tax = as.data.frame(t(aggdata))

# Combine tables
otu_tax = cbind(otu_tax, taxonomy)

# Fix names...
# No Blast Hit
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ; ", replacement = "")
# Separator
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ; ", replacement = "; ")
# Remove accession nr. identity and evalue
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = "~.*", replacement = "")

# Collapse phylogeny according selected level
# Defined above (Order in this case)
taxa_level = aggregate(otu_tax[,1:dim(otu_tax)[2]-1], by = list(otu_tax$taxonomy), FUN = sum)

# Add row names
rownames(taxa_level) = taxa_level[,1]
taxa_level = taxa_level[,2:dim(taxa_level)[2]]

# Transform to relative abundance
taxa_level = scale(taxa_level, center=F, scale=colSums(taxa_level))
taxa_level = taxa_level*100
colSums(taxa_level)

# Group low abundant taxon as rare taxa 
# Pick only phyla >= 2%
percentage = 2  
#            ^ number to filter phylum/class/order abundance

abundant.phyla = apply(taxa_level, 1, function(x){perc = sum(x>percentage)/length(x); perc>0})
rare.phyla = apply(taxa_level, 1, function(x){perc = sum(x<percentage)/length(x); perc>((100-percentage)/100)})
filtered_table = taxa_level[abundant.phyla,]
filtered_rare_table = taxa_level[rare.phyla,]

# Sustract >= phyla from 100% to get percentage of rare_taxa
rare_taxa = 100-colSums(filtered_table)

# Add rare_taxa as column to table
filtered_table = rbind(filtered_table, rare_taxa)

# Melt for ggplot
filtered_table_melt = melt(filtered_table)

colnames(filtered_table_melt)[1] = "Taxon"
colnames(filtered_table_melt)[2] = "Month_year"
colnames(filtered_table_melt)[3] = "Abundance"

# Stacked bar charts (ggplot)
barplot = ggplot(filtered_table_melt, aes(x=Month_year, y = Abundance, fill = Taxon))
barplot = barplot + geom_bar(stat = "identity", position = "stack")

# Define colors for each Phylum #
colors = rownames(filtered_table)

# Write taxa table for easy concatenation of function and colors in Excel or Libre calc 
write.table(colors, file="colors_bar.txt", col.names=F, row.names=F, quote=FALSE)

# ABUNDANT taxa color list #
#Color palette with 2% rare taxa

Taxa_color=rep('black',nrow(filtered_table))
Taxa_color[rownames(filtered_table)=="Alveolata; Ciliophora; Litostomatea"]='#d15c85'
Taxa_color[rownames(filtered_table)=="Archaeplastida; Streptophyta; Embryophyceae"]='#936aca'
Taxa_color[rownames(filtered_table)=="Excavata; Metamonada; Trichomonadea"]='#6faa41'
Taxa_color[rownames(filtered_table)=="k__Unclassified"]='#d44258'
Taxa_color[rownames(filtered_table)=="Opisthokonta; Fungi; Ascomycota"]='#c49c32'
Taxa_color[rownames(filtered_table)=="Opisthokonta; Fungi; Basidiomycota"]='#7c8bc0'
Taxa_color[rownames(filtered_table)=="Opisthokonta; Metazoa; Arthropoda"]='#5ea164'
Taxa_color[rownames(filtered_table)=="Opisthokonta; Metazoa; Craniata"]='#c8643b'
Taxa_color[rownames(filtered_table)=="Opisthokonta; Metazoa; Platyhelminthes"]='#892863' 
Taxa_color[rownames(filtered_table)=="Opisthokonta; Metazoa; Nematoda"]='#47ae9b'
Taxa_color[rownames(filtered_table)=="rare_taxa"]='#c27271'

#Barplot
barplot + 
  theme_light() +
  scale_x_discrete(limits=c(order_x_by = c("05_May18", "06_Jun18", "07_July18", "08_Aug18", "09_Sep18", 
                                           "10_Oct18", "11_Nov18", "12_Dec18", "01_Jan19", "02_Feb19",
                                           "03_Mar19", "04_Apr19")))+
  ggtitle("Group J")+
  xlab("Month")+
  scale_y_continuous(name="Relative abundance (%)", breaks=seq(0,100, 10)) +
  theme(axis.text.x = element_text(size = 14 , angle = 90, vjust = 0.4, hjust = 1), 
        axis.text.y = element_text(size = 14), panel.grid.major = element_line(size = 0.1),
        legend.title = element_text(size = 18), legend.text = element_text(size = 14),
        axis.title.x = element_text(size =18), axis.title.y = element_text(size =18),
        title = element_text(size = 18)) +
  scale_fill_manual(values=Taxa_color) + geom_bar(size=0.25, stat = "identity", position = "stack", colour="black")
