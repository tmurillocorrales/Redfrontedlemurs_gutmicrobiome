#### Bacteria 16S rRNA gene amplicon data ####

# Dominik Schneider, Dirk Berkelmann and Avril von Hoyningen-Huene contributed to these scripts

# Load libraries ####
library(ampvis2)
library(ape)
library(phangorn)
library(stringr)
library(reshape2)
library(picante)
library(rayshader)
library(ggdark)
library(tibble)
library(viridis)
library(scales)
library(dplyr)
library(tidyverse)

#
# Load data ####

# Load 16S rRNA gene count table
myotutable = read.table("16S_ASV_table.txt", sep = "\t", header=T, check.names=F, comment.char="", skip = 1, quote = "")

# Load sample metadata
mymetadata = read.table("16S_metadata.csv", sep=",", header=T, comment.char = "", quote = "\"")

#Metadata with behavioral data
mymetadata2 = read.table("metadata_alphadiv_socint_para_feed.txt", sep="\t", header=T, comment.char = "", quote = "\"")

# Modify data for ampvis2
# Get taxonomy from OTU table
row.names(myotutable)=myotutable$`#OTU ID`
tax_info = data.matrix(myotutable$taxonomy)
OTUID = row.names(myotutable)

# Construct a data frame from taxonomy info and separate by semicolon
tax_info = str_split_fixed(tax_info, '\\~',5)
tax_info = gsub(x = tax_info, pattern = "^$", replacement = "0")
tax_info = data.frame(tax_info, stringsAsFactors = F)
tax_info$X3 = as.numeric(tax_info$X3)
tax_info$X4 = as.numeric(tax_info$X4)
#Check names
str(tax_info)

# Correct for coverage (at least 90% coverage to query)
tax_info$X1[ tax_info$X4 < 90 ] = "No blast hit"
tax_info = data.frame(cbind(tax_info$X3,tax_info$X1), stringsAsFactors = F)

tax_info = cbind(tax_info$X1, str_split_fixed(tax_info$X2, '\\; ',7))
tax_info = gsub(x = tax_info, pattern = "~.*", replacement = "")
tax_info = gsub(x = tax_info, pattern = "^$", replacement = "unclassified")
tax_info = gsub(x = tax_info, pattern = "_", replacement = " ")
tax_info = data.frame(tax_info, stringsAsFactors = F)
tax_info$X1 = as.numeric(tax_info$X1)

# Apply Yarza recommended thresholds for taxonomical classification from 16S amplicons 
#(https://www.nature.com/articles/nrmicro3330)
#(<98.7 Species, <94.5 Genus,<86.5 Family, <82.0 Order, <78.5 Class, <75 Phylum)
#
# Species
tax_info$X8[ tax_info$X1 < 98.7 ] = paste0("Unclassified (o__", tax_info$X5[ tax_info$X1 < 98.7 ], ";f__", tax_info$X6[ tax_info$X1 < 98.7 ], ";g__", tax_info$X7[ tax_info$X1 < 98.7 ],")")
# Genus
tax_info$X7[ tax_info$X1 < 94.5 ] = paste0("Unclassified (o__", tax_info$X5[ tax_info$X1 < 94.5 ],";f__", tax_info$X6[ tax_info$X1 < 94.5 ],")")
# Family
tax_info$X6[ tax_info$X1 < 86.5 ] = paste0("Unclassified (p__", tax_info$X3[ tax_info$X1 < 86.5 ],";o__", tax_info$X5[ tax_info$X1 < 86.5 ],")")
# Order
tax_info$X5[ tax_info$X1 < 82.0 ] = paste0("Unclassified (p__", tax_info$X3[ tax_info$X1 < 82.0 ],";c__", tax_info$X4[ tax_info$X1 < 82.0 ],")")
# Class
tax_info$X4[ tax_info$X1 < 78.5 ] = paste0("Unclassified (p__", tax_info$X5[ tax_info$X1 < 78.5 ],")")
# Phylum
tax_info$X3[ tax_info$X1 < 75.0 ] = paste0("No blast hit")

# Add column names
colnames(tax_info) = c("identity", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_info$identity = NULL
rownames(tax_info) = OTUID

# Remove old taxonomy from OTU table
myotutable$taxonomy = NULL

# Add new taxonomy by merging both tables
myotutable = cbind(myotutable, tax_info)
myotutable$tax_info = NULL

# Rename first column
colnames(myotutable)[1] = "OTU"
i = sapply(myotutable, is.factor)
myotutable[i] = lapply(myotutable[i], as.character)

# Combine the data to ampvis2 object
dataset_yarza = amp_load(otutable = myotutable,
                         metadata = mymetadata)
dataset_yarza
#909 samples and 7285 OTUS

# Remove domains which are not Bacteria
dataset_yarza2 = amp_subset_taxa(dataset_yarza, tax_vector=c("Chloroplast", "Eukaryota", "Mitochondria", "Archaea"), remove = T, normalise = F)
#After: 7235 OTUs

#Remove samples with read counts
dataset_yarza2 = amp_subset_samples(dataset_yarza2, !nsample %in% c("759", "758", "784", "721", "878", "653", "672", "742", 
                                                                  "878", "653", "672", "742", "680", "730", "693", "701", 
                                                                  "520", "715", "696", "714", "547", "612", "588", "910", "767"))
#After: 888 samples and 7215 OTUs

#Remove samples from RNA generated amplicons, as they belong to a previous study (Murillo et al, 2022; https://www.nature.com/articles/s43705-021-00086-0).
dataset_yarza2 <- amp_subset_samples(dataset_yarza2, !nucleicacid %in% c("RNA"))
#After: 826 samples and 7213 OTUs

#REmove data collected on April 2018 as they are not part of this study.
dataset_yarza2 <- amp_subset_samples(dataset_yarza2, !month_year %in% c("04_Apr18"))
#After: 799 samples and 7213 OTUs

# Filter reads <0.25% relative abundance ####
# Remove spurious reads as recommended in https://www.nature.com/articles/s43705-021-00033-z.
# Metadata and OTU tables should be sorted by sample name
#### Table with counts but without rare taxa (<0.25%)
dataset_subset = dataset_yarza2

# Export data
amp_export_otutable(dataset_subset, "table_filtersamples", sep = "\t", extension = "tsv")
write.table(dataset_subset$metadata, file = "metadata_filtersamples.tsv", sep = "\t")

# Load data
otu = read.table("table_filtersamples.tsv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata_filtersamples.tsv", row.names = 1, header = T)

# Check if rows and collumns match (very important, otherwise everything will be mixed in the aggregation step, must be TRUE)
all(rownames(metadata) %in% colnames(otu))

# Collapse/Aggregate taxonomy
taxonomy = paste0(otu$Kingdom,'; ',otu$Phylum,'; ',otu$Class,'; ',otu$Order,'; ',otu$Family,'; ',otu$Genus,'; ',otu$Species,'; ',row.names(otu))

# Combine tables while removing old taxonomy string
otu_tax = cbind(otu[,1:(dim(otu)[2]-7)], taxonomy)

# Collapse phylogeny according selected level
taxa_level = aggregate(otu_tax[,1:dim(otu_tax)[2]-1], by = list(otu_tax$taxonomy), FUN = sum) #different

# Sanity check
colSums(taxa_level[,2:dim(taxa_level)[2]])

# Add row names
rownames(taxa_level) = taxa_level[,1]
taxa_level = taxa_level[,2:dim(taxa_level)[2]]

# Transform to relative abundance
taxa_level_RA = (scale(taxa_level, center=F, scale=colSums(taxa_level)))*100
#Check sum
colSums(taxa_level_RA)

# Remove ASVs with abundances < 0.25%
percentage = 0.25
#            ^ number to filter ASV abundance 
#(0.25 = everything below 0.25% will be summarized as "rare taxa" in this case)

abundant.phyla = apply(taxa_level_RA, 1, function(x){perc = sum(x>percentage)/length(x); perc>0})
taxa_level_filtered = taxa_level[abundant.phyla,]

#Prepare table for ampvis2
#Bind filtered table with ASV info
table_counts = data.frame(taxa_level_filtered, check.names = F)
table_counts$taxonomy = row.names(table_counts)
myotutable = rownames_to_column(table_counts, "OTU")

# Modify data for ampvis2
# Get taxonomy from OTU table
row.names(myotutable)=myotutable$OTU
tax_info = data.frame(myotutable$taxonomy)
OTUID = row.names(myotutable)

tax_info = cbind(str_split_fixed(tax_info$myotutable.taxonomy, '\\; ',8))
tax_info = data.frame(tax_info, stringsAsFactors = F)

# Add column names
colnames(tax_info) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species","ASV")
rownames(tax_info) = OTUID

# Remove old taxonomy from OTU table
myotutable$taxonomy = NULL

# Add new taxonomy (merge both tables)
myotutable = cbind(myotutable, tax_info)
myotutable$tax_info = NULL

# Rename first column
i = sapply(myotutable, is.factor)
myotutable[i] = lapply(myotutable[i], as.character)
myotutable$OTU = myotutable$ASV
myotutable$ASV = NULL

# Load data to ampvis2 ####
#Load phylogenetic tree created aligning all sequences with MAFFT v7.407-1 at 100 iterations, 
#calculated using FastTreeMP v2.1.7 and midpoint-rooted using FigTreev 1.4.4. 
tree = read.tree("tree_16S")

# Combine the data to ampvis2 object
dataset2 = amp_load(otutable = myotutable,
                    metadata = mymetadata2,
                    tree = tree)
dataset2
#After filtering and removing spuriours reads -> 799 samples, 1028 OTUs and minimum reads 8236.

### Rarefy ####
#Rarefy to minimum amount of reads (8236)
subset_rarefy2 <- amp_subset_samples(dataset2, minreads = 8236, rarefy = 8236, normalise = FALSE, removeAbsents = TRUE)
#Check object
subset_rarefy2
#799 samples

### Calculate PCOA ####
# All samples ####

#PCOA showed in figure 3a -> sample color by month_year and 3b -> sample color by sex.
PCOA1 = amp_ordinate(subset_rarefy2,
                     num_threads = 4,
                     type = "PCOA",
                     distmeasure = "wunifrac",
                     sample_color_by = "month_year",
                     #sample_color_by = "sex",
                     sample_shape_by = "group",
                     transform = "none",
                     sample_colorframe = FALSE,
                     sample_point_size = 4,
                     detailed_output = T)


PCOA1

#Permanova with behavioral data #####
#Remove samples without feeding and social interaction data points (first dataset).
dataset_perma <- amp_subset_samples(dataset2, !nsample %in% c("102","105", "106", "110", "117", "118", "123", "137", "150", "741",
                                                              "748", "750", "773", "774", "776", "777", "781", "832", "833", "835",
                                                              "89", "90", "92",  "93", "95", "98"))

dataset_perma
#After: 773 samples, 1027 OTUs and min reads 8236

dataset_perma_rarefy <- amp_subset_samples(dataset_perma, minreads = 8236, rarefy = 8236, normalise = FALSE, removeAbsents = TRUE)

### Calculate PCOA ####
#This figure was not included in the manuscript but PCOA was calculated to obtain weighted unifrac matrix.
PCOA2 = amp_ordinate(dataset_perma_rarefy,
                     num_threads = 4,
                     type = "PCOA",
                     distmeasure = "wunifrac",
                     sample_color_by = "group",   
                     sample_shape_by = "group",
                     transform = "none",
                     sample_colorframe = FALSE,
                     sample_point_size = 4,
                     detailed_output = T)


PCOA2

### Permanova #######

## Extract metadata from PCOA ###
metadata2_all = PCOA2$dsites

### Extract weighted unifrac matrix from PCOA ##
wunifrac2_all = as.matrix(PCOA2$plot$plot_env$distmatrix)

### Calculate permanova #
permax = adonis(wunifrac2_all ~ group + soc_int + sex + age_months + rain + fle_prop + ffl_prop + ffr_prop,
                metadata2_all, permutations = 10000, strata = metadata2_all$individual)

permax
#Results:
#           Df SumsOfSqs  MeanSqs F.Model      R2     Pr(>F)    
#group        3    0.3046 0.101534  9.6410 0.03459 0.00009999 ***
#soc_int      1    0.0374 0.037416  3.5528 0.00425    0.01250 *  
#sex          1    0.0468 0.046784  4.4423 0.00531 0.00009999 ***
#age_months   1    0.0166 0.016598  1.5760 0.00188    0.07439 .  
#rain         1    0.1864 0.186360 17.6955 0.02116 0.00009999 ***
#fle_prop     1    0.0407 0.040716  3.8661 0.00462    0.00220 ** 
#ffl_prop     1    0.0730 0.073022  6.9337 0.00829 0.00009999 ***
#ffr_prop     1    0.0757 0.075739  7.1916 0.00860 0.00009999 ***
#Residuals  762    8.0250 0.010532         0.91129               
#Total      772    8.8062                  1.00000               

#Create dataframe with results
permanova_soc_feed = as.data.frame(permax$aov.tab)

#### Multiple testing correction ##
#Order pvalues from min to max
pvals = permax$aov.tab$`Pr(>F)`
pvals
#Set pvalues to numeric
pvals = as.numeric(pvals)

#Calculate BH correction ##
BH  = p.adjust(pvals, method = "BH")
#Add corrected pvalues to data frame
permanova_soc_feed$BH = BH
#Export results which are shown in supplementary table S10.
write.table(permanova_soc_feed, file = "permanova_soc_feed.csv", sep = ",", row.names = TRUE, col.names = TRUE)

## Permanova for parasite richness #####
#Remove samples without parasite richness data (second dataset).
dataset_para <- amp_subset_samples(dataset2, !nsample %in% c("89", "90", "92", "93", "95", "96", "98", "102", "105", "106", "110", "113", 
                                                             "117", "118","123","127","131","132", "137",
                                                             "150", "151", "153","155","169","179","188",
                                                             "205","230", "240", "264","476","505", "509","528","529","534",
                                                             "550", "569","571","582","593","596","622","625","630","635",
                                                             "645", "646","648","649","655","661","663","664","669","673",
                                                             "675","678","681","682","685","688","691","695","708","716","723",
                                                             "724","725","728","729","732", "737", "738", "739", "741", "748", "750","771","773", 
                                                             "774","776", "777", "781", "807", "828", "832", "833", "835", "839","846", "852","853",
                                                             "856","857","866","867","868",
                                                             "879","880","881","882","884", "887","894","900","912","916","920",
                                                             "924","926", "929","931","932","933","939", "942"))

dataset_para
#After: 682 samples, 1027 OTUs and min reads 9855

#Rarefy to minimum amount of reads
dataset_para_rarefy <- amp_subset_samples(dataset_para, minreads = 9855, rarefy = 9855, normalise = FALSE, removeAbsents = TRUE)

### PCOA ####
#This figure was not included in the manuscript but PCOA was calculated to obtain weighted unifrac matrix.
PCOA3 = amp_ordinate(dataset_para_rarefy,
                     num_threads = 4,
                     type = "PCOA",
                     distmeasure = "wunifrac",
                     sample_color_by = "month_year",
                     sample_shape_by = "group",
                     transform = "none",
                     sample_colorframe = FALSE,
                     sample_point_size = 3,
                     detailed_output = T)


PCOA3

### Permanova #######
## Extract metadata from PCOA ###
metadata3_all = PCOA3$dsites
## Extract wunifrac matrix ##
wunifrac3_all = as.matrix(PCOA3$plot$plot_env$distmatrix)
## Calculate PCOA ##
perma3 = adonis(wunifrac3_all ~ group + soc_int + sex + age_months + richness_para + rain + fle_prop + ffl_prop + ffr_prop,
                metadata3_all, permutations = 10000, strata = metadata3_all$individual)

perma3

#Results:
#               Df SumsOfSqs  MeanSqs F.Model      R2     Pr(>F)    
#group           3    0.2949 0.098315 10.1607 0.04077 0.00009999 ***
#soc_int         1    0.0402 0.040180  4.1525 0.00555     0.0038 ** 
#sex             1    0.0363 0.036325  3.7541 0.00502 0.00009999 ***
#age_months      1    0.0139 0.013912  1.4378 0.00192     0.0045 ** 
#richness_para   1    0.0456 0.045620  4.7147 0.00631 0.00009999 ***
#rain            1    0.1771 0.177127 18.3057 0.02449 0.00009999 ***
#fle_prop        1    0.0367 0.036673  3.7901 0.00507     0.0018 ** 
#ffl_prop        1    0.0399 0.039893  4.1228 0.00551     0.0003 ***
#ffr_prop        1    0.0660 0.066032  6.8243 0.00913     0.0002 ***
#Residuals     670    6.4829 0.009676         0.89622               
#Total         681    7.2337                  1.00000  

#Create data frame with results
permanova_para = as.data.frame(perma3$aov.tab)

#### Multiple testing correction ###
#Order pvalues from min to max
pvals3 = perma3$aov.tab$`Pr(>F)`
pvals3
#Set pvalues to numeric
pvals3 = as.numeric(pvals3)

#Calculate BH correction ###
BH3  = p.adjust(pvals3, method = "BH")
#Add corrected pvalues to data frame.
permanova_para$BH = BH3
#Export results which are shown in supplementary table S11.
write.table(permanova_para, file = "permanova_para.csv", sep = ",", row.names = TRUE, col.names = TRUE)

# Alpha diversity ####
# Calculate alpha diversity indices
alphadiv = amp_alphadiv(subset_rarefy2, measure= c("observed","shannon","simpson","invsimpson"),richness = TRUE)

# Calculate Faiths Phylogenetic diversity
# load phylogenetic tree
OTU_tree = read.tree("tree")

# Calculate dataframe with PD and species richness
#Currently this script works with unrooted trees, if yours is rooted you may want to include the root
faithsPD = pd(t(subset_rarefy2$abund), OTU_tree, include.root = T)

# combine tables for measures
alpha = cbind(alphadiv, faithsPD)

### Organize month data according to sampling timeline
alpha$month_year <- factor(alpha$month_year, 
                           levels = c("05_May18", "06_Jun18","07_July18",
                                      "08_Aug18","09_Sep18", "10_Oct18", "11_Nov18", "12_Dec18", 
                                      "01_Jan19", "02_Feb19", "03_Mar19","04_Apr19"))


# Write metadata file with diversity indices
write.table(alpha, file = "metadata_incl_alphadiversity_0.25.txt", sep = "\t", row.names = F)

#Stats ###
stats = alpha %>% group_by(month_year, group) %>% summarise(min.PD = min(PD), max.PD = max(PD), mean.PD = mean(PD), sd.PD = sd(PD))
write.table(stats, file = "stats_alpha_PD.csv", sep = ",", row.names = F)

#Box plot to prepare figure 1c ####
d2=ggplot(alpha, aes(x= month_year, y=PD, fill=month_year)) +
  geom_boxplot(aes()) +
  theme_gray()+
  labs(y = "PD") + 
  scale_fill_manual(values = c("05_May18" = "#bc7369", "06_Jun18" = "#bc7369", "07_July18" = "#bc7369",
                               "08_Aug18" = "#bc7369", "09_Sep18" = "#bc7369", "10_Oct18" = "#bc7369",
                               "11_Nov18" = "#5ea07c", "12_Dec18" = "#5ea07c", "01_Jan19" = "#5ea07c",
                               "02_Feb19" = "#5ea07c", "03_Mar19" = "#5ea07c", "04_Apr19" = "#bc7369")) +
  scale_y_continuous(limits = c(15, 60)) +
  theme_gray(base_size = 14)+ 
  facet_wrap( ~ group, ncol=4) +
  theme(axis.title.x = element_blank(),
        legend.position = "none", axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 11), 
        axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5))
d2

# Tests on samples with glucocorticoid data ####
#Subset glucocorticoid data
dataset_fgc <- amp_subset_samples(dataset2, !n11oxo_CM_wet_feces %in% c("no"))
#After: 641 samples and 1027 OTUs

#Remove samples with no feeding, social interactions and parasite data (third dataset) 
dataset_fgc <- amp_subset_samples(dataset_fgc, !nsample %in% c("89", "90", "92", "93", "95", "96", "98", "102", "105", "106", "110", "113", 
                                                               "117", "118","123","127","131","132", "137",
                                                               "150", "151", "153","155","169","179","188",
                                                               "205","230", "240", "264","476","505", "509","528","529","534",
                                                               "550", "569","571","582","593","596","622","625","630","635",
                                                               "645", "646","648","649","655","661","663","664","669","673",
                                                               "675","678","681","682","685","688","691","695","708","716","723",
                                                               "724","725","728","729","732", "737", "738", "739", "741", "748", "750","771","773", 
                                                               "774", "776", "777", "781", "807", "828", "832", "833", "835", "839","846", "852","853",
                                                               "856","857","866","867","868",
                                                               "879","880","881","882","884", "887","894","900","912","916","920",
                                                               "924","926", "929","931","932","933","939", "942"))

#After 547 samples, 1026 OTUs and 9855 min reads

#Rarefy to minimum amount of reads
dataset_fgc_rarefy <- amp_subset_samples(dataset_fgc, minreads = 9855, rarefy = 9855, normalise = FALSE, removeAbsents = TRUE)
#Convert column with glucocorticoid data to numeric
dataset_fgc_rarefy$metadata$n11oxo_CM_wet_feces = as.numeric(dataset_fgc_rarefy$metadata$n11oxo_CM_wet_feces)

### PCOA ####
#This figure was not included in the manuscript but PCOA was calculated to obtain weighted unifrac matrix.
PCOA4 = amp_ordinate(dataset_fgc_rarefy,
                     num_threads = 4,
                     type = "PCOA",
                     distmeasure = "wunifrac",
                     sample_color_by = "month_year",
                     sample_shape_by = "group",
                     transform = "none",
                     sample_colorframe = FALSE,
                     sample_point_size = 3,
                     detailed_output = T)


PCOA4

### Permanova #######

##Extract metadata from PCOA ##
metadata4_all = PCOA4$dsites
### Extract weighted unifrac matrix from PCOA##
wunifrac4_all = as.matrix(PCOA4$plot$plot_env$distmatrix)
#Calculate PCOA
perma4 = adonis(wunifrac4_all ~ group + soc_int + sex + age_months + n11oxo_CM_wet_feces + rain + rain + richness_para + fle_prop + ffl_prop + ffr_prop,
                metadata4_all, permutations = 10000, strata = metadata4_all$individual)
perma4
#Results:
#                     Df SumsOfSqs  MeanSqs F.Model      R2     Pr(>F)    
#group                 3    0.3096 0.103212 10.8534 0.05212 0.00009999 ***
#soc_int               1    0.0308 0.030793  3.2381 0.00518  0.0207979 *  
#sex                   1    0.0399 0.039935  4.1995 0.00672 0.00009999 ***
#age_months            1    0.0178 0.017838  1.8758 0.00300  0.0008999 ***
#n11oxo_CM_wet_feces   1    0.1665 0.166542 17.5129 0.02803 0.00009999 ***
#rain                  1    0.1336 0.133553 14.0439 0.02248 0.00009999 ***
#richness_para         1    0.0328 0.032803  3.4494 0.00552  0.0003000 ***
#fle_prop              1    0.0237 0.023717  2.4940 0.00399  0.0120988 *  
#ffl_prop              1    0.0358 0.035774  3.7619 0.00602  0.0008999 ***
#ffr_prop              1    0.0719 0.071874  7.5580 0.01210 0.00009999 ***
#Residuals           534    5.0782 0.009510         0.85482               
#Total               546    5.9406                  1.00000  

#Create data frame with results
permanova_fgc  = as.data.frame(perma4$aov.tab)

# Multiple testing correction ##
#Order pvalues from min to max
pvals4 = perma4$aov.tab$`Pr(>F)`
pvals4
#Convert pvalues to numeric
pvals4 = as.numeric(pvals4)

#Calculate BH correction #
BH4  = p.adjust(pvals4, method = "BH")
#Add BH pvalue to data frame
permanova_fgc$BH = BH4
#Export results which are shown in supplementary table S12.
write.table(permanova_fgc, file = "permanova_fgc.csv", sep = ",", row.names = TRUE, col.names = TRUE)

## Boxplots of monthly fluctuations on fecal glucorticoids levels ####

### Organize months according to sampling timeline
metadata4_all$month_year <- factor(metadata4_all$month_year, 
                                   levels = c("05_May18", "06_Jun18","07_July18",
                                              "08_Aug18","09_Sep18", "10_Oct18", "11_Nov18", "12_Dec18", 
                                              "01_Jan19", "02_Feb19", "03_Mar19","04_Apr19"))

#Boxplot to create figure 1e####
fgc=ggplot(metadata4_all, aes(x= month_year, y=n11oxo_CM_wet_feces, fill=month_year)) +
  geom_boxplot(aes()) +
  theme_gray()+
  labs(y = "fGCM") + 
  scale_fill_manual(values = c("05_May18" = "#bc7369", "06_Jun18" = "#bc7369", "07_July18" = "#bc7369",
                               "08_Aug18" = "#bc7369", "09_Sep18" = "#bc7369", "10_Oct18" = "#bc7369",
                               "11_Nov18" = "#5ea07c", "12_Dec18" = "#5ea07c", "01_Jan19" = "#5ea07c",
                               "02_Feb19" = "#5ea07c", "03_Mar19" = "#5ea07c", "04_Apr19" = "#bc7369")) +
  theme_gray(base_size = 14)+ 
  facet_wrap( ~ group, ncol=4) +
  ylim(0, 2000) +
  theme(axis.title.x = element_blank(),
        legend.position = "none", axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), 
        axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5))
fgc

stats_fgc = metadata4_all %>% group_by(month_year, group) %>% summarise(mean.fGC = mean(n11oxo_CM_wet_feces), sd.fGC = sd(n11oxo_CM_wet_feces))

# MaAslin2 ####
# Prepare abundance table
# Export data

#All samples####
maaslin2_data = dataset2

amp_export_otutable(maaslin2_data, "table", sep = "\t", extension = "tsv")
write.table(maaslin2_data$metadata, file = "metadata.tsv", sep = "\t")

# Load data back in
otu = read.table("table.tsv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata.tsv", row.names = 1, header = T)

# Check if rows and collumns match (very important, otherwise everything will be mixed in the aggregation step, must be TRUE)
all(rownames(metadata) %in% colnames(otu))

# Collapse/Aggregate taxonomy
taxonomy = paste0(otu$Kingdom,'; ',otu$Phylum,'; ',otu$Class,'; ',otu$Order,'; ',otu$Family,'; ',otu$Genus)

# Combine tables while removing old taxonomy string
otu_tax = cbind(otu[,1:(dim(otu)[2]-7)], taxonomy)

# Collapse phylogeny according selected level
taxa_level = aggregate(otu_tax[,1:dim(otu_tax)[2]-1], by = list(otu_tax$taxonomy), FUN = sum)

# Add row names
rownames(taxa_level) = taxa_level[,1]
taxa_level = taxa_level[,2:dim(taxa_level)[2]]

#Transpose otu table
otu_data = as.data.frame(t(taxa_level))

# prepare metadata
maaslin2_metadata = metadata

library(Maaslin2)

#Run MaAsLin
fit_data_random = Maaslin2(
  input_data = otu_data, 
  input_metadata = maaslin2_metadata, 
  output = "results_maaslin_soc.int_diet_new", 
  fixed_effects = c("group", "soc_int", "richness_para", "age_months", "sex", "fle_prop", "ffl_prop", "ffr_prop", "rain"),
  random_effects = c("individual"),
  reference = c("group,A", "sex,male"),
  normalization = "CLR",
  cores = 5, 
  max_significance = 0.05)

#Results are shown as a heatmap in figure 4a and in supplementary table S13.

## Samples with fecal glucocorticoid measurements ####

# Filter samples in ampvis2
#
maaslin2_data <- amp_subset_samples(dataset2, !n11oxo_CM_wet_feces %in% c("no"))
# After 641 samples 1027 OTUs

# Not the nicest solution, but hey it works
amp_export_otutable(maaslin2_data, "table", sep = "\t", extension = "tsv")
write.table(maaslin2_data$metadata, file = "metadata.tsv", sep = "\t")

# Load data back in
otu = read.table("table.tsv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata.tsv", row.names = 1, header = T)

# Check if rows and collumns match (very important, otherwise everything will be mixed in the aggregation step, must be TRUE)
all(rownames(metadata) %in% colnames(otu))

# Collapse/Aggregate taxonomy
taxonomy = paste0(otu$Kingdom,'; ',otu$Phylum,'; ',otu$Class,'; ',otu$Order,'; ',otu$Family,'; ',otu$Genus)

# Combine tables while removing old taxonomy string
otu_tax = cbind(otu[,1:(dim(otu)[2]-7)], taxonomy)

# Collapse phylogeny according selected level
taxa_level = aggregate(otu_tax[,1:dim(otu_tax)[2]-1], by = list(otu_tax$taxonomy), FUN = sum)

# Add row names
rownames(taxa_level) = taxa_level[,1]
taxa_level = taxa_level[,2:dim(taxa_level)[2]]

#Transpose otu table
otu_data = as.data.frame(t(taxa_level))

# prepare metadata
maaslin2_metadata = metadata

#Run MaAsLin
fit_data_random = Maaslin2(
  input_data = otu_data, 
  input_metadata = maaslin2_metadata, 
  output = "results_maaslin_soc.int_diet_fGCM", 
  fixed_effects = c("group", "soc_int", "richness_para", "age_months", "sex", "n11oxo_CM_wet_feces", "fle_prop", "ffl_prop", "ffr_prop", "rain"),
  random_effects = c("individual"),
  reference = c("group,A", "sex,male"),
  normalization = "CLR",
  cores = 5, 
  max_significance = 0.05)

#Results are shown as a heatmap in figure 4b and in supplementary table S14.

#Normalized data #####
#Remove samples with low read counts

dataset_yarza
#909 samples

#Remove samples with low read counts
dataset = amp_subset_samples(dataset_yarza, !nsample %in% c("759", "758", "784", "721", "878", "653", "672", "742", 
                                                      "878", "653", "672", "742", "680", "730", "693", "701", "767",
                                                      "520", "715", "696", "714", "547", "612", "588", "910"))

#After: 888 samples and 7215 OTUs

#Normalization with GMPR #
library(GMPR)
# Add GMPR
gmpr.size.factor = GMPR(t(dataset$abund), min_ct = 2, intersect_no = 4) # ignore error message in Windows, 
#it seems to lead to the same results as Linux (without error there)
otu.tab.norm = data.frame(t(t(dataset$abund) / gmpr.size.factor), check.names = FALSE)
dataset$GMPR = otu.tab.norm
dataset$READ_COUNTS = dataset$abund

# Replace read counts with GMPR
dataset$abund = dataset$GMPR
dataset$abund = dataset$READ_COUNTS

# Remove extrinsic domains etc.
dataset = amp_subset_taxa(dataset, tax_vector=c("Chloroplast", "Eukaryota", "Mitochondria", "Archaea"), 
                          remove=T, normalise = F)

#50 OTUs were removed

dataset
#After 889 samples, 7235 OTUs and min reads 8838

#
# Filtering and subsetting ##

#Remove samples from RNA generated amplicons, as they belong to a previous study (Murillo et al, 2022; https://www.nature.com/articles/s43705-021-00086-0).
dataset_yarza2 <- amp_subset_samples(dataset_yarza2, !nucleicacid %in% c("RNA"))
#After: 826 samples and 7213 OTUs

#Remove data collected on April 2018 as they are not part of this study.
dataset_yarza2 <- amp_subset_samples(dataset_yarza2, !month_year %in% c("04_Apr18"))
#After: 799 samples and 7213 OTUs

#Bar charts####
#Chunk of script to produce figure 1a and supplementary table S5. 
#Bar charts were created separately for each group and panel figures were prepared in inkscape.

#### Group A ####
#Subset samples in ampvis2
group.a <- amp_subset_samples(dataset_subset, group %in% c("A"))
#After: 132 samples and 6973 OTUs

# Export data
amp_export_otutable(group.a, "otu_normalized_counts_groupa", sep = "\t")
write.table(group.a$metadata, file = "metadata_normalized_counts_groupa.csv", sep = "\t")

# Load data with normalized counts from GMPR
otu = read.table("otu_normalized_counts_groupa.csv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata_normalized_counts_groupa.csv", row.names = 1, header = T)

# Check if rows and collumns match (very important, otherwise everything will be mixed in the aggregation step, must be TRUE)
all(rownames(metadata) %in% colnames(otu))

# Collapse/Aggregate samples by metadata column
taxonomy = paste0(otu$Kingdom,"; ",otu$Phylum)
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

# Sanity check
colSums(taxa_level[,2:dim(taxa_level)[2]])

# Add row names
rownames(taxa_level) = taxa_level[,1]
taxa_level = taxa_level[,2:dim(taxa_level)[2]]

# Transform to relative abundance
#taxa_level_RA = scale(taxa_level, center=F, scale=colSums(taxa_level))
taxa_level = scale(taxa_level, center=F, scale=colSums(taxa_level))
taxa_level = taxa_level*100
colSums(taxa_level)

# Group low abundant phyla as rare taxa
# Pick only phyla >= 1%
percentage = 1
#            ^ number to filter phylum/class/order abundance (1 = everything below 1% will be summarized as "rare taxa" in this case) 

abundant.phyla = apply(taxa_level, 1, function(x){perc = sum(x>percentage)/length(x); perc>0})
rare.phyla = apply(taxa_level, 1, function(x){perc = sum(x<percentage)/length(x); perc>((100-percentage)/100)})
filtered_table = taxa_level[abundant.phyla,]
filtered_rare_table = taxa_level[rare.phyla,]

# Sustract >= phyla from 100% to get percentage of rare_taxa
rare_taxa = 100-colSums(filtered_table)

# Add rare_taxa as column to table
filtered_table = rbind(filtered_table, rare_taxa)

# Melt for barplots
filtered_table_melt = melt(filtered_table)

colnames(filtered_table_melt)[1] = "Taxon"
colnames(filtered_table_melt)[2] = "Month_year"
colnames(filtered_table_melt)[3] = "Abundance"

# Stacked barplots (ggplot)
barplot = ggplot(filtered_table_melt, aes(x=Month_year, y = Abundance, fill = Taxon))
barplot = barplot + geom_bar(stat = "identity", position = "stack")

# Define colors for each Phylum ####
colors = rownames(filtered_table)

#Color palette with 1% rare taxa
Taxa_color=rep('black',nrow(filtered_table))
Taxa_color[rownames(filtered_table)=="Bacteria; Actinobacteriota"]='#c452bb'
Taxa_color[rownames(filtered_table)=="Bacteria; Bacteroidota"]='#51b48c' 
Taxa_color[rownames(filtered_table)=="Bacteria; Campylobacterota"]="#c75d6c"
Taxa_color[rownames(filtered_table)=="Bacteria; Cyanobacteria"]='#7664cd'
Taxa_color[rownames(filtered_table)=="Bacteria; Firmicutes"]='#c99b41'
Taxa_color[rownames(filtered_table)=="Bacteria; Fusobacteriota"]='#7ab743'
Taxa_color[rownames(filtered_table)=="Bacteria; Patescibacteria"]='#618dcc'
Taxa_color[rownames(filtered_table)=="Bacteria; Proteobacteria"]='#d15133'
Taxa_color[rownames(filtered_table)=="Bacteria; Spirochaetota"]='#6a8039'
Taxa_color[rownames(filtered_table)=="Bacteria; Synergistota"]='#d14a73'
Taxa_color[rownames(filtered_table)=="Bacteria; Verrucomicrobiota"]='#d38dcb'
Taxa_color[rownames(filtered_table)=="Bacteria; Fibrobacterota"]='#bb6f54'
Taxa_color[rownames(filtered_table)=="rare_taxa"]='#4aa0d3'

# Custom colors
barplot + 
  theme_light() +
  scale_x_discrete(limits=c(order_x_by = c("05_May18", "06_June18", "07_Jul18", "08_Aug18", "09_Sep18", 
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
group.b <- amp_subset_samples(dataset_subset, group %in% c("B"))
#216 samples and 7108 OTUs

# Export data
amp_export_otutable(group.b, "otu_normalized_counts_groupb", sep = "\t")
write.table(group.b$metadata, file = "metadata_normalized_counts_groupb.csv", sep = "\t")

# Load data with normalized counts from GMPR
otu = read.table("otu_normalized_counts_groupb.csv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata_normalized_counts_groupb.csv", row.names = 1, header = T)

# Check if rows and collumns match (very important, otherwise everything will be mixed in the aggregation step, must be TRUE)
all(rownames(metadata) %in% colnames(otu))

# Collapse/Aggregate samples by metadata column (also used for sorting samples)
taxonomy = paste0(otu$Kingdom,"; ",otu$Phylum)
unique(taxonomy)

otu_transposed = t(otu[,1:(dim(otu)[2]-7)])

aggdata = aggregate(otu_transposed, by=list(metadata$month_year), FUN=sum)
#                                                     ^ Argument to aggregate and sort

# Renaming
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

# Sanity check
colSums(taxa_level[,2:dim(taxa_level)[2]])

# Add row names
rownames(taxa_level) = taxa_level[,1]
taxa_level = taxa_level[,2:dim(taxa_level)[2]]

# Transform to relative abundance
taxa_level = scale(taxa_level, center=F, scale=colSums(taxa_level))
taxa_level = taxa_level*100
colSums(taxa_level)

# Group low abundant phyla as rare taxa
# Pick only phyla >= 1%
percentage = 1  
#            ^ number to filter phylum/class/order abundance (1 = everything below 1% will be summarized as "rare taxa" in this case) #

abundant.phyla = apply(taxa_level, 1, function(x){perc = sum(x>percentage)/length(x); perc>0})
rare.phyla = apply(taxa_level, 1, function(x){perc = sum(x<percentage)/length(x); perc>((100-percentage)/100)})
filtered_table = taxa_level[abundant.phyla,]
filtered_rare_table = taxa_level[rare.phyla,]

# Substract >= phyla from 100% to get percentage of rare_taxa
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

# Define colors for each phylum
colors = rownames(filtered_table)

#Color palette with 1% rare taxa

Taxa_color=rep('black',nrow(filtered_table))
Taxa_color[rownames(filtered_table)=="Bacteria; Actinobacteriota"]='#c452bb'
Taxa_color[rownames(filtered_table)=="Bacteria; Bacteroidota"]='#51b48c' 
Taxa_color[rownames(filtered_table)=="Bacteria; Campylobacterota"]="#c75d6c"
Taxa_color[rownames(filtered_table)=="Bacteria; Cyanobacteria"]='#7664cd'
Taxa_color[rownames(filtered_table)=="Bacteria; Firmicutes"]='#c99b41'
Taxa_color[rownames(filtered_table)=="Bacteria; Fusobacteriota"]='#7ab743'
Taxa_color[rownames(filtered_table)=="Bacteria; Patescibacteria"]='#618dcc'
Taxa_color[rownames(filtered_table)=="Bacteria; Proteobacteria"]='#d15133'
Taxa_color[rownames(filtered_table)=="Bacteria; Spirochaetota"]='#6a8039'
Taxa_color[rownames(filtered_table)=="Bacteria; Synergistota"]='#d14a73'
Taxa_color[rownames(filtered_table)=="Bacteria; Verrucomicrobiota"]='#d38dcb'
Taxa_color[rownames(filtered_table)=="Bacteria; Fibrobacterota"]='#bb6f54'
Taxa_color[rownames(filtered_table)=="rare_taxa"]='#4aa0d3'

# Custom colors
barplot + 
  theme_light() +
  scale_x_discrete(limits=c(order_x_by = c("05_May18", "06_June18", "07_Jul18", "08_Aug18", "09_Sep18", 
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
group.f <- amp_subset_samples(dataset_subset, group %in% c("F"))
#After: 172 samples and 7091 OTUs

# Export data
amp_export_otutable(group.f, "otu_normalized_counts_groupf", sep = "\t")
write.table(group.f$metadata, file = "metadata_normalized_counts_groupf.csv", sep = "\t")

# Load data with normalized counts from GMPR
otu = read.table("otu_normalized_counts_groupf.csv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata_normalized_counts_groupf.csv", row.names = 1, header = T)

# Check if rows and collumns match (very important, otherwise everything will be mixed in the aggregation step, must be TRUE)
all(rownames(metadata) %in% colnames(otu))

# Collapse/Aggregate samples by metadata column (also used for sorting samples)
taxonomy = paste0(otu$Kingdom,"; ",otu$Phylum)
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

# Sanity check
colSums(taxa_level[,2:dim(taxa_level)[2]])

# Add row names
rownames(taxa_level) = taxa_level[,1]
taxa_level = taxa_level[,2:dim(taxa_level)[2]]

# Transform to relative abundance
taxa_level = scale(taxa_level, center=F, scale=colSums(taxa_level))
taxa_level = taxa_level*100
colSums(taxa_level)

#  Group low abundant phyla as rare taxa
# Pick only phyla >= 1%
percentage = 1  
#            ^ number to filter phylum/class/order abundance (1 = everything below 1% will be summarized as "rare taxa" in this case) #

abundant.phyla = apply(taxa_level, 1, function(x){perc = sum(x>percentage)/length(x); perc>0})
rare.phyla = apply(taxa_level, 1, function(x){perc = sum(x<percentage)/length(x); perc>((100-percentage)/100)})
filtered_table = taxa_level[abundant.phyla,]
filtered_rare_table = taxa_level[rare.phyla,]

# Substract >= phyla from 100% to get percentage of rare_taxa
rare_taxa = 100-colSums(filtered_table)

# Add rare_taxa as column to table
filtered_table = rbind(filtered_table, rare_taxa)

# Melt for barplots
filtered_table_melt = melt(filtered_table)

colnames(filtered_table_melt)[1] = "Taxon"
colnames(filtered_table_melt)[2] = "Month_year"
colnames(filtered_table_melt)[3] = "Abundance"

# Stacked bar charts (ggplot)
barplot = ggplot(filtered_table_melt, aes(x=Month_year, y = Abundance, fill = Taxon))
barplot = barplot + geom_bar(stat = "identity", position = "stack")

# Define colors for each phylum
colors = rownames(filtered_table)

#Color palette with 1% rare taxa

Taxa_color=rep('black',nrow(filtered_table))
Taxa_color[rownames(filtered_table)=="Bacteria; Actinobacteriota"]='#c452bb'
Taxa_color[rownames(filtered_table)=="Bacteria; Bacteroidota"]='#51b48c' 
Taxa_color[rownames(filtered_table)=="Bacteria; Campylobacterota"]="#c75d6c"
Taxa_color[rownames(filtered_table)=="Bacteria; Cyanobacteria"]='#7664cd'
Taxa_color[rownames(filtered_table)=="Bacteria; Firmicutes"]='#c99b41'
Taxa_color[rownames(filtered_table)=="Bacteria; Fusobacteriota"]='#7ab743'
Taxa_color[rownames(filtered_table)=="Bacteria; Patescibacteria"]='#618dcc'
Taxa_color[rownames(filtered_table)=="Bacteria; Proteobacteria"]='#d15133'
Taxa_color[rownames(filtered_table)=="Bacteria; Spirochaetota"]='#6a8039'
Taxa_color[rownames(filtered_table)=="Bacteria; Synergistota"]='#d14a73'
Taxa_color[rownames(filtered_table)=="Bacteria; Verrucomicrobiota"]='#d38dcb'
Taxa_color[rownames(filtered_table)=="Bacteria; Fibrobacterota"]='#bb6f54'
Taxa_color[rownames(filtered_table)=="rare_taxa"]='#4aa0d3'

# Custom colors
barplot + 
  theme_light() +
  scale_x_discrete(limits=c(order_x_by = c("05_May18", "06_June18", "07_Jul18", "08_Aug18", "09_Sep18", 
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
group.j <- amp_subset_samples(dataset_subset, group %in% c("J"))
#After: 279 samples and 7148 OTUs

# Export data
amp_export_otutable(group.j, "otu_normalized_counts_groupj", sep = "\t")
write.table(group.j$metadata, file = "metadata_normalized_counts_groupj.csv", sep = "\t")

# Load data with normalized counts from GMPR
otu = read.table("otu_normalized_counts_groupj.csv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata_normalized_counts_groupj.csv", row.names = 1, header = T)

# Check if rows and collumns match (very important, otherwise everything will be mixed in the aggregation step, must be TRUE)
all(rownames(metadata) %in% colnames(otu))

# Collapse/Aggregate samples by metadata column (also used for sorting samples)
taxonomy = paste0(otu$Kingdom,"; ",otu$Phylum)
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

# Sanity check
colSums(taxa_level[,2:dim(taxa_level)[2]])

# Add row names
rownames(taxa_level) = taxa_level[,1]
taxa_level = taxa_level[,2:dim(taxa_level)[2]]

# Transform to relative abundance
taxa_level = scale(taxa_level, center=F, scale=colSums(taxa_level))
taxa_level = taxa_level*100
colSums(taxa_level)

# Group low abundant phyla as rare taxa
# Pick only phyla >= 1%
percentage = 1  
#            ^ number to filter phylum/class/order abundance (1 = everything below 1% will be summarized as "rare taxa" in this case) ####

abundant.phyla = apply(taxa_level, 1, function(x){perc = sum(x>percentage)/length(x); perc>0})
rare.phyla = apply(taxa_level, 1, function(x){perc = sum(x<percentage)/length(x); perc>((100-percentage)/100)})
filtered_table = taxa_level[abundant.phyla,]
filtered_rare_table = taxa_level[rare.phyla,]

# Subtract >= phyla from 100% to get percentage of rare_taxa
rare_taxa = 100-colSums(filtered_table)

# Add rare_taxa as column to table
filtered_table = rbind(filtered_table, rare_taxa)

# Melt for bar plots
filtered_table_melt = melt(filtered_table)

colnames(filtered_table_melt)[1] = "Taxon"
colnames(filtered_table_melt)[2] = "Month_year"
colnames(filtered_table_melt)[3] = "Abundance"

# Stacked bar plots (ggplot)
barplot = ggplot(filtered_table_melt, aes(x=Month_year, y = Abundance, fill = Taxon))
barplot = barplot + geom_bar(stat = "identity", position = "stack")

# Define colors for each phylum 
colors = rownames(filtered_table)

#Color palette with 1% rare taxa

Taxa_color=rep('black',nrow(filtered_table))
Taxa_color[rownames(filtered_table)=="Bacteria; Actinobacteriota"]='#c452bb'
Taxa_color[rownames(filtered_table)=="Bacteria; Bacteroidota"]='#51b48c' 
Taxa_color[rownames(filtered_table)=="Bacteria; Campylobacterota"]="#c75d6c"
Taxa_color[rownames(filtered_table)=="Bacteria; Cyanobacteria"]='#7664cd'
Taxa_color[rownames(filtered_table)=="Bacteria; Firmicutes"]='#c99b41'
Taxa_color[rownames(filtered_table)=="Bacteria; Fusobacteriota"]='#7ab743'
Taxa_color[rownames(filtered_table)=="Bacteria; Patescibacteria"]='#618dcc'
Taxa_color[rownames(filtered_table)=="Bacteria; Proteobacteria"]='#d15133'
Taxa_color[rownames(filtered_table)=="Bacteria; Spirochaetota"]='#6a8039'
Taxa_color[rownames(filtered_table)=="Bacteria; Synergistota"]='#d14a73'
Taxa_color[rownames(filtered_table)=="Bacteria; Verrucomicrobiota"]='#d38dcb'
Taxa_color[rownames(filtered_table)=="Bacteria; Fibrobacterota"]='#bb6f54'
Taxa_color[rownames(filtered_table)=="rare_taxa"]='#4aa0d3'

# Custom colors
barplot + 
  theme_light() +
  scale_x_discrete(limits=c(order_x_by = c("05_May18", "06_June18", "07_Jul18", "08_Aug18", "09_Sep18", 
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

#Line charts to genus level ####
#Chunk of script to produce figure 1b. Panel figures were put together in inkscape.
#Prepare group A separately as is missing two sampling months (January and February)

## Groups B,F,J #####
#Subset groups b, f and j in ampvis2
group.bfj <- amp_subset_samples(dataset_subset, !group %in% c("A"))
#After: 667 samples and 7207 OTUs

# Export data
amp_export_otutable(group.bfj, "otu_linecharts_bfj", sep = "\t")
write.table(group.bfj$metadata, file = "metadata_linecharts_bfj.csv", sep = "\t")

# Load data for groups B,F,J with normalized counts
otu = read.table("otu_linecharts_bfj.csv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata_linecharts_bfj.csv", row.names = 1, header = T)

# Check if rows and columns match (very important, otherwise everything will be mixed in the aggregation step, must be TRUE)
all(rownames(metadata) %in% colnames(otu))

# Collapse/Aggregate samples by metadata column (also used for sorting samples)
taxonomy = paste0(otu$Phylum,"; ",otu$Order,"; ",otu$Class, ",", otu$Family,"; ", otu$Genus)
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

# Remove low abundant genera
# Pick only phyla >= 0.25%
percentage = 0.25
#            ^ number to filter phylum/class/order abundance (0.25 = everything below 0.25% will be summarized as "rare taxa" in this case) ####

abundant.phyla = apply(taxa_level, 1, function(x){perc = sum(x>percentage)/length(x); perc>0})
rare.phyla = apply(taxa_level, 1, function(x){perc = sum(x<percentage)/length(x); perc>((100-percentage)/100)})
filtered_table = taxa_level[abundant.phyla,]
filtered_rare_table = taxa_level[rare.phyla,]

# Sustract >= phyla from 100% to get percentage of rare_taxa
rare_taxa = 100-colSums(filtered_table)

# Add rare_taxa as column to table
filtered_table = rbind(filtered_table, rare_taxa)

# Plot only top 5 most abundant genera
TOP = 5 # change to top X
sub = filtered_table[ order(rowMeans(filtered_table), decreasing = T), ]
sub = sub[1:TOP,]

# Add metadata to combine samples by
table_final = rbind(year = as.character(metadata$year), month = as.character(metadata$month), 
                    individual = as.character(metadata$individual), group = as.character(metadata$group), 
                    month_year = paste0(metadata$month,"_",metadata$year),sub)

table_final = data.frame(t(table_final))

#Melt table for ploting 
filtered_table_melt = melt(table_final, id = c("year","month","individual", "month_year", "group"))

#Check for presence of NAs
which(is.na.data.frame(filtered_table_melt))

#Rename columns
names(filtered_table_melt)[names(filtered_table_melt) == "variable"] = "Taxon"
names(filtered_table_melt)[names(filtered_table_melt) == "value"] = "Abundance"

#Set abundance to numeric
filtered_table_melt$Abundance = as.numeric(filtered_table_melt$Abundance)

#Order the month variable according to the sampling timeline
filtered_table_melt$month_year <- factor(filtered_table_melt$month_year, 
                                         levels = c("5_2018", "6_2018","7_2018",
                                                    "8_2018","9_2018", "10_2018", 
                                                    "11_2018", "12_2018", 
                                                    "1_2019", "2_2019", "3_2019","4_2019"))

#Order the month variable according to the sampling timeline
filtered_table_melt$month <- factor(filtered_table_melt$month, 
                                    levels = c("5", "6","7",
                                               "8","9", "10", 
                                               "11", "12", 
                                               "1", "2", "3","4"))

#Set month to numeric
filtered_table_melt$month= as.numeric(filtered_table_melt$month)

# Create lineplot
lineplot = ggplot(filtered_table_melt, aes(x = month, y = Abundance, fill = Taxon, colour = Taxon, group = month))

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
        axis.text.y = element_text(size = 12))

## Group A ####
#Group A was plotted separately as it was not sampled in January and February
#Subset data in ampvis2
group.a <- amp_subset_samples(dataset_subset, group %in% c("A"))
#After: 132 samples and 6973 OTUs

# Export data
amp_export_otutable(group.a, "otu_normalized_counts_group", sep = "\t")
write.table(group.a$metadata, file = "metadata_normalized_counts_groupa.csv", sep = "\t")

# Load data normalized counts for group A.
otu = read.table("otu_normalized_counts_groupa.csv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata_normalized_counts_groupa.csv", row.names = 1, header = T)

# Check if rows and collumns match (very important, otherwise everything will be mixed in the aggregation step, must be TRUE)
all(rownames(metadata) %in% colnames(otu))

# Collapse/Aggregate samples by metadata column (also used for sorting samples)
taxonomy = paste0(otu$Phylum,"; ",otu$Order,"; ",otu$Class, ",", otu$Family,"; ", otu$Genus)
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

# Remove, Group low abundant genera
# Pick only phyla >= 0.25%
percentage = 0.25
#            ^ number to filter phylum/class/order abundance (0.25 = everything below 1% will be summarized as "rare taxa" in this case) ####

abundant.phyla = apply(taxa_level, 1, function(x){perc = sum(x>percentage)/length(x); perc>0})
rare.phyla = apply(taxa_level, 1, function(x){perc = sum(x<percentage)/length(x); perc>((100-percentage)/100)})
filtered_table = taxa_level[abundant.phyla,]
filtered_rare_table = taxa_level[rare.phyla,]

# Subtract >= phyla from 100% to get percentage of rare_taxa
rare_taxa = 100-colSums(filtered_table)

# Add rare_taxa as column to table
filtered_table = rbind(filtered_table, rare_taxa)

# Plot only top 5 most abundant genera
TOP = 5 # change to top X
sub = filtered_table[ order(rowMeans(filtered_table), decreasing = T), ]
sub = sub[1:TOP,]

# Add metadata to combine samples by
table_final = rbind(year = as.character(metadata$year), month = as.character(metadata$month), 
                    individual = as.character(metadata$individual), group = as.character(metadata$group), 
                    month_year = paste0(metadata$month,"_",metadata$year),sub)


table_final = data.frame(t(table_final))

#Melt table for plotting
filtered_table_melt = melt(table_final, id = c("year","month","individual", "month_year", "group"))

#Check for presence of NAs
which(is.na.data.frame(filtered_table_melt))

#Rename columns
names(filtered_table_melt)[names(filtered_table_melt) == "variable"] = "Taxon"
names(filtered_table_melt)[names(filtered_table_melt) == "value"] = "Abundance"
#Set Abundance to numeric
filtered_table_melt$Abundance = as.numeric(filtered_table_melt$Abundance)

#Order the month variable according to the sampling timeline
filtered_table_melt$month_year <- factor(filtered_table_melt$month_year, 
                                         levels = c("5_2018", "6_2018","7_2018",
                                                    "8_2018","9_2018", "10_2018", 
                                                    "11_2018", "12_2018", 
                                                    "1_2019", "2_2019", "3_2019","4_2019"))


filtered_table_melt$month <- factor(filtered_table_melt$month, 
                                    levels = c("5", "6","7",
                                               "8","9", "10", 
                                               "11", "12", 
                                               "3","4"))


#Set month to numeric
filtered_table_melt$month= as.numeric(filtered_table_melt$month)

# Create lineplot
lineplot = ggplot(filtered_table_melt, aes(x = month, y = Abundance, fill = Taxon, colour = Taxon, group = month))

#Plot
lineplot +
  geom_smooth(method="loess", span = 0.3, fullrange = F, se = T, aes(group = Taxon, fill = Taxon), 
              level = 0.95, size = 1) +
  theme_gray(base_size = 14)  +  
  scale_x_continuous(breaks = 1:length(unique(filtered_table_melt$month)), labels = levels(unique(filtered_table_melt$month))) +
  scale_fill_viridis_d(option = "C", begin= 0, end = 0.8) +
  scale_color_viridis_d(option = "C", begin= 0, end = 0.8) +
  labs(y = "Relative abundance")+
  theme(axis.text.x = element_text(size = 10 , angle = 90, vjust = 0.4, hjust = 1,face = "bold"),
        axis.text.y = element_text(size = 14), panel.grid.major = element_line(size = 0.1),
        axis.title.x = element_blank()) 

# Indicator networks to ASV level ####
# To create indicator networks the packages indicspecies was used.
# The manual can be found here -> https://cran.microsoft.com/snapshot/2016-10-29/web/packages/indicspecies/indicspecies.pdf

library(GMPR)
library(indicspecies)
library(data.table)

dataset = dataset2

#Normalize data ####

# Add GMPR
gmpr.size.factor = GMPR(t(dataset$abund), min_ct = 2, intersect_no = 4) # ignore error message in Windows, 
#it seems to lead to the same results as Linux (without error there)
otu.tab.norm = data.frame(t(t(dataset$abund) / gmpr.size.factor), check.names = FALSE)
dataset$GMPR = otu.tab.norm
dataset$READ_COUNTS = dataset$abund

# Replace read counts with GMPR
dataset$abund = dataset$GMPR
dataset$abund = dataset$READ_COUNTS

#Group A ####
#Subset samples in ampvis2
dataset_a <- amp_subset_samples(dataset, group %in% c("A"))
#After: 132 samples and 1016 OTUs

# Export data
amp_export_otutable(dataset_a, "table_network_a", sep = "\t", extension = "tsv")
write.table(dataset_a$metadata, file = "metadata_network_a.tsv", sep = "\t")

# Load data
samples = read.table("table_network_a.tsv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata_network_a.tsv", row.names = 1, header = T)

# Remove taxonomy
otus = samples [,1:(dim(samples)[2]-7)]
# Collapse phylogeny (Genus example)
taxonomy = paste0(samples$Phylum,'; ',samples$Class,'; ',samples$Order,'; ',samples$Family,'; ',samples$Genus,'; ',samples$Species,'; ',row.names(samples))

# Combine tables
otu_tax = cbind(otus, taxonomy)

# Fix names...
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ;  ;  ;  ; ", replacement = "")
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ;  ;  ; ", replacement = "")
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ;  ; ", replacement = "")
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ; ", replacement = "")
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ; ", replacement = "; ")
# Remove accession nr. identity and evalue
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = "~.*", replacement = "")
# Collapse phylogeny according selected level to ASV level
taxa_level = aggregate(otu_tax[,1:dim(otu_tax)[2]-1], by = list(otu_tax$taxonomy), FUN = sum)

# Add row names
rownames(taxa_level) = taxa_level[,1]
taxa_level = taxa_level[,2:dim(taxa_level)[2]]

# Transform to relative abundance
taxa_level = (scale(taxa_level, center=F, scale=colSums(taxa_level)))*100
colSums(taxa_level)

# Remove, Group low abundant phyla -> spurious reads were previously removed
percentage = 0
#            ^ number to filter phylum/class/order abundance

abundant.phyla = apply(taxa_level, 1, function(x){perc = sum(x>percentage)/length(x); perc>0})
rare.phyla = apply(taxa_level, 1, function(x){perc = sum(x<percentage)/length(x); perc>((100-percentage)/100)})
filtered_table = taxa_level[abundant.phyla,]
filtered_rare_table = taxa_level[!abundant.phyla,]

colSums(filtered_table)
colSums(filtered_rare_table)

# Transpose configured OTU table for next steps
samples = data.frame(t(filtered_table))

# Define sample groups from metadata for indicspecies test
groups = metadata$individual 
#                   ^ compare what?

# Associations of ASVs to individuals
# Otherwise use abundance data to get 'point biserial correlation coefficient'
# r is Pearson's phi coefficient of association (can be negative)
phi = multipatt(samples, groups, control = how(nperm=999), func = "r.g")

# here without treatment combinations (used later for association strength to target)
phiori = multipatt(samples, groups, duleg = TRUE, control = how(nperm=999), func = "r.g")

# Show association strength (just the first 6 lines and rounded to 2 decimals)
# With combinations
round(head(phi$str),2)

# Without combinations
round(head(phiori$str),2)

# List taxa with associations (with and without)
summary(phi)

# For all species (default alpha = 0.05 = only significant indicator species are reported, alpha = 1 report all)
summary(phi, alpha = 1)

# Display all species
# NA = appear in all groups/samples
phi$sign

# Save as data.frame
phi_significance = data.frame(phi$sign, phiori$str, check.names = F)

# Correct column names
setnames(phi_significance, old = colnames(phi_significance[(dim(phi_significance)[2]-dim(data.frame(unique(groups)))[1]+1):dim(phi_significance)[2]]), new = c(paste0("X",colnames(phi_significance[(dim(phi_significance)[2]-dim(data.frame(unique(groups)))[1]+1):dim(phi_significance)[2]]))))
colnames(phi_significance) = gsub("s[.]", "", colnames(phi_significance))

# Only include significant indicator species (p < 0.05)
phi_significance_p_0.05 = phi_significance[which(phi_significance$p.value <= 0.05),]

# Export data frame with ASVs with significant associations and their association strength to each individual
# These results are reported in supplementary table S15
write.table(phi_significance_p_0.05, file = "phi_significance_p_0.05_a.tsv", sep="\t", dec=".",col.names=NA,row.names=T,  quote = F)

#Cytoscape ####
# Transform the indicspecies output for Cytoscape
# Create Cytoscape compatible network table from indicspecies table to create network shown in figure 5a.
Cytoscape = NULL

for(c in 1:length(unique(groups)))
{
  for(r in 1:dim(phi_significance_p_0.05)[1])
  {
    if(phi_significance_p_0.05[r,c] == 0) next
    Cytoscape = rbind(Cytoscape, c(colnames(phi_significance_p_0.05)[c],
                                   row.names(phi_significance_p_0.05)[r],
                                   phi_significance_p_0.05$stat[r],
                                   phi_significance_p_0.05[r, (paste0("X",(colnames(phi_significance_p_0.05)[c])))]))
  }
}

# Convert to data frame
Cytoscape = as.data.frame(Cytoscape)

# Clean names
Cytoscape$V2 = gsub("[.][.]", "; ", Cytoscape$V2)
Cytoscape$V2 = gsub("[.]", " ", Cytoscape$V2)

# Rename columns (Weight=stat from the association table before) and the association strength of groups
names(Cytoscape) = c("Source", "Target", "Weight_(stat)", "Association_strength_to_target")
Cytoscape$Target = gsub(" $","", Cytoscape$Target)
Cytoscape$Target = gsub(";$","", Cytoscape$Target)

# Export table for Cytoscape
write.table(Cytoscape, file = "associated_taxa_cytoscape_a.tsv", append = F, quote = F, sep = "\t", row.names = F, col.names = T)

# Add relative abundance data to Cytoscape networks
# Collapse/Aggregate samples by metadata column (also used for sorting samples!!!)
tax_info = data.matrix(row.names(samples))
colSums(t(samples))

aggdata = aggregate(samples, by=list(groups), FUN=mean)

# Some renaming and removals
rownames(aggdata) = aggdata$Group.1
aggdata$Group.1 = NULL

# "Restore" OTU table
agg_samples = as.data.frame(t(aggdata))
colSums(agg_samples)

# Add averages as additional column
average_abundance = rowMeans(agg_samples)
agg_samples = cbind(agg_samples, average_abundance)

# Clean names
row.names(agg_samples) = gsub("[.][.]", "; ", row.names(agg_samples))
row.names(agg_samples) = gsub("[.]", " ", row.names(agg_samples))
row.names(agg_samples) = gsub(" $","", row.names(agg_samples))
row.names(agg_samples) = gsub(";$","", row.names(agg_samples))

# Export table with abundances for Cytoscape
write.table(agg_samples, file="associated_taxa_cytoscape_abundance_a.tsv", col.names=NA, row.names=T, sep = "\t", quote=FALSE)

#Estimate number of shared significant ASVs between individuals for Mantel tests
phi_sig = phi_significance_p_0.05
phi_sig = rownames_to_column(phi_sig, "row_names")
colnames(phi_sig)[colnames(phi_sig) == "row_names"] = "taxon"

phi_sig = phi_sig[,c("Amorgos", "Isabella", "Kea", "Lefkada", "Luzon", "Paros", "Thassos", "Tilos")]

#Amorgos -> AAmoM
amor = phi_sig %>% filter(Amorgos %in% c("1")) 
xx = reshape2::melt(amor, id = c("Amorgos"))
names(xx)[names(xx) == "variable"] = "individual"
names(xx)[names(xx) == "value"] = "significant"
xx = xx %>% filter(significant %in% c("1"))
amor_counts = rename(count(xx, individual, significant), Freq = n)
colnames(amor_counts)[colnames(amor_counts) == "individual"] = "Amorgos"
write.table(amor_counts, file="Amorgos_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Isabella -> AIsaF
isab = phi_sig %>% filter(Isabella %in% c("1")) 
isab = reshape2::melt(isab, id = c("Isabella"))
names(isab)[names(isab) == "variable"] = "individual"
names(isab)[names(isab) == "value"] = "significant"
isab = isab %>% filter(significant %in% c("1"))
isab_counts = rename(count(isab, individual, significant), Freq = n)
colnames(isab_counts)[colnames(isab_counts) == "individual"] = "Isabella"
write.table(isab_counts, file="Isabella_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Kea -> AKeaF
kea = phi_sig %>% filter(Kea %in% c("1")) 
kea = reshape2::melt(kea, id = c("Kea"))
names(kea)[names(kea) == "variable"] = "individual"
names(kea)[names(kea) == "value"] = "significant"
kea = kea %>% filter(significant %in% c("1"))
kea_counts = rename(count(kea, individual, significant), Freq = n)
colnames(kea_counts)[colnames(kea_counts) == "individual"] = "Kea"
write.table(kea_counts, file="Kea_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Lefkada -> ALefF
lefka = phi_sig %>% filter(Lefkada %in% c("1")) 
lefka = reshape2::melt(lefka, id = c("Lefkada"))
names(lefka)[names(lefka) == "variable"] = "individual"
names(lefka)[names(lefka) == "value"] = "significant"
lefka = lefka %>% filter(significant %in% c("1"))
lefka_counts = rename(count(lefka, individual, significant), Freq = n)
colnames(lefka_counts)[colnames(lefka_counts) == "individual"] = "Lefkada"
write.table(lefka_counts, file="Lefkada_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Luzon -> ALuzM
luz = phi_sig %>% filter(Luzon %in% c("1")) 
luz = reshape2::melt(luz, id = c("Luzon"))
names(luz)[names(luz) == "variable"] = "individual"
names(luz)[names(luz) == "value"] = "significant"
luz = luz %>% filter(significant %in% c("1"))
luz_counts = rename(count(luz, individual, significant), Freq = n)
colnames(luz_counts)[colnames(luz_counts) == "individual"] = "Luzon"
write.table(luz_counts, file="Luzon_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Paros -> AParM
paros = phi_sig %>% filter(Paros %in% c("1")) 
paros = reshape2::melt(paros, id = c("Paros"))
names(paros)[names(paros) == "variable"] = "individual"
names(paros)[names(paros) == "value"] = "significant"
paros = paros %>% filter(significant %in% c("1"))
paros_counts = rename(count(paros, individual, significant), Freq = n)
colnames(paros_counts)[colnames(paros_counts) == "individual"] = "Paros"
write.table(paros_counts, file="Paros_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Thassos -> AThaM
thass = phi_sig %>% filter(Thassos %in% c("1")) 
thass = reshape2::melt(thass, id = c("Thassos"))
names(thass)[names(thass) == "variable"] = "individual"
names(thass)[names(thass) == "value"] = "significant"
thass = thass %>% filter(significant %in% c("1"))
thass_counts = rename(count(thass, individual, significant), Freq = n)
colnames(thass_counts)[colnames(thass_counts) == "individual"] = "Thassos"
write.table(thass_counts, file="Thassos_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Tilos -> ATilM
tilo = phi_sig %>% filter(Tilos %in% c("1")) 
tilo = reshape2::melt(tilo, id = c("Tilos"))
names(tilo)[names(tilo) == "variable"] = "individual"
names(tilo)[names(tilo) == "value"] = "significant"
tilo = tilo %>% filter(significant %in% c("1"))
tilo_counts = rename(count(tilo, individual, significant), Freq = n)
colnames(tilo_counts)[colnames(tilo_counts) == "individual"] = "Tilos"
write.table(tilo_counts, file="Tilos_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Group B ####
#Subset samples in ampvis2
dataset_b <- amp_subset_samples(dataset, group %in% c("B"))
#After: 216 samples and 1021 OTUs

# Export data
amp_export_otutable(dataset_b, "table_network_b", sep = "\t", extension = "tsv")
write.table(dataset_b$metadata, file = "metadata_network_b.tsv", sep = "\t")

# Load data
samples = read.table("table_network_b.tsv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata_network_b.tsv", row.names = 1, header = T)

# Remove taxonomy
otus = samples [,1:(dim(samples)[2]-7)]
# Collapse phylogeny to ASV
taxonomy = paste0(samples$Phylum,'; ',samples$Class,'; ',samples$Order,'; ',samples$Family,'; ',samples$Genus,'; ',samples$Species,'; ',row.names(samples))

# Combine tables
otu_tax = cbind(otus, taxonomy)
# Fix names...
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ;  ;  ;  ; ", replacement = "")
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ;  ;  ; ", replacement = "")
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ;  ; ", replacement = "")
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ; ", replacement = "")
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ; ", replacement = "; ")
# Remove accession nr. identity and evalue
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = "~.*", replacement = "")

# Collapse phylogeny according selected level as defined above
taxa_level = aggregate(otu_tax[,1:dim(otu_tax)[2]-1], by = list(otu_tax$taxonomy), FUN = sum)

# Add row names
rownames(taxa_level) = taxa_level[,1]
taxa_level = taxa_level[,2:dim(taxa_level)[2]]

# Transform to relative abundance
taxa_level = (scale(taxa_level, center=F, scale=colSums(taxa_level)))*100
colSums(taxa_level)

# Remove, Group low abundant phyla set to zero as spurious reads were removed previously.
percentage = 0
#            ^ number to filter phylum/class/order abundance 

abundant.phyla = apply(taxa_level, 1, function(x){perc = sum(x>percentage)/length(x); perc>0})
rare.phyla = apply(taxa_level, 1, function(x){perc = sum(x<percentage)/length(x); perc>((100-percentage)/100)})
filtered_table = taxa_level[abundant.phyla,]
filtered_rare_table = taxa_level[!abundant.phyla,]

colSums(filtered_table)
colSums(filtered_rare_table)

# Transpose configured OTU table for next steps
samples = data.frame(t(filtered_table))

# Define sample groups from metadata for indicspecies test
groups = metadata$individual #individual color by group
#                   ^ compare what?

# Associations of ASVs to individuals
# Use abundance data to get 'point biserial correlation coefficient'
# r is Pearson's phi coefficient of association 
phi = multipatt(samples, groups, control = how(nperm=999), func = "r.g")

# here without treatment combinations (used later for association strength to target)
phiori = multipatt(samples, groups, duleg = TRUE, control = how(nperm=999), func = "r.g")

# Show association strength (just the first 6 lines and rounded to 2 decimals)
# With combinations
round(head(phi$str),2)

# Without combinations
round(head(phiori$str),2)

# List taxa with associations (with and without)
summary(phi)

# For all species (default alpha = 0.05 = only significant indicator species are reported, aplha = 1 report all)
summary(phi, alpha = 1)

# Display all species
# NA = appear in all groups/samples
phi$sign

# Save as data.frame
phi_significance = data.frame(phi$sign, phiori$str, check.names = F)

# Correct column names 
setnames(phi_significance, old = colnames(phi_significance[(dim(phi_significance)[2]-dim(data.frame(unique(groups)))[1]+1):dim(phi_significance)[2]]), new = c(paste0("X",colnames(phi_significance[(dim(phi_significance)[2]-dim(data.frame(unique(groups)))[1]+1):dim(phi_significance)[2]]))))
colnames(phi_significance) = gsub("s[.]", "", colnames(phi_significance))

# Only include significant indicator species (p < 0.05)
phi_significance_p_0.05 = phi_significance[which(phi_significance$p.value <= 0.05),]

# Export data frame with ASVs with significant associations and their association strength to each individual
# These results are reported in supplementary table S17
write.table(phi_significance_p_0.05, file = "phi_significance_p_0.05_b.tsv", sep="\t", dec=".",col.names=NA,row.names=T,  quote = F)

#Cytoscape ####
# Transform the indicspecies output for Cytoscape
# Create Cytoscape compatible network table from indicspecies table to create network shown in figure 5c.
Cytoscape = NULL

for(c in 1:length(unique(groups)))
{
  for(r in 1:dim(phi_significance_p_0.05)[1])
  {
    if(phi_significance_p_0.05[r,c] == 0) next
    Cytoscape = rbind(Cytoscape, c(colnames(phi_significance_p_0.05)[c],
                                   row.names(phi_significance_p_0.05)[r],
                                   phi_significance_p_0.05$stat[r],
                                   phi_significance_p_0.05[r, (paste0("X",(colnames(phi_significance_p_0.05)[c])))]))
  }
}

# Convert to data frame
Cytoscape = as.data.frame(Cytoscape)

# Clean names
Cytoscape$V2 = gsub("[.][.]", "; ", Cytoscape$V2)
Cytoscape$V2 = gsub("[.]", " ", Cytoscape$V2)

# Rename columns (Weight=stat from the association table before) and the association strength of groups
names(Cytoscape) = c("Source", "Target", "Weight_(stat)", "Association_strength_to_target")
Cytoscape$Target = gsub(" $","", Cytoscape$Target)
Cytoscape$Target = gsub(";$","", Cytoscape$Target)

# Export table for Cytoscape
write.table(Cytoscape, file = "associated_taxa_cytoscape_b.tsv", append = F, quote = F, sep = "\t", row.names = F, col.names = T)

# Add relative abundance data to Cytoscape networks
# Collapse/Aggregate samples by metadata column (also used for sorting samples)
tax_info = data.matrix(row.names(samples))
colSums(t(samples))

aggdata = aggregate(samples, by=list(groups), FUN=mean)

# Some nasty renaming and removals
rownames(aggdata) = aggdata$Group.1
aggdata$Group.1 = NULL

# "Restore" OTU table
agg_samples = as.data.frame(t(aggdata))
colSums(agg_samples)

# Add averages as additional column
average_abundance = rowMeans(agg_samples)
agg_samples = cbind(agg_samples, average_abundance)

# Clean names
row.names(agg_samples) = gsub("[.][.]", "; ", row.names(agg_samples))
row.names(agg_samples) = gsub("[.]", " ", row.names(agg_samples))
row.names(agg_samples) = gsub(" $","", row.names(agg_samples))
row.names(agg_samples) = gsub(";$","", row.names(agg_samples))

# Export new table
write.table(agg_samples, file="associated_taxa_cytoscape_abundance_b.tsv", col.names=NA, row.names=T, sep = "\t", quote=FALSE)

#Estimate number of shared significant ASVs between individuals for Mantel tests
phi_sig = phi_significance_p_0.05
phi_sig = rownames_to_column(phi_sig, "row_names")
colnames(phi_sig)[colnames(phi_sig) == "row_names"] = "taxon"
phi_sig = phi_sig[,c("Adonara", "Aloha", "Bangladesh", "Bora", "Buru", "Jaco", "Latalata", "Oman", "Rinca",
                     "Tilos")]

#Adonara -> BAdoF
adona = phi_sig %>% filter(Adonara %in% c("1")) 
xx = reshape2::melt(adona, id = c("Adonara"))
names(xx)[names(xx) == "variable"] = "individual"
names(xx)[names(xx) == "value"] = "significant"
xx = xx %>% filter(significant %in% c("1"))
adona_counts = rename(count(xx, individual, significant), Freq = n)
colnames(adona_counts)[colnames(adona_counts) == "individual"] = "Adonara"
write.table(adona_counts, file="Adonara_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Aloha -> BAloF
aloha = phi_sig %>% filter(Aloha %in% c("1")) 
aloha = reshape2::melt(aloha, id = c("Aloha"))
names(aloha)[names(aloha) == "variable"] = "individual"
names(aloha)[names(aloha) == "value"] = "significant"
aloha = aloha %>% filter(significant %in% c("1"))
aloha_counts = rename(count(aloha, individual, significant), Freq = n)
colnames(aloha_counts)[colnames(aloha_counts) == "individual"] = "Aloha"
write.table(aloha_counts, file="Aloha_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Bangladesh -> BBanM
bangla = phi_sig %>% filter(Bangladesh %in% c("1")) 
bangla = reshape2::melt(bangla, id = c("Bangladesh"))
names(bangla)[names(bangla) == "variable"] = "individual"
names(bangla)[names(bangla) == "value"] = "significant"
bangla = bangla %>% filter(significant %in% c("1"))
bangla_counts = rename(count(bangla, individual, significant), Freq = n)
colnames(bangla_counts)[colnames(bangla_counts) == "individual"] = "Bangladesh"
write.table(bangla_counts, file="Bangladesh_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Bora -> BBorF
bora = phi_sig %>% filter(Bora %in% c("1")) 
bora = reshape2::melt(bora, id = c("Bora"))
names(bora)[names(bora) == "variable"] = "individual"
names(bora)[names(bora) == "value"] = "significant"
bora = bora %>% filter(significant %in% c("1"))
bora_counts = rename(count(bora, individual, significant), Freq = n)
colnames(bora_counts)[colnames(bora_counts) == "individual"] = "Bora"
write.table(bora_counts, file="Bora_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Buru -> BBurM
buru = phi_sig %>% filter(Buru %in% c("1")) 
buru = reshape2::melt(buru, id = c("Buru"))
names(buru)[names(buru) == "variable"] = "individual"
names(buru)[names(buru) == "value"] = "significant"
buru = buru %>% filter(significant %in% c("1"))
buru_counts = rename(count(buru, individual, significant), Freq = n)
colnames(buru_counts)[colnames(buru_counts) == "individual"] = "Buru"
write.table(buru_counts, file="Buru_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Jaco -> BJacM
jaco = phi_sig %>% filter(Jaco %in% c("1")) 
jaco = reshape2::melt(jaco, id = c("Jaco"))
names(jaco)[names(jaco) == "variable"] = "individual"
names(jaco)[names(jaco) == "value"] = "significant"
jaco = jaco %>% filter(significant %in% c("1"))
jaco_counts = rename(count(jaco, individual, significant), Freq = n)
colnames(jaco_counts)[colnames(jaco_counts) == "individual"] = "Jaco"
write.table(jaco_counts, file="Jaco_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Latalata -> BLatF
lata = phi_sig %>% filter(Latalata %in% c("1")) 
lata = reshape2::melt(lata, id = c("Latalata"))
names(lata)[names(lata) == "variable"] = "individual"
names(lata)[names(lata) == "value"] = "significant"
lata = lata %>% filter(significant %in% c("1"))
lata_counts = rename(count(lata, individual, significant), Freq = n)
colnames(lata_counts)[colnames(lata_counts) == "individual"] = "Latalata"
write.table(lata_counts, file="Latalata_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Oman -> BOmaM
oman = phi_sig %>% filter(Oman %in% c("1")) 
oman = reshape2::melt(oman, id = c("Oman"))
names(oman)[names(oman) == "variable"] = "individual"
names(oman)[names(oman) == "value"] = "significant"
oman = oman %>% filter(significant %in% c("1"))
oman_counts = rename(count(oman, individual, significant), Freq = n)
colnames(oman_counts)[colnames(oman_counts) == "individual"] = "Oman"
write.table(oman_counts, file="Oman_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Rinca -> BRinF
rinca = phi_sig %>% filter(Rinca %in% c("1")) 
rinca = reshape2::melt(rinca, id = c("Rinca"))
names(rinca)[names(rinca) == "variable"] = "individual"
names(rinca)[names(rinca) == "value"] = "significant"
rinca = rinca %>% filter(significant %in% c("1"))
rinca_counts = rename(count(rinca, individual, significant), Freq = n)
colnames(rinca_counts)[colnames(rinca_counts) == "individual"] = "Rinca"
write.table(rinca_counts, file="Rinca_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Tilos -> BTilM
tilos = phi_sig %>% filter(Tilos %in% c("1")) 
tilos = reshape2::melt(tilos, id = c("Tilos"))
names(tilos)[names(tilos) == "variable"] = "individual"
names(tilos)[names(tilos) == "value"] = "significant"
tilos = tilos %>% filter(significant %in% c("1"))
tilos_counts = rename(count(tilos, individual, significant), Freq = n)
colnames(tilos_counts)[colnames(tilos_counts) == "individual"] = "Tilos"
write.table(tilos_counts, file="Tilos_b_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Group F ####
#Subset samples in ampvis2
dataset_f <- amp_subset_samples(dataset, group %in% c("F"))
#After: 172 samples and 1026 OTUs

# Export data
amp_export_otutable(dataset_f, "table_network_f", sep = "\t", extension = "tsv")
write.table(dataset_f$metadata, file = "metadata_network_f.tsv", sep = "\t")

# Load data
samples = read.table("table_network_f.tsv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata_network_f.tsv", row.names = 1, header = T)

# Remove taxonomy
otus = samples [,1:(dim(samples)[2]-7)]
# Collapse phylogeny to ASV level
taxonomy = paste0(samples$Phylum,'; ',samples$Class,'; ',samples$Order,'; ',samples$Family,'; ',samples$Genus,'; ',samples$Species,'; ',row.names(samples))

# Combine tables
otu_tax = cbind(otus, taxonomy)
# Fix names...
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ;  ;  ;  ; ", replacement = "")
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ;  ;  ; ", replacement = "")
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ;  ; ", replacement = "")
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ; ", replacement = "")
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ; ", replacement = "; ")
# Remove accession nr. identity and evalue
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = "~.*", replacement = "")
# Collapse phylogeny according selected level
# Defined above
taxa_level = aggregate(otu_tax[,1:dim(otu_tax)[2]-1], by = list(otu_tax$taxonomy), FUN = sum)

# Add row names
rownames(taxa_level) = taxa_level[,1]
taxa_level = taxa_level[,2:dim(taxa_level)[2]]

# Transform to relative abundance
taxa_level = (scale(taxa_level, center=F, scale=colSums(taxa_level)))*100
colSums(taxa_level)

# Spurious reads were removed before
# Set to 0 
percentage = 0
#            ^ number to filter phylum/class/order abundance (0.5= everything below 1% will be summarized as "rare taxa" in this case)

abundant.phyla = apply(taxa_level, 1, function(x){perc = sum(x>percentage)/length(x); perc>0})
rare.phyla = apply(taxa_level, 1, function(x){perc = sum(x<percentage)/length(x); perc>((100-percentage)/100)})
filtered_table = taxa_level[abundant.phyla,]
filtered_rare_table = taxa_level[!abundant.phyla,]

colSums(filtered_table)
colSums(filtered_rare_table)

# Transpose configured OTU table for next steps
samples = data.frame(t(filtered_table))

# Define sample groups from metadata for indicspecies test
groups = metadata$individual #individual color by group
#                   ^ compare what?

# Associations of ASVs to individuals
# Otherwise use abundance data to get 'point biserial correlation coefficient'
# r is Pearson's phi coefficient of association (can be negative)
phi = multipatt(samples, groups, control = how(nperm=999), func = "r.g")

# here without treatment combinations (used later for association strength to target)
phiori = multipatt(samples, groups, duleg = TRUE, control = how(nperm=999), func = "r.g")

# Show association strength (just the first 6 lines and rounded to 2 decimals)
# With combinations
round(head(phi$str),2)

# Without combinations
round(head(phiori$str),2)

# List taxa with associations (with and without)
summary(phi)

# For all species (default alpha = 0.05 = only significant indicator species are reported, alpha = 1 report all)
summary(phi, alpha = 1)

# Display all species
# NA = appear in all groups/samples
phi$sign

# Save as data.frame
phi_significance = data.frame(phi$sign, phiori$str, check.names = F)

# Correct column names
setnames(phi_significance, old = colnames(phi_significance[(dim(phi_significance)[2]-dim(data.frame(unique(groups)))[1]+1):dim(phi_significance)[2]]), new = c(paste0("X",colnames(phi_significance[(dim(phi_significance)[2]-dim(data.frame(unique(groups)))[1]+1):dim(phi_significance)[2]]))))
colnames(phi_significance) = gsub("s[.]", "", colnames(phi_significance))

# Only include significant indicator species (p < 0.05)
phi_significance_p_0.05 = phi_significance[which(phi_significance$p.value <= 0.05),]

# Export data frame with ASVs with significant associations and their association strength to each individual
# These results are reported in supplementary table S19
write.table(phi_significance_p_0.05, file = "phi_significance_p_0.05_f.tsv", sep="\t", dec=".",col.names=NA,row.names=T,  quote = F)

#Cytoscape####
# Transform the indicspecies output for Cytoscape
# Create Cytoscape compatible network table from indicspecies table to create network shown in figure 5e.
Cytoscape = NULL

for(c in 1:length(unique(groups)))
{
  for(r in 1:dim(phi_significance_p_0.05)[1])
  {
    if(phi_significance_p_0.05[r,c] == 0) next
    Cytoscape = rbind(Cytoscape, c(colnames(phi_significance_p_0.05)[c],
                                   row.names(phi_significance_p_0.05)[r],
                                   phi_significance_p_0.05$stat[r],
                                   phi_significance_p_0.05[r, (paste0("X",(colnames(phi_significance_p_0.05)[c])))]))
  }
}

# Convert to data frame
Cytoscape = as.data.frame(Cytoscape)

# Clean names
Cytoscape$V2 = gsub("[.][.]", "; ", Cytoscape$V2)
Cytoscape$V2 = gsub("[.]", " ", Cytoscape$V2)

# Rename columns (Weight=stat from the association table before) and the association strength of groups
names(Cytoscape) = c("Source", "Target", "Weight_(stat)", "Association_strength_to_target")
Cytoscape$Target = gsub(" $","", Cytoscape$Target)
Cytoscape$Target = gsub(";$","", Cytoscape$Target)

# Export table for Cytoscape
write.table(Cytoscape, file = "associated_taxa_cytoscape_f.tsv", append = F, quote = F, sep = "\t", row.names = F, col.names = T)

# Add relative abundance data to Cytoscape networks
# Collapse/Aggregate samples by metadata column (also used for sorting samples)
tax_info = data.matrix(row.names(samples))
colSums(t(samples))

aggdata = aggregate(samples, by=list(groups), FUN=mean)

# Some nasty renaming and removals
rownames(aggdata) = aggdata$Group.1
aggdata$Group.1 = NULL

# "Restore" OTU table
agg_samples = as.data.frame(t(aggdata))
colSums(agg_samples)

# Add averages as additional column
average_abundance = rowMeans(agg_samples)
agg_samples = cbind(agg_samples, average_abundance)

# Clean names
row.names(agg_samples) = gsub("[.][.]", "; ", row.names(agg_samples))
row.names(agg_samples) = gsub("[.]", " ", row.names(agg_samples))
row.names(agg_samples) = gsub(" $","", row.names(agg_samples))
row.names(agg_samples) = gsub(";$","", row.names(agg_samples))

# Export new table
write.table(agg_samples, file="associated_taxa_cytoscape_abundance_f.tsv", col.names=NA, row.names=T, sep = "\t", quote=FALSE)

#Estimate number of shared significant ASVs between individuals for Mantel tests
phi_sig = phi_significance_p_0.05
phi_sig = rownames_to_column(phi_sig, "row_names")
colnames(phi_sig)[colnames(phi_sig) == "row_names"] = "taxon"
phi_sig = phi_sig[,c("Bonacca", "Caicos", "Gozo", "Lucia", "Mayaguana", "Pinos", "Tortuga")]

#Bonacca -> FBonF
bona = phi_sig %>% filter(Bonacca %in% c("1")) 
xx = reshape2::melt(bona, id = c("Bonacca"))
names(xx)[names(xx) == "variable"] = "individual"
names(xx)[names(xx) == "value"] = "significant"
xx = xx %>% filter(significant %in% c("1"))
bona_counts = rename(count(xx, individual, significant), Freq = n)
colnames(bona_counts)[colnames(bona_counts) == "individual"] = "Bonacca"
write.table(bona_counts, file="Bonacca_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Caicos -> FCaiM
Caicos = phi_sig %>% filter(Caicos %in% c("1"))
Caicos = reshape2::melt(Caicos, id = c("Caicos"))
names(Caicos)[names(Caicos) == "variable"] = "individual"
names(Caicos)[names(Caicos) == "value"] = "significant"
Caicos = Caicos %>% filter(significant %in% c("1"))
Caicos_counts = rename(count(Caicos, individual, significant), Freq = n)
colnames(Caicos_counts)[colnames(Caicos_counts) == "individual"] = "Caicos"
write.table(Caicos_counts, file="Caicos_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Gozo -> FGozM
gozo = phi_sig %>% filter(Gozo %in% c("1")) 
gozo = reshape2::melt(gozo, id = c("Gozo"))
names(gozo)[names(gozo) == "variable"] = "individual"
names(gozo)[names(gozo) == "value"] = "significant"
gozo = gozo %>% filter(significant %in% c("1"))
gozo_counts = rename(count(gozo, individual, significant), Freq = n)
colnames(gozo_counts)[colnames(gozo_counts) == "individual"] = "Gozo"
write.table(gozo_counts, file="Gozo_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Lucia -> FLuF
Lucia = phi_sig %>% filter(Lucia %in% c("1")) 
Lucia = reshape2::melt(Lucia, id = c("Lucia"))
names(Lucia)[names(Lucia) == "variable"] = "individual"
names(Lucia)[names(Lucia) == "value"] = "significant"
Lucia = Lucia %>% filter(significant %in% c("1"))
Lucia_counts = rename(count(Lucia, individual, significant), Freq = n)
colnames(Lucia_counts)[colnames(Lucia_counts) == "individual"] = "Lucia"
write.table(Lucia_counts, file="Lucia_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Mayaguana -> FMayF
maya = phi_sig %>% filter(Mayaguana %in% c("1"))
maya = reshape2::melt(maya, id = c("Mayaguana"))
names(maya)[names(maya) == "variable"] = "individual"
names(maya)[names(maya) == "value"] = "significant"
maya = maya %>% filter(significant %in% c("1"))
maya_counts = rename(count(maya, individual, significant), Freq = n)
colnames(maya_counts)[colnames(maya_counts) == "individual"] = "Mayaguana"
write.table(maya_counts, file="Mayaguana_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Pinos -> FPinM     
Pinos = phi_sig %>% filter(Pinos %in% c("1")) 
Pinos = reshape2::melt(Pinos, id = c("Pinos"))
names(Pinos)[names(Pinos) == "variable"] = "individual"
names(Pinos)[names(Pinos) == "value"] = "significant"
Pinos = Pinos %>% filter(significant %in% c("1"))
Pinos_counts = rename(count(Pinos, individual, significant), Freq = n)
colnames(Pinos_counts)[colnames(Pinos_counts) == "individual"] = "Pinos"
write.table(Pinos_counts, file="Pinos_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Tortuga ->FTorF
tortu = phi_sig %>% filter(Tortuga %in% c("1")) 
tortu = reshape2::melt(tortu, id = c("Tortuga"))
names(tortu)[names(tortu) == "variable"] = "individual"
names(tortu)[names(tortu) == "value"] = "significant"
tortu = tortu %>% filter(significant %in% c("1"))
tortu_counts = rename(count(tortu, individual, significant), Freq = n)
colnames(tortu_counts)[colnames(tortu_counts) == "individual"] = "Tortuga"
write.table(tortu_counts, file="Tortuga_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Group J ####
#Subset samples in ampvis2
dataset_j <- amp_subset_samples(dataset, group %in% c("J"))
#After: 279 samples and 1020 OTUs

# Export data
amp_export_otutable(dataset_j, "table_network_j", sep = "\t", extension = "tsv")
write.table(dataset_j$metadata, file = "metadata_network_j.tsv", sep = "\t")

# Load data
samples = read.table("table_network_j.tsv", sep = "\t", row.names = 1, header = T, check.names = F, comment.char = "", skip = 0, quote = "")
metadata = read.delim ("metadata_network_j.tsv", row.names = 1, header = T)

# Remove taxonomy
otus = samples [,1:(dim(samples)[2]-7)]
# Collapse phylogeny to ASV
taxonomy = paste0(samples$Phylum,'; ',samples$Class,'; ',samples$Order,'; ',samples$Family,'; ',samples$Genus,'; ',samples$Species,'; ',row.names(samples))

# Combine tables
otu_tax = cbind(otus, taxonomy)
# Fix names...
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ;  ;  ;  ; ", replacement = "")
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ;  ;  ; ", replacement = "")
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ;  ; ", replacement = "")
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ;  ; ", replacement = "")
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = " ; ", replacement = "; ")
# Remove accession nr. identity and evalue
otu_tax$taxonomy = gsub(x = otu_tax$taxonomy, pattern = "~.*", replacement = "")
# Collapse phylogeny according selected level
# Defined above
taxa_level = aggregate(otu_tax[,1:dim(otu_tax)[2]-1], by = list(otu_tax$taxonomy), FUN = sum)

# Add row names
rownames(taxa_level) = taxa_level[,1]
taxa_level = taxa_level[,2:dim(taxa_level)[2]]

# Transform to relative abundance
taxa_level = (scale(taxa_level, center=F, scale=colSums(taxa_level)))*100
colSums(taxa_level)

# Spurious reads were removed before
# Set to 0
percentage = 0
#            ^ number to filter phylum/class/order abundance (0.5= everything below 1% will be summarized as "rare taxa" in this case)

abundant.phyla = apply(taxa_level, 1, function(x){perc = sum(x>percentage)/length(x); perc>0})
rare.phyla = apply(taxa_level, 1, function(x){perc = sum(x<percentage)/length(x); perc>((100-percentage)/100)})
filtered_table = taxa_level[abundant.phyla,]
filtered_rare_table = taxa_level[!abundant.phyla,]

colSums(filtered_table)
colSums(filtered_rare_table)

# Transpose configured OTU table for next steps
samples = data.frame(t(filtered_table))

# Define sample groups from metadata for indicspecies test
groups = metadata$individual #individual color by group
#                   ^ compare what?

# Associations of ASVs to individuals
# Otherwise use abundance data to get 'point biserial correlation coefficient'
# r is Pearson's phi coefficient of association (can be negative)
phi = multipatt(samples, groups, control = how(nperm=999), func = "r.g")

# here without treatment combinations (used later for association strength to target)
phiori = multipatt(samples, groups, duleg = TRUE, control = how(nperm=999), func = "r.g")

# Show association strength (just the first 6 lines and rounded to 2 decimals)
# With combinations
round(head(phi$str),2)

# Without combinations
round(head(phiori$str),2)

# List taxa with associations (with and without)
summary(phi)

# For all species (default alpha = 0.05 = only significant indicator species are reported, aplha = 1 report all)
summary(phi, alpha = 1)

# Display all species
# NA = appear in all groups/samples
phi$sign

# Save as data.frame
phi_significance = data.frame(phi$sign, phiori$str, check.names = F)

# Correct column names
setnames(phi_significance, old = colnames(phi_significance[(dim(phi_significance)[2]-dim(data.frame(unique(groups)))[1]+1):dim(phi_significance)[2]]), new = c(paste0("X",colnames(phi_significance[(dim(phi_significance)[2]-dim(data.frame(unique(groups)))[1]+1):dim(phi_significance)[2]]))))
colnames(phi_significance) = gsub("s[.]", "", colnames(phi_significance))

# Only include significant indicator species (p < 0.05)
phi_significance_p_0.05 = phi_significance[which(phi_significance$p.value <= 0.05),]

# Export data frame with ASVs with significant associations and their association strength to each individual
# These results are reported in supplementary table S21
write.table(phi_significance_p_0.05, file = "phi_significance_p_0.05_j.tsv", sep="\t", dec=".",col.names=NA,row.names=T,  quote = F)

# Transform the indicspecies output for Cytoscape
# Cystoscape ####
# Create Cytoscape compatible network table from indicspecies table to create network shown in figure 5g.
Cytoscape = NULL

for(c in 1:length(unique(groups)))
{
  for(r in 1:dim(phi_significance_p_0.05)[1])
  {
    if(phi_significance_p_0.05[r,c] == 0) next
    Cytoscape = rbind(Cytoscape, c(colnames(phi_significance_p_0.05)[c],
                                   row.names(phi_significance_p_0.05)[r],
                                   phi_significance_p_0.05$stat[r],
                                   phi_significance_p_0.05[r, (paste0("X",(colnames(phi_significance_p_0.05)[c])))]))
  }
}

# Convert to data frame
Cytoscape = as.data.frame(Cytoscape)

# Clean names
Cytoscape$V2 = gsub("[.][.]", "; ", Cytoscape$V2)
Cytoscape$V2 = gsub("[.]", " ", Cytoscape$V2)

# Rename columns (Weight=stat from the association table before) and the association strength of groups
names(Cytoscape) = c("Source", "Target", "Weight_(stat)", "Association_strength_to_target")
Cytoscape$Target = gsub(" $","", Cytoscape$Target)
Cytoscape$Target = gsub(";$","", Cytoscape$Target)

# Export table for Cytoscape
write.table(Cytoscape, file = "associated_taxa_cytoscape_j.tsv", append = F, quote = F, sep = "\t", row.names = F, col.names = T)

# Add relative abundance data to Cytoscape networks
# Collapse/Aggregate samples by metadata column (also used for sorting samples!!!)
tax_info = data.matrix(row.names(samples))
colSums(t(samples))

aggdata = aggregate(samples, by=list(groups), FUN=mean)

# Some nasty renaming and removals
rownames(aggdata) = aggdata$Group.1
aggdata$Group.1 = NULL

# "Restore" OTU table
agg_samples = as.data.frame(t(aggdata))
colSums(agg_samples)

# Add averages as additional column
average_abundance = rowMeans(agg_samples)
agg_samples = cbind(agg_samples, average_abundance)

# Clean names
row.names(agg_samples) = gsub("[.][.]", "; ", row.names(agg_samples))
row.names(agg_samples) = gsub("[.]", " ", row.names(agg_samples))
row.names(agg_samples) = gsub(" $","", row.names(agg_samples))
row.names(agg_samples) = gsub(";$","", row.names(agg_samples))

# Export table for Cytoscape with abundances
write.table(agg_samples, file="associated_taxa_cytoscape_abundance_j.tsv", col.names=NA, row.names=T, sep = "\t", quote=FALSE)

#Estimate number of shared significant ASVs between individuals for Mantel tests
phi_sig = phi_significance_p_0.05
phi_sig = rownames_to_column(phi_sig, "row_names")
colnames(phi_sig)[colnames(phi_sig) == "row_names"] = "taxon"
phi_sig = phi_sig[,c("Afganistan", "Armenia", "Bahrain", "Cambodia", "Colanta", "Kasachstan", "Kuwait", "Mongolei", "Pakistan",
                     "Syria", "Taji")]

#Afganistan -> JAfgM
afga = phi_sig %>% filter(Afganistan %in% c("1"))
xx = reshape2::melt(afga, id = c("Afganistan"))
names(xx)[names(xx) == "variable"] = "individual"
names(xx)[names(xx) == "value"] = "significant"
xx = xx %>% filter(significant %in% c("1"))
afga_counts = rename(count(xx, individual, significant), Freq = n)
colnames(afga_counts)[colnames(afga_counts) == "individual"] = "afganistan"
write.table(afga_counts, file="afganistan_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Armenia -> JArmF
arme = phi_sig %>% filter(Armenia %in% c("1")) 
arme = reshape2::melt(arme, id = c("Armenia"))
names(arme)[names(arme) == "variable"] = "individual"
names(arme)[names(arme) == "value"] = "significant"
arme = arme %>% filter(significant %in% c("1"))
arme_counts = rename(count(arme, individual, significant), Freq = n)
colnames(arme_counts)[colnames(arme_counts) == "individual"] = "armenia"
write.table(arme_counts, file="armenia_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Bahrain -> JBahM
bah = phi_sig %>% filter(Bahrain %in% c("1")) 
bah = reshape2::melt(bah, id = c("Bahrain"))
names(bah)[names(bah) == "variable"] = "individual"
names(bah)[names(bah) == "value"] = "significant"
bah = bah %>% filter(significant %in% c("1"))
bah_counts = rename(count(bah, individual, significant), Freq = n)
colnames(bah_counts)[colnames(bah_counts) == "individual"] = "bahrain"
write.table(bah_counts, file="bahrain_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Cambodia -> JCamF
cam = phi_sig %>% filter(Cambodia %in% c("1")) 
cam = reshape2::melt(cam, id = c("Cambodia"))
names(cam)[names(cam) == "variable"] = "individual"
names(cam)[names(cam) == "value"] = "significant"
cam = cam %>% filter(significant %in% c("1"))
cam_counts = rename(count(cam, individual, significant), Freq = n)
colnames(cam_counts)[colnames(cam_counts) == "individual"] = "Cambodia"
write.table(cam_counts, file="Cambodia_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Colanta -> JColF
col = phi_sig %>% filter(Colanta %in% c("1"))
col = reshape2::melt(col, id = c("Colanta"))
names(col)[names(col) == "variable"] = "individual"
names(col)[names(col) == "value"] = "significant"
col = col %>% filter(significant %in% c("1"))
col_counts = rename(count(col, individual, significant), Freq = n)
colnames(col_counts)[colnames(col_counts) == "individual"] = "Colanta"
write.table(col_counts, file="Colanta_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Kasachstan -> JKasM
kasa = phi_sig %>% filter(Kasachstan %in% c("1")) 
kasa = reshape2::melt(kasa, id = c("Kasachstan"))
names(kasa)[names(kasa) == "variable"] = "individual"
names(kasa)[names(kasa) == "value"] = "significant"
kasa = kasa %>% filter(significant %in% c("1"))
kasa_counts = rename(count(kasa, individual, significant), Freq = n)
colnames(kasa_counts)[colnames(kasa_counts) == "individual"] = "Kasachstan"
write.table(kasa_counts, file="Kasachstan_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Kuwait -> JKuwM
kuwa = phi_sig %>% filter(Kuwait %in% c("1")) 
kuwa = reshape2::melt(kuwa, id = c("Kuwait"))
names(kuwa)[names(kuwa) == "variable"] = "individual"
names(kuwa)[names(kuwa) == "value"] = "significant"
kuwa = kuwa %>% filter(significant %in% c("1"))
kuwa_counts = rename(count(kuwa, individual, significant), Freq = n)
colnames(kuwa_counts)[colnames(kuwa_counts) == "individual"] = "Kuwait"
write.table(kuwa_counts, file="Kuwait_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Mongolei ->JMonM
mongo = phi_sig %>% filter(Mongolei %in% c("1")) 
mongo = reshape2::melt(mongo, id = c("Mongolei"))
names(mongo)[names(mongo) == "variable"] = "individual"
names(mongo)[names(mongo) == "value"] = "significant"
mongo = mongo %>% filter(significant %in% c("1"))
mongo_counts = rename(count(mongo, individual, significant), Freq = n)
colnames(mongo_counts)[colnames(mongo_counts) == "individual"] = "Mongolei"
write.table(mongo_counts, file="Mongolei_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Pakistan -> JPakM
paki = phi_sig %>% filter(Pakistan %in% c("1")) 
paki = reshape2::melt(paki, id = c("Pakistan"))
names(paki)[names(paki) == "variable"] = "individual"
names(paki)[names(paki) == "value"] = "significant"
paki = paki %>% filter(significant %in% c("1"))
paki_counts = rename(count(paki, individual, significant), Freq = n)
colnames(paki_counts)[colnames(paki_counts) == "individual"] = "Pakistan"
write.table(paki_counts, file="Pakistan_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Syria -> JSyrF
syri = phi_sig %>% filter(Syria %in% c("1")) 
syri = reshape2::melt(syri, id = c("Syria"))
names(syri)[names(syri) == "variable"] = "individual"
names(syri)[names(syri) == "value"] = "significant"
syri = syri %>% filter(significant %in% c("1"))
syri_counts = rename(count(syri, individual, significant), Freq = n)
colnames(syri_counts)[colnames(syri_counts) == "individual"] = "Syria"
write.table(syri_counts, file="Syria_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)

#Taji -> JTajM
taji = phi_sig %>% filter(Taji %in% c("1")) 
taji = reshape2::melt(taji, id = c("Taji"))
names(taji)[names(taji) == "variable"] = "individual"
names(taji)[names(taji) == "value"] = "significant"
taji = taji %>% filter(significant %in% c("1"))
taji_counts = rename(count(taji, individual, significant), Freq = n)
colnames(taji_counts)[colnames(taji_counts) == "individual"] = "Taji"
write.table(taji_counts, file="Taji_sharedASVcounts.csv", col.names=NA, row.names=T, sep = ",", quote=FALSE)
