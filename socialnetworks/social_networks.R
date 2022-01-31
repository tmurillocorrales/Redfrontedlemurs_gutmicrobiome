setwd("~/PhD Thesis/Manuscript_Socialtransmission/Submission_microbiome/scripts/test/Socialnetworks")


######### Calculation of social networks #####
#These scripts were developed with the support of Lea Prox.

#Load libraries ####
library(stringr)
library(vegan)
library(tidyverse)
library(igraph)
library(ape)
library(phangorn)
library(dplyr)

#load behavioral data ####
xdata <- read.table(file="focals.duration.txt", sep="\t", header=T)
focal.length <- read.table(file="length_focal.txt", sep="\t", header=T)

#Convert to character
xdata$individual <- as.character(xdata$individual)
#Convert to numeric
xdata$duration2 <- as.numeric(xdata$duration2)
focal.length$length_focalsec = as.numeric(focal.length$length_focalsec)

#From now on prepare networks for each group separately

### Group A ####
metadata = read.table(file="groupa_metadata.txt", sep="\t", header=T)

#Unite partner in one column
xdata = unite(xdata, "With", c("partner1", "partner2", "partner3", "partner4", "partner5"), sep = "; ", remove = F)
xdata$With <- as.character(xdata$With)
str(xdata)

#Subset columns of interest to bind to xdata
focal.length = focal.length[, c("code", "length_focalsec")]

#Bind to xdata
xdata = merge(xdata, focal.length, by = "code")

#Subset group A
xdata.a = xdata %>%  filter(group %in% c("A"))

#Remove stops because it has NAs
xdata.a = xdata.a%>% filter(!behavior %in% c("stop"))

#Calculate observation time ####
#create matrix to store total observation time
names <- unique(xdata.a$individual) #keep same order of names throughout network analysis
totaltime <- as.vector(NA) #create empty vector for total time 
focaltime <- as.data.frame(cbind(names,totaltime)) #bind names and empty vector together
focaltime$totaltime <- as.numeric(focaltime$totaltime)

#loop trhough all individuals and calculate observation time
for(i in 1:length(names)){
  focaltime$totaltime[focaltime$names==names[i]] <- sum(xdata.a$length_focalsec[xdata.a$individual==names[i]])
}

#Calculate grooming aspect ####
#create matrix to store groomtime
groomtime=matrix(0,length(names),length(names))
groomtime=as.data.frame(groomtime)
rownames(groomtime)=names
colnames(groomtime)=names

#Include all grooming behaviors
groom.behaviour <- c("g","mg","bg")
xdata.grooming <- xdata.a[xdata.a$behavior %in% groom.behaviour,]

#loop through all focals and loop through all partners to calculate times spent in grooming
for(i in 1:length(names)){
  for(j in 1:length(names)){
    xx <- xdata.grooming[xdata.grooming$individual==names[i],]
    xx <- xx[str_detect(xx$With,names[j]),]
    groomtime[i,j] <-  groomtime[i,j]+sum(xx$duration2,na.rm=T)
    groomtime[j,i] <-  groomtime[j,i]+sum(xx$duration2,na.rm=T)#create undirected network
  }
}


#calculate percentage of time spent in grooming
for(i in 1:length(names)){
  for(j in 1:length(names)){
    groomtime[i,j] <-  groomtime[i,j]*2/(focaltime$totaltime[i]+focaltime$totaltime[j])
  }
}

diag(groomtime) =NA #make diagonal NA

#calculate median grooming of group
groomtime <- as.matrix(groomtime)
median.groom <- median(groomtime,na.rm=T)
#calculate the mean
mean.groom <- mean(groomtime,na.rm=T)
groom.csi <- groomtime/mean.groom

groomtime.scaled = groomtime/max(groomtime,na.rm=T)

#Calculate body aspect ####
#create matrix to store time in body contact
bctime=matrix(0,length(names),length(names))
bctime=as.data.frame(bctime)
rownames(bctime)=names
colnames(bctime)=names

#Include body contacts and huddling
bc.behaviour <- c("bc","hu", "bbc", "bhu")
xdata.ac <- xdata.a[xdata.a$behavior %in% bc.behaviour,]

#loop trhough all focals and loop through all partners to calculate times spent in body contact
for(i in 1:length(names)){
  for(j in 1:length(names)){
    xx <- xdata.ac[xdata.grooming$individual==names[i],]
    xx <- xx[str_detect(xx$With,names[j]),]
    bctime[i,j] <-  bctime[i,j]+sum(xx$duration2,na.rm=T)
    bctime[j,i] <-  bctime[j,i]+sum(xx$duration2,na.rm=T)#create undirected network
  }
}

#calculate percentage of time spent in bc
for(i in 1:length(names)){
  for(j in 1:length(names)){
    bctime[i,j] <-  bctime[i,j]*2/(focaltime$totaltime[i]+focaltime$totaltime[j])
  }
}

diag(bctime) =NA #make diagonal NA

#calculate median grooming of group
bctime <- as.matrix(bctime)
mean.bc <- mean(bctime,na.rm=T)
bc.csi <- bctime/mean.bc

#Calculate DSI ####
#first check if both matrices are correlated
mantel(groom.csi,bc.csi, method="spearman")

#Mantel statistic r: 0.4954 
#Significance: 0.012 

#Upper quantiles of permutations (null model):
#90%   95% 97.5%   99% 
#0.288 0.378 0.432 0.496 
#Permutation: free
#Number of permutations: 999

#matrices are correlated -> Behaviors can be sum in the index

dsi <- (groom.csi+bc.csi)/2
dsi.scaled <- dsi/max(dsi, na.rm = T)

dsi.scaled <- dsi.scaled[order(row.names(dsi.scaled)),]
dsi.scaled <- dsi.scaled[,order(colnames(dsi.scaled))]

write.table(dsi.scaled, file = "groupa.dsi.csv", sep = ",",row.names = TRUE, col.names = TRUE)

# Create network ####

library(igraph)

#format data
grooming <- as.matrix(dsi.scaled)
net<-graph.adjacency(grooming, mode="undirected", weighted=TRUE,diag=FALSE)
plot(net)

# make labels nicer
plot(net, vertex.size=20, vertex.label.family="Arial Black", vertex.label.cex=0.60, vertex.label.color="black")

#Node atributes
str(metadata)
#Be careful individuals should be in the same order as in dsi matrix (dsi.scaled)
# link metadata data to vertices and edges
V(net)$sex <- metadata$sex
V(net)$age <- metadata$age_category

# change node shape according to age class
V(net)$shape[which(V(net)$age=="juvenile")]  <- "square"
V(net)$shape[which(V(net)$age=="adult")]  <- "circle"
V(net)$shape[which(V(net)$age=="infant")]  <- "square"

plot(net)

# change node colour according to sex
V(net)$color[which(V(net)$sex=="male")]  <- "lightblue"
V(net)$color[which(V(net)$sex=="female")]  <- "orange"
plot(net)

#Edge atributes
E(net)$color <- "black"
#let's colour strong relationships in blue
quantile(grooming,na.rm=T,probs = c(0.9)) #0.6177596
###
E(net)$color[which(E(net)$weight>0.43)] <- "blue"
plot(net, vertex.size=20, vertex.label.family="Arial Black", vertex.label.cex=0.60, vertex.label.color="black")

E(net)$width <- 5*(E(net)$weight)
plot(net, vertex.size=20, vertex.label.family="Arial Black", vertex.label.cex=0.60, vertex.label.color="black")

### Compare indicator ASVs and social network ###################

indicspeca = read.table("indic_a.csv", sep = ",", header=T, check.names = F, row.names = 1, stringsAsFactors = F) 
snaa = read.table("groupa.dsi.csv", sep = ",", header=T, check.names = F, row.names = 1, stringsAsFactors = F) 

indicspeca = as.matrix(indicspeca)
snaa = as.matrix(snaa)

### Mantel correlation test ####

mantel_a = mantel(indicspeca, snaa, method="spearman", permutations=1000)
mantel_a

#Mantel statistic r: 0.5364 
#Significance: 0.001998 

#Upper quantiles of permutations (null model):
#90%   95% 97.5%   99% 
#0.255 0.326 0.385 0.435 
#Permutation: free
#Number of permutations: 1000

#Both matrices correlate

# Group B ####
#load data
xdata <- read.table(file="focals.duration.txt", sep="\t", header=T)
focal.length <- read.table(file="length_focal.txt", sep="\t", header=T)
#convert to character
xdata$individual <- as.character(xdata$individual)
#convert to numeric
xdata$duration2 <- as.numeric(xdata$duration2)
focal.length$length_focalsec = as.numeric(focal.length$length_focalsec)
#Load metadata for networks
metadata = read.table(file="groupb_metadata.txt", sep="\t", header=T)

#Unite partner in one column
xdata = unite(xdata, "With", c("partner1", "partner2", "partner3", "partner4", "partner5"), sep = "; ", remove = F)
xdata$With <- as.character(xdata$With)
str(xdata)

#Subset columns of interest to bind to xdata
focal.length = focal.length[, c("code", "length_focalsec")]

#Bind to xdata
xdata = merge(xdata, focal.length, by = "code")

#Subset group B
xdata.b = xdata %>%  filter(group %in% c("B"))

#Remove stops because it has NAs
xdata.b = xdata.b %>% filter(!behavior %in% c("stop"))

#Calculate observation time ####
#create matrix to store total observation time
names <- unique(xdata.b$individual) #keep same order of names throughout network analysis
totaltime <- as.vector(NA) #create empty vector for total time 
focaltime <- as.data.frame(cbind(names,totaltime)) #bind names and empty vector together
focaltime$totaltime <- as.numeric(focaltime$totaltime)

#loop trhough all individuals and calculate observation time
for(i in 1:length(names)){
  focaltime$totaltime[focaltime$names==names[i]] <- sum(xdata.b$length_focalsec[xdata.b$individual==names[i]])
}

#Calculate grooming aspect ####
#create matrix to store groomtime
groomtime=matrix(0,length(names),length(names))
groomtime=as.data.frame(groomtime)
rownames(groomtime)=names
colnames(groomtime)=names

groom.behaviour <- c("g","mg","bg")
xdata.grooming <- xdata.b[xdata.b$behavior %in% groom.behaviour,]

#loop trhough all focals and loop through all partners to calculate times spent in grooming
for(i in 1:length(names)){
  for(j in 1:length(names)){
    xx <- xdata.grooming[xdata.grooming$individual==names[i],]
    xx <- xx[str_detect(xx$With,names[j]),]
    groomtime[i,j] <-  groomtime[i,j]+sum(xx$duration2,na.rm=T)
    groomtime[j,i] <-  groomtime[j,i]+sum(xx$duration2,na.rm=T)#create undirected network
  }
}


#calculate percentage of time spent in grooming
for(i in 1:length(names)){
  for(j in 1:length(names)){
    groomtime[i,j] <-  groomtime[i,j]*2/(focaltime$totaltime[i]+focaltime$totaltime[j])
  }
}

diag(groomtime) =NA #make diagonal NA

#calculate median grooming of group
groomtime <- as.matrix(groomtime)
median.groom <- median(groomtime,na.rm=T)
#Calculate the mean 
mean.groom <- mean(groomtime,na.rm=T)
groom.csi <- groomtime/mean.groom

groomtime.scaled = groomtime/max(groomtime,na.rm=T)

#Calculate body contact aspect ####
#create matrix to store bodycontacttime
bctime=matrix(0,length(names),length(names))
bctime=as.data.frame(bctime)
rownames(bctime)=names
colnames(bctime)=names

bc.behaviour <- c("bc","hu", "bbc", "bhu")
xdata.bc <- xdata.b[xdata.b$behavior %in% bc.behaviour,]

#loop trhough all focals and loop through all partners to calculate times spent in body contact
for(i in 1:length(names)){
  for(j in 1:length(names)){
    xx <- xdata.bc[xdata.grooming$individual==names[i],]
    xx <- xx[str_detect(xx$With,names[j]),]
    bctime[i,j] <-  bctime[i,j]+sum(xx$duration2,na.rm=T)
    bctime[j,i] <-  bctime[j,i]+sum(xx$duration2,na.rm=T)#create undirected network
  }
}

#calculate percentage of time spent in bc
for(i in 1:length(names)){
  for(j in 1:length(names)){
    bctime[i,j] <-  bctime[i,j]*2/(focaltime$totaltime[i]+focaltime$totaltime[j])
  }
}

diag(bctime) =NA #make diagonal NA

#calculate median grooming of group
bctime <- as.matrix(bctime)
mean.bc <- mean(bctime,na.rm=T)
bc.csi <- bctime/mean.bc

#Calculate DSI ####
#It says csi but it is dsi
#first check if both matrices are correlated
mantel(groom.csi,bc.csi, method = "spearman")

#Mantel statistic r: 0.2875 
#Significance: 0.047 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.232 0.283 0.330 0.386 
#Permutation: free
#Number of permutations: 999

#Behaviors correlate so it is possible to sum them

dsi <- (groom.csi+bc.csi)/2
dsi.scaled <- dsi/max(dsi, na.rm = T)

dsi.scaled <- dsi.scaled[order(row.names(dsi.scaled)),]
dsi.scaled <- dsi.scaled[,order(colnames(dsi.scaled))]

write.table(dsi.scaled, file = "groupb.dsi.csv", sep = ",",row.names = TRUE, col.names = TRUE)

#Create network ####

#format data
grooming <- as.matrix(dsi.scaled)
net<-graph.adjacency(grooming, mode="undirected", weighted=TRUE,diag=FALSE)
plot(net)

# make labels nicer
plot(net, vertex.size=20, vertex.label.family="Arial Black", vertex.label.cex=0.60, vertex.label.color="black")

#Node atributes
#Be careful individuals should be in the same order as in dsi matrix (dsi.scaled)
# link metadata data to vertices and edges
V(net)$sex <- metadata$sex
V(net)$age <- metadata$age_category

# change node shape according to age class
V(net)$shape[which(V(net)$age=="juvenile")]  <- "square"
V(net)$shape[which(V(net)$age=="adult")]  <- "circle"

plot(net)

# change node colour according to sex
V(net)$color[which(V(net)$sex=="male")]  <- "lightblue"
V(net)$color[which(V(net)$sex=="female")]  <- "orange"
plot(net)

#Edge atributes
E(net)$color <- "black"
# strong relationships in blue
quantile(grooming,na.rm=T,probs = c(0.9)) #0.5212719
###
E(net)$color[which(E(net)$weight>0.43)] <- "blue"
plot(net, vertex.size=20, vertex.label.family="Arial Black", vertex.label.cex=0.60, vertex.label.color="black")

E(net)$width <- 5*(E(net)$weight)
plot(net, vertex.size=20, vertex.label.family="Arial Black", vertex.label.cex=0.60, vertex.label.color="black")

### Compare indicator ASVs and social network ####

#Load matrix with shared indicator ASVs
indicspecb = read.table("indic_b.csv", sep = ",", header=T, check.names = F, row.names = 1, stringsAsFactors = F)
#Matrix with DSI measures
snab = read.table("groupb.dsi.csv", sep = ",", header=T, check.names = F, row.names = 1, stringsAsFactors = F) 

#Convert to matrices
indicspecb = as.matrix(indicspecb)
snab = as.matrix(snab)

### Mantel correlation test ####

mantel_b = mantel(indicspecb, snab, method="spearman", permutations=1000)
mantel_b

#Mantel statistic r: 0.3989 
#Significance: 0.012987 

#Upper quantiles of permutations (null model):
#90%   95% 97.5%   99% 
#0.245 0.298 0.349 0.408 
#Permutation: free
#Number of permutations: 1000
 
# Group F ####

metadata <- read.table(file="groupf_metadata.txt", sep="\t", header=T)

#Unite partner in one column
xdata = unite(xdata, "With", c("partner1", "partner2", "partner3", "partner4", "partner5"), sep = "; ", remove = F)
xdata$With <- as.character(xdata$With)

#Subset columns of interest to bind to xdata
focal.length = focal.length[, c("code", "length_focalsec")]

#Bind to xdata
xdata = merge(xdata, focal.length, by = "code")

#Subset group A
xdata.f = xdata %>%  filter(group %in% c("F"))

#Remove stops because it has NAs
xdata.f = xdata.f%>% filter(!behavior %in% c("stop"))

#Calculate observation time ####
#create matrix to store total observation time
names <- unique(xdata.f$individual) #keep same order of names throughout network analysis
totaltime <- as.vector(NA) #create empty vector for total time 
focaltime <- as.data.frame(cbind(names,totaltime)) #bind names and empty vector together
focaltime$totaltime <- as.numeric(focaltime$totaltime)

#loop through all individuals and calculate observation time
for(i in 1:length(names)){
  focaltime$totaltime[focaltime$names==names[i]] <- sum(xdata.f$length_focalsec[xdata.f$individual==names[i]])
}

#Calculate grooming aspect ####
#create matrix to store groomtime
groomtime=matrix(0,length(names),length(names))
groomtime=as.data.frame(groomtime)
rownames(groomtime)=names
colnames(groomtime)=names

#Include all grooming behaviors
groom.behaviour <- c("g","mg","bg")
xdata.grooming <- xdata.f[xdata.f$behavior %in% groom.behaviour,]

#loop through all focals and loop through all partners to calculate times spent in grooming
for(i in 1:length(names)){
  for(j in 1:length(names)){
    xx <- xdata.grooming[xdata.grooming$individual==names[i],]
    xx <- xx[str_detect(xx$With,names[j]),]
    groomtime[i,j] <-  groomtime[i,j]+sum(xx$duration2,na.rm=T)
    groomtime[j,i] <-  groomtime[j,i]+sum(xx$duration2,na.rm=T)#create undirected network
  }
}


#calculate percentage of time spent in grooming
for(i in 1:length(names)){
  for(j in 1:length(names)){
    groomtime[i,j] <-  groomtime[i,j]*2/(focaltime$totaltime[i]+focaltime$totaltime[j])
  }
}

diag(groomtime) =NA #make diagonal NA

#calculate median grooming of group
groomtime <- as.matrix(groomtime)
median.groom <- median(groomtime,na.rm=T)
#Calculate the mean
mean.groom <- mean(groomtime,na.rm=T)
groom.csi <- groomtime/mean.groom

groomtime.scaled = groomtime/max(groomtime,na.rm=T)

#Calculate body contact aspect ####
#create matrix to store bodycontact time
bctime=matrix(0,length(names),length(names))
bctime=as.data.frame(bctime)
rownames(bctime)=names
colnames(bctime)=names

#Include body contacts and huddling
bc.behaviour <- c("bc","hu", "bbc", "bhu")
xdata.fc <- xdata.f[xdata.f$behavior %in% bc.behaviour,]

#loop through all focals and loop through all partners to calculate times spent in body contact
for(i in 1:length(names)){
  for(j in 1:length(names)){
    xx <- xdata.fc[xdata.grooming$individual==names[i],]
    xx <- xx[str_detect(xx$With,names[j]),]
    bctime[i,j] <-  bctime[i,j]+sum(xx$duration2,na.rm=T)
    bctime[j,i] <-  bctime[j,i]+sum(xx$duration2,na.rm=T)#create undirected network
  }
}

#calculate percentage of time spent in bc
for(i in 1:length(names)){
  for(j in 1:length(names)){
    bctime[i,j] <-  bctime[i,j]*2/(focaltime$totaltime[i]+focaltime$totaltime[j])
  }
}

diag(bctime) =NA #make diagonal NA

#calculate median grooming of group
bctime <- as.matrix(bctime)
mean.bc <- mean(bctime,na.rm=T)
bc.csi <- bctime/mean.bc

#Calculate DSI
#first check if both matrices are correlated
mantel(groom.csi,bc.csi, method = "spearman")
#Matrices do not correlate but for comparison purposes the behaviors were summed up for calculating the DSI

dsi <- (groom.csi+bc.csi)/2
dsi.scaled <- dsi/max(dsi, na.rm = T)

dsi.scaled <- dsi.scaled[order(row.names(dsi.scaled)),]
dsi.scaled <- dsi.scaled[,order(colnames(dsi.scaled))]

write.table(dsi.scaled, file = "groupf.dsi.csv", sep = ",",row.names = TRUE, col.names = TRUE)

## Create network ####

#load libraries
library(igraph)

metadata = read.table(file="groupf_metadata.txt", sep="\t", header=T)

#format data
grooming <- as.matrix(dsi.scaled)
net<-graph.adjacency(grooming, mode="undirected", weighted=TRUE,diag=FALSE)
plot(net)

# make labels nicer
plot(net, vertex.size=20, vertex.label.family="Arial Black", vertex.label.cex=0.60, vertex.label.color="black")

#Node attributes 
#Be careful individuals should be in the same order as in dsi matrix (dsi.scaled)
# link metadata data to vertices and edges
V(net)$sex <- metadata$sex
V(net)$age <- metadata$age_cat

# change node shape according to age class
V(net)$shape[which(V(net)$age=="juvenile")]  <- "square"
V(net)$shape[which(V(net)$age=="adult")]  <- "circle"

plot(net)

# change node colour according to sex
V(net)$color[which(V(net)$sex=="male")]  <- "lightblue"
V(net)$color[which(V(net)$sex=="female")]  <- "orange"
plot(net)

#Edge atributes
E(net)$color <- "black"
#let's colour strong relationships in blue
quantile(grooming,na.rm=T,probs = c(0.9)) #0.814853 
###
E(net)$color[which(E(net)$weight>0.43)] <- "blue"
plot(net, vertex.size=20, vertex.label.family="Arial Black", vertex.label.cex=0.60, vertex.label.color="black")

E(net)$width <- 5*(E(net)$weight)
plot(net, vertex.size=20, vertex.label.family="Arial Black", vertex.label.cex=0.60, vertex.label.color="black")

### Compare indicator ASVs and social network ####

#Load matrix with shared indicator ASVs
indicspecf = read.table("indic_f.csv", sep = ",", header=T, check.names = F, row.names = 1, stringsAsFactors = F)
#Matrix with DSI measures
snaf = read.table("groupf.dsi.csv", sep = ",", header=T, check.names = F, row.names = 1, stringsAsFactors = F) 

#Convert dataframes to matrices
indicspecf = as.matrix(indicspecf)
snaf = as.matrix(snaf)

### Mantel correlation test ####
mantel_f = mantel(indicspecf, snaf, method="spearman", permutations=1000)

mantel_f

#Mantel statistic r: 0.5033 
#Significance: 0.083916 

#Upper quantiles of permutations (null model):
#90%   95% 97.5%   99% 
#0.479 0.579 0.628 0.708 
#Permutation: free
#Number of permutations: 5039

#Matrices do not correlate.

# Group J #####

#Read metadata
metadata = read.table(file="groupj_metadata.txt", sep="\t", header=T)

#Unite partner in one column
xdata = unite(xdata, "With", c("partner1", "partner2", "partner3", "partner4", "partner5"), sep = "; ", remove = F)
xdata$With <- as.character(xdata$With)
str(xdata)

#Subset columns of interest to bind to xdata
focal.length = focal.length[, c("code", "length_focalsec")]

#Bind to xdata
xdata = merge(xdata, focal.length, by = "code")

#Subset group J
xdata.j = xdata %>%  filter(group %in% c("J"))

#Remove stops because it has NAs
xdata.j = xdata.j%>% filter(!behavior %in% c("stop"))

#Calculate observation time ##
#create matrix to store total observation time
names <- unique(xdata.j$individual) #keep same order of names throughout network analysis
totaltime <- as.vector(NA) #create empty vector for total time 
focaltime <- as.data.frame(cbind(names,totaltime)) #bind names and empty vector together
focaltime$totaltime <- as.numeric(focaltime$totaltime)

#loop through all individuals and calculate observation time
for(i in 1:length(names)){
  focaltime$totaltime[focaltime$names==names[i]] <- sum(xdata.j$length_focalsec[xdata.j$individual==names[i]])
}

#Calculate grooming aspect ######
#create matrix to store groomtime
groomtime=matrix(0,length(names),length(names))
groomtime=as.data.frame(groomtime)
rownames(groomtime)=names
colnames(groomtime)=names

#Include all grooming behaviors
groom.behaviour <- c("g","mg","bg")
xdata.grooming <- xdata.j[xdata.j$behavior %in% groom.behaviour,]

#loop through all focals and loop through all partners to calculate times spent in grooming
for(i in 1:length(names)){
  for(j in 1:length(names)){
    xx <- xdata.grooming[xdata.grooming$individual==names[i],]
    xx <- xx[str_detect(xx$With,names[j]),]
    groomtime[i,j] <-  groomtime[i,j]+sum(xx$duration2,na.rm=T)
    groomtime[j,i] <-  groomtime[j,i]+sum(xx$duration2,na.rm=T)#create undirected network
  }
}


#calculate percentage of time spent in grooming
for(i in 1:length(names)){
  for(j in 1:length(names)){
    groomtime[i,j] <-  groomtime[i,j]*2/(focaltime$totaltime[i]+focaltime$totaltime[j])
  }
}

diag(groomtime) =NA #make diagonal NA

#calculate median grooming of group
groomtime <- as.matrix(groomtime)
median.groom <- median(groomtime,na.rm=T)
#calculate the mean
mean.groom <- mean(groomtime,na.rm=T)
groom.csi <- groomtime/mean.groom

groomtime.scaled = groomtime/max(groomtime,na.rm=T)

#Calculate body contact aspect #####
#create matrix to store time in body contact
bctime=matrix(0,length(names),length(names))
bctime=as.data.frame(bctime)
rownames(bctime)=names
colnames(bctime)=names

#Include all body contact and huddling
bc.behaviour <- c("bc","hu", "bbc", "bhu")
xdata.jc <- xdata.j[xdata.j$behavior %in% bc.behaviour,]

#loop through all focals and loop through all partners to calculate times spent in body contact
for(i in 1:length(names)){
  for(j in 1:length(names)){
    xx <- xdata.jc[xdata.grooming$individual==names[i],]
    xx <- xx[str_detect(xx$With,names[j]),]
    bctime[i,j] <-  bctime[i,j]+sum(xx$duration2,na.rm=T)
    bctime[j,i] <-  bctime[j,i]+sum(xx$duration2,na.rm=T)#create undirected network
  }
}

#calculate percentage of time spent in bc
for(i in 1:length(names)){
  for(j in 1:length(names)){
    bctime[i,j] <-  bctime[i,j]*2/(focaltime$totaltime[i]+focaltime$totaltime[j])
  }
}

diag(bctime) =NA #make diagonal NA

#calculate median grooming of group
bctime <- as.matrix(bctime)
mean.bc <- mean(bctime,na.rm=T)
bc.csi <- bctime/mean.bc

#Calculate DSI ####

#first check if both matrices are correlated
mantel(groom.csi,bc.csi, method = "spearman")
#matrices do  not correlate, 0.291, r= 0.07104 -> not significant however sum behaviors for keeping it standard

dsi <- (groom.csi+bc.csi)/2
dsi.scaled <- dsi/max(dsi, na.rm = T)

dsi.scaled <- dsi.scaled[order(row.names(dsi.scaled)),]
dsi.scaled <- dsi.scaled[,order(colnames(dsi.scaled))]

write.table(dsi.scaled, file = "groupj.dsi.csv", sep = ",",row.names = TRUE, col.names = TRUE)

# Create networks ####

#load libraries
library(igraph)

#format data
grooming <- as.matrix(dsi.scaled)
net<-graph.adjacency(grooming, mode="undirected", weighted=TRUE,diag=FALSE)
plot(net)

# make labels nicer
plot(net, vertex.size=20, vertex.label.family="Arial Black", vertex.label.cex=0.60, vertex.label.color="black")

#Node atributes
#Be careful individuals should be in the same order as in dsi matrix (dsi.scaled)
# link metadata data to vertices and edges
V(net)$sex <- metadata$sex
V(net)$age <- metadata$age_cat

# change node shape according to age class
V(net)$shape[which(V(net)$age=="juvenile")]  <- "square"
V(net)$shape[which(V(net)$age=="adult")]  <- "circle"

plot(net)

# change node colour according to sex
V(net)$color[which(V(net)$sex=="male")]  <- "lightblue"
V(net)$color[which(V(net)$sex=="female")]  <- "orange"
plot(net)

#Edge atributes
E(net)$color <- "black"
#let's colour strong relationships in blue
quantile(grooming,na.rm=T,probs = c(0.9)) #0.5393311 
###
E(net)$color[which(E(net)$weight>0.43)] <- "blue"
plot(net, vertex.size=20, vertex.label.family="Arial Black", vertex.label.cex=0.60, vertex.label.color="black")

E(net)$width <- 5*(E(net)$weight)
plot(net, vertex.size=20, vertex.label.family="Arial Black", vertex.label.cex=0.60, vertex.label.color="black")

#### Compare indicator ASVs and social network ####

indicspecj = read.table("indic_j.csv", sep = ",", header=T, check.names = F, row.names = 1, stringsAsFactors = F) 
snaj = read.table("groupj.dsi.csv", sep = ",", header=T, check.names = F, row.names = 1, stringsAsFactors = F) 

indicspecj = as.matrix(indicspecj)
snaj = as.matrix(snaj)

### Mantel correlation test ####

mantel_j = mantel(indicspecj, snaj, method="spearman", permutations=1000)
mantel_j

#Mantel statistic r: 0.2347 
#Significance: 0.05994 

#Upper quantiles of permutations (null model):
#90%   95% 97.5%   99% 
#0.192 0.246 0.296 0.364 
#Permutation: free
#Number of permutations: 1000


