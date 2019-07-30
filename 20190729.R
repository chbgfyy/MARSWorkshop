##### install needed packages #####

# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("vegan")
# install.packages("indicspecies")
# install.packages("tidyr")
# install.packages("tidyverse")
library(ggplot2)
library(dplyr)
library(vegan)
library(indicspecies)
library(tidyr)
library(tidyverse)

# create function for reading mothur lower triangle matrices
parseDistanceDF = function(phylip_file) {
  
  # Read the first line of the phylip file to find out how many sequences/samples it contains
  temp_connection = file(phylip_file, 'r')
  len = readLines(temp_connection, n=1)
  len = as.numeric(len)
  len = len +1
  close(temp_connection)
  
  
  phylip_data = read.table(phylip_file, fill=T, row.names=1, skip=1, col.names=1:len)
  colnames(phylip_data) <- row.names(phylip_data)
  return(phylip_data)
}

#### Read in data ####
# otu
otu <- read.table(file = "../workshop2019.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.0.03.subsample.shared", 
                  header = TRUE, stringsAsFactors = FALSE)

# we do not need "label" and "numOtus"
otu <- select(otu, -label, -numOtus)

# making list of max value for each OTU
maxab <- apply(otu, 2, max)
# make list of sample names that are smaller than cutoff
n1 <- names(which(maxab < 50))
# remove all OTU aren't in n1
otu.ab <- otu[ , -which(names(otu) %in% n1)]
# remove all OTU that aren't in n1 from taxa
taxa.ab <- taxa[-which(taxa$OTU %in% n1),]

# taxonomy
taxa <- read.table(textConnection(gsub("\\(.+?\\);", "\t", readLines("../workshop2019.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy"))), 
                   col.names=c("OTU", "Size", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), skip=1)

# alpha diversity

alpha <- read.table(file="../workshop2019.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.groups.ave-std.summary", 
                    header = T, stringsAsFactors = FALSE)


# Our alpha diversity has 26, we only need 13
### need to use filter, cann't remember exactly how
str(iris)
alpha <- filter(alpha, method == "ave")

# beta diversity

jc <- parseDistanceDF("../workshop2019.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.jest.0.03.lt.ave.dist")
bc <- parseDistanceDF("../workshop2019.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.ave.dist")
tyc <-parseDistanceDF("../workshop2019.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.thetayc.0.03.lt.ave.dist")     

# environmental

envdata <- read.table(file="../workshop2019.env.txt", header = T)

# sample group link
samp <- read.table(file = "../workshop2019.sample.txt", header = TRUE)

# join sample sanme to env data

expdata <- left_join(envdata, samp, on ="Sample" )

alpha.expdata <- left_join(alpha, expdata, on = "group")



#### Alpha Diversity ####

# alpha diversity is diversity within a single sample

# sobs == species observed = upweighting importance of rare species
ggplot(data = alpha.expdata, aes(x=Type, y = sobs)) +
  geom_boxplot() + 
  geom_jitter(width = 0.2) +
  theme_bw() +
  ggtitle("Bacterial Richness by Sample Type")

# Shannon diversity index, balance richness and evenness
ggplot(data = alpha.expdata, aes(x=Type, y = shannon)) +
  geom_boxplot() + 
  geom_jitter(width = 0.2) +
  theme_bw() +
  ggtitle("Bacterial Shannon Diversity by Sample Type")

# Simpson index (invsimpson) downweight impact of rare species
ggplot(data = alpha.expdata, aes(x=Type, y = invsimpson)) +
  geom_boxplot() + 
  geom_jitter(width = 0.2) +
  theme_bw() +
  ggtitle("Bacterial Simpson Index by Sample Type")

# Simpson by Site
ggplot(data = alpha.expdata, aes(x=Type, y = invsimpson)) +
  geom_boxplot() + 
  geom_jitter(width = 0.2) +
  theme_bw() +
  ggtitle("Bacterial Simpson Index by Sample Type and Site") + 
  facet_grid(.~Site)

# sobs by Site
ggplot(data = alpha.expdata, aes(x=Type, y = sobs)) +
  geom_boxplot() + 
  geom_jitter(width = 0.2) +
  theme_bw() +
  ggtitle("Bacterial Richness by Sample Type and Site") + 
  facet_grid(.~Site)

# Shannon diversity index by Site
ggplot(data = alpha.expdata, aes(x=Type, y = shannon)) +
  geom_boxplot() + 
  geom_jitter(width = 0.2) +
  theme_bw() +
  ggtitle("Bacterial Shannon Diversity by Sample Type and Site") + 
  facet_grid(.~Site)

# is there sig difference between sample type (sobs)
summary(aov(alpha.expdata$sobs~alpha.expdata$Type))
TukeyHSD(aov(alpha.expdata$sobs~alpha.expdata$Type))

# is there sig difference between sample type (Shannon)
summary(aov(alpha.expdata$shannon~alpha.expdata$Type))
TukeyHSD(aov(alpha.expdata$shannon~alpha.expdata$Type))

# is there sig difference between sample type (invsimpson)
summary(aov(alpha.expdata$invsimpson~alpha.expdata$Type))
TukeyHSD(aov(alpha.expdata$invsimpson~alpha.expdata$Type))

# save TukeyHSD
tuk <- TukeyHSD(aov(alpha.expdata$shannon~alpha.expdata$Type))
str(tuk)

# access the list
tuk$`alpha.expdata$Type`

# shannoneven, more even is better?

# compare alpha indices to each other
ggplot(data=alpha.expdata, aes(x=sobs, y=shannoneven, color= as.factor(Type))) +
  geom_point(size=3)

# July 25 2019 #
#### Beta diversity  ####

# diversity between samples, pairwise dissimilarity

# visualize beta diversity by ordination, 
# NMS, non-metricmultidimensional scaling


# NMS jc, Jaccard dissimilarity, only community membership, presence/absence = upweights rare species
jc.nms <- metaMDS(as.dist(jc), k=2, trymin = 50, trymax = 500, wascores = F)

jc.points <- data.frame(jc.nms$points)

x <- max(jc.points$MDS1)/1.5

y <- min(jc.points$MDS2)

jc.plot <- ggplot(data = jc.points, aes(x=MDS1, y=MDS2, label = rownames(jc)))

jc.plot + geom_point(aes(color=alpha.expdata$Type, shape=alpha.expdata$Site), size = 3) +
  # geom_text() +
  annotate("text", x, y, label = paste("stress = ", round(bc.nms$stress, digits = 3))) +
  ggtitle("Bacterial Community Jaccard")

# make nms plot for bray-curtis, incoorporates abundance, 
# tries to balance community memebership and community structure
bc.nms <- metaMDS(as.dist(bc), k=2, trymin = 50, trymax = 500, wascores = F)

bc.points <- data.frame(bc.nms$points)

x <- max(bc.points$MDS1)/1.5

y <- min(bc.points$MDS2)

bc.plot <- ggplot(data = bc.points, aes(x=MDS1, y=MDS2, label = rownames(bc)))

bc.plot + geom_point(aes(color=alpha.expdata$Type, shape=alpha.expdata$Site), size = 3) +
  annotate("text", x, y, label = paste("stress = ", round(bc.nms$stress, digits = 3))) +
  ggtitle("Bacterial Community Bray-Curtis")

# make nms plot for Theta YC, tries to punish abundanent unshared species   
tyc.nms <- metaMDS(as.dist(tyc), k=2, trymin = 50, trymax = 500, wascores = F)

tyc.points <- data.frame(tyc.nms$points)

x <- max(tyc.points$MDS1)/1.5

y <- min(tyc.points$MDS2)

tyc.plot <- ggplot(data = tyc.points, aes(x=MDS1, y=MDS2, label = rownames(tyc)))

tyc.plot + geom_point(aes(color=alpha.expdata$Type, shape=alpha.expdata$Site), size = 3) +
  annotate("text", x, y, label = paste("stress = ", round(tyc.nms$stress, digits = 3))) +
  ggtitle("Bacterial Community Theta YC")

# Hypothesis testing on beta diverstiy

# permanova (adonis in vegan)

# Jaccard
adonis(as.dist(jc)~alpha.expdata$Type, permutations = 99, rm.na=TRUE)
adonis(as.dist(jc)~alpha.expdata$Type*alpha.expdata$Site, permutations = 99, rm.na=TRUE)

# Bray-Curtis
adonis(as.dist(bc)~alpha.expdata$Type, permutations = 99, rm.na=TRUE)
adonis(as.dist(bc)~alpha.expdata$Type*alpha.expdata$Site, permutations = 99, rm.na=TRUE)

# Theta YC
adonis(as.dist(tyc)~alpha.expdata$Type, permutations = 9999, rm.na=TRUE)
adonis(as.dist(tyc)~alpha.expdata$Type*alpha.expdata$Site, permutations = 999, rm.na=TRUE)


#### indicator species analysis ####

indic <- multipatt(otu.ab[, -1], alpha.expdata$Type, control=how(nperm = 999))

summary(indic)

write.csv(file="indicator.species.cvs", indic$sign %>% 
            rownames_to_column(var = "OTU") %>%
            mutate(p.fdr = round(p.adjust(p.value,"fdr"), 3)) %>%
            right_join(taxa, by = "OTU") %>%
            filter(p.fdr < 0.1) %>%
            arrange (index)
)

### Here I have some problem with the code
