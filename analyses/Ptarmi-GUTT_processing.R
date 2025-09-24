# 
# GUTT-P Mothur Analysis
# Amanda Stromecki
# Muscarella Lab
# Sept 2025
#

# Set working directory
setwd("/Users/amandastromecki/GitHub/Ptarmi-GUTT/data")

# Set Source R Tools
source("../analyses/DiversityFunctions.R")
source("../analyses/MothurTools.R")


# Load Req'd Packages
library(vegan)
library(RColorBrewer)
library(tidyr)
library(readr)
library(Hmisc)
library(reshape2)
library(dplyr)
library(tidyverse)
library(otuSummary)# cite
library(corrplot)
library(naniar)
library(car)

# Set Std Err Functions
se <- function(x, ...) {
  sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))
}

ci <- function(x, ...) {
  1.96 * sd(x, na.rm = TRUE)
}

cv <- function(x, ...) {
  sd(x) / mean(x)
}



# Loading data files ----
metadata <- read.csv("Only_ceca_subsampling_2022.csv", row.names = 1, header = TRUE, na.strings = "")
diet_data <- read.csv("IcelandPtarmiganID_GoodPoorFoodPercent_BCI.csv", header = TRUE, row.names = 1, check.names = TRUE)
health_data <- read.csv("PtarmiganHealthData_Dec_2023.csv", header = TRUE, row.names = 1)
diet_diversity <- read.csv("Diet_diversity.csv", header = TRUE, row.names = 1, na.strings = "")

# load asv file
asv_raw <- read.table("GUTT.trim.contigs.good.unique.good.filter.good.unique.precluster.denovo.vsearch.pick.ASV.abund.asv.shared", header = T)
rownames(asv_raw) <- asv_raw$Group
rownames(asv_raw) <- gsub("22_LM_C_", "", rownames(asv_raw))
asv_raw <- asv_raw[,-c(1:3)]
asv_raw <- as.matrix(asv_raw)

# load taxonomy file
tax_raw <- read.table("GUTT.trim.contigs.good.unique.good.filter.good.unique.precluster.denovo.vsearch.pick.ASV.abund.nr_v138_1.wang.taxonomy", header = F)

## Cleaning up ----
# Site-by-Species Matrix Curation
dim(asv_raw)
summary(rowSums(asv_raw))
summary(colSums(asv_raw))

# remove rare ASVs, low coverage samples
rare <- which(colSums(asv_raw) < 4)
length(rare)

hist(log10(rowSums(asv_raw)),main="Coverage Histogram",xlab="Coverage", breaks = 20)
low_coverage <- which(rowSums(asv_raw) < 10000)
asv_trim <- asv_raw[-low_coverage,-rare]
dim(asv_trim)
# 1052 55261
hist(log10(rowSums(asv_trim)),main="Coverage Histogram",xlab="Coverage", breaks = 20)

### Replicates ----
# pull out all replicate samples for QC
# remove replicate samples (retain those with highest number of reads)
rowSums(asv_trim[grep("^0104",rownames(asv_trim)),])
# # remove 0104, 0104_1078, 0104_1098
rowSums(asv_trim[grep("^0105",rownames(asv_trim)),])
# remove 0105, 0105_1059, 0105_1099, 0105_119
rowSums(asv_trim[grep("^0106",rownames(asv_trim)),])
# remove 0106, 0106_1060, 0106_1100, 0106_1120
rowSums(asv_trim[grep("^0262",rownames(asv_trim)),])
# remove 0262, 0262_1064, 0262_1084, 0262_1104
rowSums(asv_trim[grep("^0263",rownames(asv_trim)),])
# remove 0263, 0263_1065, 0263_1105, 0263_1125
rowSums(asv_trim[grep("^0264",rownames(asv_trim)),])
# remove 0264, 0264_1066, 0264_1086, 0264_1106
rowSums(asv_trim[grep("^0370",rownames(asv_trim)),])
# remove 0370, 0370_1061, 0370_1101, 0370_1121
rowSums(asv_trim[grep("^0372",rownames(asv_trim)),])
# remove 0372, 0372_1063, 0372_1083, 0372_1123
rowSums(asv_trim[grep("^0508",rownames(asv_trim)),])
# remove 0508, 0508_1062, 0508_1102, 0508_1122
rowSums(asv_trim[grep("^0614",rownames(asv_trim)),])
# remove 0614, 0614_1087, 0614_1107, 0614_1127
rowSums(asv_trim[grep("^0615",rownames(asv_trim)),])
# remove 0615, 0615_1088, 0615_1108, 0615_1128
rowSums(asv_trim[grep("^0616",rownames(asv_trim)),])
# remove 0616, 0616_1089, 0616_1109, 0616_1129
rowSums(asv_trim[grep("^0745",rownames(asv_trim)),])
# remove 0745, 0745_1070, 0745_1090, 0745_1110
rowSums(asv_trim[grep("^0746",rownames(asv_trim)),])
# remove 0746, 0746_1071, 0746_1091, 0746_1111
rowSums(asv_trim[grep("^0749",rownames(asv_trim)),])
# remove 0749, 0749_1072, 0749_1092, 0749_1112
rowSums(asv_trim[grep("^0751",rownames(asv_trim)),])
# remove 0751, 0751_1073, 0751_1093, 0751_1113
rowSums(asv_trim[grep("^0917",rownames(asv_trim)),])
# remove 0917, 0917_1074, 0917_1094, 0917_1134
rowSums(asv_trim[grep("^0918",rownames(asv_trim)),])
# remove 0918, 0918_1075, 0918_1095, 0918_1115
rowSums(asv_trim[grep("^0920",rownames(asv_trim)),])
# remove 0920, 0920_1076, 0920_1096, 0920_1116
rowSums(asv_trim[grep("^0921",rownames(asv_trim)),])
# remove 0921, 0921_1097, 0921_1137
replicates_all <- c("0104", "0104_1078", "0104_1098", "0104_1118",
                    "0105", "0105_1059", "05105_1079", "0105_1099", "0105_1119",
                    "0106", "0106_1060", "0105_1080", "0106_1100", "0106_1120",
                    "0262", "0262_1064", "0262_1084", "0262_1104", "0262_1124",
                    "0263", "0263_1065", "0263_1105", "0263_1125", "0263_1085",
                    "0264", "0264_1066", "0264_1086", "0264_1106", "0264_1126",
                    "0370", "0370_1061", "0370_1101", "0370_1121", "0370_1081",
                    "0372", "0372_1063", "0372_1083", "0372_1123", "0372_1103",
                    "0508", "0508_1062", "0508_1102", "0508_1122", "0508_1082",
                    "0614", "0614_1087", "0614_1107", "0614_1127", "0614_1067",
                    "0615", "0615_1088", "0615_1108", "0615_1128", "0615_1068",
                    "0616", "0616_1089", "0616_1109", "0616_1129", "0616_1069",
                    "0745", "0745_1070", "0745_1090", "0745_1110", "0745_1130",
                    "0746", "0746_1071", "0746_1091", "0746_1111", "0746_1131",
                    "0749", "0749_1072", "0749_1092", "0749_1112", "0749_1132",
                    "0751", "0751_1073", "0751_1093", "0751_1113", "0751_1133",
                    "0917", "0917_1074", "0917_1094", "0917_1134", "0917_1114",
                    "0918", "0918_1075", "0918_1095", "0918_1115", "0918_1135",
                    "0920", "0920_1076", "0920_1096", "0920_1116", "0920_1136",
                    "0921", "0921_1097", "0921_1137", "0921_1117")
 rep_index <- which(rownames(asv_trim) %in% replicates_all)
 asv_reps <- asv_trim[c(rep_index), ]
 dim(asv_reps)
# 96 55261
 
 asv_reps_rel <- asv_reps
 for (i in 1:dim(asv_reps)[1]){
   asv_reps_rel[i, ] <- asv_reps[i, ]/sum(asv_reps[i,])
 }
 rowSums(asv_reps_rel)
 all(row.names(asv_reps_rel) == row.names(asv_reps))

 # Remove leading zeros only from row names that don't contain "_"
 row.names(asv_reps) <- ifelse(
   grepl("_", row.names(asv_reps)),
   row.names(asv_reps),  # keep as-is if it contains "_"
   sub("^0+", "", row.names(asv_reps))  # remove leading zeros otherwise
 )
 
 # Remove leading zeros only from row names that don't contain "_"
 row.names(asv_reps_rel) <- ifelse(
   grepl("_", row.names(asv_reps_rel)),
   row.names(asv_reps_rel),  # keep as-is if it contains "_"
   sub("^0+", "", row.names(asv_reps_rel))  # remove leading zeros otherwise
 )
 
 all(row.names(asv_reps_rel) == row.names(asv_reps))

#### asv_reps_rel.csv ----
write.csv(asv_reps_rel, "../data/asv_reps_rel.csv")

# subset metadata for replicates only
# replace hyphens with _ in row.names
metadata2 <- metadata
rownames(metadata2) <- sub("-", "_", rownames(metadata2)) 
metadata_reps <- metadata2[which(row.names(metadata2) %in% row.names(asv_reps_rel)), ]

# reorder metadata_reps to match asv_reps_rel order
asv_reps_rel_order <- row.names(asv_reps_rel)
metadata_reps <- metadata_reps[asv_reps_rel_order, ]
all(row.names(asv_reps_rel) == row.names(metadata_reps))

#### metadata_reps.csv ----
write.csv(metadata_reps, "../data/metadata_reps.csv")


# remove replicate samples from asv table (retain those with highest number of reads)
replicates <- c("0104", "0104_1078", "0104_1098", "0105", "0105_1059",
                "0105_1099", "0105_119",
                "0106", "0106_1060", "0106_1100", "0106_1120",
                "0262", "0262_1064", "0262_1084", "0262_1104",
                "0263", "0263_1065", "0263_1105", "0263_1125",
                "0264", "0264_1066", "0264_1086", "0264_1106",
                "0370", "0370_1061", "0370_1101", "0370_1121",
                "0372", "0372_1063", "0372_1083", "0372_1123",
                "0508", "0508_1062", "0508_1102", "0508_1122",
                "0614", "0614_1087", "0614_1107", "0614_1127",
                "0615", "0615_1088", "0615_1108", "0615_1128",
                "0616", "0616_1089", "0616_1109", "0616_1129",
                "0745", "0745_1070", "0745_1090", "0745_1110",
                "0746", "0746_1071", "0746_1091", "0746_1111",
                "0749", "0749_1072", "0749_1092", "0749_1112",
                "0751", "0751_1073", "0751_1093", "0751_1113",
                "0917", "0917_1074", "0917_1094", "0917_1134",
                "0918", "0918_1075", "0918_1095", "0918_1115",
                "0920", "0920_1076", "0920_1096", "0920_1116",
                "0921", "0921_1097", "0921_1137")
drop_index <- which(rownames(asv_trim) %in% replicates)
asv_trim <- asv_trim[-c(drop_index), ]
dim(asv_trim)
# 975 55261
# remove sample 0467 (no metadata)
asv_trim <- asv_trim[-which(rownames(asv_trim)=="0467"),]
dim(asv_trim)
# 974 55261

### Controls ----
# remove mock communities, negative controls
mock <- rownames(asv_trim[grep("^Mock",rownames(asv_trim)),])
negative <- rownames(asv_trim[grep("^Negative",rownames(asv_trim)),])
asv_trim <- asv_trim[-which(rownames(asv_trim) %in% c(mock, negative)), ]
dim(asv_trim)
# 967 55261

# remove leading 0s from sample names 
rownames(asv_trim) <- sub("^0+", "", rownames(asv_trim)) 

# update metadata for asv_trim
# remove 105_1079
asv_trim <- asv_trim[-which(row.names(asv_trim)=="105_1079"),]
dim(asv_trim)
# 966 55261

# # rename replicates so IDs match metadata
rename <- row.names(asv_trim[grep("\\_",row.names(asv_trim)),])
newnames <- c("25", "26", "28", "29", "30", "31", "32", "33", "34", "35", "36", 
              "37", "38", "39", "41", "42", "43", "44", "45", "46", "47", "48",
              "49", "50", "51", "52", "53", "54", "55", "56", "61", "62", "63",
              "66", "67", "68", "104", "105", "106", "262", "263", "264", "370",
              "372", "508", "614", "615", "616", "745", "746", "749", "751",
              "917", "918", "920", "921", "945", "1444", "1447", "1448", "1490",
              "1491", "1492", "1493", "1494", "1495", "1496", "1497", "1498",
              "1499", "1500", "407")
 rownames(asv_trim)[which(rownames(asv_trim) %in% rename)] <- newnames 

# remove "late phase" samples
nolate <- metadata[metadata$Phase == "early",]
asv_trim <- asv_trim[which(row.names(asv_trim) %in% row.names(nolate)),]
dim(asv_trim)
# 908 55261

metadata_trim <- metadata[which(row.names(metadata) %in% rownames(asv_trim)),]
dim(metadata_trim)
# 908 16

# create presence/absence and relative abundance matrices
asv_pa <- (asv_trim > 0) * 1
rowSums(asv_pa)

asv_rel <- asv_trim
for (i in 1:dim(asv_trim)[1]){
  asv_rel[i, ] <- asv_trim[i, ]/sum(asv_trim[i,])
}
rowSums(asv_rel)
all(row.names(asv_rel) == row.names(asv_trim))
all(row.names(asv_pa) == row.names(asv_trim))

## asv_rel.csv ----
write.csv(asv_rel, "../data/asv_rel.csv")

# Reorder based on asv table order
asv_order <- rownames(asv_trim)
metadata_trim <- metadata_trim[asv_order, ]
all(row.names(metadata_trim) == asv_order)

# update diet data for samples in asv_trim, metadata
diet_data <- diet_data[which(row.names(diet_data) %in% metadata_trim$BirdID),]

# combine metadata and diet data
metadata_diet <- metadata_trim
metadata_diet <- metadata_diet[which((metadata_diet$BirdID) %in% row.names(diet_data)),]
diet_data_order <- metadata_diet$BirdID
diet_data <- diet_data[diet_data_order, ]
all(row.names(diet_data) == metadata_diet$BirdID)
# TRUE
metadata_diet <- cbind(metadata_diet, diet_data)

# update asv table for only birds with diet data
asv_rel_diet <- asv_rel[which(row.names(asv_rel) %in% row.names(metadata_diet)),]

# update diet diversity data for samples in asv_trim, metadata
diet_diversity <- diet_diversity[which(row.names(diet_diversity) %in% metadata_trim$BirdID),]

# combine metadata and diet diversity data
metadata_diet_diversity <- metadata_trim
metadata_diet_diversity <- metadata_diet_diversity[which((metadata_diet_diversity$BirdID) %in% row.names(diet_diversity)),]
diet_div_order <- metadata_diet_diversity$BirdID
diet_diversity <- diet_diversity[diet_div_order, ]
all(row.names(diet_div_order) == metadata_diet_diversity$BirdID)
# TRUE
metadata_diet_diversity <- cbind(metadata_diet_diversity, diet_diversity)

# update asv table for only birds with diet diversity data
asv_rel_diet_diversity <- asv_rel[which(row.names(asv_rel) %in% row.names(metadata_diet_diversity)),]

# update health data with sample IDs
health_data2 <- health_data[which(row.names(health_data) %in% metadata_trim$BirdID),]

# combine metadata and diet diversity data
metadata_order <- metadata_trim$BirdID
health_data3 <- health_data2[metadata_order, ]
all(row.names(health_data3)==metadata_trim$BirdID)
row.names(health_data3) <- row.names(metadata_trim)
all(row.names(health_data3)==row.names(metadata_trim))
# Replace GizzardFFDM -99 with NA
health_data3[health_data3 == -99] <- NA
metadata_health <- cbind(health_data3,metadata_trim)

## metadata_health.csv ----
# Save file
write.csv(metadata_health, "../data/metadata_health.csv")


## Add taxonomy ----
# rename, subset tax table
tax_raw2 <- tax_raw
rownames(tax_raw2) <- colnames(asv_raw)
all(row.names(tax_raw2) == colnames(asv_raw))
tax_trim <- tax_raw2[which(rownames(tax_raw2) %in% colnames(asv_trim)),]
colnames(tax_trim) <- c("V1","taxonomy")
tax_trim <- tax_trim %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement=""))

# combine taxonomy into genera, put ASVs into column, sum asv table for those ASVs to get genera abundance 
tax_trim2 <- tax_trim %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep=";")

# add taxonomy row to asv table
all(row.names(tax_trim)==colnames(asv_trim))
new_row <- tax_trim$taxonomy
asv_trim2 <- rbind(asv_trim,new_row)
row.names(asv_trim2)[909] <- "taxonomy"
asv_trim3 <- as.data.frame(asv_trim2)


# Alpha Diversity ----
summary(rowSums(asv_trim))

S <- specnumber(asv_trim) # observed number of species
raremax <- min(rowSums(asv_trim))

asv_rare <- rrarefy(asv_trim, raremax)

Srare <- rarefy(asv_trim, raremax)

# Species Richness
richness <- rowSums((asv_rare >= 1))

# Shannon Diversity
shannon <- diversity(asv_rare, "shannon")

# Simpson's Evenness
simp.even <- apply(asv_rare, 1, simp_even)

# Diversity: combined richness, diversity, evenness
diversity <- cbind(metadata_trim,richness,shannon,simp.even)

## diversity.csv ----
write.csv(diversity, "../data/diversity.csv")



# Table 1 ----
# top 30 ASVs - taxonomy, % Prevalence, CV
tax_top30 <- tax_trim2[c(1:30),]
# remove V1 column
tax_top30 <- tax_top30[,-which(colnames(tax_top30)=="V1")]
# rename
tax_top30$genus[tax_top30$genus=="Clostridia_UCG-014_ge"] <- "Clostridia UCG-014"
tax_top30$order[tax_top30$order=="Clostridia_UCG-014"] <- "Clostridia UCG-014"
tax_top30$family[tax_top30$family=="Clostridia_UCG-014_fa"] <- "Clostridia UCG-014"
tax_top30$family[tax_top30$family=="Coriobacteriales_unclassified"] <- "Unclassified"
tax_top30$genus[tax_top30$genus=="Coriobacteriales_unclassified"] <- "Unclassified"
tax_top30$family[tax_top30$family=="uncultured"] <- "Uncultured"
tax_top30$genus[tax_top30$genus=="uncultured_ge"] <- "Uncultured"
tax_top30$family[tax_top30$family=="Oscillospirales_fa"] <- "Unclassified"
tax_top30$genus[tax_top30$genus=="Oscillospiraceae_unclassified"] <- "Unclassified"
tax_top30$genus[tax_top30$genus=="Lachnospiraceae_unclassified"] <- "Unclassified"
# capitalize columns
colnames(tax_top30) <- c("Kingdom","Phylum","Class","Order","Family","Genus")

# add % prevalence
asv_pc <- as.data.frame(colSums(asv_pa)/nrow(asv_pa))
asv_pc_top30 <- as.data.frame(asv_pc[c(1:30),])
tax_top30$"Prevalence (%)" <- asv_pc_top30$`asv_pc[c(1:30), ]`
# calculate CV, add
asv_top30 <- asv_rel[,which(colnames(asv_rel) %in% row.names(tax_top30))]
asv_top30_cv <- apply(asv_top30, 2, FUN = cv)
tax_top30$"Coefficient of Variation" <- asv_top30_cv
## tax_top30.csv ----
# save
write.csv(tax_top30, "tax_top30.csv")



# Beta Diversity ----
# Bray Curtis Distance
REL_dist <- vegdist(asv_rel, method="bray")

# Principal Coordinates Analysis 
# Classical (Metric) Multidimensional Scaling; returns PCoA coordinates
# eig=TRUE returns eigenvalues; k = # of dimensions to calculate
PCoA_BC <- cmdscale(REL_dist, k=3, eig=TRUE, add=FALSE)

## Statistics  ----
# age
set.seed(38)
age <- as.factor(metadata_trim$Age)
asv.age<-adonis2(REL_dist~age, permutations=999)

#sex
set.seed(38)
asv.sex<-adonis2(REL_dist~as.factor(metadata_trim$Sex), data=metadata_trim, permutations=999)

asv.age.sex<-adonis2(REL_dist~ Sex * Age, data=metadata_trim, permutations=999)

# year
asv.year<-adonis2(REL_dist ~ CollYear, data=metadata_trim, permutations=999)

## Indicator taxa ----
# Which ASVs are most strongly correlated with beta diversity?
env_top30 = envfit(PCoA_BC, asv_top30, permutations = 999, strata = NULL, choices=c(1,2), display = "sites", w = weights(PCoA_BC, display), na.rm = TRUE)

# three indicator taxa (r2 > 0.6)
# "ASV000001","ASV000002", "ASV000005"


# Diet ----
metadata_diet_annual <- metadata_diet %>%
  group_by(CollYear) %>%
  summarize(Fruits_mean=mean(Fruits_pcbiomass),
            Infructescence_mean=mean(Infructescence_pcbiomass),
            Catkins_mean=mean(Catkins_pcbiomass),
            Leaves_mean=mean(Leaves_pcbiomass),
            `Dry oct leaf_mean`=mean(`Dry.oct.leaf_pcbiomass`))
row.names(metadata_diet_annual) <- metadata_diet_annual$CollYear

metadata_health$CollYear <- as.factor(metadata_health$CollYear)
metadata_diet_annual$CollYear <- as.factor(metadata_diet_annual$CollYear)

metadata_health_diet_annual <- metadata_health %>%
  left_join(metadata_diet_annual, by = "CollYear")
row.names(metadata_health_diet_annual) <- row.names(metadata_health)
diet_means <- c("Fruits_mean","Infructescence_mean", "Leaves_mean")
metadata_diet_annual_means <- metadata_health_diet_annual[,which(colnames(metadata_health_diet_annual) %in% diet_means
)]

## metadata_diet.csv ----
write.csv(metadata_diet_annual_means, "metadata_diet.csv")

# diet data - what % of birds had either fruits, infructscence, or leaves 
# detected in the crop?
non_zero_any <- with(diet_data, 
                     +                      Fruits_pcbiomass != 0 | 
                       +                          Infructescence_pcbiomass != 0 | 
                       +                          Leaves_pcbiomass != 0)
