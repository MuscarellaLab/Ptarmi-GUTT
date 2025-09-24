# 
# GUTT-P Figure Generation
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
library(otuSummary)
library(corrplot)
library(naniar)
library(car)
library(randomForest)
library(ggsignif)

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

# Load data
metadata_health <- read.csv("metadata_health.csv", row.names = 1)
diversity <- read.csv("diversity.csv", row.names = 1)
asv_rel <- read.csv("asv_rel.csv", row.names = 1)
metadata_diet <- read.csv("metadata_diet.csv", row.names = 1)
gizzard_new <- read.csv("2018GizzardFFDM.csv", row.names = 1)
asv_reps_rel <- read.csv("asv_reps_rel.csv", row.names = 1)
metadata_reps <- read.csv("metadata_reps.csv", row.names = 1)


# Diet digestibility ----

## gizzard mass ----
metadata_health$Age <- as.factor(metadata_health$Age)
metadata_health$Sex <- as.factor(metadata_health$Sex)
metadata_health$CollYear <- as.factor(metadata_health$CollYear)

colnames(metadata_health)[colnames(metadata_health) == 'GizzardEmptyFFDM'] <- 'Gizzard Mass'
# Update 2018 gizzard mass
# Ensure BirdID is character
metadata_health$BirdID <- as.character(metadata_health$BirdID)

# Filter indices of birds collected in 2018
rows_2018 <- which(metadata_health$CollYear == 2018)

# BirdIDs for 2018 birds
bird_ids_2018 <- metadata_health$BirdID[rows_2018]

# Match BirdIDs to gizzard_new rownames
match_indices <- match(bird_ids_2018, rownames(gizzard_new))

# Valid matches (non-NA)
valid <- !is.na(match_indices)

# Replace Gizzard Mass in metadata_health with GizzardEmptyFFDM values from gizzard_new
metadata_health$`Gizzard Mass`[rows_2018[valid]] <- gizzard_new$GizzardEmptyFFDM[match_indices[valid]]

colnames(metadata_health)[colnames(metadata_health) == 'TotalFat'] <- 'Total Fat'

# sex
var.test(metadata_health$`Gizzard Mass` ~ metadata_health$Sex, alternative = "two.sided")
t.test(metadata_health$`Gizzard Mass` ~ metadata_health$Sex, var.equal = FALSE)
# standard deviation
library(dplyr)

metadata_health %>%
  group_by(Sex) %>%
  summarise(
    sd_gizzard_mass = sd(`Gizzard Mass`, na.rm = TRUE)
  )


# age
var.test(metadata_health$`Gizzard Mass` ~ metadata_health$Age, alternative = "two.sided")
t.test(metadata_health$`Gizzard Mass` ~ metadata_health$Age, var.equal = TRUE)
# standard deviation 
metadata_health %>%
  group_by(Age) %>%
  summarise(
    sd_gizzard_mass = sd(`Gizzard Mass`, na.rm = TRUE)
  )

# sex*age
summary(aov(metadata_health$`Gizzard Mass` ~ metadata_health$Age*metadata_health$Sex))
# year
summary(aov(metadata_health$`Gizzard Mass` ~ metadata_health$CollYear))
# sex*age*year
summary(aov(metadata_health$`Gizzard Mass` ~ metadata_health$Sex*metadata_health$Age*metadata_health$CollYear))


# regress gizzard mass and body size, save residuals
gizzard_body_lm <- lm(metadata_health$`Gizzard Mass` ~ metadata_health$BodySize)
# keep NAs for rows without gizzard mass
rows_with_gizzard_mass <- !is.na(metadata_health$`Gizzard Mass`)
metadata_health$gizzard_body_residuals[rows_with_gizzard_mass] <- gizzard_body_lm$residuals

# rerun t.test with gizzard mass adjusted for body size
var.test(metadata_health$gizzard_body_residuals ~ metadata_health$Age, alternative = "two.sided")
t.test(metadata_health$gizzard_body_residuals ~ metadata_health$Sex, var.equal = TRUE)
t.test(metadata_health$gizzard_body_residuals ~ metadata_health$Age, var.equal = TRUE)

# rerun anova with gizzard mass adjusted for body size
summary(aov(metadata_health$gizzard_body_residuals ~ metadata_health$Age*metadata_health$Sex))
summary(aov(metadata_health$gizzard_body_residuals ~ metadata_health$CollYear))
summary(aov(metadata_health$gizzard_body_residuals ~ metadata_health$CollYear*metadata_health$Age*metadata_health$Sex))

# Body condition ----
## total fat ----
# sex
var.test(metadata_health$`Total Fat` ~ metadata_health$Sex, alternative = "two.sided")
t.test(metadata_health$`Total Fat` ~ metadata_health$Sex, var.equal = FALSE)
# standard deviation 
metadata_health %>%
  group_by(Sex) %>%
  summarise(
    sd_total_fat = sd(`Total Fat`, na.rm = TRUE)
  )

# age
var.test(metadata_health$`Total Fat` ~ metadata_health$Age, alternative = "two.sided")
t.test(metadata_health$`Total Fat` ~ metadata_health$Age, var.equal = TRUE)

# year
summary(aov(metadata_health$`Total Fat` ~ metadata_health$CollYear))
# sex*age*year
summary(aov(metadata_health$`Total Fat` ~ metadata_health$Sex*metadata_health$Age*metadata_health$CollYear))


# Microbiome ----
## Alpha diversity ----
diversity$Age <- as.factor(diversity$Age)
diversity$Sex <- as.factor(diversity$Sex)
diversity$CollYear <- as.factor(diversity$CollYear)

# age
var.test(diversity$shannon ~ diversity$Age, alternative = "two.sided")
t.test(diversity$shannon ~ diversity$Age, var.equal = TRUE)
# standard deviation 
diversity %>%
  group_by(Age) %>%
  summarise(
    sd_diversity_age = sd(shannon, na.rm = TRUE)
  )

# sex
var.test(diversity$shannon ~ diversity$Sex, alternative = "two.sided")
t.test(diversity$shannon ~ diversity$Sex, var.equal = TRUE)
# standard deviation
diversity %>%
  group_by(Sex) %>%
  summarise(
    sd_diversity_sex = sd(shannon, na.rm = TRUE)
  )

# year
shannon_year <- aov(diversity$shannon ~ diversity$CollYear)
TukeyHSD(shannon_year, ordered = TRUE)
summary(aov(diversity$shannon ~ diversity$Sex*diversity$Age*diversity$CollYear))


# Fig 1A ----
c("#88CCEE","#332288","#CC6677","#882255","#DDCC77","#117733","#44AA99","#AA4499")


# Set up dimensions
dpi <- 600
width_cm <- 8.5
aspect_ratio <- 0.75
height_cm <- width_cm * aspect_ratio

# Convert to inches
width_in <- width_cm / 2.54
height_in <- height_cm / 2.54

# Font sizes
label_pt <- 8    # Axis label font size
tick_pt <- 6     # Tick label font size

# Convert to cex
cex.lab <- label_pt / 12
cex.axis <- tick_pt / 12

# Open JPEG device
jpeg("../figures/Fig1A.jpg", width = width_in, height = height_in, units = "in", res = dpi)

# Base plotting parameters
par(mfrow = c(1, 1),
    mar = c(3.2, 4.2, 1, 1),      # Add a bit more space on the left (y-axis)
    cex.lab = cex.lab,
    cex.axis = cex.axis,
    mgp = c(2.5, 0.6, 0))         # More generous spacing for y-axis label & ticks

# Suppress default axes
boxplot(diversity$shannon ~ diversity$CollYear,
        xlab = "", ylab = "Diversity",
        las = 1,
        ylim = c(1.5, 4.25),
        col = "#882255",
        axes = FALSE)

# Draw y-axis (uses global mgp for better spacing)
axis(2, las = 1, cex.axis = cex.axis)

# Draw x-axis with tighter spacing (custom mgp)
axis(1, at = 1:length(unique(diversity$CollYear)),
     labels = sort(unique(diversity$CollYear)),
     tick = TRUE,
     cex.axis = cex.axis,
     mgp = c(1.6, 0.2, 0))  # Tighter tick-label distance for x-axis

# Manually add x-axis label
mtext("Year", side = 1, line = 1.6, cex = cex.lab)

# Draw box
box(lwd = 2)

# Add 'A' to the top-left corner outside the plot
mtext("(a)", side = 3, adj = -0.25, line = 0, cex = 1.2, font = 2)

# Close the device
dev.off()


# Rename factor levels for Age and Sex
diversity$Age <- factor(diversity$Age, levels = c("A", "J"), labels = c("Adult", "Juvenile"))
diversity$Sex <- factor(diversity$Sex, levels = c("F", "M"), labels = c("Female", "Male"))



# Fig 1B ----
# Create side-by-side boxplots

# === Setup ===
dpi <- 600
total_width_cm <- 8.5
aspect_ratio <- 0.75     # same as Fig 1A for consistent height
total_height_cm <- total_width_cm * aspect_ratio

width_in <- total_width_cm / 2.54
height_in <- total_height_cm / 2.54

label_pt <- 8
tick_pt <- 6
cex.lab <- label_pt / 12
cex.axis <- tick_pt / 12

jpeg("./figures/Fig1B.jpg", width = width_in, height = height_in, units = "in", res = dpi)

# === Set up plotting area ===
par(mfrow = c(1, 2),
    mar = c(4.2, 4.2, 1, 0.2),  # more bottom space, less right margin
    cex.lab = cex.lab,
    cex.axis = cex.axis,
    mgp = c(2.2, 0.6, 0),       # slight upward shift of labels
    xaxs = "i",
    yaxs = "i")

# === Boxplot for Sex ===
boxplot(diversity$shannon ~ diversity$Sex,
        col = c("#332288", "#117733"),
        axes = FALSE,
        xlab = "",
        ylab = "Diversity",
        ylim = c(1.5, 4.5),
        boxwex = 0.9,          # make boxplots wider
        at = c(1, 2))
box(lwd = 2)

# Axes for Sex plot
axis(2, las = 1, cex.axis = cex.axis)
axis(1, at = c(1, 2), labels = c("Female", "Male"), cex.axis = cex.axis)
mtext("Sex", side = 1, line = 2.4, cex = cex.lab)  # move label down a bit

# === Boxplot for Age ===
boxplot(diversity$shannon ~ diversity$Age,
        col = c("#88CCEE", "#44AA99"),
        axes = FALSE,
        xlab = "",
        ylab = "",
        ylim = c(1.5, 4.5),
        boxwex = 0.9,
        at = c(1, 2))
box(lwd = 2)

# Axes for Age plot
axis(2, las = 1, cex.axis = cex.axis)
axis(1, at = c(1, 2), labels = c("Adult", "Juvenile"), cex.axis = cex.axis)
mtext("Age", side = 1, line = 2.4, cex = cex.lab)

# Add 'B' to the top-left corner outside the plot
mtext("(b)", side = 3, adj = -4.15, line = 0, cex = 1.2, font = 2)

dev.off()



# Fig 2A ----
## Bray Curtis Distance
REL_dist <- vegdist(asv_rel, method="bray")

# Principal Coordinates Analysis 
# Classical (Metric) Multidimensional Scaling; returns PCoA coordinates
# eig=TRUE returns eigenvalues; k = # of dimensions to calculate
PCoA_BC <- cmdscale(REL_dist, k=3, eig=TRUE, add=FALSE)

explainvar1 <- round(PCoA_BC$eig[1] / sum(PCoA_BC$eig), 3) * 100
explainvar2 <- round(PCoA_BC$eig[2] / sum(PCoA_BC$eig), 3) * 100
explainvar3 <- round(PCoA_BC$eig[3] / sum(PCoA_BC$eig), 3) * 100
sum.eig <- sum(explainvar1, explainvar2)

explainvar1 #29.2
explainvar2 #10.2

# envfit with health data (2024)
health <- c("Wing_mm", "Head_mm", "SternL_mm", "BodySize", "Gizzard Mass", "Total Fat", "LDBM", "Oxalat", "OxalatIntens", "DigestivTrLength..SI.LI.2xceca.")
metadata_health_env <- metadata_health[which(colnames(metadata_health) %in% health)]

names(metadata_health_env)[names(metadata_health_env) == 'DigestivTrLength..SI.LI.2xceca.'] <- 'Gut Length'


env_health = envfit(PCoA_BC, metadata_health_env, permutations = 999, strata = NULL, choices=c(1,2), display = "sites", w = weights(PCoA_BC, display), na.rm = TRUE)


# plotting Fig 2A


# Open JPEG device
jpeg("../figures/Fig2A.jpg", width = 18, height = 12, units = "cm", res = dpi)
par(mar = c(5, 6, 4, 2)) 
plot.new()
## 
shannon <- diversity$shannon
# Define breaks for the Shannon values (6 break points = 5 categories)
breaks <- seq(min(shannon), max(shannon), length.out = 6)  # 6 break points = 5 categories

# Use cut() with right = TRUE to create 5 categories (intervals)
shannon_cut <- cut(shannon, breaks = breaks, include.lowest = TRUE, labels = FALSE, right = TRUE)

# Define colors for each of the 5 categories
colors <- c('#88CCEE', '#44AA99', '#117733', "#DDCC77", "#CC6677")  # 5 colors for 5 categories

# Plot the points with colors based on the Shannon categories
plot(PCoA_BC$points[, 1], PCoA_BC$points[, 2], 
     xlab = "", ylab = "",
     xlim = c(-0.65, 0.58), ylim = c(-0.45, 0.55),
     pch = 21, 
     bg = colors[shannon_cut],  # Color the points based on the Shannon categories
     cex = 1.5, xaxt = "n", yaxt = "n", cex.lab = 1.5, cex.axis = 1.2)

# Add axes and lines
axis(side = 1, las = 1, label = TRUE, lwd.ticks = 2)
axis(side = 2, las = 1, label = TRUE, lwd.ticks = 2)
abline(h = 0, lty = "dotted")
abline(v = 0, lty = "dotted")
box(lwd = 2)

# Add health variable (assuming env_health exists)
plot(env_health, p.max = 0.001, col = "red")

# Add axis labels
mtext(paste("PCoA 1 (", explainvar1, "%)", sep = ""), side = 1, cex = 1.25, outer = FALSE, line = 2.5)
mtext(paste("PCoA 2 (", explainvar2, "%)", sep = ""), side = 2, cex = 1.25, outer = FALSE, line = 3.5)

# Create the legend to display only the 5 intervals (categories), not the break points
# Use dash ("-") instead of "to" for the intervals
# Custom circles in the legend
legend("bottomleft", 
       legend = paste(format(signif(breaks[-length(breaks)], digits = 3)), 
                      "-", 
                      format(signif(breaks[-1], digits = 3))),  # Use dash "-" between intervals
       pch = 21,  # Set pch = 21 to make points circles
       pt.bg = colors,  # Assign the color for each circle in the legend
       title = "Shannon Diversity", 
       bty = "n", pt.cex = 1.5, 
       title.cex = 0.8,  # Size of the legend title
       cex = 0.8)  # Size of the legend text


# Add 'A' to the top-left corner outside the plot
mtext("(a)", side = 3, adj = 0, line = 1, cex = 1.2, font = 2)


## END

dev.off()

# Count the number of samples in each category of shannon_cut
category_counts <- table(shannon_cut)

# View the result
print(category_counts)


# Fig 2B ----
# envfit
top3 <- c("ASV000001","ASV000002", "ASV000005")
top3_network_rel <- asv_rel[,which(colnames(asv_rel) %in% top3)]
colnames(top3_network_rel) <- c("Clostridia_UCG-014", "Actinomyces", "Coriobacteriales")
env_top3 = envfit(PCoA_BC, top3_network_rel, permutations = 999, strata = NULL, choices=c(1,2), display = "sites", w = weights(PCoA_BC, display), na.rm = TRUE)


# Define breaks for the Shannon values (6 break points = 5 categories)
breaks <- seq(min(shannon), max(shannon), length.out = 6)  # 6 break points = 5 categories

# Use cut() with right = TRUE to create 5 categories (intervals)
shannon_cut <- cut(shannon, breaks = breaks, include.lowest = TRUE, labels = FALSE, right = TRUE)

# Define colors for each of the 5 categories
colors <- c('#88CCEE', '#44AA99', '#117733', "#DDCC77", "#CC6677")  # 5 colors for 5 categories

# Open JPEG device
jpeg("../figures/Fig2B.jpg", width = 18, height = 12, units = "cm", res = dpi)

# Increase left margin to 6 lines (default is ~4)
 par(mar = c(5, 6, 4, 2))  # bottom, left, top, right

plot.new()

# Plot points
plot(PCoA_BC$points[, 1], PCoA_BC$points[, 2], 
     xlab = "", ylab = "",
     xlim = c(-0.65, 0.58), ylim = c(-0.45, 0.55),
     pch = 21, 
     bg = colors[shannon_cut],
     cex = 1.5, xaxt = "n", yaxt = "n", cex.lab = 1.5, cex.axis = 1.2)

# Add axes and lines
axis(side = 1, las = 1, lwd.ticks = 2)
axis(side = 2, las = 1, lwd.ticks = 2)
abline(h = 0, lty = "dotted")
abline(v = 0, lty = "dotted")
box(lwd = 2)

# Add health variable (assuming env_health exists)
plot(env_top3, p.max = 0.05, col = "red", cex = 0.8)

# Add axis labels
mtext(paste("PCoA 1 (", explainvar1, "%)", sep = ""), side = 1, cex = 1.25, line = 2.5)
mtext(paste("PCoA 2 (", explainvar2, "%)", sep = ""), side = 2, cex = 1.25, line = 3.5)  # Increase line for more space

# Legend and label
legend("bottomleft", 
       legend = paste(format(signif(breaks[-length(breaks)], digits = 3)), "-", format(signif(breaks[-1], digits = 3))),
       pch = 21,
       pt.bg = colors,
       title = "Shannon Diversity",
       bty = "n",
       pt.cex = 1.5,
       title.cex = 0.8,
       cex = 0.8)

# Add 'B' label top-left
mtext("(b)", side = 3, adj = 0, line = 1, cex = 1.2, font = 2)

dev.off()
## END



# Fig 3A ----
top3_rel_logit <- logit(top3_network_rel)
all(row.names(top3_rel_logit)==row.names(diversity))
diversity_top3 <- cbind(top3_rel_logit,diversity)


plot.new()
# Open JPEG device
jpeg("../figures/Fig3A.jpg", width = 8.5, height = 6.4, units = "cm", res = dpi, pointsize = 10)

# Adjust margins to maximize plot area and prevent label cutoff
par(mar = c(4.2, 4.5, 2.5, 1))  # Bottom, left, top, right

# Set font scaling relative to 10pt base size
par(cex.lab = 1)     # Axis labels = 10pt
par(cex.axis = 0.8)  # Tick labels ~8pt
par(cex.main = 1)    # (if you had a title)

# Create the initial plot
plot(diversity_top3$shannon ~ diversity_top3$Coriobacteriales,
     xlab = "Coriobacteriales",
     ylab = "Diversity",
     bg = "#AA4499",
     pch = 21,
     lwd.ticks = 2,
     cex = 1.5,
     cex.axis = 0.8,
     cex.lab = 1)

box(lwd = 2)

# Fit a smoothing spline to the data
spline_fit <- smooth.spline(diversity_top3$Coriobacteriales, diversity_top3$shannon)

# Add the spline curve to the plot
lines(spline_fit, col = "blue", lwd = 2)

# Spearman's
corr_corio <- cor.test(x=diversity_top3$Coriobacteriales, y=diversity_top3$shannon, method = 'spearman')

# Add 'A' to the top-left corner outside the plot
mtext("(a)", side = 3, adj = 0, line = 1, cex = 1.2, font = 2)

dev.off()



# Fig 3B ----

plot.new()
# Open JPEG device (in cm)
jpeg("../figures/Fig3B.jpg", width = 8.5, height = 6.4, units = "cm", res = dpi, pointsize = 10)

# Adjust margins to maximize plot area and prevent label cutoff
par(mar = c(4.2, 4.5, 2.5, 1))  # Bottom, left, top, right

# Set font scaling relative to 10pt base size
par(cex.lab = 1)     # Axis labels = 10pt
par(cex.axis = 0.8)  # Tick labels ~8pt
par(cex.main = 1)    # (if you had a title)

# Create the plot
plot(diversity_top3$shannon ~ diversity_top3$`Clostridia_UCG-014`,
     xlab = "Clostridia UCG-014",
     ylab = "Diversity",
     bg = "#DDCC77",
     pch = 21,
     lwd.ticks = 2,
     cex = 1.5,
     cex.axis = 0.8,
     cex.lab = 1)

# Add box around plot
box(lwd = 2)

# Fit smoothing spline
spline_fit <- smooth.spline(diversity_top3$`Clostridia_UCG-014`, diversity_top3$shannon)
lines(spline_fit, col = "blue", lwd = 2)

# Spearman correlation (optional)
corr_clost <- cor.test(
  x = diversity_top3$`Clostridia_UCG-014`,
  y = diversity_top3$shannon,
  method = 'spearman'
)

# Add 'B' to top-left corner, slightly lowered
mtext("(b)", side = 3, adj = 0, line = 1, cex = 1.2, font = 2)

# Close device
dev.off()




# Fig 3C ----

# Open JPEG device (in cm with 10pt font cap)
jpeg("../figures/Fig3C.jpg", width = 8.5, height = 6.4, units = "cm", res = dpi, pointsize = 10)

# Adjust margins to maximize plot space and prevent cutoff
par(mar = c(4.2, 4.5, 2.5, 1))  # Bottom, left, top, right

# Font scaling (relative to pointsize = 10)
par(cex.lab = 1)     # Axis labels at 10pt
par(cex.axis = 0.8)  # Tick labels ~8pt
par(cex.main = 1)    # Main title at 10pt (if used)

# Create the plot
plot(diversity_top3$shannon ~ diversity_top3$Actinomyces,
     xlab = "Actinomyces",
     ylab = "Diversity",
     bg = "#117733",
     pch = 21,
     lwd.ticks = 2,
     cex = 1.5,
     cex.axis = 0.8,
     cex.lab = 1)

# Draw box around plot
box(lwd = 2)

# Fit and draw smoothing spline
spline_fit <- smooth.spline(diversity_top3$Actinomyces, diversity_top3$shannon)
lines(spline_fit, col = "blue", lwd = 2)

# Spearman correlation
corr_actino <- cor.test(
  x = diversity_top3$Actinomyces,
  y = diversity_top3$shannon,
  method = 'spearman'
)

# Add 'C' to top-left corner, slightly lower
mtext("(c)", side = 3, adj = 0, line = 1, cex = 1.2, font = 2)

# Close device
dev.off()




# Fig 4 ----

# Wrangle data

# # remove NAs from diet data
metadata_diet_clean <- na.omit(metadata_diet)

# combine diet, health, and top3 data
health_diet <- metadata_health[which(row.names(metadata_health) %in% row.names(metadata_diet_clean)),]
top3_diet <- top3_rel_logit[which(row.names(top3_rel_logit) %in% row.names(metadata_diet_clean)),]
top3_diet <- as.data.frame(top3_diet)
all(row.names(health_diet)==row.names(top3_diet))
all(row.names(health_diet)==row.names(metadata_diet_clean))
health_diet_top3 <- cbind(metadata_diet_clean,health_diet,top3_diet)
diversity_diet <- diversity[which(row.names(diversity) %in% row.names(health_diet_top3)),]
all(row.names(health_diet_top3)==row.names(diversity_diet))
health_diet_top3_diversity <- cbind(health_diet_top3, diversity_diet)
health_diet_top3_diversity <- health_diet_top3_diversity[!is.na(health_diet_top3_diversity$`Gizzard Mass`), ]

# Random forest of total fat
rf <- randomForest(health_diet_top3_diversity$`Total Fat` ~ 
                     health_diet_top3_diversity$Coriobacteriales + 
                     health_diet_top3_diversity$shannon + 
                     health_diet_top3_diversity$Fruits_mean +
                     health_diet_top3_diversity$Infructescence_mean +
                     health_diet_top3_diversity$Leaves_mean + 
                     health_diet_top3_diversity$`Gizzard Mass` +
                     health_diet_top3_diversity$Actinomyces +
                     health_diet_top3_diversity$`Clostridia_UCG-014`, ntree = 10000, do.trace = T, importance = T)
print(rf) # view results
inc_mse <- as.data.frame(importance(rf, type = 1))
inc_mse <- inc_mse[order(inc_mse$`%IncMSE`, decreasing = F), ]


# rf 2 
plot.new()
# Open JPEG device
jpeg("../figures/Fig4.jpg", width = 8.5, height = 6.4, units = "cm", res = dpi, pointsize = 10)

# Adjust margins to maximize plot space and prevent cutoff
par(mar = c(4.5, 7.5, 2.5, 1))  # Bottom, left, top, right

# Font scaling (relative to pointsize = 10)
par(cex.lab = 1)     # Axis labels at 10pt
par(cex.axis = 0.7)  # Tick labels ~8pt
par(cex.main = 1)    # Main title at 10pt (if used)

barplot(inc_mse, horiz = T, 
        names.arg = rev(c("Leaves", "Fruits", "Infructescence", "Coriobacteriales", "Gizzard Mass", "Clostridia UCG-014", "Microbiome Diversity", "Actinomyces")), las = 1, xlim = c(0, 150),
        xlab = "% Increase Mean Square Error",
        col = "#88CCEE",
        lwd.ticks = 2)
box(lwd = 2)
dev.off()



# Fig S1 ----
# test for effect of seq run on BC 
# PCoA of replicate samples 

## Bray Curtis Distance
REL_dist_reps <- vegdist(asv_reps_rel, method="bray")

# Principal Coordinates Analysis 
# Classical (Metric) Multidimensional Scaling; returns PCoA coordinates
# eig=TRUE returns eigenvalues; k = # of dimensions to calculate
PCoA_BC_reps <- cmdscale(REL_dist_reps, k=3, eig=TRUE, add=FALSE)

explainvar1.reps <- round(PCoA_BC_reps$eig[1] / sum(PCoA_BC_reps$eig), 3) * 100
explainvar2.reps <- round(PCoA_BC_reps$eig[2] / sum(PCoA_BC_reps$eig), 3) * 100
explainvar3.reps <- round(PCoA_BC_reps$eig[3] / sum(PCoA_BC_reps$eig), 3) * 100
sum.eig.reps <- sum(explainvar1.reps, explainvar2.reps)

explainvar1.reps #34.4
explainvar2.reps #18

# add metadata column for replicate ID
# Start from the rownames
sample_names <- rownames(metadata_reps)

# Extract base sample ID: take everything before "_" if present, or the full name
base_ids <- sub("_.*", "", sample_names)

# Remove leading zeros
replicate_ids <- sub("^0+", "", base_ids)

# Add as a new column
metadata_reps$SampleID <- replicate_ids

# envfit with sequencing run 
metadata_reps_seqrun <- metadata_reps[which(colnames(metadata_reps)=="SeqRun")]
metadata_reps_seqrun$SeqRun <- as.factor(metadata_reps_seqrun$SeqRun)

env_seqrun_reps = envfit(PCoA_BC_reps, metadata_reps_seqrun, permutations = 999, strata = NULL, choices=c(1,2), display = "sites", w = weights(PCoA_BC_reps, display), na.rm = TRUE)


# compare original versus repeat run 
# is variation between pairs smaller than between samples 
# shape by seq run or color by sample


dev.off()
plot.new()

jpeg("../figures/FigS1.jpg")

# Plot the points with colors based on replicate ID and shapes based on sequencing run
  # Ensure correct metadata order
  metadata_ordered <- metadata_reps[rownames(PCoA_BC_reps$points), ]

# Factorize grouping variables
sample_factor <- factor(metadata_ordered$SampleID)
seqrun_factor <- factor(metadata_ordered$SeqRun)

# Load colors
library(RColorBrewer)

# Choose a good qualitative palette with strong contrast
n_reps <- length(levels(sample_factor))
color_palette <- brewer.pal(min(n_reps, 8), "Set1")  # Set1 is high-contrast
if (n_reps > 8) {
  color_palette <- colorRampPalette(brewer.pal(8, "Set1"))(n_reps)
}
colors <- color_palette[sample_factor]

# Use filled shapes (21–25), wrap if needed
shape_palette <- c(21, 22, 23, 24, 25, 21, 22)
shapes <- shape_palette[(as.integer(seqrun_factor) - 1) %% length(shape_palette) + 1]

# Setup plot area with room for legends
par(mar = c(5, 4, 4, 8), xpd = TRUE)  # expand right margin for legends

# Plot PCoA
plot(PCoA_BC_reps$points[, 1], PCoA_BC_reps$points[, 2], 
     xlab = "PCoA 1", ylab = "PCoA 2",
     pch = shapes,
     bg = colors,
     col = "black",
     cex = 1.5, lwd = 1.2,
     xaxt = "n", yaxt = "n")

# Axes and reference lines
axis(1, las = 1, lwd.ticks = 2)
axis(2, las = 1, lwd.ticks = 2)
box(lwd = 2)

# Legend for Replicate ID (colors), placed outside plot to the right
legend("topright", inset = c(-0.25, 0), legend = levels(sample_factor),
       pt.bg = color_palette, pch = 21, title = "Sample ID", 
       cex = 0.7, pt.cex = 1.2, bty = "n", x.intersp = 0.8)

# Legend for Sequencing Run (shapes), placed below it
legend("topright", inset = c(-0.1, 0), legend = levels(seqrun_factor),
       pch = shape_palette[1:length(levels(seqrun_factor))],
       pt.bg = "gray90", col = "black", title = "SeqRun", 
       cex = 0.7, pt.cex = 1.2, bty = "n", x.intersp = 0.8)

dev.off()


# PERMANOVA
rep.adonis <- adonis2(REL_dist_reps ~ sample_factor*seqrun_factor, data = metadata_ordered, permutations = 9999)



# Fig S2 ----
corr_keep <- c("shannon", "Gizzard Mass", "Total Fat", "Clostridia_UCG-014", "Actinomyces", "Coriobacteriales", "Fruits_mean", "Infructescence_mean", "Leaves_mean")
corr_data <- health_diet_top3_diversity[,which(colnames(health_diet_top3_diversity) %in% corr_keep)]

order <- c("Actinomyces", "Clostridia_UCG-014", "Coriobacteriales", "shannon", "Gizzard Mass", "Total Fat", "Fruits_mean", "Infructescence_mean", "Leaves_mean")
corr_data <- corr_data[, order]
colnames(corr_data) <- c("Actinomyces", "Clostridia UCG-014", "Coriobacteriales", "Microbiome Diversity", "Gizzard Mass", "Total Fat", "Fruits", "Infructescence", "Leaves")

plot.new()

jpeg("../figures/FigS2.jpg")
corr.data.corr.matrix <- cor(corr_data, method = "spearman", use = "pairwise.complete.obs")
corrplot(corr.data.corr.matrix)
dev.off()

