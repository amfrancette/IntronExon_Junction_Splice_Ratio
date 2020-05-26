library(tidyverse)
library(here)
library(dplyr)
library(ggplot2)
here()
here::set_here()
print(getwd())
here::dr_here()

sample_names <- c(as.matrix(read.csv("res/sample_names.txt", header = F, sep = "\n")))
totalCounts <- read.csv("res/totalSpliceJunctionCount.tab", header = F, sep = "\t")
intronCounts <- read.csv("res/unsplicedSpliceJunctionCount.tab", header = F, sep = "\t")

colnames(totalCounts) <- c("feature","locus", sample_names)
colnames(intronCounts) <- c("feature","locus", sample_names)
number_of_samples <- ncol(totalCounts)-2

# to only consider junctions with reasonable sequencing depth, this command removes any loci where the minimum 
# less than a specified number (minimumDepthCutoff) of reads crossing the 5' splice junction 
depthCutoff <- 5
message("Depth Cutoff = ", depthCutoff)
passesCutoff <- (apply(totalCounts[,3:(2+number_of_samples)], 1, FUN=median)) >= depthCutoff
# can change cutoff to be median to only exclude loci with poor coverage in general
# passesCutoff <- (apply(totalCounts[,3:(2+number_of_samples)], 1, FUN=median)) >= depthCutoff

totalCounts_filt <- totalCounts[passesCutoff,]
intronCounts_filt <- intronCounts[passesCutoff,]
message ("Number of Introns Pre-Filtering: ", dim(intronCounts)[1])
message ("Number of Introns Passing Filter: ", dim(intronCounts_filt)[1])

IEJR <- cbind(totalCounts_filt[,1:2],
               intronCounts_filt[,3:(2+number_of_samples)]/totalCounts_filt[,3:(2+number_of_samples)])

# replacing Na with 0 
#IEJR[is.na(IEJR)] <- 0

write.csv(IEJR, file = "IEJR.csv", row.names = FALSE)

# >----- Example IEJR Analysis of Barass et al. Data -------<
# averaging replicates
IEJR_avg <- cbind(totalCounts_filt[,1:2], rowMeans(IEJR[,3:5]), 
                  rowMeans(IEJR[,6:8]), rowMeans(IEJR[,9:11]), IEJR[,12])
# appending correct column names to averaged matrix
colnames(IEJR_avg)[3:6] <- c("1.5min", "2.5min", "5min", "Total")

# make biplots comparing replicates
ggplot(IEJR, aes(x = `1point5min_Rep1`, y = `1point5min_Rep2`)) + geom_point(color = "black") +
  scale_x_continuous(name="1point5min_Rep1", limits = c(0,1)) +
  scale_y_continuous(name="1point5min_Rep2", limits = c(0,1)) +
  theme_grey(base_size = 15) +
  ggtitle("Comparing Replicates w/ Biplot") +
  theme(legend.position = "none") +
  geom_abline(slope=1, intercept = 0)

# make biplots comparing conditions
ggplot(IEJR_avg, aes(x = `1.5min`, y = `2.5min`)) + geom_point(color = "black") +
  scale_x_continuous(name="1.5min", limits = c(0,2)) +
  scale_y_continuous(name="2.5min", limits = c(0,2)) +
  theme_grey(base_size = 10) +
  ggtitle("Comparing Conditions w/ Biplot") +
  theme(legend.position = "none") +
  geom_abline(slope=1, intercept = 0)

# calculating log2 difference of conditions relative to 5 min labeling
methodDif <- cbind(IEJR_avg[,1:2], log2(IEJR_avg[,3]/IEJR_avg[,5]), 
                   log2(IEJR_avg[,4]/IEJR_avg[,5]), 
                   log2(IEJR_avg[,6]/IEJR_avg[,5]))

# appending column names to matrix
colnames(methodDif)[3:5] <- c("1.5vs5", "2.5vs5", "Totalvs5")

# replacing Na with 0 
methodDif[is.na(methodDif)] <- 0

# Generates Long-form table for easier plot building
IEJR_long <- pivot_longer(IEJR, -feature:-locus, names_to = "Condition", values_to = "IEJR" )

# # This should plot individual IEJR distributuions of everything used as input
# ggplot(IEJR_long[IEJR_long$feature=="gene",], aes(x = `IEJR`, fill = `Condition`)) + geom_density(color="black", alpha = 0.3) +
#   labs(y="Density of Observations", x = "IEJR Distributions") +
#   ggtitle("IEJR Distributions For All Samples") +
#   theme_grey(base_size = 15) +
#   geom_vline(xintercept=0) + xlim(0, 1.5) + coord_cartesian(ylim=c(0, 6)) 

# plots individual IEJR distribtutions
ggplot(IEJR[IEJR$feature=="gene",], aes(x = `Total_Rep2`)) + geom_density(color="black", fill="red", alpha = 0.3) +
  geom_density(data = IEJR[IEJR$feature=="gene",], aes(x = `Total_Rep3`) ,color="black", fill="red", alpha = 0.3) +
  geom_density(data = IEJR[IEJR$feature=="gene",], aes(x = `5min_Rep1`) ,color="black", fill="orange", alpha = 0.3) +
  geom_density(data = IEJR[IEJR$feature=="gene",], aes(x = `5min_Rep2`) ,color="black", fill="orange", alpha = 0.3) +
  geom_density(data = IEJR[IEJR$feature=="gene",], aes(x = `5min_Rep3`) ,color="black", fill="orange", alpha = 0.3) +
  geom_density(data = IEJR[IEJR$feature=="gene",], aes(x = `2point5min_Rep1`) ,color="black", fill="yellow", alpha = 0.3) +
  geom_density(data = IEJR[IEJR$feature=="gene",], aes(x = `2point5min_Rep2`) ,color="black", fill="yellow", alpha = 0.3) +
  geom_density(data = IEJR[IEJR$feature=="gene",], aes(x = `2point5min_Rep3`) ,color="black", fill="yellow", alpha = 0.3) +
  geom_density(data = IEJR[IEJR$feature=="gene",], aes(x = `1point5min_Rep1`) ,color="black", fill="green", alpha = 0.3) +
  geom_density(data = IEJR[IEJR$feature=="gene",], aes(x = `1point5min_Rep2`) ,color="black", fill="green", alpha = 0.3) +
  geom_density(data = IEJR[IEJR$feature=="gene",], aes(x = `1point5min_Rep3`) ,color="black", fill="green", alpha = 0.3) +
  labs(y="Density of Observations", x = "IEJR Distributions") +
  ggtitle("IEJR Distributions For All Samples") +
  theme_grey(base_size = 15) +
  geom_vline(xintercept=0) + xlim(0, 1.5) + coord_cartesian(ylim=c(0, 6)) 


# Generates Long-form table for easier plot building
IEJR_avg_long <- pivot_longer(IEJR_avg, -feature:-locus, names_to = "Condition", values_to = "IEJR" )

png("imgs/Barass_IEJR_density.png")
ggplot(IEJR_avg_long[IEJR_avg_long$feature=="gene",], aes(x = `IEJR`, fill = `Condition`)) + 
  geom_density(color="black", alpha = 0.3) +
  labs(y="Density of Observations", x = "IEJR Distributions by Condition") +
  theme_grey(base_size = 15) +
  theme(legend.background = element_rect()) +
  geom_vline(xintercept=0) + xlim(0, 1.2) + coord_cartesian(ylim=c(0, 4)) +
  scale_fill_manual(values=c("green", "yellow", "orange", "red"), labels = c("1.5 min 4tU", "2.5 min 4tU", "5 min 4tU", "Total RNA")) 
dev.off()

# Generates Long-form table for easier plot building
methodDiff_long <- pivot_longer(methodDif, -feature:-locus, names_to = "Condition", values_to = "log2fc" )

pdf("imgs/Barass_methodDif_density.pdf")
ggplot(methodDiff_long[methodDiff_long$feature=="gene",], aes(x = `log2fc`, fill = `Condition`)) + 
  geom_density(color="black", alpha = 0.3) +
  labs(y="Density of Observations", x = "log2fc") +
  scale_fill_manual(values=c("green", "orange", "red"), labels = c("1.5 min / 5 min", "2.5 min / 5 min", "Total RNA / 5 min")) +
  theme_grey(base_size = 15) +
  geom_vline(xintercept=0)  + xlim(-2, 4)
dev.off()

message("column median IEJR values: ")
apply(IEJR[3:(2 + number_of_samples)], 2, median, na.rm = T)


