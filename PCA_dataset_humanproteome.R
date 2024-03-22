#library(readxl)
library(ggfortify)
library(ggplot2)
library(factoextra)
library(MetBrewer)
library(tidyverse)
library(data.table)

#setwd("/Users/ayyamibrahim/Documents/Rhcalc/PCA_dataset.txt")
#PCAdata <- read.csv("PCA_dataset_proteome.csv")

PCAdata<- fread("/Users/ayyamibrahim/Documents/Rhcalc/PCA_dataset.txt")

#PCAdata <- fread("PCA_dataset.txt")

# Get PCAs
 names(PCAdata)[12] <- "sheet.3"
pca <- prcomp(PCAdata[,3:36], center = TRUE, scale = TRUE)
dim(pca)
# this command gives a Scree plot
#fviz_eig(pca, barcolor = "black", barfill = "white")

# get the data
scree_data <- summary(pca)$importance[,1:10]
colnames(scree_data) <- as.character(1:10)
scree_data <- as.data.frame(scree_data)
scree_data$type <- rownames(scree_data)
long_scree_data <- scree_data %>%
  pivot_longer(cols = 1:10, 
               names_to = "PC",
               values_to = "value") %>%
  filter(type == "Proportion of Variance")

long_scree_data$PC <- factor(long_scree_data$PC,
                             levels = long_scree_data$PC,
                             ordered = T)

# Plot
ggplot(long_scree_data, aes(x = PC, y = value)) +
  geom_col(fill = "light grey", color = "black") +
  geom_boxplot(width = 0.15) +
  geom_line(group = 1, size = 2.5) +
  geom_point(color = "black", size = 4) +
  theme_classic() + theme_minimal() +
  ylab("Explained Variances") + xlab("Dimensions") +
  ggtitle("Scree Plot of Human Proteome") + 
  theme(plot.title = element_text(hjust = 0.5,size=25)) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22))

ggsave("Scree_plot_humanproteome.png")
# 
# 
# PCA biplot
names(pca$center)
thick <- c(2,2,rep(0.5,26),2,rep(0.4,5))
grp <- c("black",rep("#4000FF",15),rep("#00A300",18))
#names(PCAdata)[10] <- "sheet.3"
fviz_pca_var(pca, col.var = grp, 
             labelsize=4,
             arrowsize=thick,
             palette=c("#00A300", "#4000FF", "black"),
             #gradient.cols = c("#2271B2", "darkgreen", "black"), 
             repel = TRUE, geom = "arrow") +
  ggtitle("PCA of Human Proteome") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        text = element_text(size=20))
  #annotate("text",x=-0.5,y=-0.9,label="Nu-Model ")
  
ggsave("PCA_dataset_humanproteome.png")
ggsave("test.png", width=5, height=7, units="cm", dpi=300)

# Check which are the categories for the thick arrows
thick

thick_arrows <- thick == 2

rownames(pca$rotation)[thick_arrows]
             
