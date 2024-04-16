# -------------------------------------------
# Using PCA to Analyze Gene Differences
# -------------------------------------------

library(tidyverse)
library(gridExtra)

# import data
cpg_content <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Ananas_comosus_CpG_content_for_all_list_genes.txt", header = F, sep = "\t")
type1_list <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Ananas_comosus_CpG_content_for_high_GC_list_genes.txt",header = F, sep = "\t")
type2_list <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Ananas_comosus_CpG_content_for_low_GC_list_genes.txt",header = F, sep = "\t")

# adjust column and row names
rownames(cpg_content) <- cpg_content[,1]
cpg_content <- cpg_content[,c(2:61)]
colnames(cpg_content) <- c(1:60)

# add type numbers into data frame
type1_labels <- type1_list[,1]
type2_labels <- type2_list[,1]

type1_labels <- as.data.frame(type1_labels) %>%
  mutate(type = "type1")
type2_labels <- as.data.frame(type2_labels) %>%
  mutate(type = "type2")

# -------------------------------------------

# make pca and compile data
pca <- prcomp(cpg_content)
pca.data <- data.frame(X=pca$x[,1], Y=pca$x[,2])

# combine data with types
colnames(type1_labels) <- c("gene", "type")
colnames(type2_labels) <- c("gene", "type")
pca.data <- pca.data %>% 
  rownames_to_column("gene") %>%
  full_join(type1_labels, by = "gene") %>%
  full_join(type2_labels, by = "gene") %>%
  mutate(type = ifelse(!is.na(type.x), "type_1",
                             ifelse(!is.na(type.y), "type_2", NA)))

pca.data <- pca.data[,c(1:3,6)] # remove redundant columns
pca.data <- pca.data %>%
  filter(!is.na(type))

## make a scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) # percentage of variation that each pc accounts for

# -------------------------------------------
# Plots

# both
pca_both <- ggplot(pca.data, aes(X, Y, color = type))+
  geom_point(alpha = 0.6) +
  ggtitle("PCA of type 1 and type 2 genes in Ananas comosus")+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  scale_color_manual(values = c("#00203F", "#82c293"))+
  theme_bw()
pca_both

# using smoothScatter
smoothScatter(pca.data$X, pca.data$Y, col = ifelse(pca.data$type == "type_1", "#00203F", "#82c293"), pch = 16, cex = 0.6,
     xlab = paste("PC1 - ", pca.var.per[1], "%", sep = ""), 
     ylab = paste("PC2 - ", pca.var.per[2], "%", sep = ""),
     main = "PCA of type 1 and type 2 genes in Ananas comosus")

# type 1 only
pca_type1 <- ggplot(pca.data, aes(X, Y, color = type))+
  geom_point(data = filter(pca.data, type == "type_1"), alpha = 0.4) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA of type 1 genes in Ananas comosus")+
  scale_color_manual(values = c("#00203F"))+
  theme_bw()

# type 2 only
pca_type2 <- ggplot(pca.data, aes(X, Y, color = type))+
  geom_point(data = filter(pca.data, type == "type_2"), alpha = 0.4) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA of type 2 genes in Ananas comosus")+
  scale_color_manual(values = c("#82c293"))+
  theme_bw()

grid.arrange(pca_both, pca_type1, pca_type2, ncol=2)







# -------------------------------------------
# Removing NA types before making PCA 

all_type1_2 <- rbind(type1_labels, type2_labels)
rows_to_keep <- (rownames(cpg_content) %in% all_type1_2$gene)
cpg_content_nona <- cpg_content[rows_to_keep, ]

# make pca
pca_nona <- prcomp(cpg_content_nona)
pca_data_nona <- data.frame(X=pca_nona$x[,1], Y=pca_nona$x[,2])

# organize
pca_data_nona <- pca_data_nona %>% 
  rownames_to_column("gene") %>%
  full_join(type1_labels, by = "gene") %>%
  full_join(type2_labels, by = "gene") %>%
  mutate(type = ifelse(!is.na(type.x), "type_1",
                       ifelse(!is.na(type.y), "type_2", NA)))

pca_data_nona <- pca_data_nona[,c(1:3,6)] # remove redundant columns

#plot
pca_nona_both <- ggplot(pca_data_nona, aes(X, Y, color = type))+
  geom_point(alpha = 0.6) +
  ggtitle("PCA of type 1 and type 2 genes in Ananas comosus")+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  scale_color_manual(values = c("#00203F", "#82c293"))+
  theme_bw()

# -------------------------------------------
# Other species
# -------------------------------------------

# Carex cristatella -------------------------

# import data
cpg_content <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Carex_cristatella_CpG_content_for_all_list_genes.txt", header = F, sep = "\t")
type1_list <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Carex_cristatella_CpG_content_for_high_GC_list_genes.txt",header = F, sep = "\t")
type2_list <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Carex_cristatella_CpG_content_for_low_GC_list_genes.txt",header = F, sep = "\t")

# adjust column and row names
rownames(cpg_content) <- cpg_content[,1]
cpg_content <- cpg_content[,c(2:61)]
colnames(cpg_content) <- c(1:60)

# add type numbers into data frame
type1_labels <- type1_list[,1]
type2_labels <- type2_list[,1]

type1_labels <- as.data.frame(type1_labels) %>%
  mutate(type = "type1")
type2_labels <- as.data.frame(type2_labels) %>%
  mutate(type = "type2")

# -------------------------------------------

# make pca and compile data
pca <- prcomp(cpg_content)
pca.data <- data.frame(X=pca$x[,1], Y=pca$x[,2])

# combine data with types
colnames(type1_labels) <- c("gene", "type")
colnames(type2_labels) <- c("gene", "type")
pca.data <- pca.data %>% 
  rownames_to_column("gene") %>%
  full_join(type1_labels, by = "gene") %>%
  full_join(type2_labels, by = "gene") %>%
  mutate(type = ifelse(!is.na(type.x), "type_1",
                       ifelse(!is.na(type.y), "type_2", NA)))

pca.data <- pca.data[,c(1:3,6)] # remove redundant columns
pca.data <- pca.data %>%
  filter(!is.na(type))

## make a scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) # percentage of variation that each pc accounts for

# -------------------------------------------
# Plots

# both
pca_both <- ggplot(pca.data, aes(X, Y, color = type))+
  geom_point(alpha = 0.6) +
  ggtitle("PCA of type 1 and type 2 genes in Carex cristatella")+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  scale_color_manual(values = c("#00203F", "#82c293"))+
  theme_bw()

# type 1 only
pca_type1 <- ggplot(pca.data, aes(X, Y, color = type))+
  geom_point(data = filter(pca.data, type == "type_1"), alpha = 0.4) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA of type 1 genes in Carex cristatella")+
  scale_color_manual(values = c("#00203F"))+
  theme_bw()

# type 2 only
pca_type2 <- ggplot(pca.data, aes(X, Y, color = type))+
  geom_point(data = filter(pca.data, type == "type_2"), alpha = 0.4) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA of type 2 genes in Carex cristatella")+
  scale_color_manual(values = c("#82c293"))+
  theme_bw()

grid.arrange(pca_both, pca_type1, pca_type2, ncol=2)

# Juncus effusus -------------------------

# import data
cpg_content <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Juncus_effusus_CpG_content_for_all_list_genes.txt", header = F, sep = "\t")
type1_list <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Juncus_effusus_CpG_content_for_high_GC_list_genes.txt",header = F, sep = "\t")
type2_list <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Juncus_effusus_CpG_content_for_low_GC_list_genes.txt",header = F, sep = "\t")

# adjust column and row names
rownames(cpg_content) <- cpg_content[,1]
cpg_content <- cpg_content[,c(2:61)]
colnames(cpg_content) <- c(1:60)

# add type numbers into data frame
type1_labels <- type1_list[,1]
type2_labels <- type2_list[,1]

type1_labels <- as.data.frame(type1_labels) %>%
  mutate(type = "type1")
type2_labels <- as.data.frame(type2_labels) %>%
  mutate(type = "type2")

# -------------------------------------------

# make pca and compile data
pca <- prcomp(cpg_content)
pca.data <- data.frame(X=pca$x[,1], Y=pca$x[,2])

# combine data with types
colnames(type1_labels) <- c("gene", "type")
colnames(type2_labels) <- c("gene", "type")
pca.data <- pca.data %>% 
  rownames_to_column("gene") %>%
  full_join(type1_labels, by = "gene") %>%
  full_join(type2_labels, by = "gene") %>%
  mutate(type = ifelse(!is.na(type.x), "type_1",
                       ifelse(!is.na(type.y), "type_2", NA)))

pca.data <- pca.data[,c(1:3,6)] # remove redundant columns
pca.data <- pca.data %>%
  filter(!is.na(type))

## make a scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) # percentage of variation that each pc accounts for

# -------------------------------------------
# Plots

# both
pca_both <- ggplot(pca.data, aes(X, Y, color = type))+
  geom_point(alpha = 0.6) +
  ggtitle("PCA of type 1 and type 2 genes in Juncus effusus")+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  scale_color_manual(values = c("#00203F", "#82c293"))+
  theme_bw()

# type 1 only
pca_type1 <- ggplot(pca.data, aes(X, Y, color = type))+
  geom_point(data = filter(pca.data, type == "type_1"), alpha = 0.4) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA of type 1 genes in Juncus effusus")+
  scale_color_manual(values = c("#00203F"))+
  theme_bw()

# type 2 only
pca_type2 <- ggplot(pca.data, aes(X, Y, color = type))+
  geom_point(data = filter(pca.data, type == "type_2"), alpha = 0.4) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA of type 2 genes in Juncus effusus")+
  scale_color_manual(values = c("#82c293"))+
  theme_bw()

grid.arrange(pca_both, pca_type1, pca_type2, ncol=2)


# Musa balbisiana -------------------------

# import data
cpg_content <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Musa_balbisiana_CpG_content_for_all_list_genes.txt", header = F, sep = "\t")
type1_list <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Musa_balbisiana_CpG_content_for_high_GC_list_genes.txt",header = F, sep = "\t")
type2_list <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Musa_balbisiana_CpG_content_for_low_GC_list_genes.txt",header = F, sep = "\t")

# adjust column and row names
rownames(cpg_content) <- cpg_content[,1]
cpg_content <- cpg_content[,c(2:61)]
colnames(cpg_content) <- c(1:60)

# add type numbers into data frame
type1_labels <- type1_list[,1]
type2_labels <- type2_list[,1]

type1_labels <- as.data.frame(type1_labels) %>%
  mutate(type = "type1")
type2_labels <- as.data.frame(type2_labels) %>%
  mutate(type = "type2")

# -------------------------------------------

# make pca and compile data
pca <- prcomp(cpg_content)
pca.data <- data.frame(X=pca$x[,1], Y=pca$x[,2])

# combine data with types
colnames(type1_labels) <- c("gene", "type")
colnames(type2_labels) <- c("gene", "type")
pca.data <- pca.data %>% 
  rownames_to_column("gene") %>%
  full_join(type1_labels, by = "gene") %>%
  full_join(type2_labels, by = "gene") %>%
  mutate(type = ifelse(!is.na(type.x), "type_1",
                       ifelse(!is.na(type.y), "type_2", NA)))

pca.data <- pca.data[,c(1:3,6)] # remove redundant columns
pca.data <- pca.data %>%
  filter(!is.na(type))

## make a scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) # percentage of variation that each pc accounts for

# -------------------------------------------
# Plots

# both
pca_both <- ggplot(pca.data, aes(X, Y, color = type))+
  geom_point(alpha = 0.6) +
  ggtitle("PCA of type 1 and type 2 genes in Musa_balbisiana")+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  scale_color_manual(values = c("#00203F", "#82c293"))+
  theme_bw()

# type 1 only
pca_type1 <- ggplot(pca.data, aes(X, Y, color = type))+
  geom_point(data = filter(pca.data, type == "type_1"), alpha = 0.4) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA of type 1 genes in Musa_balbisiana")+
  scale_color_manual(values = c("#00203F"))+
  theme_bw()

# type 2 only
pca_type2 <- ggplot(pca.data, aes(X, Y, color = type))+
  geom_point(data = filter(pca.data, type == "type_2"), alpha = 0.4) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA of type 2 genes in Musa_balbisiana")+
  scale_color_manual(values = c("#82c293"))+
  theme_bw()

grid.arrange(pca_both, pca_type1, pca_type2, ncol=2)

# Oryza sative -------------------------

# import data
cpg_content <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Oryza_sativa_CpG_content_for_all_list_genes.txt", header = F, sep = "\t")
type1_list <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Oryza_sativa_CpG_content_for_high_GC_list_genes.txt",header = F, sep = "\t")
type2_list <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Oryza_sativa_CpG_content_for_low_GC_list_genes.txt",header = F, sep = "\t")

# adjust column and row names
rownames(cpg_content) <- cpg_content[,1]
cpg_content <- cpg_content[,c(2:61)]
colnames(cpg_content) <- c(1:60)

# add type numbers into data frame
type1_labels <- type1_list[,1]
type2_labels <- type2_list[,1]

type1_labels <- as.data.frame(type1_labels) %>%
  mutate(type = "type1")
type2_labels <- as.data.frame(type2_labels) %>%
  mutate(type = "type2")

# -------------------------------------------

# make pca and compile data
pca <- prcomp(cpg_content)
pca.data <- data.frame(X=pca$x[,1], Y=pca$x[,2])

# combine data with types
colnames(type1_labels) <- c("gene", "type")
colnames(type2_labels) <- c("gene", "type")
pca.data <- pca.data %>% 
  rownames_to_column("gene") %>%
  full_join(type1_labels, by = "gene") %>%
  full_join(type2_labels, by = "gene") %>%
  mutate(type = ifelse(!is.na(type.x), "type_1",
                       ifelse(!is.na(type.y), "type_2", NA)))

pca.data <- pca.data[,c(1:3,6)] # remove redundant columns
pca.data <- pca.data %>%
  filter(!is.na(type))

## make a scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) # percentage of variation that each pc accounts for

# -------------------------------------------
# Plots

# both
pca_both <- ggplot(pca.data, aes(X, Y, color = type))+
  geom_point(alpha = 0.6) +
  ggtitle("PCA of type 1 and type 2 genes in Oryza_sativa")+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  scale_color_manual(values = c("#00203F", "#82c293"))+
  theme_bw()

# type 1 only
pca_type1 <- ggplot(pca.data, aes(X, Y, color = type))+
  geom_point(data = filter(pca.data, type == "type_1"), alpha = 0.4) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA of type 1 genes in Oryza_sativa")+
  scale_color_manual(values = c("#00203F"))+
  theme_bw()

# type 2 only
pca_type2 <- ggplot(pca.data, aes(X, Y, color = type))+
  geom_point(data = filter(pca.data, type == "type_2"), alpha = 0.4) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA of type 2 genes in Oryza_sativa")+
  scale_color_manual(values = c("#82c293"))+
  theme_bw()

grid.arrange(pca_both, pca_type1, pca_type2, ncol=2)

