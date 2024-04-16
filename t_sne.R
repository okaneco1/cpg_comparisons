#------------------------------------------------
# t-SNE plots

library(tidyverse)
library(Rtsne)

#------------------------------------------------
# data and cleaning

# import data (Musa_balbisiana)
cpg_content <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Musa_balbisiana_CpG_content_for_all_list_genes.txt", header = F, sep = "\t")
type1_list <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Musa_balbisiana_CpG_content_for_high_GC_list_genes.txt",header = F, sep = "\t")
type2_list <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Musa_balbisiana_CpG_content_for_low_GC_list_genes.txt",header = F, sep = "\t")

# import data (Juncus_effusus)
cpg_content <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Juncus_effusus_CpG_content_for_all_list_genes.txt", header = F, sep = "\t")
type1_list <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Juncus_effusus_CpG_content_for_high_GC_list_genes.txt",header = F, sep = "\t")
type2_list <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Juncus_effusus_CpG_content_for_low_GC_list_genes.txt",header = F, sep = "\t")

# adjust column and row names
rownames(cpg_content) <- cpg_content[,1]
cpg_content <- cpg_content[,c(2:61)]
colnames(cpg_content) <- paste("B", 1:60, sep = "")

# create lists with labels
type1_labels <- type1_list[,1]
type2_labels <- type2_list[,1]


#------------------------------------------------
# t-SNE setup

# assigning type labels (or NA for neither)
cpg_content <- cpg_content %>%
  mutate(type = case_when(
    rownames(cpg_content) %in% type1_labels ~ "type_1",
    rownames(cpg_content) %in% type2_labels ~ "type_2",
    TRUE ~ NA_character_
  )) 

# sample of 1000
small_subset <- cpg_content %>%
  filter(!is.na(type)) %>%
  sample_n(1000)  # samples a number of rows from a table



#------------------------------------------------
# run t-SNE 

tsne_result <- Rtsne(small_subset, dims = 2, perplexity = 30, verbose = TRUE)
tsne_df <- data.frame(X = tsne_result$Y[, 1], Y = tsne_result$Y[, 2])

# add the type column back
tsne_df$type <- small_subset$type

# (save as specific species)
plot_mb <- ggplot(tsne_df, aes(x = X, y = Y, color = type)) +
  geom_point(size = 1.5) +
  labs(x = "Dimension 1", y = "Dimension 2") +
  ggtitle("t-SNE of type 1 and type 2 genes in Musa balbisiana")+
  theme_minimal()


#------------------------------------------------
# putting plots together

grid.arrange(plot_ac, plot_cc, plot_je, plot_mb, plot_os, ncol = 3) 


#------------------------------------------------
# looking at some statistics

type1_hist <- filter(cpg_content, type == "type_1")
type2_hist <- filter(cpg_content, type == "type_2")

ggplot(cpg_content, aes(x = type, fill = type)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.5) +
  labs(title = "Histogram of type_1 and type_2", x = "X", y = "Frequency") +
  scale_fill_manual(values = c("type_1" = "blue", "type_2" = "red"))

hist(cpg_content$type)
