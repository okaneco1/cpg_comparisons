##################### General Script for Command Line ##########################

#!/usr/bin/env Rscript

# Conor O'Kane
# Michigan State University
# 08.28.23
# calculating euclidean distance between type 1 and 2 genes

# Define input files from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
full_cpg_content <- args[1]  # CpG content file for all genes of desired species
cpg_avg_type1 <- args[2]  # CpG average across type 1 genes of desired species
cpg_avg_type2 <- args[3]  # CpG average across type 2 genes of desired species
type1_list_labels <- args[4]  # list of type 1 (low GC) genes for label assignment
type2_list_labels <- args[5]  # list of type 2 (high GC) genes for label assignment

suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"
suppressPackageStartupMessages(require(tidyverse))

# Create list of options on whether to use type 
#option_list = list(
#make_option(c("-a"), action="store_true",help="Use cpg_avg_type1"),
#make_option(c("-b"), action="store_true", help="Use cpg_avg_type2")
#)

# Read in args
full <- read.delim(full_cpg_content, header = F, sep = "\t")
type1_avg <- read.delim(cpg_avg_type1, header = F, sep = "\t")
type2_avg <- read.delim(cpg_avg_type2, header = F, sep = "\t")
type1_list <- read.delim(highGC_list_labels, header = F, sep = "\t")
type2_list <- read.delim(lowGC_list_labels, header = F, sep = "\t")

# Clean and organize data for analysis
full_labels <- full[,1]  # extract the gene labels
full_labels <- as.data.frame(full_labels) 
full <- full[,2:61]  # trim data from gene labels

type1_avg <- type1_avg[,1:60]  # removes empty final line
type2_avg <- type2_avg[,1:60]

type1_list <- type1_list[,1]  # gets the labels for type 1 genes
type1_list <- as.data.frame(type1_list) %>%
  mutate(GC_content = "type1")  # creates data frame indicating those genes as type 1

type2_list <- type2_list[,1]  # gets the labels for type 2 genes
type2_list <- as.data.frame(type2_list) %>%
  mutate(GC_content = "type2")

colnames(full_labels)[1] <- "label"  # align column names to merge
colnames(type1_list)[1] <- "label"
colnames(type2_list)[1] <- "label"

full_labels_df <- full_labels %>%
  left_join(type1_list, by = "label") %>%
  left_join(type2_list, by = "label") %>%
  # Combine the data from both type 1 and type 2 GC content columns
  mutate(GC_content = ifelse(!is.na(GC_content.x), "type 1",
                             ifelse(!is.na(GC_content.y), "type 2", NA)))
full_labels_df <- full_labels_df[, c(1,4)]  # trim out redundant columns

# Euclidean distance function
e_dist <- function(cpg_content, cpg_avg) {
  # cpg_content is the CpG content for a species for all genes
  # cpg_avg is the aveage CpG content per gene for either type 1 or 2
  
  # check to make sure column number is the same
  if (!identical(ncol(cpg_content), ncol(cpg_avg))) {
    stop("Must have the same number of columns.")
  }
  
  # assign the output variable to store distances
  distances <- numeric(nrow(cpg_content))
  
  # calculate euclidean distance on a per-row basis, storing each value in "distances"
  for (i in 1:nrow(cpg_content)) {
    squared_dif <- (cpg_content[i, ] - cpg_avg1[1, ])^2
    sum_squared_dif <- sum(squared_dif)
    distances[i] <- sqrt(sum_squared_dif) 
  }
  
  # output
  return(distances)
}

# Generate distances for type 1 and type 2 genes
type1_dist <- e_dist(full_cpg_content, cpg_avg_type1)
type2_dist <- e_dist(full_cpg_content, cpg_avg_type2)

# Combine distances in data frame
dist_df <- data.frame(type1_dist, type2_dist)

# Combine labels with distances and set column names
dist_df <- cbind(dist_df, full_labels_df)
dist_df <- dist_df[, c(3, 1:2, 6)]
colnames(dist_df)[2:3] <- c("type1_distance", "type2_distance")

#---- Plot

# Filter out rows with NA
filtered_dist_df <- dist_df %>%
  filter(!is.na(GC_content))

# Plot without NA values
dist_plot <- ggplot(filtered_dist_df, aes(type1_distance, type2_distance, color = GC_content))+
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "#706f6f") +
  scale_color_manual(values = c("#ffc240", "#56B4E9"))+
  theme_light()

# Save plot in current working directory
ggsave("euclidean_distance_type1_type2_genes_plot.jpg", gg, width = 8, height = 6, dpi = 300)

