# ------------------------------------------
# MDS plots for gene differences
# ------------------------------------------

library(tidyverse)
library(gridExtra)

# import data
cpg_content <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Ananas_comosus_CpG_content_for_all_list_genes.txt", header = F, sep = "\t")
type1_list <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Ananas_comosus_CpG_content_for_high_GC_list_genes.txt",header = F, sep = "\t")
type2_list <- read.delim("/Users/conorokane/Desktop/MSU/GC\ Project/distance/Ananas_comosus_CpG_content_for_low_GC_list_genes.txt",header = F, sep = "\t")

# adjust column and row names
rownames(cpg_content) <- cpg_content[,1]
cpg_content <- cpg_content[,c(2:61)]
colnames(cpg_content) <- paste("B", 1:60, sep = "")


# add type numbers into data frame
type1_labels <- type1_list[,1]
type2_labels <- type2_list[,1]

type1_labels <- as.data.frame(type1_labels) %>%
  mutate(type = "type1")
type2_labels <- as.data.frame(type2_labels) %>%
  mutate(type = "type2")


## first, calculate the distance matrix using the Euclidian distance.
random_subset <- cpg_content[sample(nrow(cpg_content), 1000), ]
distance.matrix <- dist(scale(random_subset), method="maximum")


## do the MDS math (this is basically eigen value decomposition)
mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)
mds.var.per

## now make a fancy looking plot that shows the MDS axes and the variation:
mds.values <- mds.stuff$points
mds.data <- data.frame(X=mds.values[,1],
                       Y=mds.values[,2])
mds.data

# combine data with types
colnames(type1_labels) <- c("gene", "type")
colnames(type2_labels) <- c("gene", "type")
mds.data <- mds.data %>% 
  rownames_to_column("gene") %>%
  full_join(type1_labels, by = "gene") %>%
  full_join(type2_labels, by = "gene") %>%
  mutate(type = ifelse(!is.na(type.x), "type_1",
                       ifelse(!is.na(type.y), "type_2", NA)))

mds.data <- mds.data[,c(1:3,6)] # remove redundant columns

mds.data <- mds.data %>%
  filter(!is.na(type))

# both
ggplot(mds.data, aes(X, Y, color = type))+
  geom_point(alpha = 0.6) +
  ggtitle("MDS Plot (Euclidean distance) of type 1 and type 2 genes in Ananas comosus")+
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
  scale_color_manual(values = c("#00203F", "#82c293"))+
  theme_bw()





# ---------------------
# some data exploration
#----------------------
small_subset <- cpg_content[sample(nrow(cpg_content), 1000), ]

ggplot(small_subset, aes(B53, B24))+
  geom_point()



# ---------------------
# NMDS through vegan
#----------------------
install.packages("vegan")
library("vegan")

#########################
#subset the dataframe on which to base the ordination (dataframe 1)
data_1 <- cpg_content[sample(nrow(cpg_content), 1000), ]

#Identify the columns that contains the descriptive/environmental data (dataframe 2)
# combine data with types
colnames(type1_labels) <- c("gene", "type")
colnames(type2_labels) <- c("gene", "type")
data_2_all <- data_1 %>% 
  rownames_to_column("gene") %>%
  full_join(type1_labels, by = "gene") %>%
  full_join(type2_labels, by = "gene") %>%
  mutate(type = ifelse(!is.na(type.x), "type_1",
                       ifelse(!is.na(type.y), "type_2", NA)))

data_2 <- data_2_all[63:64]# remove redundant columns
#pca.data <- pca.data %>%
  #filter(!is.na(type))


#ordination by NMDS
NMDS <- metaMDS(data_1, distance = "bray", k = 2)



#-----------------------------
# Data visualisation via ggplot

install.packages("ggplot2")
library("ggplot2")

#Extract the axes scores
datascores <- as.data.frame(scores(NMDS))  #extract the site scores

#Add/calculate spider diagram
scores <- cbind(as.data.frame(datascores), gene_type = data_2$type)
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ type, data = scores, FUN = mean)
seg <- merge(scores, setNames(centroids, c('gene type','oNMDS1','oNMDS2')),
             by = 'Habitat', sort = FALSE)

#plot
ggplot(scores, aes(x = NMDS1, y = NMDS2, colour = type)) +
  geom_segment(data = seg,
               mapping = aes(xend = oNMDS1, yend = oNMDS2)) + # add spiders
  geom_point(data = centroids, size = 4) +                    # add centroids
  geom_point() +                                              
  coord_fixed()+                                              
  theme_bw()+ 
  theme(legend.position="right",legend.text=element_text(size=10),legend.direction='vertical')


#----------------------------------------------------------------
# t-SNE

library(Rtsne)

# first setting up the data with types labelled, then removing the NA types before the subsample
cpg_content_no_na <- cpg_content %>%
  mutate(type = case_when(
    rownames(cpg_content) %in% type1_labels ~ "type_1",
    rownames(cpg_content) %in% type2_labels ~ "type_2",
    TRUE ~ NA_character_
  ))

small_subset <- cpg_content[sample(nrow(!is.na(cpg_content)), 1000), ]

# removing duplicates
duplicates <- duplicated(small_subset)
clean_subset <- small_subset[!duplicates, ]

# add a type column for the subsample
clean_subset <- mutate(clean_subset, type = NA)
for (i in 1:nrow(clean_subset)) {
  if (rownames(clean_subset[i, ]) %in% type1_labels) {
    clean_subset[i, "type"] = "type_1"
  } else if (rownames(clean_subset[i, ]) %in% type2_labels){
    clean_subset[i, "type"] = "type_2"
  } else {
    clean_subset[i, "type"] = NA
  } 
}

"Aco004101.1.v3 " %in% type2_labels

# run PCA for tsne
tsne_result <- Rtsne(clean_subset, dims = 2, perplexity = 30, verbose = TRUE)
tsne_df <- data.frame(X = tsne_result$Y[, 1], Y = tsne_result$Y[, 2])

ggplot(tsne_df, aes(x = X, y = Y)) +
  geom_point(size = 3, color = "blue") +
  labs(x = "Dimension 1", y = "Dimension 2") +
  theme_minimal()

