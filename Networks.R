################Network############

###### Libraries #####

library(igraph)
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(SpiecEasi)

#### Other packages needed, but installed on the HPC #### 
# NetCoMi, SPRING

# set options for scientific numbers to not be displayed
options(scipen=10000) 

# color blind pallet used throughout 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#### Read in RDS file #### 
# using the non-normalized reads since spieceasi has its own normalizaiton methods
bac_sperm <- readRDS(file = "Bacteria_spermosphere_nonnorm_112922.rds")
fungi_sperm <- readRDS(file = "Fungi_spermosphere_unedited_083022.rds")

# Filter the datasets 

# have to remove the samples that are not the same between the two datasets
network.bacteria <- bac_sperm %>%  
  subset_samples(Crop %in% c("Cotton ", "Soybean") & !Fungal_Code %in% c("SB.6.5", "C.12.6")) %>%
  filter_taxa(function(x) sum(x > 100) > (0.10*length(x)), TRUE)

network.fungi <- fungi_sperm %>%  
  subset_samples(Crop %in% c("Cotton ", "Soybean") & Code != "C.18.6") %>%
  filter_taxa(function(x) sum(x > 100) > (0.10*length(x)), TRUE)

# Store count matrices (taxa are columns)
counts_network.bacteria <- as.matrix(t(phyloseq::otu_table(network.bacteria)@.Data))
counts_network.fungi <- as.matrix(t(phyloseq::otu_table(network.fungi)@.Data))

                                            # Print data frame subset


# setting the rownames of the count matrices to the same codes between the two sample sources for cross-domain interactions in spieceasi
rownames(counts_network.bacteria) <- network.bacteria@sam_data$Fungal_Code

# Sanity check - has to be true to proceed. 
all.equal(rownames(counts_network.bacteria), rownames(counts_network.fungi))

# Subsetting based on cotton versus soybean

# Bacteria
rows_keep_cotton_bacteria <- grep("C", network.bacteria@sam_data$Fungal_Code) 
data_cotton_bacteria <- counts_network.bacteria[rows_keep_cotton_bacteria, ]  # Extract rows from data
data_cotton_bacteria    
dim(data_cotton_bacteria)

rows_keep_soybean_bacteria <- grep("S", network.bacteria@sam_data$Fungal_Code) 
data_soybean_bacteria <- counts_network.bacteria[rows_keep_soybean_bacteria, ]  # Extract rows from data
data_soybean_bacteria  
dim(data_soybean_bacteria)
# Fungi
rows_keep_cotton_fungi <- grep("C", network.fungi@sam_data$Code) 
data_cotton_fungi <- counts_network.fungi[rows_keep_cotton_fungi, ]  # Extract rows from data
data_cotton_fungi    
dim(data_cotton_fungi)

rows_keep_soybean_fungi <- grep("S", network.fungi@sam_data$Code) 
data_soybean_fungi <- counts_network.fungi[rows_keep_soybean_fungi, ]  # Extract rows from data
data_soybean_fungi 
dim(data_soybean_fungi)

set.seed(123456)

# Run SpiecEasi and create association matrix for group 1
spiec_result_gr1 <- multi.spiec.easi(list(data_cotton_bacteria, data_cotton_fungi), 
                                     method='mb', nlambda=40, 
                                     lambda.min.ratio=1e-2, 
                                     pulsar.params = list(thresh = 0.05))


assoMat1 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiec_result_gr1), mode = "ave")

# Run SpiecEasi and create association matrix for group 2
spiec_result_gr2 <- multi.spiec.easi(list(data_soybean_bacteria, data_soybean_fungi), 
                                     method='mb', nlambda=40, 
                                     lambda.min.ratio=1e-2, 
                                     pulsar.params = list(thresh = 0.05))

assoMat2 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiec_result_gr2), mode = "ave")

taxnames <- c(taxa_names(network.bacteria), taxa_names(network.fungi))

assoMat1 <- as.matrix(assoMat1)
colnames(assoMat1) <- rownames(assoMat1) <- taxnames
diag(assoMat1) <- 1

colnames(assoMat2) <- rownames(assoMat2) <- taxnames
assoMat2 <- as.matrix(assoMat2)
diag(assoMat2) <- 1

# Network construction (pass association matrices to netConstruct)
# - sparsMethod must be set to "none" because sparsification is already included in SpiecEasi
sperm.crossdomain.cottonvsoybean <- netConstruct(data = assoMat1, data2 = assoMat2, 
                                                 dataType = "condDependence", 
                                                 sparsMethod = "none")

# Network analysis
netprops_sperm.crossdomain.cottonvsoybean <- netAnalyze(sperm.crossdomain.cottonvsoybean, hubPar = "eigenvector")


# Network comparison
# - Permutation tests cannot be performed because the association matrices are
#   used for network construction. For permutation tests, however, the count 
#   data are needed.
netcomp_soybean_cotton <- netCompare(netprops_sperm.crossdomain.cottonvsoybean, permTest = FALSE)

summary(netcomp_soybean_cotton, groupNames = c("Soybean", "Cotton"))

saveRDS(sperm.crossdomain.cottonvsoybean, "netConstruct.rds")
saveRDS(netprops_sperm.crossdomain.cottonvsoybean, "netanalyse.rds")
saveRDS(netcomp_soybean_cotton, "netCompare.rds")

# Generate network plots
nodeCols <- c(rep(cbbPalette[[2]], ntaxa(network.bacteria)), rep(cbbPalette[[8]], ntaxa(network.fungi)))
names(nodeCols) <- taxnames

NetCoMi::plot(net.analyse, 
     sameLayout = TRUE, 
     layoutGroup = "union",
     nodeColor = "colorVec", 
     colorVec = nodeCols,
     nodeSize = "eigen", 
     nodeSizeSpread = 2,
     labelScale = FALSE,
     cexNodes = 2, 
     cexLabels = 2,
     cexHubLabels = 2.5,
     cexTitle = 3.8,
     groupNames = c("Soybean", "Cotton"))

legend(-0.2, 1.2, cex = 3, pt.cex = 4, 
       legend = c("Soybean Spermosphere", "Cotton Spermosphere"), col = c(cbbPalette[[2]], cbbPalette[[8]]), 
       bty = "n", pch = 16) 

tax.bac <- data.frame(network.bacteria@tax_table)
soybean.hubs <- tax.bac %>%
  subset(OTU %in% net.analyse$hubs$hubs1)
soybean.hubs
cotton.hubs <- tax.bac %>%
  subset(OTU %in% net.analyse$hubs$hubs2)
cotton.hubs


