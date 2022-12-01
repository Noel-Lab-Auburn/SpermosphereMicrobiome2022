################Network############

###### Libraries #####
library(igraph)
library(phyloseq)
library(vegan)
library(tidyverse)
library(metagenomeSeq)
library(ggplot2)
library(ggpubr)
library(NetCoMi)
library(SPRING)

# set options for scientific numbers to not be displayed
options(scipen=10000) 

# color blind pallet used throughout 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#### Read in RDS file #### 
# using the non-normalized reads since spieceasi has its own normalizaiton methods
bac_sperm <- readRDS(file = "Bacteria/Bacteria_spermosphere_nonnorm_112922.rds")
fungi_sperm <- readRDS(file = "Fungi/Fungi_spermosphere_unedited_083022.rds")

# Network 1 with cross domain interactions: soybean

# Filter the datasets 

# have to remove the samples that are not the same between the two datasets
network.soybean.bacteria <- bac_sperm %>%  
  subset_samples(Crop == "Soybean" & Fungal_Code != "SB.6.5") %>%
  filter_taxa(function(x) sum(x > 100) > (0.10*length(x)), TRUE)

network.soybean.fungi <- fungi_sperm %>%  
  subset_samples(Crop == "Soybean") %>%
  filter_taxa(function(x) sum(x > 100) > (0.10*length(x)), TRUE)

# Store count matrices (taxa are columns)
counts_network.soybean.bacteria <- as.matrix(t(phyloseq::otu_table(network.soybean.bacteria)@.Data))
counts_network.soybean.fungi <- as.matrix(t(phyloseq::otu_table(network.soybean.fungi)@.Data))

# setting the rownames of the count matrices to the same codes between the two sample sources for cross-domain interactions in spieceasi
rownames(counts_network.soybean.bacteria) <- network.soybean.bacteria@sam_data$Fungal_Code

# Sanity check - has to be true to proceed. 
all.equal(rownames(counts_network.soybean.bacteria), rownames(counts_network.soybean.fungi))

set.seed(123456)

# Run SpiecEasi and create association matrix for group 1
spiec_result_gr1 <- multi.spiec.easi(list(counts_network.soybean.bacteria, counts_network.soybean.fungi), 
                                     method='mb', nlambda=40, 
                                     lambda.min.ratio=1e-2, 
                                     pulsar.params = list(thresh = 0.05))

assoMat1 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiec_result_gr1), mode = "ave")

assoMat1 <- as.matrix(assoMat1)

# Network 2 with cross domain interactions: soybean

network.cotton.bacteria <- bac_sperm %>%  
  subset_samples(Crop == "Cotton "  & Fungal_Code != "C.12.6") %>%
  filter_taxa(function(x) sum(x > 100) > (0.10*length(x)), TRUE)

network.cotton.fungi <- fungi_sperm %>%  
  subset_samples(Crop == "Cotton " & Code != "C.18.6") %>%
  filter_taxa(function(x) sum(x > 100) > (0.10*length(x)), TRUE)

# Store count matrices (taxa are columns)
counts_network.cotton.bacteria <- as.matrix(t(phyloseq::otu_table(network.cotton.bacteria)@.Data))
counts_network.cotton.fungi <- as.matrix(t(phyloseq::otu_table(network.cotton.fungi)@.Data))

# setting the rownames of the count matrices to the same codes between the two sample sources for cross-domain interactions in spieceasi
rownames(counts_network.cotton.bacteria) <- network.cotton.bacteria@sam_data$Fungal_Code

# Sanity check - has to be true to proceed. 
all.equal(rownames(counts_network.cotton.bacteria), rownames(counts_network.cotton.fungi))

# Run SpiecEasi and create association matrix for group 1
spiec_result_gr2 <- multi.spiec.easi(list(counts_network.cotton.bacteria, counts_network.cotton.fungi), 
                                     method='mb', nlambda=40, 
                                     lambda.min.ratio=1e-2, 
                                     pulsar.params = list(thresh = 0.05))

assoMat2 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiec_result_gr2), mode = "ave")

assoMat2 <- as.matrix(assoMat1)

#### Not sure how to make sure this is correct yet.... 
# Get taxa names
taxnames <- c(taxa_names(hmp216S), taxa_names(hmp2prot))

colnames(assoMat1) <- rownames(assoMat1) <- taxnames


colnames(assoMat2) <- rownames(assoMat2) <- taxnames

# setting the diagnol to 1 because its between the same taxa. 
diag(assoMat1) <- 1
diag(assoMat2) <- 1

# Network construction (pass association matrices to netConstruct)
# - sparsMethod must be set to "none" because sparsification is already included in SpiecEasi
sperm.crossdomain.cottonvsoybean <- netConstruct(data = assoMat1, data2 = assoMat2, 
                                 dataType = "condDependence", 
                                 sparsMethod = "none")

# Network analysis
netprops_sperm.crossdomain.cottonvsoybean <- netAnalyze(sperm.crossdomain.cottonvsoybean, hubPar = "eigenvector")



