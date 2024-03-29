################Transmission Analysis############

###### Libraries #####
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(microbiome)
library(ggVennDiagram)
library(VennDiagram)

##### Set global options #####

# no scientific notation
options(scipen=10000) 

# color blind pallet used throughout 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ibm.cbb <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")
tol.cbb <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")

### Transmission analysis ##### 

# The purpose of this analysis is to determine which OTUs are only present in the spermosphere samples compared to bulk soil. Which may indicate that they came from the seed.
# We will transform the dataset to presence absence, then compared what was preseent in only spermospheres. 

##### Bacteria ####
#### Soybean #####
soybean.sub <- bac_sperm %>% 
  subset_samples(Crop == "Soybean")
otu.soybean <- soybean.sub@otu_table %>%
  as("matrix")
map.soybean <- soybean.sub@sam_data %>%
  as("data.frame")
otu_PA_soybean <- 1*((otu.soybean>0)==1)   
otu_occ_soybean <- data.frame(rowSums(otu_PA_soybean)); colnames(otu_occ_soybean) <- "Abundance"

otu_occ_soybean$OTU <- rownames(otu_occ_soybean)
otu_occ_soybean$Crop <- "Soybean"

##### Cotton ##### 
cotton.sub <- bac_sperm %>% 
  subset_samples(Crop == "Cotton ")
otu.cotton <- cotton.sub@otu_table %>%
  as("matrix")
map.cotton <- cotton.sub@sam_data %>%
  as("data.frame")
otu_PA_cotton <- 1*((otu.cotton>0)==1)   
otu_occ_cotton <- data.frame(rowSums(otu_PA_cotton)); colnames(otu_occ_cotton) <- "Abundance"

otu_occ_cotton$OTU <- rownames(otu_occ_cotton)
otu_occ_cotton$Crop <- "Cotton"

##### Bulk Soil #####
bulk.sub <- bac_sperm %>% 
  subset_samples(Crop == "Bulk Soil")
otu.bulk <- bulk.sub@otu_table %>%
  as("matrix")
map.bulk <- bulk.sub@sam_data %>%
  as("data.frame")
otu_PA_bulk <- 1*((otu.bulk>0)==1)   
otu_occ_bulk <- data.frame(rowSums(otu_PA_bulk)); colnames(otu_occ_bulk) <- "Abundance"

otu_occ_bulk$OTU <- rownames(otu_occ_bulk)
otu_occ_bulk$Crop <- "Bulk Soil"

#### Combine ######
combined.occ <- rbind.data.frame(otu_occ_bulk, otu_occ_cotton, otu_occ_soybean)
combined.occ$pres.abs <- ifelse(combined.occ$Abundance > 0, 1, 0)
combined.occ$otu.pres <- interaction(combined.occ$Crop, combined.occ$pres.abs)

BulkSoil.absent <- combined.occ$OTU[combined.occ$otu.pres == "Bulk Soil.0"]
Soybean.absent <- combined.occ$OTU[combined.occ$otu.pres == "Soybean.0"]
Cotton.absent <- combined.occ$OTU[combined.occ$otu.pres == "Cotton.0"]

BulkSoil.present <- combined.occ$OTU[combined.occ$otu.pres == "Bulk Soil.1"]
Soybean.present <- combined.occ$OTU[combined.occ$otu.pres == "Soybean.1"]
Cotton.present <- combined.occ$OTU[combined.occ$otu.pres == "Cotton.1"]

x <- list(
  "Cotton" = Cotton.present,
  "Soybean" = Soybean.present,
  "Bulk Soil" = BulkSoil.present
)
# ggVennDiagram return a ggplot object, the fill/edge colors can be further modified with ggplot functions.

ven <- ggVennDiagram(x, color = 1, lwd = 0.8, lty = 1) +
  #scale_fill_gradient(low="blue",high = "red") +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  scale_color_manual(values = c("black", "black", "black"))
ven

overlap <- calculate.overlap(
  x = list(
    "Cotton" = Cotton.present,
    "Soybean" = Soybean.present,
    "Bulk Soil" = BulkSoil.present));

length(overlap[[1]]) # intersect of all
length(overlap[[2]]) # coton and soybean spermosphere 
length(overlap[[3]]) # bulk soil and cotton
length(overlap[[4]]) # soybean and bulk soil
length(overlap[[5]]) # Cotton unique
length(overlap[[6]]) # soybean unique
length(overlap[[7]]) # Bulk soil unique

all <- bac_sperm@tax_table %>%
  as.matrix() %>%
  as.data.frame() %>%
  subset(OTU %in% overlap[[1]]) %>%
  mutate(Crop = "All")

cotton.soybean <- bac_sperm@tax_table %>%
  as.matrix() %>%
  as.data.frame() %>%
  subset(OTU %in% overlap[[2]]) %>%
  mutate(Crop = "Cotton and Soybean")

cotton.bulk <- bac_sperm@tax_table %>%
  as.matrix() %>%
  as.data.frame() %>%
  subset(OTU %in% overlap[[3]]) %>%
  mutate(Crop = "Cotton and Bulk Soil")

soybean.bulk <- bac_sperm@tax_table %>%
  as.matrix() %>%
  as.data.frame() %>%
  subset(OTU %in% overlap[[4]]) %>%
  mutate(Crop = "Soybean and Bulk Soil")

cotton.unique <- bac_sperm@tax_table %>%
  as.matrix() %>%
  as.data.frame() %>%
  subset(OTU %in% overlap[[5]]) %>%
  mutate(Crop = "Cotton")

soybean.unique <- bac_sperm@tax_table %>%
  as.matrix() %>%
  as.data.frame() %>%
  subset(OTU %in% overlap[[6]]) %>%
  mutate(Crop = "Soybean")

bulk.unique <- bac_sperm@tax_table %>%
  as.matrix() %>%
  as.data.frame() %>%
  subset(OTU %in% overlap[[7]]) %>%
  mutate(Crop = "Bulk Soil")

combined.unique <- rbind.data.frame(soybean.unique, cotton.unique, bulk.unique)
combined.all <- rbind.data.frame(soybean.unique, 
                                 cotton.unique, 
                                 bulk.unique, 
                                 all,
                                 cotton.soybean,
                                 cotton.bulk,
                                 soybean.bulk)
hubs <- c("BOTU_49", "BOTU_11", "BOTU_1559", "BOTU_132", "BOTU_36",
          "BOTU_12", "BOTU_29", "BOTU_46", "BOTU_119", "BOTU_1009", "BOTU_19")

subset(combined.all, OTU %in% hubs)


combined.unique %>%
  group_by(Crop, Phylum) %>%
  count() %>%
  arrange(-n) %>%
  mutate(Genus.label = ifelse(n > 5, Phylum, "Other")) %>%
  ggplot(aes(x = Crop, y = n, fill = Genus.label)) +
  geom_bar(position="fill", stat = "identity") + 
  theme_classic() + 
  xlab("")+
  ylab("Proportion")+
  scale_fill_manual(values=c(cbbPalette, ibm.cbb, tol.cbb, "grey", "green", "blue", "orange", "purple"))

combined.unique %>%
  subset(Phylum == "Firmicutes" & Crop == "Soybean") %>%
  group_by(Crop, Label) %>%
  count() %>%
  arrange(-n)
  


soy.cot[soy.cot$Genus == "Paenibacillus",]
soybean.unique[soybean.unique$Genus == "Paenibacillus",]
cotton.unique[cotton.unique$Genus == "Paenibacillus",]




##### Fungi ####
fungi.filt <- core(fungi_sperm, detection = 0.01, prevalence = 0.1)
#### Soybean #####
soybean.sub <- fungi.filt %>% 
  subset_samples(Crop == "Soybean")
otu.soybean <- soybean.sub@otu_table %>%
  as("matrix")
map.soybean <- soybean.sub@sam_data %>%
  as("data.frame")
otu_PA_soybean <- 1*((otu.soybean>0)==1)   
otu_occ_soybean <- data.frame(rowSums(otu_PA_soybean)); colnames(otu_occ_soybean) <- "Abundance"

otu_occ_soybean$OTU <- rownames(otu_occ_soybean)
otu_occ_soybean$Crop <- "Soybean"

##### Cotton ##### 
cotton.sub <- fungi.filt %>% 
  subset_samples(Crop == "Cotton ")
otu.cotton <- cotton.sub@otu_table %>%
  as("matrix")
map.cotton <- cotton.sub@sam_data %>%
  as("data.frame")
otu_PA_cotton <- 1*((otu.cotton>0)==1)   
otu_occ_cotton <- data.frame(rowSums(otu_PA_cotton)); colnames(otu_occ_cotton) <- "Abundance"

otu_occ_cotton$OTU <- rownames(otu_occ_cotton)
otu_occ_cotton$Crop <- "Cotton"

##### Bulk Soil #####
bulk.sub <- fungi.filt %>% 
  subset_samples(Crop == "Bulk Soil")
otu.bulk <- bulk.sub@otu_table %>%
  as("matrix")
map.bulk <- bulk.sub@sam_data %>%
  as("data.frame")
otu_PA_bulk <- 1*((otu.bulk>0)==1)   
otu_occ_bulk <- data.frame(rowSums(otu_PA_bulk)); colnames(otu_occ_bulk) <- "Abundance"

otu_occ_bulk$OTU <- rownames(otu_occ_bulk)
otu_occ_bulk$Crop <- "Bulk Soil"

#### Combine ######
combined.occ <- rbind.data.frame(otu_occ_bulk, otu_occ_cotton, otu_occ_soybean)
combined.occ$pres.abs <- ifelse(combined.occ$Abundance > 5, 1, 0)
combined.occ$otu.pres <- interaction(combined.occ$Crop, combined.occ$pres.abs)

BulkSoil.absent <- combined.occ$OTU[combined.occ$otu.pres == "Bulk Soil.0"]
Soybean.absent <- combined.occ$OTU[combined.occ$otu.pres == "Soybean.0"]
Cotton.absent <- combined.occ$OTU[combined.occ$otu.pres == "Cotton.0"]

BulkSoil.present <- combined.occ$OTU[combined.occ$otu.pres == "Bulk Soil.1"]
Soybean.present <- combined.occ$OTU[combined.occ$otu.pres == "Soybean.1"]
Cotton.present <- combined.occ$OTU[combined.occ$otu.pres == "Cotton.1"]

x <- list(
  "Cotton" = Cotton.present,
  "Soybean" = Soybean.present,
  "Bulk Soil" = BulkSoil.present
)
# ggVennDiagram return a ggplot object, the fill/edge colors can be further modified with ggplot functions.

ven <- ggVennDiagram(x, color = 1, lwd = 0.8, lty = 1) +
  #scale_fill_gradient(low="blue",high = "red") +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  scale_color_manual(values = c("black", "black", "black"))
ven

overlap <- calculate.overlap(
  x = list(
    "Cotton" = Cotton.present,
    "Soybean" = Soybean.present,
    "Bulk Soil" = BulkSoil.present));

length(overlap[[1]]) # intersect of all
length(overlap[[2]]) # coton and soybean spermosphere 
length(overlap[[3]]) # bulk soil and cotton
length(overlap[[4]]) # soybean and bulk soil
length(overlap[[5]]) # Cotton unique
length(overlap[[6]]) # soybean unique
length(overlap[[7]]) # Bulk soil unique

soybean.unique <- fungi_sperm@tax_table %>%
  as.matrix() %>%
  as.data.frame() %>%
  subset(OTU %in% overlap[[5]]) %>%
  mutate(Crop = "Soybean")

cotton.unique <- fungi_sperm@tax_table %>%
  as.matrix() %>%
  as.data.frame() %>%
  subset(OTU %in% overlap[[6]]) %>%
  mutate(Crop = "Cotton")

bulk.unique <- fungi_sperm@tax_table %>%
  as.matrix() %>%
  as.data.frame() %>%
  subset(OTU %in% overlap[[7]]) %>%
  mutate(Crop = "Bulk Soil")

soy.cot <- fungi_sperm@tax_table %>%
  as.matrix() %>%
  as.data.frame() %>%
  subset(OTU %in% overlap[[2]]) %>%
  mutate(Crop = "Soybean and Cotton")

combined.unique <- rbind.data.frame(soy.cot, soybean.unique, cotton.unique, bulk.unique)

combined.unique %>%
  group_by(Crop, Phylum) %>%
  count() %>%
  mutate(Genus.label = ifelse(n >= 1, Phylum, "Other")) %>%
  ggplot(aes(x = Crop, y = n, fill = Genus.label)) +
  geom_bar(position="fill", stat = "identity") + 
  theme_classic() + 
  xlab("")+
  ylab("Proportion")+
  scale_fill_manual(values=c(cbbPalette, ibm.cbb, tol.cbb, "grey", "green"))
