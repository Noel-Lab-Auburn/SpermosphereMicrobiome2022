################Data Preprocessing############

###### Libraries #####
library(phyloseq)
library(decontam)
library(vegan)
library(tidyverse)
library(metagenomeSeq)
library(ggplot2)
library(ggpubr)
library(Biostrings)
library(microbiome)

##### Set global options #####

# no scientific notation
options(scipen=10000) 

# color blind pallets used throughout 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ibm.cbb <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")
tol.cbb <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")


########### Bacteria #####

###### Read in data ####

# Metadata #
samp_dat_bac <- read.csv("Bacteria/spermospheremetadata.csv", na.strings = "NA")

rownames(samp_dat_bac) <- samp_dat_bac$Code #row names must match OTU table headers
SAMP.bac <- phyloseq::sample_data(samp_dat_bac)

# OTU table #
otu_bac <- read.csv("Bacteria/otu_table_16s.csv")
rownames(otu_bac) <- otu_bac$OTU_ID
otu_bac <- otu_bac[,-1]
OTU.bac <- phyloseq::otu_table(otu_bac, taxa_are_rows = TRUE)

any(is.na(otu_bac)) # no NA in the OTU table

# Taxonomy #
taxonomy.bac <- read.csv("Bacteria/16s_taxonomy.csv")
rownames(taxonomy.bac) <- taxonomy.bac$OTU
TAX.bac <- phyloseq::tax_table(as.matrix(taxonomy.bac))

all.equal(rownames(samp_dat_bac), colnames(otu_bac))

# Fasta #
FASTA.bac <- readDNAStringSet("Bacteria/otus_16s.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

# Phylogentic tree #
tree <- phyloseq::read_tree("Bacteria/otus_16s_midpoint.tre")

###### Create Initial Phyloseq object #####
# Merge reads into Phyloseq object #
bac.unedited <- phyloseq::phyloseq(OTU.bac, TAX.bac, FASTA.bac, SAMP.bac, tree)

###### Decontaminate #####
bac.unedited@sam_data$Sample_or_Control <- ifelse(bac.unedited@sam_data$Crop %in% c("NEC", "Water"), "Control Sample", "True Sample")
sample_data(bac.unedited)$is.neg <- sample_data(bac.unedited)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(bac.unedited, method="prevalence", neg="is.neg", threshold = 0.1, normalize = TRUE)
badTaxa <- rownames(contamdf.prev[contamdf.prev$contaminant == TRUE,])

print(badTaxa)

ps.pa <- transform_sample_counts(bac.unedited, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
decontaminate.bac <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
  geom_point() +
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (True Samples)") + 
  scale_color_manual(values = cbbPalette)+ 
  ggtitle("Prokaryote") +
  theme_classic()

goodTaxa <- setdiff(taxa_names(bac.unedited), badTaxa)
bac_sub_no_bad <- prune_taxa(goodTaxa, bac.unedited)


###### Taxonomy filtering #####
# remove OTUs that are mitochondria, chloroplast, or unidentified at the kingdom level 
bac_no_chloro <- bac_sub_no_bad %>% 
  phyloseq::subset_taxa(Order != "Chloroplast") %>%
  phyloseq::subset_taxa(Family != "Mitochondria") %>%
  phyloseq::subset_taxa(Kingdom != "unidentified")

###### Mock Community analysis ##### 
# positive controls
bac_mock <- bac_no_chloro %>% 
  subset_samples(Crop == "MOCK") %>%
  phyloseq::filter_taxa(function(x) sum(x) > 2, TRUE) # filter OTUs to have more than 1 read in mock samples

mock2 <- microbiome::transform(bac_mock, "compositional") # relative abundance transform

sequenced.mock.bac <- mock2 %>%
  psmelt() %>% 
  ggplot(aes(Sample, Abundance, fill = Label)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values= c(cbbPalette, ibm.cbb, tol.cbb)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Relative abundance (%)",
       title = "Prokaryote") + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.text = element_text(face = "italic", size = 5),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, 'cm')) 
sequenced.mock.bac

# Adding in theoretical distribution - the last two are fungi and are not expected to be amplified with 16S
Label <- c("Pseudomonas aeruginosa", 
             "Escherichia coli",
             "Salmonella enterica", 
             "Lactobacillus fermentum", 
             "Enterococcus faecalis", 
             "Staphylococcus aureus", 
             "Listeria monocytogenes", 
             "Bacillus subtilis")

# theoretical species composition in the mock community
Abundance <- c(rep(0.125, 8))

th.mock <- data.frame(Label, Abundance)
th.mock$Sample <- "Theoretical"

th.mock$Label <- factor(th.mock$Label, levels = c("Lactobacillus fermentum", 
                                                  "Staphylococcus aureus", 
                                                  "Bacillus subtilis",
                                                  "Escherichia coli",
                                                  "Listeria monocytogenes",
                                                  "Enterococcus faecalis",
                                                  "Salmonella enterica",
                                                  "Pseudomonas aeruginosa"))


theory.mock <- ggplot(th.mock, aes(Sample, Abundance, fill = Label)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values= c(cbbPalette[[1]], 
                              cbbPalette[[2]], 
                              cbbPalette[[3]], 
                              cbbPalette[[4]], 
                              cbbPalette[[5]],
                              cbbPalette[[6]],
                              cbbPalette[[8]],
                              "violet", "pink", "grey", "black", "blue")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Relative abundance (%)",
       title = "Theoretical composition") + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank())

# I think maybe the theoretical mock community can also be mentioned in the figure legend. 

mock.composition <- mock2 %>%
  psmelt() %>%
  group_by(Label) %>%
  summarise(MeanRelAbund = mean(Abundance)) %>%
  arrange(-MeanRelAbund)

# these 8 OTUs made up 99.9% of the mock composition. These OTUs also match the 8 supposed to be in the mock
sum(mock.composition[1:8,]$MeanRelAbund)

###### Data filtering #####
# remove samples with less than 10,000 reads
# removed emergence data since not all samples emerged and their inclusion was questionable from the begining. 
# Removed one outlier-S185_86 with very dominant community structure; infected?)

bac_sperm_noE <- bac_no_chloro %>% 
  subset_samples(Time.Point %in% c("0", "6", "12", "18")) %>%
  prune_samples(sample_sums(.) > 10000, .) %>% # remove samples below 10,000 reads
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) # remove taxa with less than 1 reads (i.e., those not present in objective 1)

bac_sperm <- subset_samples(bac_sperm_noE, Code != "S185_86")

###### RDS of Non-normalized Prokaryote data ######
# Save an object to a file
saveRDS(bac_sperm, file = "Bacteria/Bacteria_spermosphere_nonnorm_112922.rds")
# Restore the object
bac_sperm <- readRDS(file = "Bacteria/Bacteria_spermosphere_nonnorm_112922.rds")

###### READS PER SAMPLE ######
sample.sums <- sort.DataFrame(data.frame(sample_sums(bac_sperm)))

read.dist.bac <- ggplot(sample.sums, aes(x = sample_sums.bac_sperm.)) +
  geom_histogram(color = "black", fill = cbbPalette[[4]]) + 
  theme_classic() +
  xlab("Read Depth") + 
  ggtitle("Prokaryote")

sum(sample_sums(bac_sperm)) # total reads = 2,090,814
median(sample_sums(bac_sperm)) # 29237.5

###### Rarefaction analysis #####
sam.data <- data.frame(bac_sperm@sam_data)
bOTU.table <- otu_table(bac_sperm) %>%
  as.data.frame() %>%
  as.matrix()
  
raremax <- min(rowSums(t(bOTU.table)))
rare.fun <- rarecurve(t(bOTU.table), step = 1000, sample = raremax, tidy = T)

bac.rare.curve.extract2 <- left_join(sam.data, rare.fun, by = c("Code" = "Site"))

bac.rare <- ggplot(bac.rare.curve.extract2, aes(x = Sample, y = Species, group = Code, color = Crop)) + 
  #geom_point() +
  geom_line() + 
  xlab("Reads") + 
  ylab("Number of OTUs") +
  ggtitle("Prokaryote") +
  theme_classic() + 
  geom_vline(xintercept = median(sample_sums(bac_sperm)), linetype = "dashed") +
  scale_color_manual(values = cbbPalette)

###### Metagenome CSS normalization ######
MGS <- phyloseq_to_metagenomeSeq(bac_sperm) #converts to metagenomeseq format
p <- metagenomeSeq::cumNormStatFast(MGS)
MGS <- metagenomeSeq::cumNorm(MGS, p =p)
metagenomeSeq::normFactors(MGS) # exports the normalized factors for each sample
norm.bac <- metagenomeSeq::MRcounts(MGS, norm = T) 
norm.bac.OTU <- phyloseq::otu_table(norm.bac, taxa_are_rows = TRUE) #exports the new otu table
bac.css.norm <- phyloseq::phyloseq(norm.bac.OTU, FASTA.bac, SAMP.bac, TAX.bac, tree) #new otu table phyloseq object

saveRDS(bac.css.norm, file = "Bacteria/Bacteria_spermosphere_CSS_112922.rds")
# Restore the object
bac.css.norm <- readRDS(file = "Bacteria/Bacteria_spermosphere_CSS_112922.rds")


########### Fungi #####

###### Read in data ####

#Loading the mapping file
samp_dat <- read.csv("Fungi/METADATA.csv", na.strings = "NA")

rownames(samp_dat) <- samp_dat$Code #row names must match OTU table headers
SAMP.fungi <- phyloseq::sample_data(samp_dat)

# OTU table 
otu <- read.csv("Fungi/OTU_Table.csv")
rownames(otu) <- otu$OTU
otu <- otu[,-1]
OTU.fungi <- phyloseq::otu_table(otu, taxa_are_rows = TRUE)

any(is.na(otu)) # no NA in the OTU table

# Taxonomy
taxonomy.fungi <- read.csv("Fungi/fungal_taxonomy_DADA2_NBC.csv")
rownames(taxonomy.fungi) <- taxonomy.fungi$OTU
taxonomy.fungi2 <- taxonomy.fungi %>%
  subset(Kingdom != "unidentified")

TAX.fungi <- phyloseq::tax_table(as.matrix(taxonomy.fungi2))

# Fasta
FASTA.fungi <- readDNAStringSet("Fungi/otus_R1.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

###### Create Initial Phyloseq object ######
fungi.unedited <- phyloseq::phyloseq(OTU.fungi, TAX.fungi, FASTA.fungi, SAMP.fungi)

###### Decontaminate  ######
fungi.unedited@sam_data$Sample_or_Control <- ifelse(fungi.unedited@sam_data$Crop == "NEC", "Control Sample", "True Sample")
sample_data(fungi.unedited)$is.neg <- sample_data(fungi.unedited)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(fungi.unedited, method="prevalence", neg="is.neg", threshold = 0.1, normalize = TRUE)
badTaxa <- rownames(contamdf.prev[contamdf.prev$contaminant == TRUE,])

print(badTaxa)

ps.pa <- transform_sample_counts(fungi.unedited, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
decontaminate <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
  geom_point() +
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (True Samples)") + 
  ggtitle("Fungi") +
  theme_classic() + 
  scale_color_manual(values = cbbPalette)

goodTaxa <- setdiff(taxa_names(fungi.unedited), badTaxa)
fungi_sub_no_bad <- prune_taxa(goodTaxa, fungi.unedited)

###### Mock Community #######
fungi_mock <- fungi_sub_no_bad %>% 
  subset_samples(Crop == "MOCK") %>%
  phyloseq::filter_taxa(function(x) sum(x) > 2, TRUE)

mock2 <- microbiome::transform(fungi_mock, "compositional") # relative abundance transform

sequenced.mock.fungi <- mock2 %>%
  psmelt() %>% 
  ggplot(aes(Sample, Abundance, fill = Species)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values= c(cbbPalette, ibm.cbb, tol.cbb, "violet", "pink", "grey", "black", "blue", "green")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Relative abundance (%)",
       title = "Fungi") + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.text = element_text(face = "italic", size = 5),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, 'cm')) 
sequenced.mock.fungi

mock.composition <- mock2 %>%
  psmelt() %>%
  group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species, Label) %>%
  summarise(MeanRelAbund = mean(Abundance)) %>%
  arrange(-MeanRelAbund) %>%
  ungroup() %>%
  group_by(Kingdom) %>%
  summarise(SumAbund = sum(MeanRelAbund)) 
mock.composition
 
# OTUs classified into the mock kingdom made up 99.9% of the reads in mock samples

fungi_not_mock <- fungi_sub_no_bad %>% 
  subset_samples(Crop != "MOCK") %>%
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

not_mock <- microbiome::transform(fungi_not_mock, "compositional") # relative abundance transform

not.mock.composition <- not_mock %>%
  psmelt() %>%
  group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species, Label) %>%
  summarise(MeanRelAbund = mean(Abundance)) %>%
  arrange(-MeanRelAbund) %>%
  ungroup() %>%
  group_by(Kingdom) %>%
  summarise(SumAbund = sum(MeanRelAbund)) %>%
  arrange(-SumAbund)
not.mock.composition
sum(not.mock.composition$SumAbund[c(1,3:10)])

# OTUs classified into the Fungal kingdom made up 98.8% of the reads in real samples

###### Taxonomy and Sample filtering #####
# remove OTUs that are not fungi, or unidentified at the kingdom level 
fungi_sperm <- fungi_sub_no_bad %>% 
  phyloseq::subset_taxa(Kingdom %in% c("Fungi")) %>%
  subset_samples(Crop %in% c("Bulk Soil", "Cotton ", "Soybean") & Time.Point != "E") %>%
  prune_samples(sample_sums(.) > 1000, .) %>% # remove samples below 1,000 reads
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) # remove taxa not present in samples cut out

###### RDS of Non-normalized Fungi data #####
# Save an object to a file
saveRDS(fungi_sperm, file = "Fungi/Fungi_spermosphere_unedited_083022.rds")
# Restore the object
fungi.obj1.unedited <- readRDS(file = "Fungi/Fungi_spermosphere_unedited_083022.rds")

###### READS PER SAMPLE #######
sample.sums <- sort.DataFrame(data.frame(sample_sums(fungi.obj1.unedited)))

read.dist <- ggplot(sample.sums, aes(x = sample_sums.fungi.obj1.unedited.)) +
  geom_histogram(color = "black", fill = cbbPalette[[4]]) + 
  theme_classic() +
  xlab("Read Depth") +
  ylab("Number of samples") + 
  ggtitle("Fungi")

sum(taxa_sums(fungi.obj1.unedited)) # total reads = 2534301

mean(sample_sums(fungi.obj1.unedited)) # 36729
median(sample_sums(fungi.obj1.unedited)) # 37933 reads

######## Rarefaction anlaysis ######## 
sam.data <- data.frame(fungi.obj1.unedited@sam_data)
fOTU.table <- otu_table(fungi.obj1.unedited) %>%
  as.data.frame() %>%
  as.matrix()
rare.fun <- rarecurve(t(fOTU.table), step = 1000, sample = raremax, tidy = TRUE)

fungi.rare.curve.extract2 <- left_join(rare.fun, sam.data, by = c("Site" = "Code"))

fungi.rare <- ggplot(fungi.rare.curve.extract2, aes(x = Sample, y = Species, group = Site, color = Crop)) + 
  #geom_point() +
  scale_color_manual(values = cbbPalette)+
  geom_line() + 
  xlab("Reads") + 
  ylab("Number of OTUs") + 
  ggtitle("Fungi")+
  theme_classic() + 
  geom_vline(xintercept = 40212, linetype = "dashed") +
  ggtitle("") 

######### Metagenome CSS normalization #########
MGS <- phyloseq_to_metagenomeSeq(fungi.obj1.unedited)
p <- metagenomeSeq::cumNormStatFast(MGS)

MGS <- metagenomeSeq::cumNorm(MGS, p =p)

metagenomeSeq::normFactors(MGS) # exports the normalized factors for each sample

norm.fungi <- metagenomeSeq::MRcounts(MGS, norm = T)

norm.fungi.OTU <- phyloseq::otu_table(norm.fungi, taxa_are_rows = TRUE)

fungi.css.norm <- phyloseq::phyloseq(norm.fungi.OTU, TAX.fungi, FASTA.fungi, SAMP.fungi)

######## Save CSS object to a file ########
saveRDS(fungi.css.norm, file = "Fungi/Fungi_CSSNorm_083022.rds")
# Restore the object
fungi.css.norm <- readRDS(file = "Fungi/Fungi_CSSNorm_083022.rds")


###### Supplemental Figure 1 ######

mock.combined <- ggpubr::ggarrange(sequenced.mock.bac, sequenced.mock.fungi, nrow = 1)
rare.combined <- ggpubr::ggarrange(bac.rare, fungi.rare, nrow = 1, common.legend = T)
decontaminate.combined <- ggpubr::ggarrange(decontaminate.bac, decontaminate, nrow = 1, common.legend = T)
read.dist.combined <- ggpubr::ggarrange(read.dist.bac, read.dist, nrow = 1)

supp.fig.1 <- ggpubr::ggarrange(mock.combined, 
                                rare.combined, 
                                decontaminate.combined, 
                                read.dist.combined, nrow = 4, ncol = 1, labels = c("a", "b", "c", "d"))

