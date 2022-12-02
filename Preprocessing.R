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

##### Set global options #####

# no scientific notation
options(scipen=10000) 

# color blind pallet used throughout 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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

###### Create a phyloseq object #####
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
decontaminate <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") + 
  scale_color_manual(values = cbbPalette)+ 
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
  phyloseq::filter_taxa(function(x) sum(x) > 1, TRUE) # filter OTUs to have more than 1 read in mock samples

mock2 <- microbiome::transform(bac_mock, "compositional") # relative abundance transform

sequenced.mock <- mock2 %>%
  psmelt() %>% 
  ggplot(aes(Sample, Abundance, fill = Label)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values= c(cbbPalette, "violet", "pink", "grey", "black", "blue")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Relative abundance (%)",
       title = "Sequenced") + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank()) 

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

read.dist <- ggplot(sample.sums, aes(x = sample_sums.bac_sperm.)) +
  geom_histogram(color = "black") + 
  theme_classic() +
  xlab("Read Depth")

sum(sample_sums(bac_sperm)) # total reads = 2,090,814
median(sample_sums(bac_sperm)) # 29237.5

###### Rarefaction analysis #####
sam.data <- data.frame(bac_sperm@sam_data)
bOTU.table <- bac_sperm@otu_table
S <- specnumber(t(bOTU.table)) # observed number of species
raremax <- min(rowSums(t(bOTU.table)))
Srare <- rarefy(t(bOTU.table), raremax)
rare.fun <- rarecurve(t(bOTU.table), step = 1000, sample = raremax, cex = 0.6)

bac.rare.curve.extract <- NULL
for(i in 1:length(rare.fun)){
  sample.200 <- data.frame(rare.spec = rare.fun[[i]])
  sample.200$read_depth <- attr(rare.fun[[i]], "Subsample")
  sample.200$Code <- rownames(t(bOTU.table[,i]))
  bac.rare.curve.extract <- rbind.data.frame(bac.rare.curve.extract, sample.200)
}
bac.rare.curve.extract2 <- left_join(sam.data, bac.rare.curve.extract, by = "Code")

bac.rare <- ggplot(bac.rare.curve.extract2, aes(x = read_depth, y = rare.spec, group = Code, color = Crop)) + 
  #geom_point() +
  geom_line() + 
  xlab("Reads") + 
  ylab("Number of OTUs") + 
  theme_classic() + 
  geom_vline(xintercept = median(sample_sums(bac_sperm)), linetype = "dashed") +
  scale_color_manual(values = cbbPalette)

###### Supplemental Figure 1 ######
ggpubr::ggarrange(bac.rare, decontaminate, read.dist, ggpubr::ggarrange(theory.mock, sequenced.mock), nrow = 2, ncol = 2, labels = c("a", "b", "c", "d"))

####### Metagenome CSS normalization ######
MGS <- phyloseq_to_metagenomeSeq(bac_sperm) #converts to metagenomeseq format
p <- metagenomeSeq::cumNormStatFast(MGS)
MGS <- metagenomeSeq::cumNorm(MGS, p =p)
metagenomeSeq::normFactors(MGS) # exports the normalized factors for each sample
norm.bac <- metagenomeSeq::MRcounts(MGS, norm = T) 
norm.bac.OTU <- phyloseq::otu_table(norm.bac, taxa_are_rows = TRUE) #exports the new otu table
bac.css.norm <- phyloseq::phyloseq(norm.bac.OTU, FASTA.bac, SAMP.bac, TAX.bac, tree) #new otu table phyloseq object

