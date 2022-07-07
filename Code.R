################BACTERIOME############

library(phyloseq)
library(decontam)
library(vegan)
library(tidyverse)
library(metagenomeSeq)
library(ggplot2)
library(ggpubr)
library(Biostrings)

source('functions_themes.R')
options(scipen=10000) 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##### Bacteria #####

samp_dat_bac <- read.csv("~/Kemi/My Research/Spermosphere/Spermosphere Microbiome/Bacteria_Merged/spermospheremetadata.csv", na.strings = "NA")

rownames(samp_dat_bac) <- samp_dat_bac$Code #row names must match OTU table headers
SAMP.bac <- phyloseq::sample_data(samp_dat_bac)

# OTU table 
otu_bac <- read.csv("~/Kemi/My Research/Spermosphere/Spermosphere Microbiome/Bacteria_Merged/otu_table_16s.csv")
rownames(otu_bac) <- otu_bac$OTU_ID
otu_bac <- otu_bac[,-1]
OTU.bac <- phyloseq::otu_table(otu_bac, taxa_are_rows = TRUE)

any(is.na(otu_bac)) # no NA in the OTU table

# Taxonomy
taxonomy.bac <- read.csv("~/Kemi/My Research/Spermosphere/Spermosphere Microbiome/Bacteria_Merged/16s_taxonomy.csv")
rownames(taxonomy.bac) <- taxonomy.bac$OTU
TAX.bac <- phyloseq::tax_table(as.matrix(taxonomy.bac))

all.equal(rownames(samp_dat_bac), colnames(otu_bac))

# Fasta
FASTA.bac <- readDNAStringSet("~/Kemi/My Research/Spermosphere/Spermosphere Microbiome/Bacteria_Merged/otus_16s.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

#Merge reads into Phyloseq object
bac.unedited <- phyloseq::phyloseq(OTU.bac, TAX.bac, FASTA.bac, SAMP.bac)

# Decontaminate
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

## remove OTUs that are mitochondria or chloroplast
bac_no_chloro <- bac_sub_no_bad %>% 
  phyloseq::subset_taxa(Order != "Chloroplast") %>%
  phyloseq::subset_taxa(Family != "Mitochondria") %>%
  phyloseq::subset_taxa(Kingdom != "unidentified")

# positive controls
bac_mock <- bac_no_chloro %>% 
  subset_samples(Crop == "MOCK") %>%
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

#unique(bac_mock@tax_table@.Data[,2])

mock = filter_taxa(bac_mock, function(x) sum(x > 1) > (0.1*length(x)), TRUE)

mock2 <- microbiome::transform(mock, "compositional")


mock.community <- microbiome::plot_composition(mock2,
                                               taxonomic.level = "Genus",
                                               x.label = "Sample_ID") +
  theme_classic()+
  #scale_fill_manual(values= cbbPalette) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Relative abundance (%)",
       title = "Relative abundance data") + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.text = element_text(face = "italic")) 

Mock3 <- mock.community$data
# at the 0.00015 and 0.00015 prevelence and detection cutoffs 99.9% of the positive control samples were classified as mock community. 

Mock4 <- left_join(Mock3, taxonomy.bac, by = c("Tax" = "OTU"))

#Plot for mock communities 
mock.community1 <- ggplot(Mock4, aes(Sample, Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  #scale_fill_manual(values= cbbPalette) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Relative abundance (%)",
       title = "Relative abundance data") + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.text = element_text(face = "italic")) 

#Removed Emergence data and outlier-S185_86)
bac_sperm_noE <- bac_no_chloro %>% 
  subset_samples(Time.Point %in% c("0", "6", "12", "18")) %>%
  prune_samples(sample_sums(.) > 10000, .) %>% # remove samples below 10,000 reads
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) # remove taxa with zero reads (i.e., those not present in objective 1)

bac_sperm <- subset_samples(bac_sperm_noE, Code != "S185_86")

# Save an object to a file
saveRDS(bac_sperm, file = "Bacteria_spermosphere_nonnorm_041922.rds")
# Restore the object
bac_sperm <- readRDS(file = "Bacteria_spermosphere_nonnorm_041922.rds")

# READS PER SAMPLE
sample.sums <- sort.DataFrame(data.frame(sample_sums(bac_sperm)))

read.dist <- ggplot(sample.sums, aes(x = sample_sums.bac_sperm.)) +
  geom_histogram(color = "black") + 
  theme_classic() +
  xlab("Read Depth")

sum(sample_sums(bac_sperm)) # total reads = 2,456,177
#after removing emergence data: 2,090,814

# minimum number of reads per sample is 14,568, which is really good. 

mean(sample_sums(bac_sperm)) # 28896.2 ##29868.77
median(sample_sums(bac_sperm)) # 28151 reads ##29237.5

# Rarefaction analysis 
sam.data <- data.frame(bac_sperm@sam_data)
bOTU.table <- bac_sperm@otu_table
S <- specnumber(t(bOTU.table)) # observed number of species
raremax <- min(rowSums(t(bOTU.table)))
Srare <- rarefy(t(bOTU.table), raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rare.fun <- rarecurve(t(bOTU.table), step = 1000, sample = raremax, col = "blue", cex = 0.6)

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
  geom_vline(xintercept = 28773, linetype = "dashed") +
  scale_color_manual(values = cbbPalette)

# Supplemental Figure 1
ggpubr::ggarrange(bac.rare, decontaminate, read.dist, mock.community1, nrow = 2, ncol = 2, labels = c("a", "b", "c", "d"))

bac_sperm2 <- filter_taxa(bac_sperm, function(x) sum(x > 1) > (0.05*length(x)), TRUE)

# Metagenome CSS normalization
MGS <- phyloseq_to_metagenomeSeq(bac_sperm)
p <- metagenomeSeq::cumNormStatFast(MGS)

MGS <- metagenomeSeq::cumNorm(MGS, p =p)

metagenomeSeq::normFactors(MGS) # exports the normalized factors for each sample

norm.bac <- metagenomeSeq::MRcounts(MGS, norm = T)

norm.bac.OTU <- phyloseq::otu_table(norm.bac, taxa_are_rows = TRUE)

bac.css.norm <- phyloseq::phyloseq(norm.bac.OTU, FASTA.bac, SAMP.bac, TAX.bac)

bac.bray = phyloseq::distance(bac.css.norm, "bray") # create bray-curtis distance matrix

# PERMANOVA - testing for differences in centroids
set.seed(12325)
adonis2(bac.bray~Crop*Time.Point, as(sample_data(bac.css.norm), "data.frame"))

GP.ord <- phyloseq::ordinate(bac_sperm, "MDS", "bray")
p2 = phyloseq::plot_ordination(bac_sperm, GP.ord, type="samples", color="Crop") 
p2

bac.ggplot.data <- p2$data

global <- ggplot(bac.ggplot.data, aes(x = Axis.1, y = Axis.2, fill = Crop, shape = Time.Point)) +
  geom_point(size=4, alpha = 0.7)+
  scale_shape_manual(values=c(21,22, 23, 24,25,26,27), name = "")+
  scale_fill_manual(values=cbbPalette, name = "")+
  xlab("PcoA1 (30.3%)") + 
  ylab("PcoA2 (15.3%)") +
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(fill=guide_legend(override.aes=list(shape= 21)),
         shape=guide_legend(override.aes=list(fill= "black")))+
  ggtitle("") +
  theme_bw()



#Separating based on time-points (re-run for each time-point)
bac_sperm_18 <- bac.css.norm %>% 
  subset_samples(Time.Point == "18") %>% 
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

bac.bray.18 = phyloseq::distance(bac_sperm_18, "bray") # create bray-curtis distance matrix

set.seed(12325)
adonis.output  <- adonis2(bac.bray.18~Crop, as(sample_data(bac_sperm_18), "data.frame"))

GP.ord <- phyloseq::ordinate(bac_sperm_18, "MDS", "bray")
p2 = phyloseq::plot_ordination(bac_sperm_18, GP.ord, type="samples", color="Crop") 
p2

bac.ggplot.data <- p2$data

bac.sprm.18.pcoa <- ggplot(bac.ggplot.data, aes(x = Axis.1, y = Axis.2, color = Crop)) +
  #geom_label()+
  geom_point(size=4, alpha = 0.7)+
  scale_color_manual(values=cbbPalette, name = "")+
  xlab("PcoA1") + 
  ylab("PcoA2") +
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  ggtitle("Prokaryote - 18 hrs") +
  theme_bw() 



ggpubr::ggarrange(bac.sprm.0.pcoa, 
                  bac.sprm.6.pcoa,
                  bac.sprm.12.pcoa, 
                  bac.sprm.18.pcoa, nrow = 1, ncol = 4, labels = c("a", "b", "c", "d"), common.legend = T)
# No difference after 0 hours P = 0.452
# 6 hours P < 0.001
# 12 hours P < 0.001
# 18 hours P < 0.002
# Emergence P = 0.727

# DISPERSION - testing for differences in variation due to management
beta.disp <- betadisper(bac.bray, bac.css.norm@sam_data$Time.Point)
permutest(beta.disp)
TukeyHSD(beta.disp)
plot(beta.disp)
boxplot(beta.disp)
# Which groups are different from others?

#Dispersion for each Time-Point
beta.disp <- betadisper(bac.bray.18, bac_sperm_18@sam_data$Crop)
permutest(beta.disp)
TukeyHSD(beta.disp)
plot(beta.disp)
boxplot(beta.disp)

# ANOSIM - Rank based differences in centroids
anosim(prok.dist.bray, physeq.css@sam_data$Time.Point) #Are there significant changes?
# Is differences in management due to differences in centroids (means) or differences in dispersion?

# Perform indicator species analysis just considering the four original groups
indicator.dist.crop <- indicspecies::multipatt(as.data.frame(t(bac.css.norm@otu_table)), cluster = bac.css.norm@sam_data$Crop, max.order = 1)
# summary of results
summary(indicator.dist.crop, indvalcomp = TRUE)

# get all data for each OTU
all.OTUs <- indicator.dist.crop$sign
all.OTUs$OTU <- rownames(all.OTUs)
sig.otus <- na.omit(all.OTUs[all.OTUs$p.value <= 0.01,])

frequency.groups <- data.frame(table(all.OTUs$index))
frequency.groups$categories <- colnames(indicator.dist.crop$str)

sig.crop <- left_join(sig.otus, taxonomy.bac, by = "OTU")

sig.crop$cropaffected <- ifelse(sig.crop$index == 1, "Soybean", 
                                ifelse(sig.crop$index == 2, "Cotton", "Bulk Soil"))

##Community Composition acccording to Microbiome Package##
library(microbiome)
# Make sure we use functions from correct package
transform <- microbiome::transform

pseq <- transform(bac.css.norm, "compositional")
pseq <- aggregate_rare(pseq, level = "Genus", detection = 1/100, prevalence = 2/100)

# Pick sample subset
library(phyloseq)

install.packages("hrbrthemes", repos = c("https://cinc.rud.is", "https://cloud.r-project.org/"))
library(hrbrthemes)
install.packages("gcookbook")
library(gcookbook)
library(tidyverse)
p <- plot_composition(pseq,
                      taxonomic.level = "Genus",
                      sample.sort = "Crop",
                      x.label = "Crop") +
  scale_fill_brewer("Genera", palette = "Paired") +
  guides(fill = guide_legend(ncol = 1)) +
  scale_y_percent() +
  labs(x = "Samples", y = "Relative abundance (%)",
       title = "Relative abundance data") + 
  theme_ipsum(grid="Y") +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.text = element_text(face = "italic"))
print(p)  

##COmmunity Composition with psmelt

#Remove OTUs not greater than 10 reads
keepTaxa = apply(X = as(otu_table(bac.css.norm), "matrix") > 10, MARGIN = 1, FUN = sum) >= 10
bac.prop = prune_taxa(keepTaxa, bac.css.norm)


glom <- tax_glom(bac.prop, taxrank = 'Class')


genus_colors <-c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", 
                 "#FB9A99", "#E31A1C", "#FDBF6F", "#00FFFF", 
                 "#FF7F00", "#CAB2D6","#8A7C64","#652926",
                 "#6A3D9A", "#B15928", "#FFC000")

# agglomerate taxa at Genus level
all_genus <- tax_glom(bac.prop, "Genus", NArm = TRUE)
# Get top 15 genera
top15_genera <- names(sort(taxa_sums(all_genus), decreasing=TRUE))[1:20]
# Transform Taxa counts to relative abundance
all_genus_relabun <- transform_sample_counts(all_genus, function(OTU) OTU/sum(OTU) * 100)
# Extract the top 15 taxa 
all_genus_top15 <- prune_taxa(top15_genera, all_genus_relabun)
# Convert into dataframe
taxa_abundance_table_genus <- psmelt(all_genus_top15)
#Plot top 15 genus
T= StackedBarPlot_genus_endo <- taxa_abundance_table_genus %>% 
  ggplot(aes(x =Crop, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Sample",
       y = "Relative Abundance",
       title = "Genus Relative Abundance (Spermosphere)") +
  facet_grid(~ Time.Point, scales = "free") +
  theme(
    axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 12))+
  scale_fill_manual(values=genus_colors)

#To put the samples in a specific order (according to time-point)

T$data$Time.Point <- factor(T$data$Time.Point, c("0","6","12","18"))

print(T)

ggsave("taxonomy1.pdf",width = 15, height = 7, units="in", dpi=700)


# Core - abundance occupancy modeling
core.prioritizing <- function(phyloseq.object){
  
  set.seed(19)
  rare.phyloseq.object <- rarefy_even_depth(phyloseq.object, replace=TRUE)
  
  nReads=sample_sums(rare.phyloseq.object)[[1]]                                                                 # input dataset needs to be rarified and the rarifaction depth included 
  otu <- rare.phyloseq.object@otu_table %>%
    as("matrix")
  map <- rare.phyloseq.object@sam_data %>%
    as("data.frame")
  
  otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
  otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
  otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # mean relative abundance
  occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') # combining occupancy and abundance data frame
  
  # Ranking OTUs based on their occupancy
  # For caluclating raking index we included following conditions:
  #   - time-specific occupancy (sumF) = frequency of detection within time point (genotype or site)
  #   - replication consistency (sumG) = has occupancy of 1 in at least one time point (genotype or site) (1 if occupancy 1, else 0)
  
  PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>% 
    gather(Code, abun, -otu) %>%
    left_join(map, by = 'Code') %>%
    group_by(otu, Time.Point) %>%
    dplyr::summarise(time_freq=sum(abun>0)/length(abun),            # frequency of detection between time points
                     coreTime=ifelse(time_freq == 1, 1, 0)) %>%     # 1 only if occupancy 1 with specific time, 0 if not
    group_by(otu) %>%
    dplyr::summarise(sumF=sum(time_freq),
                     sumG=sum(coreTime),
                     nS=length(Time.Point)*2,           
                     Index=(sumF+sumG)/nS)                 # calculating weighting Index based on number of time points detected and 
  
  otu_ranked <- occ_abun %>%
    left_join(PresenceSum, by='otu') %>%
    transmute(otu=otu,
              rank=Index) %>%
    arrange(desc(rank))
  
  # Calculating the contribution of ranked OTUs to the BC similarity
  BCaddition <- NULL
  
  # calculating BC dissimilarity based on the 1st ranked OTU
  # with 36 samples there should be 630 combinations n!/r!
  otu_start=otu_ranked$otu[1]                   
  start_matrix <- as.matrix(otu[otu_start,])
  start_matrix <- t(start_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_s <- data.frame(x_names,x)
  df_s$rank_count <- 1
  BCaddition <- rbind(BCaddition,df_s)
  # calculating BC dissimilarity based on additon of ranked OTUs from 2nd to 500th. Can be set to the entire length of OTUs in the dataset, however it might take some time if more than 5000 OTUs are included.
  for(i in 2:500){                              
    otu_add=otu_ranked$otu[i]                       
    add_matrix <- as.matrix(otu[otu_add,])
    add_matrix <- t(add_matrix)
    start_matrix <- rbind(start_matrix, add_matrix)
    x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
    #x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
    df_a <- data.frame(x_names,x)
    df_a$rank_count <- i 
    BCaddition <- rbind.data.frame(BCaddition, df_a)
  }
  # calculating the BC dissimilarity of the whole dataset (not needed if the second loop is already including all OTUs) 
  x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))   
  x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
  df_full <- data.frame(x_names,x)
  df_full$rank_count <- length(rownames(otu))
  BCfull <- rbind.data.frame(BCaddition, df_full)
  
  BC_ranked <- BCfull %>%
    group_by(rank_count) %>%
    dplyr::summarise(MeanBC=mean(x)) %>%            # mean Bray-Curtis dissimilarity
    arrange(desc(-MeanBC)) %>%
    mutate(proportionBC=MeanBC/max(MeanBC))   # proportion of the dissimilarity explained by the n number of ranked OTUs
  Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
  increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
  increaseDF$rank <- as.numeric(increaseDF$rank)
  BC_ranked <- left_join(BC_ranked, increaseDF, by = c("rank_count" = "rank"))
  BC_ranked <- BC_ranked[-nrow(BC_ranked),]
  
  #Creating threshold for core inclusion - last call method
  
  #B) Final increase in BC similarity of equal or greater then 2% 
  lastCall <- last(as.numeric(BC_ranked$rank_count[(BC_ranked$IncreaseBC>=1.02)]))
  
  #Creating plot of Bray-Curtis similarity
  plot <- ggplot(BC_ranked[1:100,], aes(x=factor(BC_ranked$rank_count[1:100], levels=BC_ranked$rank_count[1:100]))) +
    geom_point(aes(y=proportionBC)) +
    theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
    geom_vline(xintercept=last(as.numeric(BC_ranked$rank_count[(BC_ranked$IncreaseBC>=1.02)])), lty=3, col='black', cex=.5) +
    labs(x='ranked OTUs',y='Bray-Curtis similarity') +
    annotate(geom="text", x=last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))+3, y=.5, label=paste("Last 2% increase (",last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)])),")",sep=''), color="black")
  
  core.otus.CSS.mean.T1 <- otu_ranked$otu[1:lastCall]
  return_list <- list(core.otus.CSS.mean.T1, plot, otu_ranked, occ_abun)
  return(return_list)
}

bac.no.norm.soybean <- bac_sperm %>% 
  subset_samples(Crop == "Soybean" & Code != "S185_86" & Time.Point %in% c("0","6", "12", "18")) %>% 
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)


core.rare.soybean <- core.prioritizing(bac.no.norm.soybean)

occ.abund.soybean <- core.rare.soybean[[4]]
occ.abund2.soybean <- left_join(occ.abund.soybean, taxonomy.bac, by = c("otu" = "OTU"))
occ.abund2.soybean$core <- ifelse(occ.abund2.soybean$otu %in% core.rare.soybean[[1]], "Core", "Not Core")

soybean <- ggplot() + 
  geom_point(data = occ.abund2.soybean[occ.abund2.soybean$core == "Core",], aes(x = log10(otu_rel), y = otu_occ, color = Phylum)) + 
  geom_point(data = occ.abund2.soybean[occ.abund2.soybean$core == "Not Core",], aes(x = log10(otu_rel), y = otu_occ), shape = 2, color = "grey", alpha = 0.5) + 
  theme_classic() 
#scale_color_manual(values = c(cbbPalette, "purple"))

fungicide <- ggplot() + 
  geom_point(data = occ.abund2.fungicide[occ.abund2.fungicide$core == "Core",], aes(x = log10(otu_rel), y = otu_occ, color = Genus)) + 
  geom_point(data = occ.abund2.fungicide[occ.abund2.fungicide$core == "Not Core",], aes(x = log10(otu_rel), y = otu_occ), shape = 2, color = "grey", alpha = 0.5) + 
  theme_classic() + 
  scale_color_manual(values = c(cbbPalette, "purple", "orange"))

ggpubr::ggarrange(control, fungicide, nrow = 2, ncol = 1, labels = c("a", "b"), common.legend = T)


library(tyRa)
library(minpack.lm)
library(Hmisc)
set.seed(19)
rare.phyloseq.object <- rarefy_even_depth(bac.no.norm.soybean, replace=TRUE)

nReads=sample_sums(rare.phyloseq.object)[[1]]                                                                 # input dataset needs to be rarified and the rarifaction depth included 
otu <- rare.phyloseq.object@otu_table %>%
  as("matrix")
taxa <- rownames(otu)
map <- rare.phyloseq.object@sam_data %>%
  as("data.frame")
spp.out <- tyRa::fit_sncm(spp = t(otu), pool=NULL, taxon=taxa)

predictions <- spp.out$predictions
predictions$otu <- rownames(predictions)
predictions$core <- ifelse(predictions$otu %in% core.rare.soybean[[1]], "core", "not core")
predictions2 <- left_join(predictions, taxonomy.bac, by = c("otu" = "OTU"))
predictions2[predictions2$fit_class == "Above prediction" & predictions2$core == "core",]

ggplot() +
  geom_point(data = predictions2, aes(x = log10(p), y = freq, color = fit_class, shape = core), alpha = 0.8, size = 2) +
  geom_line(color='black', data=predictions2, size=1, aes(y=predictions2$freq.pred, x=log10(predictions2$p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=predictions2, aes(y=predictions2$pred.upr, x=log10(predictions2$p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=predictions2, aes(y=predictions2$pred.lwr, x=log10(predictions2$p)), alpha=.25)+
  labs(x="log10(Mean relative abundance)", y="Occupancy") + 
  theme_classic() + 
  ylim(c(0, 1))+
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9")) +
  geom_label(data = predictions2[predictions2$core == "core" & predictions2$fit_class == "Below prediction" & log10(predictions2$p) < -2,], 
             aes(x = log10(p), y = freq, label = Label))


#####Fungi########
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("metagenomeSeq")

library(BiocManager)
BiocManager::install("microbiome")

library(devtools) # Load the devtools package
install_github("microbiome/microbiome") # Install the package

library(microbiome) 
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(Biostrings)
library(ggpubr)
library(decontam)
library(metagenomeSeq)
library(indicspecies)

# color blind pallet
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# load in our files from the HPC for input into phyloseq. 

# Metadata

samp_dat <- read.csv("~/Kemi/My Research/Spermosphere/Spermosphere Microbiome/METADATA.csv")
rownames(samp_dat) <- samp_dat$Code
#samp_dat <- samp_dat[,-1]
SAMP <- phyloseq::sample_data(samp_dat)

# Taxonomy 

tax <- read.csv("~/Kemi/My Research/Spermosphere/Spermosphere Microbiome/Taxonomy_Spermosphere_F_Rdp.csv")
rownames(tax) <- tax$OTU
TAX.fungi <- phyloseq::tax_table(as.matrix(tax))
head(tax)

#OTU table 
table <- read.csv("~/Kemi/My Research/Spermosphere/Spermosphere Microbiome/OTU_Table.csv")
rownames(table) <- table$OTU
table <- table[,-1]
OTU.fungi <- phyloseq::otu_table(table, taxa_are_rows = TRUE)
head(OTU.fungi)
#Meta data 


#fasta file
FASTA.fungi <- Biostrings::readDNAStringSet("~/Kemi/My Research/Spermosphere/Spermosphere Microbiome/otus_R1.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

#Phylogeny
# alternative; table2= read.csv(filechoose())


phyloseq.start <- phyloseq(SAMP, TAX.fungi, OTU.fungi, FASTA.fungi)

#to confirm sample names
sample.names(phyloseq.start)


phyloseq.start@otu_table # the OTU table
phyloseq.start@tax_table # the taxonomy table
phyloseq.start@sam_data # the metadata
phyloseq.start@refseq # the sequences 

# removing chloroplast or taxa not assigned at the domain level
physeq.no.chloro <- phyloseq.start %>% subset_taxa(Kingdom!= "unidentified")


## DECONTAMINATE
#Use the full dataset to call contaminants, then remove them, if they exist in the non plant OTU dataset
sample_data(physeq.no.chloro)$neg <- sample_data(physeq.no.chloro)$Crop == "NEC"
contamdf.prev <- isContaminant(physeq.no.chloro, method="prevalence", neg="neg", threshold = 0.1, normalize = TRUE)
badTaxa <- rownames(contamdf.prev[contamdf.prev$contaminant == TRUE,])

print(badTaxa)



# transform data to presence absence
ps.pa <- transform_sample_counts(physeq.no.chloro, function(abund) 1*(abund>0))

# making a dataframe for both negative and positive smaples.
ps.pa.neg <- prune_samples(sample_data(ps.pa)$neg == "TRUE", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$neg == "FALSE", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
decontaminated <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") + 
  theme_classic() + 
  scale_color_manual(values = c(cbbPalette[[1]], cbbPalette[[2]]))

print(decontaminated)

#Take out the contaminants 
goodTaxa <- setdiff(taxa_names(physeq.no.chloro), badTaxa)
physeq.clean <- prune_taxa(goodTaxa, physeq.no.chloro)

# Now that we have removed the contaminants using the negative controls lets get rid of those samples
physeq.clean.samples <- physeq.clean %>% 
  subset_samples(neg == "FALSE" & Crop != "MOCK") %>%
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

# Now lets look at the read distribution per sample and decide if we need to get rid of some samples because of low sequence depth
# a good general rule of thumb is samples below 1000 reads could be eliminated, although this isn't a hard rule, you can remove at 10,000 or more if you want
sort(sample_sums(physeq.clean.samples), decreasing = T) # read distribution

#to remove samples with low reads
fungi_spermosphere_1000reads <- prune_samples(sample_sums(physeq.clean.samples) > 1000, physeq.clean.samples) %>%
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

# how many total reads are we now working with?
sum(sample_sums(fungi_spermosphere_1000reads))


# What is our mean and median read depth per sample? 
mean(sample_sums(fungi_spermosphere_1000reads))

median(sample_sums(fungi_spermosphere_1000reads))
#median for oomycetes can be around 10,000 reads
#bacterian 50000
#fungi 20,000- 30,000

# positive controls
fungi_mock <- phyloseq.start %>%
  subset_samples(Crop == "MOCK") %>%
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)



unique(fungi_mock@tax_table@.Data[,2])



mock <- microbiome::transform(fungi_mock, "compositional")
mock2 <- microbiome::aggregate_rare(mock, level = "Genus", prevalence = 0.00015, detection = 0.00015)



mock.community <- microbiome::plot_composition(mock2,
                                               taxonomic.level = "Genus",
                                               x.label = "Sample_ID") +
  theme_classic()+
  #scale_fill_manual(values= cbbPalette) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Samples", y = "Relative abundance (%)",
       title = "Relative abundance data") +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.text = element_text(face = "italic"))

mock.community$data
#so at those prevelence and abundance cutoffs, 99.9% of the sequences in the positive controls were classified as mock community.

fungi_sperm_noE <- fungi_spermosphere_1000reads %>% 
  subset_samples(Time.Point %in% c("0", "6", "12", "18")) %>%
  prune_samples(sample_sums(.) > 10000, .) %>% # remove samples below 10,000 reads
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) # remove taxa with zero reads (i.e., those not present in objective 1)
sam.data$Code <- rownames(sam.data)
fungi_sperm <- subset_samples(fungi_sperm_noE, Code != "S.18.4")


#Lets make a histogram of the read distribution and put the median read depth
read.depths <- data.frame(sample_sums(fungi_sperm))
colnames(read.depths) <- "read.depth"
read.depth.plot <- ggplot(read.depths, aes(read.depth)) +
  geom_histogram(fill = cbbPalette[[3]], color = "black") + 
  geom_vline(xintercept = median(sample_sums(fungi_sperm)), linetype = "dashed") + 
  theme_classic() + 
  xlab("Read Depth")

# Rarefaction analysis
# how much of the species did we sample in our dataset?
sam.data <- data.frame(fungi_sperm@sam_data)
pOTU.table <- fungi_sperm@otu_table
S <- specnumber(t(pOTU.table)) # observed number of species
raremax <- min(rowSums(t(pOTU.table)))
Srare <- rarefy(t(pOTU.table), raremax) #if you sequence well enough, this step doesn't matter
#rareify to min read dept would lose samples as curve doesn't plateau at min value ~7000
rare.fun <- rarecurve(t(pOTU.table), step = 1000, sample = raremax, col = "blue", cex = 0.6)

# lets extract that data and make a nice plot in ggplot, because nobody likes base R plotting
prok.rare.curve.extract <- NULL
for(i in 1:length(rare.fun)){
  sample.200 <- data.frame(rare.spec = rare.fun[[i]])
  sample.200$read_depth <- attr(rare.fun[[i]], "Subsample")
  sample.200$Code <- rownames(t(pOTU.table[,i]))
  prok.rare.curve.extract <- rbind.data.frame(prok.rare.curve.extract, sample.200)
}
sam.data$Code <- rownames(sam.data)
prok.rare.curve.extract2 <- left_join(sam.data, prok.rare.curve.extract, by = "Code")

rare.curve <- ggplot(prok.rare.curve.extract2, aes(x = read_depth, y = rare.spec, group = Code)) + 
  #geom_point() +
  geom_line(color = "grey") + 
  xlab("Reads") + 
  ylab("Number of OTUs") + 
  theme_classic() + 
  geom_vline(xintercept = median(sample_sums(fungi_sperm)), linetype = "dashed")



# Stick many figures together for a publication ready figure
SuplementalFig1 <- ggarrange(rare.curve, read.depth.plot, decontaminated, nrow = 1, labels = c("a", "b", "c"))
ggsave("SupplementalFig1.pdf", dpi=300, width = 40, height = 10, units = "cm")


# Normalize Sampling reads based on proportions
physeq.prop <- transform_sample_counts(fungi_sperm, function(x) x/sum(x)) # normalize using proportions

# Normalize Sampling reads based on cumulative sum scaling (CSS normalization)
MGS <- phyloseq_to_metagenomeSeq(fungi_sperm)
p <- metagenomeSeq::cumNormStatFast(MGS)
MGS <- metagenomeSeq::cumNorm(MGS, p =p)
metagenomeSeq::normFactors(MGS) # exports the normalized factors for each sample
norm.bacteria <- metagenomeSeq::MRcounts(MGS, norm = T)
norm.bacteria.OTU <- phyloseq::otu_table(norm.bacteria, taxa_are_rows = TRUE)

physeq.css <- phyloseq::phyloseq(norm.bacteria.OTU, SAMP, TAX.fungi, FASTA.fungi)

# Save RDS files of each type of normalized reads for easy loading in the future and reproducibility 
# Save an object to a file
saveRDS(physeq.css, file = "Bacteria_Soil_CSSnorm_102021.rds")
saveRDS(physeq.prop, file = "Bacteria_Soil_propnorm_102021.rds")
saveRDS(physeq.clean.samples, file = "Bacteria_Soil_nonnorm_102021.rds")

# This is how to read in an RDS file
# just uncomment the line below to read it
physeq.css <- readRDS(file = "Bacteria_Soil_CSSnorm_102021.rds")
physeq.prop <- readRDS(file = "Bacteria_Soil_propnorm_102021.rds")
physeq.clean.samples <- readRDS(file = "Bacteria_Soil_nonnorm_102021.rds")


# Getting comfortable 
# most abundant OTUs
most.abundant <- data.frame(sort(rowSums(physeq.clean.samples@otu_table), decreasing = TRUE))
head(most.abundant)
physeq.clean.samples@refseq$FOTU_2 #find BOTU2 sequence
tax[tax$OTU == "FOTU_2",]

# Occupancy (i.e., number of samples the fungi were observed)
presence.absence <- transform_sample_counts(physeq.clean.samples, function(abund) 1*(abund>0))
presence.absence@otu_table
OTU.occupancy <- data.frame(sort(rowSums(presence.absence@otu_table)))
colnames(OTU.occupancy) <- "Occupancy"
OTU.occupancy$OTU <- rownames(OTU.occupancy)
OTU.occupancy$percent.occupancy <- (100*OTU.occupancy$Occupancy/max(OTU.occupancy$Occupancy))

# How many bacteria are present in every sample?
OTU.occupancy[OTU.occupancy$percent.occupancy >= 100,]

# ALPHA DIVERSITY
fungi_sperm@sam_data$shannon <- estimate_richness(fungi_sperm, measures=c("Shannon"))$Shannon
fungi_sperm@sam_data$invsimpson <- estimate_richness(fungi_sperm, measures=c("InvSimpson"))$InvSimpson
fungi_sperm@sam_data$richness <- estimate_richness(fungi_sperm, measures=c("Observed"))$Observed
fungi_sperm@sam_data$even <- fungi_sperm@sam_data$shannon/log(fungi_sperm@sam_data$richness)

sample.data.fungi <- data.frame(fungi_sperm@sam_data)

sample.data.fungi$Time.Point <- factor(sample.data.fungi$Time.Point, levels = c("0", "6", "12", "18", "E"))

# Richness over time
Soy.bull.richness <- sample.data.fungi %>%
  ggplot(aes(x = Time.Point, y = even, group = Crop, color = Crop, linetype = Crop)) +
  stat_summary(fun.y=mean,geom="line") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
  ylab("S") +
  #ggtitle("Fungi") +
  xlab("")+
  scale_color_manual(values=cbbPalette) +
  #scale_linetype_manual(values = c(rep("solid", 2), rep("dashed", 2))) +
  stat_compare_means(method = "anova") +
  theme_classic()

richness.time <- ggplot(sample.data.fungi, aes(x = Time.Point, y = even)) + 
  geom_boxplot() +
  geom_jitter() + 
  ylab("Richness") +
  stat_compare_means(method = "kruskal", hide.ns = TRUE) + 
  xlab("")+
  theme_classic() +
  facet_wrap(~Crop)

# Richness by treatment
richness.crop <- ggplot(sample.data.fungi, aes(x = Crop, y = even)) + 
  geom_boxplot() +
  geom_jitter() + 
  ylab("Richness") + 
  stat_compare_means(method = "kruskal") + 
  xlab("")+
  theme_classic() +
  facet_wrap(~Time.Point)
ggsave("Microbiome_Assignment2.png", dpi=300)

# Richness by treatment*collection interaction
richness.crop.time <- ggplot(sample.data.fungi, aes(x = Crop, y = richness, fill= Crop)) + 
  geom_boxplot() +
  geom_jitter() + 
  ylab("Richness") + 
  stat_compare_means(method = "kruskal") + 
  xlab("")+
  theme_classic() +
  facet_wrap(~Time.Point)

# Evenness diversity by management
even.Time.Point <- ggplot(sample.data.fungi, aes(x = Time.Point, y = even)) + 
  geom_boxplot() +
  geom_jitter() + 
  ylab("Evenness") + 
  stat_compare_means(method = "kruskal") + 
  xlab("")+
  theme_classic() 

# Shannon diversity by management
shannon.management <- ggplot(sample.data.fungi, aes(x = Crop, y = shannon)) + 
  geom_boxplot() +
  geom_jitter() + 
  ylab("Shannon") + 
  stat_compare_means(method = "kruskal") + 
  xlab("")+
  theme_classic() 

# BETA DIVERSITY

physeq.css.12 <- physeq.css %>%
  subset_samples(Time.Point %in% c("0", "18", "R"))

physeq.css@sam_data$Sample <- rownames(physeq.css@sam_data)

physeq.css.new <- physeq.css %>%
  subset_samples(Sample!= "S.18.4" & Time.Point != "E")


sample.data.fungi$Sample <- row.names(sample.data.fungi)



# Principle coordinates analysis with Bray-Curtis distances
ordination.pcoa <- ordinate(physeq.css, "PCoA", "bray") # calculate the resemblance and ordinate using PCoA
ordination.pcoa$vectors # positions of your points on the pcoa graph
ordination.pcoa$values #values to calculate the variance explained on each axis (dimension)
ordination.pcoa$points
pcoa <- plot_ordination(physeq.css, ordination = ordination.pcoa, type = "samples", color = "Time.Point") +
  theme_classic() + 
  scale_color_manual(values = cbbPalette)
pcoa

vectors <- data.frame(ordination.pcoa$vectors)

vectors$Sample <- row.names(vectors)
sample.data.fungi2 <- left_join(sample.data.fungi, vectors, by = "Sample")


ggplot(sample.data.fungi2[sample.data.fungi2$Time.Point == "18",], aes(x = even, y = Axis.1)) + 
  geom_point() +
  geom_smooth(method = "lm")

cor.test(sample.data.fungi2$Axis.1[sample.data.fungi2$Time.Point == "18"],
         sample.data.fungi2$even[sample.data.fungi2$Time.Point == "18"])


pcoa.data <- pcoa$data # taking the data to make a fancy plot

ggplot() + 
  geom_point(data = pcoa.data, aes(x = Axis.1, y = Axis.2, shape = Crop, fill = Time.Point), alpha = 0.8, size = 4) +
  theme_bw() +
  ylab("PCoA2 (13.7%)") + 
  xlab("PCoA1 (18.0%)") +
  scale_fill_manual(values=cbbPalette) +
  #stat_ellipse(data = global.pcoa.obj1.data, aes(x = Axis.1, y = Axis.2, group = Compartment), type = "norm", linetype = 2) +
  scale_shape_manual(values=c(21, 24, 23)) +
  theme(axis.text.x=element_text(size = 20),
        axis.text.y=element_text(size = 20), 
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15), 
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  guides(fill=guide_legend(override.aes=list(shape= 21)),
         shape=guide_legend(override.aes=list(fill= "black")))

ggsave("Microbiome_Assignment.png", dpi=300)

# PERMANOVA - testing for differences in centroids
prok.dist.bray = phyloseq::distance(physeq.css, "bray") # create bray-curtis distance matrix
adonis2(prok.dist.bray~Crop*Time.Point, as(sample_data(physeq.css), "data.frame")) #Are there significant changes?
# The R2 is the variation due to different factors 
# pvalue is the pvalue for the factor

# DISPERSION - testing for differences in variation due to management
beta.disp <- betadisper(prok.dist.bray, physeq.css@sam_data$Time.Point)
permutest(beta.disp)
TukeyHSD(beta.disp)
plot(beta.disp)
boxplot(beta.disp)
# Which groups are different from others?

# ANOSIM - Rank based differences in centroids
anosim(prok.dist.bray, physeq.css@sam_data$Time.Point) #Are there significant changes?
# Is differences in management due to differences in centroids (means) or differences in dispersion?

fungi_spermosphere_1000reads@sam_data$CropXTime.Point <- interaction(fungi_spermosphere_1000reads@sam_data$Crop, 
                                                                     fungi_spermosphere_1000reads@sam_data$Time.Point)
physeq.glom <- fungi_spermosphere_1000reads %>%
  tax_glom(taxrank = "Family") %>%  # agglomerate at phylum level
  merge_samples("CropXTime.Point") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Order)   # Sort data frame alphabetically by phylum

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "blue", "brown"
)

# Plot 
ggplot(physeq.glom [physeq.glom$Family != "unidentified", ], aes(x = Sample, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity") + #without this, all bars are plotted on each other
  scale_fill_manual(values = phylum_colors) +
  theme_classic() +
  # Remove x axis title
  theme(legend.position="bottom") + 
  ylab("Relative Abundance (Phyla > 1%)") +
  scale_y_continuous(labels = scales::percent) 

# Differential abundance analysis 
# Perform indicator species analysis just considering the original groups
#similar to differential abundance used in RNA-Seq
#unique taxa based on management
#t = transform
#css normalized reads are used in this type analysis so samples are comparative
#specify management systems to get their indicator species specifically
indicator.management <- indicspecies::multipatt(as.data.frame(t(physeq.css@otu_table)), cluster = physeq.css@sam_data$Crop)

# summary of results
indic <- summary(indicator.management, indvalcomp = TRUE)
indicators <- indicator.management$sign
indicators$OTU <- rownames(indicators)

convetional.indicators <- subset(indicators, index == 2 & p.value <= 0.01)

convetional.indicators2 <- left_join(convetional.indicators, tax, by = "OTU")

library("readr")
write_tsv(convetional.indicators2, "zotuindicators")

write_csv(convetional.indicators2, "indicators2")

tax[tax$OTU == "BOTU_4281",]
# Explore some of these taxa, what are they? What might they do?


