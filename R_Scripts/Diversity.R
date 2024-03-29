################Diversity Analysis############

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
library(picante)
library(codyn)

##### Set global options #####

# no scientific notation
options(scipen=10000) 

# color blind pallet used throughout 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ibm.cbb <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")
tol.cbb <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")

####### Read in Data #######

# using the non-normalized reads since spieceasi has its own normalizaiton methods
bac_sperm <- readRDS(file = "Bacteria/Bacteria_spermosphere_nonnorm_112922.rds")
bac.css.norm <- readRDS(file = "Bacteria/Bacteria_spermosphere_CSS_112922.rds")
fungi_sperm <- readRDS(file = "Fungi/Fungi_spermosphere_unedited_083022.rds")
fungi.css.norm <- readRDS(file = "Fungi/Fungi_CSSNorm_083022.rds")

######## PROKARYOTE ########
###### Beta-diversity using bray-curtis ########

##### All Samples; Bray #####
bac.bray = phyloseq::distance(bac.css.norm, "bray") # create bray-curtis distance matrix

# PERMANOVA - testing for differences in centroids
set.seed(12325)
adonis2(bac.bray~Crop*Time.Point, as(sample_data(bac.css.norm), "data.frame"))

global.ord.bray <- phyloseq::ordinate(bac.css.norm, "MDS", "bray")
variation.axis1.bray <- round(100*global.ord.bray$values$Relative_eig[[1]], 2) # % variation on axis one
variation.axis2.bray <- round(100*global.ord.bray$values$Relative_eig[[2]], 2) # % variation on axis two

bac.ggplot.data.bray <- data.frame(global.ord.bray$vectors) # get the point values
bac.ggplot.data.bray$Code <- rownames(bac.ggplot.data.bray) # add rownames
bac.ggplot.data.bray2 <- left_join(data.frame(bac.css.norm@sam_data), bac.ggplot.data.bray,  by = "Code") # add rest of metadata
bac.ggplot.data.bray2$Time.Point <- factor(bac.ggplot.data.bray2$Time.Point, levels = c("0", "6", "12", "18")) # order the time points in chronological order

global.bray <- ggplot(bac.ggplot.data.bray2, aes(x = Axis.1, y = Axis.2, fill = Crop, shape = Time.Point)) +
  geom_point(size=4, alpha = 0.7)+
  scale_shape_manual(values=c(21,22, 23, 24,25,26,27), name = "Hours post sowing")+
  scale_fill_manual(values=cbbPalette, name = "Environment", labels = c("Soil no seeds", "Cotton spermosphere", "Soybean spermosphere"))+
  xlab(paste("PcoA1 -", variation.axis1.bray, "%")) + 
  ylab(paste("PcoA2 -", variation.axis2.bray, "%")) +
  guides(fill=guide_legend(override.aes=list(shape= 21)),
         shape=guide_legend(override.aes=list(fill= "black")))+
  theme_bw() +
  ggtitle("Bray-Curtis")

# Separate analysis based on time-points #


#### 12 v. 18 hours
bac_sperm_12v18 <- bac.css.norm %>% 
  subset_samples(Time.Point %in% c("6", "12", "18")) %>% 
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

# bray curtis distance
bac.bray.12v18 = phyloseq::distance(bac_sperm_12v18, "bray") # create bray-curtis distance matrix

# permanova
set.seed(12325)
adonis.output.12v18  <- adonis2(bac.bray.12v18~Crop*Time.Point, as(sample_data(bac_sperm_12v18), "data.frame"))

# since no difference between 12 and 18 maybe combine for power in other analyses

### 0 hours; Bray ####

# filter to time point
bac_sperm_0 <- bac.css.norm %>% 
  subset_samples(Time.Point == "0") %>% 
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

# bray curtis distance
bac.bray.0 = phyloseq::distance(bac_sperm_0, "bray") # create bray-curtis distance matrix

# permanova
set.seed(12325)
adonis.output.0  <- adonis2(bac.bray.0~Crop, as(sample_data(bac_sperm_0), "data.frame"))

pvalue <- adonis.output.0$`Pr(>F)`[[1]]
R2 <- adonis.output.0$R2[[1]]


#Dispersion for each Time-Point
beta.disp.0 <- betadisper(bac.bray.0, bac_sperm_0@sam_data$Crop)

adonis.output.0
permutest(beta.disp.0)

sperm.0.ord <- phyloseq::ordinate(bac_sperm_0, "MDS", "bray")
sperm.0.ord.vectors <- data.frame(sperm.0.ord$vectors)
sperm.0.ord.vectors$Code <- rownames(sperm.0.ord.vectors)
sperm.0.ord.vectors2 <- left_join(data.frame(bac_sperm_0@sam_data), sperm.0.ord.vectors,  by = "Code")
sperm.0.ord.values <- data.frame(sperm.0.ord$values)


# percent variation
variation.axis1.0.bray <- round(100*sperm.0.ord.values$Relative_eig[[1]], 2)
variation.axis2.0.bray <- round(100*sperm.0.ord.values$Relative_eig[[2]], 2)

sperm.0.ord.plot <- ggplot(sperm.0.ord.vectors2, aes(x = Axis.1, y = Axis.2, color = Crop)) +
  geom_point(size=4, alpha = 0.7) +
  scale_color_manual(values=cbbPalette, name = "", labels = c("Soil no seeds", "Cotton spermosphere", "Soybean spermosphere")) +  
  xlab(paste("PcoA1 -", variation.axis1.0.bray, "%")) + 
  ylab(paste("PcoA2 -", variation.axis1.0.bray, "%")) +
  ggtitle(paste("0 hrs, P =", round(pvalue, 3))) +
  theme_bw() +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        legend.key = element_blank(),
        plot.title=element_text(size=10, face = "bold")) 

### 6 hours; Bray ####

# filter to time point
bac_sperm_6 <- bac.css.norm %>% 
  subset_samples(Time.Point == "6") %>% 
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

# bray curtis distance
bac.bray.6 = phyloseq::distance(bac_sperm_6, "bray") # create bray-curtis distance matrix

# permanova
set.seed(12325)
adonis.output.6  <- adonis2(bac.bray.6~Crop, as(sample_data(bac_sperm_6), "data.frame"))

pvalue <- adonis.output.6$`Pr(>F)`[[1]]
R2 <- adonis.output.6$R2[[1]]


#Dispersion for each Time-Point
beta.disp.6 <- betadisper(bac.bray.6, bac_sperm_6@sam_data$Crop)

adonis.output.6
permutest(beta.disp.6)

sperm.6.ord <- phyloseq::ordinate(bac_sperm_6, "MDS", "bray")
sperm.6.ord.vectors <- data.frame(sperm.6.ord$vectors)
sperm.6.ord.vectors$Code <- rownames(sperm.6.ord.vectors)
sperm.6.ord.vectors2 <- left_join(data.frame(bac_sperm_6@sam_data), sperm.6.ord.vectors,  by = "Code")
sperm.6.ord.values <- data.frame(sperm.6.ord$values)


# percent variation
variation.axis1.6.bray <- round(166*sperm.6.ord.values$Relative_eig[[1]], 2)
variation.axis2.6.bray <- round(166*sperm.6.ord.values$Relative_eig[[2]], 2)

sperm.6.ord.plot <- ggplot(sperm.6.ord.vectors2, aes(x = Axis.1, y = Axis.2, color = Crop)) +
  geom_point(size=4, alpha = 0.7) +
  scale_color_manual(values=cbbPalette, name = "", labels = c("Soil no seeds", "Cotton spermosphere", "Soybean spermosphere")) +
  xlab(paste("PcoA1 -", variation.axis1.6.bray, "%")) + 
  ylab(paste("PcoA2 -", variation.axis2.6.bray, "%")) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  ggtitle(paste("6 hrs, P <", round(pvalue, 3))) +
  theme_bw() +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        legend.key = element_blank(),
        plot.title=element_text(size=10, face = "bold"))
  


### 12 hours; Bray ####

# filter to time point
bac_sperm_12 <- bac.css.norm %>% 
  subset_samples(Time.Point == "12") %>% 
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

# bray curtis distance
bac.bray.12 = phyloseq::distance(bac_sperm_12, "bray") # create bray-curtis distance matrix

# permanova
set.seed(12325)
adonis.output.12  <- adonis2(bac.bray.12~Crop, as(sample_data(bac_sperm_12), "data.frame"))

pvalue <- adonis.output.12$`Pr(>F)`[[1]]
R2 <- adonis.output.12$R2[[1]]


#Dispersion for each Time-Point
beta.disp.12 <- betadisper(bac.bray.12, bac_sperm_12@sam_data$Crop)

adonis.output.12
permutest(beta.disp.12)

sperm.12.ord <- phyloseq::ordinate(bac_sperm_12, "MDS", "bray")
sperm.12.ord.vectors <- data.frame(sperm.12.ord$vectors)
sperm.12.ord.vectors$Code <- rownames(sperm.12.ord.vectors)
sperm.12.ord.vectors2 <- left_join(data.frame(bac_sperm_12@sam_data), sperm.12.ord.vectors,  by = "Code")
sperm.12.ord.values <- data.frame(sperm.12.ord$values)


# percent variation
variation.axis1.12.bray <- round(100*sperm.12.ord.values$Relative_eig[[1]], 2)
variation.axis2.12.bray <- round(100*sperm.12.ord.values$Relative_eig[[2]], 2)

sperm.12.ord.plot <- ggplot(sperm.12.ord.vectors2, aes(x = Axis.1, y = Axis.2, color = Crop)) +
  geom_point(size=4, alpha = 0.7) +
  scale_color_manual(values=cbbPalette, name = "", labels = c("Soil no seeds", "Cotton spermosphere", "Soybean spermosphere")) +
  xlab(paste("PcoA1 -", variation.axis1.12.bray, "%")) + 
  ylab(paste("PcoA2 -", variation.axis2.12.bray, "%")) +
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank()) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  ggtitle(paste("12 hrs, P <", round(pvalue, 3))) +
  theme_bw() +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        legend.key = element_blank(),
        plot.title=element_text(size=10, face = "bold"))


### 18 hours; Bray ####

# filter to time point
bac_sperm_18 <- bac.css.norm %>% 
  subset_samples(Time.Point == "18") %>% 
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

# bray curtis distance
bac.bray.18 = phyloseq::distance(bac_sperm_18, "bray") # create bray-curtis distance matrix

# permanova
set.seed(18325)
adonis.output.18  <- adonis2(bac.bray.18~Crop, as(sample_data(bac_sperm_18), "data.frame"))

pvalue <- adonis.output.18$`Pr(>F)`[[1]]
R2 <- adonis.output.18$R2[[1]]

#Dispersion for each Time-Point
beta.disp.18 <- betadisper(bac.bray.18, bac_sperm_18@sam_data$Crop)

adonis.output.18
permutest(beta.disp.18)

sperm.18.ord <- phyloseq::ordinate(bac_sperm_18, "MDS", "bray")
sperm.18.ord.vectors <- data.frame(sperm.18.ord$vectors)
sperm.18.ord.vectors$Code <- rownames(sperm.18.ord.vectors)
sperm.18.ord.vectors2 <- left_join(data.frame(bac_sperm_18@sam_data), sperm.18.ord.vectors,  by = "Code")
sperm.18.ord.values <- data.frame(sperm.18.ord$values)

# percent variation
variation.axis1.18.bray <- round(100*sperm.18.ord.values$Relative_eig[[1]], 2)
variation.axis2.18.bray <- round(100*sperm.18.ord.values$Relative_eig[[2]], 2)

sperm.18.ord.plot <- ggplot(sperm.18.ord.vectors2, aes(x = Axis.1, y = Axis.2, color = Crop)) +
  geom_point(size=4, alpha = 0.7) +
  scale_color_manual(values=cbbPalette, name = "", labels = c("Soil no seeds", "Cotton spermosphere", "Soybean spermosphere")) +
  xlab(paste("PcoA1 -", variation.axis1.18.bray, "%")) + 
  ylab(paste("PcoA2 -", variation.axis2.18.bray, "%")) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  ggtitle(paste(paste("18 hrs P <", round(pvalue, 3)))) +
  theme_bw() +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        legend.key = element_blank(),
        plot.title=element_text(size=10, face = "bold"))

###### Beta-Diversity using Weighted Unifrac ####

##### All Samples ########
bac.unifrac = UniFrac(bac.css.norm, weighted = T) # create weighted unifrac distance matrix

# PERMANOVA - testing for differences in centroids
set.seed(12325)
adonis2(bac.unifrac~Crop*Time.Point, as(sample_data(bac.css.norm), "data.frame"))

global.ord.uni <- phyloseq::ordinate(bac.css.norm, "MDS", "unifrac", weighted = TRUE)
variation.axis1 <- round(100*global.ord.uni$values$Relative_eig[[1]], 2)
variation.axis2 <- round(100*global.ord.uni$values$Relative_eig[[2]], 2)

bac.ggplot.data.uni <- data.frame(global.ord.uni$vectors)
bac.ggplot.data.uni$Code <- rownames(bac.ggplot.data.uni)
bac.ggplot.data.uni2 <- left_join(data.frame(bac.css.norm@sam_data), bac.ggplot.data.uni,  by = "Code")
bac.ggplot.data.uni2$Time.Point <- factor(bac.ggplot.data.uni2$Time.Point, levels = c("0", "6", "12", "18"))

global.uni <- ggplot(bac.ggplot.data.uni2, aes(x = Axis.1, y = Axis.2, fill = Crop, shape = Time.Point)) +
  geom_point(size=4, alpha = 0.7)+
  scale_shape_manual(values=c(21,22, 23, 24,25,26,27), name = "Hours post sowing")+
  scale_fill_manual(values=cbbPalette, name = "Environment", labels = c("Soil no seeds", "Cotton spermosphere", "Soybean spermosphere"))+
  xlab(paste("PcoA1 -", variation.axis1, "%")) + 
  ylab(paste("PcoA2 -", variation.axis2, "%")) +
  guides(fill=guide_legend(override.aes=list(shape= 21)),
         shape=guide_legend(override.aes=list(fill= "black")))+
  theme_bw() + 
  ggtitle("Weighted Unifrac")

### Supplemental Figure 2
Supplemental.figure2 <- ggpubr::ggarrange(global.uni, global.bray, common.legend = TRUE, labels = c("a", "b"), legend = "right")

### 0 hours; Unifrac #####
bac_sperm_0 <- bac.css.norm %>% 
  subset_samples(Time.Point == "0") %>% 
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

# PERMANOVA - testing for differences in centroids
prok.dist.uni.0hrs = UniFrac(bac_sperm_0, weighted = T) # create distance matrix
adonis.0.uni <- adonis2(prok.dist.uni.0hrs~Crop, as(sample_data(bac_sperm_0), "data.frame")) #Are there significant changes?
adonis.0.uni
# The R2 is the variation due to different factors 
# pvalue is the pvalue for the factor

pvalue <- adonis.0.uni$`Pr(>F)`[[1]]
R2 <- adonis.0.uni$R2[[1]]

# ploting
ordination.unifrac.pcoa.0hrs <- ordinate(bac_sperm_0, "PCoA", distance = prok.dist.uni.0hrs) # calculate the resemblance and ordinate using PCoA
uni.0.values <- data.frame(ordination.unifrac.pcoa.0hrs$values)

# percent variation
variation.axis1.0.uni <- round(100*uni.0.values$Relative_eig[[1]], 2)
variation.axis2.0.uni <- round(100*uni.0.values$Relative_eig[[2]], 2)

# data from ordinations
uni.0.vectors <- data.frame(ordination.unifrac.pcoa.0hrs$vectors) # positions of your points on the pcoa graph
uni.0.vectors$Code <- rownames(uni.0.vectors)
uni.0.vectors2 <- left_join(data.frame(bac_sperm_0@sam_data), uni.0.vectors,  by = "Code")


pcoa.unifrac.0hrs.new  <- ggplot(uni.0.vectors2, aes(x = Axis.1, y = Axis.2, color = Crop)) +
  geom_point(size=4, alpha = 0.7) +
  scale_color_manual(values=cbbPalette, name = "", labels = c("Soil no seeds", "Cotton spermosphere", "Soybean spermosphere")) +
  xlab(paste("PcoA1 -", variation.axis1.0.uni, "%")) + 
  ylab(paste("PcoA2 -", variation.axis2.0.uni, "%")) +
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank()) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  ggtitle(paste("0 hrs, P =", round(pvalue, 3))) +
  theme_bw() +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        legend.key = element_blank(),
        plot.title=element_text(size=10, face = "bold"))


#Dispersion for each Time-Point
beta.disp.0.uni <- betadisper(prok.dist.uni.0hrs, bac_sperm_0@sam_data$Crop)
permutest(beta.disp.0.uni)



### 6 hours; Unifrac #####
bac_sperm_6 <- bac.css.norm %>% 
  subset_samples(Time.Point == "6") %>% 
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

# PERMANOVA - testing for differences in centroids
prok.dist.uni.6hrs = UniFrac(bac_sperm_6, weighted = T) # create distance matrix
adonis.6.uni <- adonis2(prok.dist.uni.6hrs~Crop, as(sample_data(bac_sperm_6), "data.frame")) #Are there significant changes?
adonis.6.uni
# The R2 is the variation due to different factors 
# pvalue is the pvalue for the factor

pvalue <- adonis.6.uni$`Pr(>F)`[[1]]
R2 <- adonis.6.uni$R2[[1]]

# ploting
ordination.unifrac.pcoa.6hrs <- ordinate(bac_sperm_6, "PCoA", distance = prok.dist.uni.6hrs) # calculate the resemblance and ordinate using PCoA
uni.6.values <- data.frame(ordination.unifrac.pcoa.6hrs$values)

# percent variation
variation.axis1.6.uni <- round(100*uni.6.values$Relative_eig[[1]], 2)
variation.axis2.6.uni <- round(100*uni.6.values$Relative_eig[[2]], 2)

# data from ordinations
uni.6.vectors <- data.frame(ordination.unifrac.pcoa.6hrs$vectors) # positions of your points on the pcoa graph
uni.6.vectors$Code <- rownames(uni.6.vectors)
uni.6.vectors2 <- left_join(data.frame(bac_sperm_6@sam_data), uni.6.vectors,  by = "Code")


pcoa.unifrac.6hrs.new  <- ggplot(uni.6.vectors2, aes(x = Axis.1, y = Axis.2, color = Crop)) +
  geom_point(size=4, alpha = 0.7) +
  scale_color_manual(values=cbbPalette, name = "") +
  xlab(paste("PcoA1 -", variation.axis1.6.uni, "%")) + 
  ylab(paste("PcoA2 -", variation.axis2.6.uni, "%")) +
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank()) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  ggtitle(paste("6 hrs, P <", round(pvalue, 3))) +
  theme_bw() +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        legend.key = element_blank(),
        plot.title=element_text(size=10, face = "bold"))

  

#Dispersion for each Time-Point
beta.disp.6.uni <- betadisper(prok.dist.uni.6hrs, bac_sperm_6@sam_data$Crop)
permutest(beta.disp.6.uni)


### 12 hours; Unifrac #####
bac_sperm_12 <- bac.css.norm %>% 
  subset_samples(Time.Point == "12") %>% 
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

# PERMANOVA - testing for differences in centroids
prok.dist.uni.12hrs = UniFrac(bac_sperm_12, weighted = T) # create distance matrix
adonis.12.uni <- adonis2(prok.dist.uni.12hrs~Crop, as(sample_data(bac_sperm_12), "data.frame")) #Are there significant changes?
adonis.12.uni
# The R2 is the variation due to different factors 
# pvalue is the pvalue for the factor

pvalue <- adonis.12.uni$`Pr(>F)`[[1]]
R2 <- adonis.12.uni$R2[[1]]

# ploting
ordination.unifrac.pcoa.12hrs <- ordinate(bac_sperm_12, "PCoA", distance = prok.dist.uni.12hrs) # calculate the resemblance and ordinate using PCoA
uni.12.values <- data.frame(ordination.unifrac.pcoa.12hrs$values)

# percent variation
variation.axis1.12.uni <- round(100*uni.12.values$Relative_eig[[1]], 2)
variation.axis2.12.uni <- round(100*uni.12.values$Relative_eig[[2]], 2)

# data from ordinations
uni.12.vectors <- data.frame(ordination.unifrac.pcoa.12hrs$vectors) # positions of your points on the pcoa graph
uni.12.vectors$Code <- rownames(uni.12.vectors)
uni.12.vectors2 <- left_join(data.frame(bac_sperm_12@sam_data), uni.12.vectors,  by = "Code")


pcoa.unifrac.12hrs.new  <- ggplot(uni.12.vectors2, aes(x = Axis.1, y = Axis.2, color = Crop)) +
  geom_point(size=4, alpha = 0.7) +
  scale_color_manual(values=cbbPalette, name = "") +
  xlab(paste("PcoA1 -", variation.axis1.12.uni, "%")) + 
  ylab(paste("PcoA2 -", variation.axis2.12.uni, "%")) +
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank()) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  ggtitle(paste("12 hrs P <", round(pvalue, 3))) +
  theme_bw() +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        legend.key = element_blank(),
        plot.title=element_text(size=10, face = "bold"))


#Dispersion for each Time-Point
beta.disp.12.uni <- betadisper(prok.dist.uni.12hrs, bac_sperm_12@sam_data$Crop)
permutest(beta.disp.12.uni)


### 18 hours; Unifrac #####
bac_sperm_18 <- bac.css.norm %>% 
  subset_samples(Time.Point == "18") %>% 
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

# PERMANOVA - testing for differences in centroids
prok.dist.uni.18hrs = UniFrac(bac_sperm_18, weighted = T) # create distance matrix
adonis.18.uni <- adonis2(prok.dist.uni.18hrs~Crop, as(sample_data(bac_sperm_18), "data.frame")) #Are there significant changes?
adonis.18.uni
# The R2 is the variation due to different factors 
# pvalue is the pvalue for the factor

pvalue <- adonis.18.uni$`Pr(>F)`[[1]]
R2 <- adonis.18.uni$R2[[1]]

# ploting
ordination.unifrac.pcoa.18hrs <- ordinate(bac_sperm_18, "PCoA", distance = prok.dist.uni.18hrs) # calculate the resemblance and ordinate using PCoA
uni.18.values <- data.frame(ordination.unifrac.pcoa.18hrs$values)

# percent variation
variation.axis1.18.uni <- round(100*uni.18.values$Relative_eig[[1]], 2)
variation.axis2.18.uni <- round(100*uni.18.values$Relative_eig[[2]], 2)

# data from ordinations
uni.18.vectors <- data.frame(ordination.unifrac.pcoa.18hrs$vectors) # positions of your points on the pcoa graph
uni.18.vectors$Code <- rownames(uni.18.vectors)
uni.18.vectors2 <- left_join(data.frame(bac_sperm_18@sam_data), uni.18.vectors,  by = "Code")


pcoa.unifrac.18hrs.new  <- ggplot(uni.18.vectors2, aes(x = Axis.1, y = Axis.2, color = Crop)) +
  geom_point(size=4, alpha = 0.7) +
  scale_color_manual(values=cbbPalette, name = "") +
  xlab(paste("PcoA1 -", variation.axis1.18.uni, "%")) + 
  ylab(paste("PcoA2 -", variation.axis2.18.uni, "%")) +
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank()) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  ggtitle(paste("18 hrs P <", round(pvalue, 3))) +
  theme_bw() +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        legend.key = element_blank(),
        plot.title=element_text(size=10, face = "bold"))

#Dispersion for each Time-Point
beta.disp.18.uni <- betadisper(prok.dist.uni.18hrs, bac_sperm_18@sam_data$Crop)
permutest(beta.disp.18.uni)


##### Figure 1; Stick all together #####

figure1 <- ggpubr::ggarrange(pcoa.unifrac.0hrs.new, 
                                       pcoa.unifrac.6hrs.new,
                                       pcoa.unifrac.12hrs.new, 
                                       pcoa.unifrac.18hrs.new,
                                       sperm.0.ord.plot, 
                                       sperm.6.ord.plot,
                                       sperm.12.ord.plot, 
                                       sperm.18.ord.plot,
                                       labels = "auto",
                                       nrow = 2, ncol = 4, common.legend = T)

####### Prokaryote ALPHA DIVERSITY #########
bac_sperm@sam_data$shannon <- estimate_richness(bac_sperm, measures=c("Shannon"))$Shannon
bac_sperm@sam_data$invsimpson <- estimate_richness(bac_sperm, measures=c("InvSimpson"))$InvSimpson
bac_sperm@sam_data$simpson <- estimate_richness(bac_sperm, measures=c("Simpson"))$Simpson
bac_sperm@sam_data$richness <- estimate_richness(bac_sperm, measures=c("Observed"))$Observed
bac_sperm@sam_data$even <- bac_sperm@sam_data$shannon/log(bac_sperm@sam_data$richness)

sample.data.bac <- data.frame(bac_sperm@sam_data)

sample.data.bac$Time.Point <- factor(sample.data.bac$Time.Point, levels = c("0", "6", "12", "18"))

# Faith's phylogenetic diversity
otu <- bac_sperm@otu_table %>%
  as.data.frame()

phylo.diversity <- pd(t(otu), bac_sperm@phy_tree, include.root = T)
phylo.diversity$Samples <- rownames(phylo.diversity) 

phylo.diversity.join <- left_join(phylo.diversity, as.data.frame(bac_sperm@sam_data), by = c("Samples" = "Code"))
phylo.diversity.join$Time.Point <- factor(phylo.diversity.join$Time.Point, levels = c("0", "6", "12", "18"))
####### Richness over time; supplemental figure part 1 #######
compare_means(richness ~ Crop, sample.data.bac, group.by = "Time.Point")
compare_means(richness ~ Time.Point, sample.data.bac)

bac.richness <- sample.data.bac %>%
  ggplot(aes(x = Time.Point, y = richness, group = Crop, color = Crop)) +
  stat_summary(fun.y=mean,geom="line") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
  ylab("Richness") +
  xlab("Hours post sowing")+
  geom_jitter(width = 0.1, alpha = 0.8)+
  scale_color_manual(values=cbbPalette, name = "", labels = c("Soil no seeds", "Cotton spermosphere", "Soybean spermosphere")) +  
  #stat_compare_means(method = "kruskal", hide.ns = TRUE) +
  theme_classic()+
  theme(legend.position="none")

####### Faith's phylogenetic diversity over time; supplemental figure #######
compare_means(PD ~ Crop, phylo.diversity.join, group.by = "Time.Point")
compare_means(PD ~ Time.Point, phylo.diversity.join, method = "kruskal")

bac.phylo.div <- phylo.diversity.join %>%
  ggplot(aes(x = Time.Point, y = PD, group = Crop, color = Crop)) +
  stat_summary(fun.y=mean,geom="line") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
  ylab("Faith's PD") +
  xlab("Hours post sowing")+
  geom_jitter(width = 0.1, alpha = 0.8)+
  scale_color_manual(values=cbbPalette, name = "", labels = c("Soil no seeds", "Cotton spermosphere", "Soybean spermosphere")) +  
  #stat_compare_means(method = "kruskal", hide.ns = TRUE) +
  theme_classic()+
  theme(legend.position="none")

####### Bacterial evenness; Figure 2a ########
compare_means(even ~ Crop, sample.data.bac, group.by = "Time.Point")

bac.even <- sample.data.bac %>%
  ggplot(aes(x = Time.Point, y = even, color = Crop)) +
  geom_boxplot(position = position_dodge2(0.85, preserve = "single")) + 
  geom_point(position=position_jitterdodge(0.05)) +
  #stat_summary(fun.y=mean,geom="bar", position = "dodge", width = 0.5) +
  #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
  ylab("Pielou's evenness") +
  xlab("Hours post sowing")+
  scale_color_manual(values=cbbPalette, name = "", labels = c("Soil no seeds", "Cotton spermosphere", "Soybean spermosphere")) +  
  theme_classic() +
  theme(legend.position="none")
bac.even
####### Water Imbibition correlate with bacterial evenness; Figure 2b #####

water.imbibed <- ggplot(na.omit(sample.data.bac), aes(as.factor(Time.Point), 1000*Water_Imbibed, color = Crop)) + 
  geom_jitter(width = 0.5, alpha = 0.5) +
  stat_summary(fun =mean,geom="line", aes(group = Crop)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
  #stat_compare_means(aes(group = Crop, label = ..p.signif..)) + 
  xlab("Hours post sowing") +
  ylab("Water Imbibed (mg)") + 
  scale_color_manual(values=c(cbbPalette[[2]], cbbPalette[[3]]), name = "", labels = c("", "")) + 
  theme_classic() +
  theme(strip.background = element_blank(),
        legend.position="right") + 
  facet_wrap(~Crop, scales = "free")

# is there a statistical correlation? - yes
sample.data.bac.nosoil <- na.omit(sample.data.bac)
cotton.cor <- subset(sample.data.bac.nosoil, Crop == "Cotton ") 
soybean.cor <- subset(sample.data.bac.nosoil, Crop == "Soybean") 

cor.test(cotton.cor$Water_Imbibed, cotton.cor$even, method = "spearman")
cor.test(soybean.cor$Water_Imbibed, soybean.cor$even, method = "spearman")

##### Figure 2c #####
water.imbibed.cor <- ggplot(na.omit(sample.data.bac), aes(y = even, x = 1000*Water_Imbibed, color = Crop)) + 
  geom_point(aes(shape = Time.Point)) +
  geom_smooth(se = FALSE, method = lm) + 
  xlab("Water Imbibed (mg)") +
  ylab("Pielou's evenness") + 
  scale_color_manual(values=c(cbbPalette[[2]], cbbPalette[[3]]), name = "", labels = c("", "")) + 
  scale_shape_manual(values=c(15,16,17,18), name = "", labels = c("0 hrs", "6 hrs", "12 hrs", "18 hrs")) +
  theme_classic() +
  guides(color="none") + 
  theme(strip.background = element_blank(),
        legend.position="right") + 
  facet_wrap(~Crop, scales = "free")

#### Figure 2; significance levels added with Adobe or powerpoint #### 

figure2 <- ggpubr::ggarrange(water.imbibed,
                             bac.even, 
                             water.imbibed.cor, 
                             labels = "auto",
                             nrow = 3, ncol = 1, legend = F)

Supplemental.figure3 <- ggpubr::ggarrange(bac.richness, 
                  bac.phylo.div,
                  labels = "auto",
                  nrow = 1, ncol = 2, common.legend = T)

######## FUNGI ######
##### All Samples; Bray #####
fungi.bray = phyloseq::distance(fungi.css.norm, "bray") # create bray-curtis distance matrix

# PERMANOVA - testing for differences in centroids
set.seed(12325)
adonis2(fungi.bray~Crop*Time.Point, as(sample_data(fungi.css.norm), "data.frame"))

global.ord.bray <- phyloseq::ordinate(fungi.css.norm, "MDS", "bray")
variation.axis1.bray <- round(100*global.ord.bray$values$Relative_eig[[1]], 2) # % variation on axis one
variation.axis2.bray <- round(100*global.ord.bray$values$Relative_eig[[2]], 2) # % variation on axis two

fungi.ggplot.data.bray <- data.frame(global.ord.bray$vectors) # get the point values
fungi.ggplot.data.bray$Code <- rownames(fungi.ggplot.data.bray) # add rownames
fungi.ggplot.data.bray2 <- left_join(data.frame(fungi.css.norm@sam_data), fungi.ggplot.data.bray,  by = "Code") # add rest of metadata
fungi.ggplot.data.bray2$Time.Point <- factor(fungi.ggplot.data.bray2$Time.Point, levels = c("0", "6", "12", "18")) # order the time points in chronological order

global.bray.fungi <- ggplot(fungi.ggplot.data.bray2, aes(x = Axis.1, y = Axis.2, fill = Crop, shape = Time.Point)) +
  geom_point(size=4, alpha = 0.7)+
  scale_shape_manual(values=c(21,22, 23, 24,25,26,27), name = "Hours post sowing")+
  scale_fill_manual(values=cbbPalette, name = "Environment", labels = c("Soil no seeds", "Cotton spermosphere", "Soybean spermosphere"))+
  xlab(paste("PcoA1 -", variation.axis1.bray, "%")) + 
  ylab(paste("PcoA2 -", variation.axis2.bray, "%")) +
  guides(fill=guide_legend(override.aes=list(shape= 21)),
         shape=guide_legend(override.aes=list(fill= "black")))+
  theme_bw() +
  ggtitle("Bray-Curtis")


####### Fungal ALPHA DIVERSITY #########
fungi_sperm@sam_data$shannon <- estimate_richness(fungi_sperm, measures=c("Shannon"))$Shannon
fungi_sperm@sam_data$invsimpson <- estimate_richness(fungi_sperm, measures=c("InvSimpson"))$InvSimpson
fungi_sperm@sam_data$simpson <- estimate_richness(fungi_sperm, measures=c("Simpson"))$Simpson
fungi_sperm@sam_data$richness <- estimate_richness(fungi_sperm, measures=c("Observed"))$Observed
fungi_sperm@sam_data$even <- fungi_sperm@sam_data$shannon/log(fungi_sperm@sam_data$richness)

sample.data.fungi <- data.frame(fungi_sperm@sam_data)

sample.data.fungi$Time.Point <- factor(sample.data.fungi$Time.Point, levels = c("0", "6", "12", "18"))

####### Richness over time #######
fungi.richness <- sample.data.fungi %>%
  ggplot(aes(x = Time.Point, y = richness, group = Crop, color = Crop)) +
  stat_summary(fun=mean,geom="line") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
  ylab("Richness") +
  xlab("Hours post sowing")+
  geom_jitter(width = 0.1, alpha = 0.8)+
  scale_color_manual(values=cbbPalette, name = "", labels = c("Soil no seeds", "Cotton spermosphere", "Soybean spermosphere")) +  
  #stat_compare_means(method = "kruskal", hide.ns = TRUE) +
  theme_classic()+
  theme(legend.position="none")
fungi.richness

####### Bacterial evenness; Figure 2a ########
fungi.even <- sample.data.fungi %>%
  ggplot(aes(x = Time.Point, y = even, group = Crop, color = Crop)) +
  stat_summary(fun.y=mean,geom="line") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
  ylab("Pielou's evenness") +
  xlab("Hours post sowing")+
  geom_jitter(width = 0.1, alpha = 0.8)+
  scale_color_manual(values=cbbPalette, name = "", labels = c("Soil no seeds", "Cotton spermosphere", "Soybean spermosphere")) +  
  #stat_compare_means(method = "kruskal", hide.ns = TRUE) +
  theme_classic()+
  theme(legend.position="none")
fungi.even

###### Supplemental fungal alpha diversity figure ####
Supplemental.figure4 <- ggpubr::ggarrange(fungi.richness, 
                                          fungi.even,
                                          labels = "auto",
                                          nrow = 1, ncol = 2, common.legend = T)

##### Community Composition of top 20 most abundant OTUs ######

##### Prokaryote #####
topx <- top_taxa(bac_sperm, n = 20)

bac.composition <- bac_sperm %>%
  subset_taxa(OTU %in% topx) %>%
  microbiome::transform("compositional") %>%
  psmelt() %>% 
  group_by(Crop, Time.Point, Label) %>%
  mutate(Time.Point = factor(Time.Point, levels = c("0", "6", "12", "18"))) %>%
  summarise(MeanRelAbund = mean(Abundance)) %>%
  arrange(-MeanRelAbund) %>%
  left_join(as.data.frame(tax_table(bac_sperm), by = "Label")) %>% 
  ggplot(aes(Crop, MeanRelAbund, fill = Label)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_x_discrete(labels=c('Soil no seed', 'Cotton Sp.', 'Soybean sp.'))+
  scale_fill_manual(values= c(cbbPalette, ibm.cbb, tol.cbb)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Relative abundance (%)",
       title = "Prokaryote") + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.text = element_text(face = "italic", size = 5),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, 'cm')) +
  facet_wrap(~Time.Point, nrow = 1)
bac.composition



##### Fungi #####
topx.fungi <- top_taxa(fungi_sperm, n = 20)

fung.composition <- fungi_sperm %>%
  subset_taxa(OTU %in% topx.fungi) %>%
  microbiome::transform("compositional") %>%
  psmelt() %>% 
  group_by(Crop, Time.Point, Label) %>%
  mutate(Time.Point = factor(Time.Point, levels = c("0", "6", "12", "18"))) %>%
  summarise(MeanRelAbund = mean(Abundance)) %>%
  arrange(-MeanRelAbund) %>%
  left_join(as.data.frame(tax_table(fungi_sperm), by = "Label")) %>% 
  ggplot(aes(Crop, MeanRelAbund, fill = Label)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_x_discrete(labels=c('Soil no seed', 'Cotton Sp.', 'Soybean sp.'))+
  scale_fill_manual(values= c(cbbPalette, ibm.cbb, tol.cbb)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Relative abundance (%)",
       title = "Fungi") + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.text = element_text(face = "italic", size = 5),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, 'cm')) +
  facet_wrap(~Time.Point, nrow = 1)
fung.composition

##### Compositional Figure ####
compositional.fig <- ggpubr::ggarrange(bac.composition, 
                                       fung.composition,
                             labels = "auto",
                             nrow = 2, ncol = 1, align = "v")
