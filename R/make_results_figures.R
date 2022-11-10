# library(vegan)
# library(ape)
# library(iNEXT)
# library(picante)
# library(MuMIn)
# library(UpSetR)
# library(lmodel2)
# library(betapart)
# library(fuzzySim)
# library(viridis)
# library(RColorBrewer)
# library(kableExtra)
# library(VennDiagram)
# library(ComplexHeatmap)
# library(estimateR)
# library(knitr)
# library(ade4)
# library(ggfortify)
# library(ggplot2)
# library(viridis)

library(tidyverse)
library(patchwork)
library(vegan)
library(ggthemes)

theme_set(theme_bw(base_size = 12))

# -----------------------------------------------------------------------------------------------------------------
# Data import
# -----------------------------------------------------------------------------------------------------------------
species_descriptions <- read.csv("data-raw/ARCHIVED_DATA/TEMP_papers_table_1980_on_2022-11-05.csv") %>% 
      pivot_longer(cols = 6:8, names_to = "var", values_to = "number")
phyla_overview <- read.csv("data-raw/ARCHIVED_DATA/TEMP_SUMMARY_FIG2_ALL_PHYLA_2022-11-05.csv") %>% 
      group_by(Phylum) %>% 
      mutate(total_phylum = sum(Total)) %>% 
      group_by(Data) %>% 
      mutate(perc = paste0(round((Total/total_phylum)*100), "%"))
phyla_overview$ypos <- ifelse(phyla_overview$Data == "named species", phyla_overview$total_phylum, phyla_overview$total_phylum + 70)

community_matrix_CCZ <- read.table('data-processed/community_matrix_CCZ.txt')
specaccum_sites_df <- read.csv('data-processed/CCZ_specaccum_sites.csv')
CCZ_rarecurve <- read.csv("data-processed/CCZ_rarecurve.csv")

Hills_q_CCZ_df <- read.csv('data-processed/Hills_q_CCZ_df.csv') %>% 
      filter(size_based.Method != "Observed") %>% 
      mutate(size_based.Method = factor(size_based.Method, levels = c("Rarefaction", "Extrapolation")))

Hills_q0_inc_df <- read.csv('data-processed/Hills_q_0_inc_df.csv') %>% 
      filter(size_based.Method != "Observed") %>% 
      mutate(size_based.Method = factor(size_based.Method, levels = c("Rarefaction", "Extrapolation")))

taxon_rank_data <- read.csv("data-raw/ARCHIVED_DATA/temp_log_v2_2022-11-06.csv")

data_regions <- read.csv('data-raw/CCZ_ALL_SPP_DATA_V2_REGION_2022-11-08.csv') %>% 
      drop_na()
data_regions <- data_regions[nzchar(data_regions$Site), ]

test <- read.csv('data-raw/CCZ_ALL_TAXA_v2_2022-11-08.csv')
load('data-processed/CCZ_com_matrix_standardised.RData') 

# -----------------------------------------------------------------------------------------------------------------
# Figure 2: Raw data
# -----------------------------------------------------------------------------------------------------------------
descriptions_figure <- ggplot(species_descriptions, aes(x = Year, y = number, colour = var, fill = var)) +
      geom_line(size = 1) +
      labs(x = "Year", y = "Cumulative totals") +
      theme(legend.justification = c(-0.35, 1), 
            legend.position = c(0, 1),
            legend.box.margin=margin(c(20, 20, 20, 20)),
            legend.title = element_blank(),
            legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
      scale_color_manual(breaks = c("cumul_desc", "cumul_spp", "cumul_pubs"),
                         labels = c("All descriptions", "New species", "Publications"),
                         values = c('black', 'steelblue', 'coral2'))

phyla_figure <- ggplot(phyla_overview, aes(Phylum, Total, fill = Data)) + 
      geom_bar(stat="identity", position="stack", color="black") + 
      theme(legend.justification = c(-0.35, 1), 
            legend.position = c(0, 1),
            legend.box.margin=margin(c(20, 20, 20, 20)),
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
      geom_text(aes(label = perc, y = ypos, col = Data), vjust = -0.6,
                position = ggstance::position_dodgev(height = 1),
                size = 3.5, fontface='bold') +
      scale_fill_manual(values = c('steelblue', 'coral2'),
                        labels = c('Unnamed taxa', 'Named species')) +
      scale_color_manual(values = c('steelblue', 'coral2'), guide = "none") +
      ylim(0, 2000); phyla_figure

figure_2 <- descriptions_figure / phyla_figure + plot_annotation(tag_levels = 'A')

ggsave(figure_2,
       filename = 'output-figures/figure_2.tiff', 
       width = 8.5, height = 10, units = 'in', dpi = 150)
ggsave(descriptions_figure,
       filename = 'output-figures/descriptions_figure.tiff', 
       width = 8, height = 6, units = 'in', dpi = 150)
ggsave(phyla_figure,
       filename = 'output-figures/phyla_figure.tiff', 
       width = 8, height = 6, units = 'in', dpi = 150)

# -----------------------------------------------------------------------------------------------------------------
# Figure 3: Species accumulation curves for CCZ
# -----------------------------------------------------------------------------------------------------------------
A_Chao1 <- ggplot() +
      geom_line(data = Hills_q_CCZ_df, aes(x = size_based.m, y = size_based.qD, lty = size_based.Method),
                col = "coral2", cex = 1) +
      # geom_line(data = CCZ_rarecurve, aes(Individuals, Species), 
      #           cex = 1, col = "green") +
      theme(legend.justification = c(0, 1), 
            legend.position = c(.5, .4),
            legend.box.margin=margin(c(20, 20, 20, 20)),
            legend.margin = margin(1, 1, 1, 1),
            legend.spacing.x = unit(0, "mm"),
            legend.spacing.y = unit(0, "mm"),
            legend.title = element_blank(),
            legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
      xlab("Individuals") +
      ylab("Species Richness") +
      scale_y_continuous(breaks = seq(0, 5000, 1000)); A_Chao1

B_rarefaction_CCZ <- ggplot(CCZ_rarecurve) +
      geom_line(aes(Individuals, Species), cex = 1, col = "coral2") +
      xlab("Individuals") +
      ylab("Species Richness"); B_rarefaction_CCZ

C_Chao2 <- ggplot() +
      # geom_ribbon(data = specaccum_sites_df, aes(sites, ymin = richness - 1.96*sd, ymax = richness + 1.96*sd), alpha = .3, fill = 'blue') +
      # geom_line(data = specaccum_sites_df, aes(sites, richness), col = "blue") +
      geom_ribbon(data = Hills_q0_inc_df, aes(x = size_based.t, ymin = size_based.qD.LCL, ymax = size_based.qD.UCL), alpha = 0.3, show.legend = FALSE) +      
      geom_line(data = Hills_q0_inc_df, aes(x = size_based.t, y = size_based.qD, lty = size_based.Method), col = "coral2", cex = 1) +
      theme(legend.justification = c(0, 1), 
            legend.position = c(.5, .4),
            legend.box.margin=margin(c(20, 20, 20, 20)),
            legend.margin = margin(1, 1, 1, 1),
            legend.spacing.x = unit(0, "mm"),
            legend.spacing.y = unit(0, "mm"),
            legend.title = element_blank(),
            legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
      xlab("Sites") +
      ylab("Species Richnessy") +
      scale_y_continuous(breaks = seq(0, 5000, 1000)); C_Chao2

D_species_accum <- ggplot(specaccum_sites_df) +
      geom_ribbon(aes(sites, ymin = richness - 1.96*sd, ymax = richness + 1.96*sd), alpha = .3, fill = 'blue') +
      geom_line(aes(sites, richness)) +
      xlab("Sites") +
      ylab("Species Richness"); D_species_accum


figure_3 <- (A_Chao1 | B_rarefaction_CCZ) / (C_Chao2 | D_species_accum) + plot_annotation(tag_levels = 'A')

ggsave(figure_3,
       filename = 'output-figures/figure_3.tiff', 
       width = 8.5, height = 7, units = 'in', dpi = 150)



 # -----------------------------------------------------------------------------------------------------------------
# Higher taxon richness extrapolation
# -----------------------------------------------------------------------------------------------------------------
mod <- lm(log ~ order, data = taxon_rank_data)
taxon_rank_data <- rbind(taxon_rank_data, data.frame(order = 6, log = predict(mod, newdata = data.frame(order = 6), type = "response")))
taxon_rank_data$ID <- c(rep("obs", 5), "pred")
taxon_rank_data$number <- exp(taxon_rank_data$log)

C_taxon_rank <- ggplot(taxon_rank_data, aes(x = order, y = log, col = ID)) +
      geom_point(cex = 3) +
      geom_smooth(method = "lm", se = F, col = "coral2", lty = 2, alpha = 0.5) +
      ylab("Log total") +
      xlab("") +
      theme(legend.title = element_blank()) +
      scale_x_continuous(breaks = 1:6, labels = c("Phyla", "Class", "Order", "Family", "Genus", "Species")) +
      scale_color_manual(values = c("black", "coral2"), labels = c("observed", "predicted")); C_taxon_rank

ggsave(C_taxon_rank,
       filename = 'output-figures/taxa_int.tiff', 
       width = 8, height = 5, units = 'in', dpi = 150)

# -----------------------------------------------------------------------------------------------------------------
# Supplement: NMDS plot for regions
# -----------------------------------------------------------------------------------------------------------------
all_taxa_data <- read.csv('data-raw/CCZ_ALL_TAXA_v2_2022-11-08.csv') %>% 
      mutate(max_depth.start_depth = as.numeric(max_depth.start_depth),
             count = as.integer(COUNT)) %>% 
      drop_na()
depth_quantile <- quantile(all_taxa_data$max_depth.start_depth)
all_taxa_data$depth_quantile <- NA
all_taxa_data$depth_quantile[all_taxa_data$max_depth.start_depth <= depth_quantile[1]] <- "0-1667"
all_taxa_data$depth_quantile[all_taxa_data$max_depth.start_depth > depth_quantile[1] & all_taxa_data$max_depth.start_depth <= depth_quantile[2]] <- "1668-4139"
all_taxa_data$depth_quantile[all_taxa_data$max_depth.start_depth > depth_quantile[2] & all_taxa_data$max_depth.start_depth <= depth_quantile[3]] <- "4140-4449"
all_taxa_data$depth_quantile[all_taxa_data$max_depth.start_depth > depth_quantile[3] & all_taxa_data$max_depth.start_depth <= depth_quantile[4]] <- "4450-4888"
all_taxa_data$depth_quantile[all_taxa_data$max_depth.start_depth > depth_quantile[4] ] <- ">5850"

all_taxa_data <- all_taxa_data %>% 
      dplyr::select(SITE_STRING_V1, region, depth_quantile, SPECIES, count) %>% 
      rename(site = SITE_STRING_V1, species = SPECIES)
all_taxa_data <- all_taxa_data[nzchar(all_taxa_data$site), ]
all_taxa_data <- all_taxa_data[nzchar(all_taxa_data$species), ]

samplePoints.sp <- all_taxa_data %>% 
      dplyr::group_by(region, site, depth_quantile, species) %>% 
      dplyr::summarise(cover = sum(count)) %>%  
      reshape::cast(.,  region + site + depth_quantile ~ species, value = "cover")
samplePoints.sp[is.na(samplePoints.sp)] <-  0
rownames(samplePoints.sp) <- with(samplePoints.sp, paste(region, site, depth_quantile, sep="_"))


# -----------------------------------------------------------------------------------------------------------------
# Supplement: Heatmap of effort
# -----------------------------------------------------------------------------------------------------------------







# -----------------------------------------------------------------------------------------------------------------
# Mu's hodge podge code corner of doom
# -----------------------------------------------------------------------------------------------------------------


## FIG 3 SPECIES RICHNESS ESTIMATES + CURVES 

rm(list=ls())

library(vegan)
library(ape)
library(iNEXT)
library(picante)
library(MuMIn)
library(UpSetR)
library(lmodel2)
library(betapart)
library(fuzzySim)
library(viridis)
library(RColorBrewer)
library(kableExtra)
library(VennDiagram)
library(ComplexHeatmap)
library(estimateR)
library(knitr)
library(ade4)
library(ggfortify)
library(ggplot2)
library(spadeR)
library(tidyverse)
library(reshape2)
library(viridis)
library(ggthemes)

#########################
## FIG 2

## 2B CUMULATIVE SPECIES DESCRIPTIONS

data <- read.csv("CCZ_LITERATURE_PAPERS_2022-11-05.csv")
dim(data)
names(data)
## summary by year

data3 <- tapply(data[, "reference"], data[, "year_published"],
                function(x) { length(unique(x))})
write.csv(data3, "OUTPUT.csv")

data <- read.csv("TEMP_papers_table_1980_on_2022-11-05.csv")

## cumulative totals - run with and without legend
p1 <- ggplot(data = data, aes(x = Year)) + theme_bw() +
   geom_line(aes(y = cumul_desc, colour = "cumul_desc"), size = 1) +
   geom_line(aes(y = cumul_spp, colour = "cumul_spp"), size = 1) +
   geom_line(aes(y = cumul_pubs, colour = "cumul_pubs"), size = 1) +
   labs(x = "Year", y = "Cumulative totals") +
   #  ggtitle("Cumulative totals of publications and descriptions (year 2000- 2021)") +
   theme(text=element_text(size=14), #change font size of all text
         axis.text=element_text(size=14), #change font size of axis text
         axis.title=element_text(size=14), #change font size of axis titles
         plot.title=element_text(size=14), #change font size of plot title
         legend.text=element_text(size=14), #change font size of legend text
         legend.title=element_text(size=14)) + 
   theme(legend.position = "none") +
   scale_colour_manual("", 
                       breaks = c("cumul_desc", "cumul_spp", "cumul_pubs"),
                       values = c("black", "coral2", "steelblue"),
                       labels = c("all descriptions", "new species", " publications"))

#############################

## FIG 2B TOTAL SPECIES BY PHYLA

library(tidyverse)

data <- read.csv("CCZ_CHECKLIST_2022-11-06.csv")

## filter all named spp 
data2 <- data %>%
   filter(SPP_CHECKLIST == "yes")
dim(data2) # 647

## filter named spp no qual
data2 <- data %>%
   filter(short_spp_cl == "yes")
dim(data2) # 424
length(unique(data2$scientificName))

## whole checklist
data3 <- tapply(data[, "scientificName"], data[, "taxonRank"],
                function(x) { length(unique(x))})

## spp only
data3 <- tapply(data2[, "scientificName"], data2[, "Phylum"],
                function(x) { length(unique(x))})

data3 <- tapply(data2[, "scientificName"], data2[, "size_cat"],
                function(x) { length(unique(x))})

data3 <- tapply(data2[, "scientificName"], data2[, "habitat_interim"],
                function(x) { length(unique(x))})

write.csv(data3, "OUTPUT.csv")

###############################################

## MORPHOSPP
data <- read.csv("CCZ_MSPP_CHECKLIST_2022-11-04.csv")
length(unique(data$taxonConceptID))

## spp only
data3 <- tapply(data[, "taxonConceptID"], data[, "Phylum"],
                function(x) { length(unique(x))})

data3 <- tapply(data[, "taxonConceptID"], data[, "size_cat"],
                function(x) { length(unique(x))})

data3 <- tapply(data[, "taxonConceptID"], data[, "habitat_interim"],
                function(x) { length(unique(x))})


write.csv(data3, "OUTPUT.csv")

CHECK_SUM <- read.csv("TEMP_SUMMARY_CHECK-SHORT_PHYLUM_2022-09-21.csv")
CHECK_MSPP <- read.csv("TEMP_SUMMARY_CHECK-SHORT_PHYLUM_2022-09-21.csv")

## COMBINED INTO ONE FILE - ALL PHYLA VERSION
SUM <- read.csv("TEMP_SUMMARY_FIG2_ALL_PHYLA_2022-11-05.csv")

## SELECTED PHYLA
SUM <- read.csv("TEMP_SUMMARY_FIG2_2022-09-21.csv")

## border around bands
ggplot(SUM, aes(Phylum, Total, fill = Data)) + 
   geom_bar(stat="identity", position="dodge", color="black") + coord_flip()

## FINAL PLOT ## RUN WITH AND WITHOUT LEGEND
ggplot(SUM, aes(Phylum, Total, fill = Data)) + theme_bw() +
   geom_bar(stat="identity", position="dodge", color="black") + coord_flip() +
   theme(text=element_text(size=12), #change font size of all text
         axis.text=element_text(size=12), #change font size of axis text
         axis.title=element_text(size=12), #change font size of axis titles
         plot.title=element_text(size=12), #change font size of plot title
         legend.text=element_text(size=12), #change font size of legend text
         legend.title=element_text(size=12)) + coord_flip() + 
   theme(legend.position = "none") +
   scale_fill_manual('Position', values=c('steelblue', 'pink'))


## fig 2 ## stacked bar
p1 <- ggplot(SUM, aes(Phylum, Total, fill = Data)) + theme_bw() +
   geom_bar(stat="identity", position="stack", color="black") + 
   theme(text=element_text(size=14), #change font size of all text
         axis.text=element_text(size=14), #change font size of axis text
         axis.title=element_text(size=14), #change font size of axis titles
         plot.title=element_text(size=14), #change font size of plot title
         legend.text=element_text(size=14), #change font size of legend text
         legend.title=element_text(size=14)) + coord_flip() + 
   theme(legend.position = "none") +
   scale_fill_manual('Position', values=c('steelblue', 'pink'))


## fig 2 vertical

p1 <- ggplot(SUM, aes(Phylum, Total, fill = Data)) + theme_bw() +
   geom_bar(stat="identity", position="stack", color="black") + 
   theme(text=element_text(size=14), #change font size of all text
         axis.text=element_text(size=14), #change font size of axis text
         axis.title=element_text(size=14), #change font size of axis titles
         plot.title=element_text(size=14), #change font size of plot title
         legend.text=element_text(size=14), #change font size of legend text
         legend.title=element_text(size=14)) + #coord_flip() + 
   theme(legend.position = "none") +
   scale_fill_manual('Position', values=c('steelblue', 'pink')) +
   theme(axis.title.x=element_blank(),
         axis.text.x = element_text(angle = 55, hjust = 1))

#########################

## FIG 3 SPECIES RICHNESS ESTIMATES + CURVES 

rm(list=ls())


##########################
## READ IN DATA
## all data

## chao 2 in specpool
data <- read.csv("temp_spp_v1_abundance_2022-11-08.csv") ## 4410 SPP 
data <- read.csv("temp_spp_v1_presence_2022-11-08.csv") ## 4410 SPP 

data <- read.csv("temp_fam_v1_presence_2022-11-08.csv") ## 4410 SPP 


## by size class
## by region

## 8TH NOV- LIT DATA ONLY- NO DUPLIC ACROSS SAMPLING UNITS ACROSS DATASETS- all zero rows incl
#data <- (read.csv("TEMP_ALL_SPP_LIT_ONLY_NO_ZERO_2022-11-08.csv",header=T,sep=",",fileEncoding="latin1")[-2,])

head(data)
dim(data) ## 39116 ## no nas- 36 842

## MAKE SPP MATRIX

data.matrix <- sample2matrix(data)
data.matrix[1:5, 1:5]

length(unique(data$Species[data$Abundance > 0])) ## 4372

## CLEANUP

## length(names(apply(data.matrix, is.numeric))) == length(colnames(data.matrix))
## data.matrix$Abundance <- as.numeric(data.matrix$Abundance)

## SPECIES ACCUMULATION CURVE, SPEC POOOL- VEGAN

data3<- specpool(data.matrix)
write.csv(data3, "OUTPUT.csv")

## DO SPP ACCUM CURVE

data.curve <- specaccum(data.matrix, method = "random", permutations = 1000)

# Look at the results
data.curve

# Make a new dataframe
data.curve.df <- data.frame(Sites = data.curve$sites,
                            Richness = data.curve$richness,
                            sd = data.curve$sd)
head(data.curve.df)
tail(data.curve.df)

ggplot(data.curve.df, aes(x = Sites, y = Richness)) +
   # Add line
   geom_line() +
   # Add confidence intervals
   geom_ribbon(aes(ymin = Richness - 1.96*sd, ymax = Richness + 1.96*sd), 
               alpha = 0.2, colour = "lightgrey") +
   theme(text=element_text(size=14), #change font size of all text
         axis.text=element_text(size=20), #change font size of axis text
         axis.title=element_text(size=20), #change font size of axis titles
         plot.title=element_text(size=20), #change font size of plot title
         legend.text=element_text(size=20), #change font size of legend text
         legend.title=element_text(size=20)) + theme_bw() #
labs(title="Species accumulation curve- DeepData (named species and morphospecies)")


############################

## RAREFACTIION CURVE FOR CCZ ONLY

data <- read.csv("temp_spp_v1_abundance_2022-11-08.csv") %>% mutate(Site = "CCZ")
data <- read.csv("temp_fam_v1_abundance_2022-11-08.csv") %>% mutate(Site = "CCZ")

data.matrix <- sample2matrix(data)
data.matrix[1:5, 1:5]

species_no <- specnumber(data.matrix)
raremax <- min(rowSums(data.matrix))
Srare <- rarefy(data.matrix, raremax)
plot(species_no, Srare, xlab = "Observed No. of Species", ylab = "Species")
abline(0, 1)

rarecurve(data.matrix, step = 20, sample = raremax, cex = 1, col = viridis(n=4)[-4], 
          lwd = 2, lty = c("solid", rep("dashed", 2)), label = F,xlab = "No. of individuals", 
          ylab = "No. of Species", cex.axis = 1)
# legend("bottomright", legend = rownames(com_mat), col = viridis(n=4)[-4], pch=19, cex = 0.5, y.intersp = 0.8, pt.cex = 1)
rarecurve(data.matrix, step = 20)

## estimateR for CCZ only data

test <- estimateR(community_matrix_CCZ) 


############################

## CHAO1 AND CHAO2 IN INEXT- 

#### CHAO1 
data <- read.csv("temp_spp_v1_abundance_2022-11-08.csv") %>% mutate(Site = "CCZ")

data.matrix <- sample2matrix(data)
data.matrix[1:5, 1:5]

## transpose
com_matT <- t(data.matrix)

Hills_q <- iNEXT(com_matT, q=0, datatype = "abundance", nboot = 2)

test <- ChaoRichness((com_matT))

write.csv(test, "OUTPUT.csv")

ggiNEXT(Hills_q, type=1)

plot_q <- ggiNEXT(Hills_q, type=1)+theme_clean()


Hills_q <- iNEXT(com_matT, q=c(0,1,2,3,4,5), datatype = "abundance", nboot = 1)

#####################################

#### CHAO2 in INEXT - INCIDENCE

## 9th NOVEMBER
## READ IN AND CLEAN UP DATA - don't use drop.na for abundance- NA values are y for presence

data <- read.csv("temp_spp_v1_presence_2022-11-08.csv")

data <- (read.csv("temp_spp_v1_presence_2022-11-08.csv",header=T,sep=",",fileEncoding="latin1")[-2,])

data <- data[nzchar(data$Site),]
com_mat <-picante::sample2matrix(data)
com_mat[1:5, 1:5]

## com_mat$Abundance <- as.numeric(com_mat$Abundance)

empty_sites <- rownames(com_mat[rowSums(com_mat) == 0, ])
com_mat <- com_mat[!rownames(com_mat) %in% empty_sites, ]
dim(com_mat) ## [1] 1637 4440
dim(data.matrix) ## [1] 1736 4770

## INCIDENCE

com_mat_inc <- com_mat
com_mat_inc[which(com_mat>1)] <- 1
inc_freq <- as.incfreq(t(com_mat_inc))
Hills_q_0_inc <- iNEXT(inc_freq, q=0, datatype = "incidence_freq", nboot = 2)
test <- ChaoRichness((inc_freq))
write.csv(as.data.frame(Hills_q_0_inc, "OUTPUT.csv"))
#Error in as.data.frame.default(Hills_q_0_inc, "OUTPUT.csv") : cannot coerce class ‘"iNEXT"’ to a data.frame

plot <- ggiNEXT(Hills_q_0_inc, type = 1)
plot <- ggiNEXT(Hills_q_0_inc, type = 1) + theme_clean()

##################################

## OTHER VERSION INEXT INCIDENCE

m_dat <- read_csv("temp_spp_v1_presence_2022-11-08.csv")%>%
   filter(!is.na(Site))
dim(m_dat)

com_dat <- dcast(m_dat, Site~Species, fun.aggregate = sum, value.var = "Abundance", na.rm=T)
com_mat <- as.matrix(com_dat[,-1])
rownames(com_mat) <- com_dat[,1]

com_mat_inc <- com_mat
com_mat_inc[which(com_mat_inc>1)] <- 1
inc_freq <- as.incfreq(t(com_mat_inc))

Hills_q_0_inc <- iNEXT(inc_freq, q=0, datatype = "incidence_freq", nboot = 2)
write.csv(as.data.frame(Hills_q_0_inc, "OUTPUT.csv"))

test <- ChaoRichness((inc_freq))
write.csv(test, "OUTPUT.csv")

## PLOT

plot_q <- ggiNEXT(Hills_q_0_inc, type=1)+
   theme_clean() +
   scale_fill_manual(values = viridis(n=4)[-4])+
   scale_colour_manual(values = viridis(n=4)[-4])+
   facet_wrap(~order, scales = "free", ncol = 1)

##############################################
## rerun iNEXT- with all q nos and 100 bootstraps
##############################################

## REDO WITH ALL Q NOS- change to 100 bootstraps for overnight run
Hills_q <- iNEXT(com_matT, q=c(0,1,2,3,4,5), datatype = "abundance", nboot = 1)

saveRDS(Hills_q, "Hills_q-5_orders.rds")
Hills_q <- readRDS("Hills_q-5_orders.rds")
Hills_df <- bind_rows(Hills_q$iNextEst, .id="site")

## SAVE FILE
## saveRDS(Hills_q, "Hills_q-5_orders.rds")
write.table(community_matrix, file = 'data-processed/CCZ_community_matrix.txt')
write.csv(specaccum_df, file = 'CCZ_specaccum.csv')
save(i.out, file = 'data-processed/iNEXT_abundance.RData')
Hills_q

# -----------------------------------------------------------------------------------------------------------------

## S FIG HIGHER TAXON RICHNESS INTERPOLATION

# log of total number of taxa by taxonomic level- ie no of phlya- to genera - log this and use slope to get a prdicted species no


data <- read.csv("CCZ_CHECKLIST_2022-11-04.CSV")

## get total per taxon level- ie total phyla- order-calss-fam-genera recorded
data2 <- tapply(data[, "scientificName"], data[, "taxonRank"],
                function(x) { length(unique(x))})

## output summary as csv
write.csv(data2, "OUTPUT.csv")

## log of total per taxa level
test <- log(data$Total)
## append this as a column in file (could do within script- here I did in excel)
write.csv(test, "OUTPUT.csv")

## read in file with column of log value and taxon order as number- 1 - phlya, 2- class, 3 - order, 4 fam, 5 genus
data <- read.csv("temp_summary_checklist_taxonrank.csv")

## model
model1 <- lm(log ~ order, data = data)

autoplot(model1, smooth.colour = NA)
anova(model1)
summary(model1)

## PLOT
plot(data$olog, data$order, pch = 16)
## add line to basic plot -  line best fit with linear model
abline(lm(log ~ order, data = data)) ## 

## how to get log value for order = 6 (i.e. next taxon level- sp will be 6)
fitted <- predict(lm(data$log ~ data$order))

data$order[7];data$log[7]
order[7];log[7]; data = data

# -----------------------------------------------------------------------------------------------------------------


## all legit taxa records- not just spp level

## supp fig by depth - violin

data1 <- read.csv("temp_checklist_2022-11-06.csv")
data2 <- read.csv("CCZ_ALL_TAXA_2022-11-06.csv")
merga <- merge(data1, data2, by="scientificName", all.y = T)
dim(merga) ## 511   7
dim(data1) ## 268   6
dim(data2) ## 511   2
write.csv(merga, "OUTPUT.csv")

data <- read.csv("CCZ_ALL_TAXA_v1_2022-11-06.csv")
## >4000 replaced with 4000 for now

## remove nas in depth
data2 <- data %>%
   filter(Depth != "NA")
dim(data2) # [1] 62908    17
dim(data) [1] 64541    17

## too many phyla- weed out- where <100 records

data3 <- tapply(data2[, "all_taxa_rec"], data2[, "Phylum"],
                function(x) { length(unique(x))})
write.csv(data3, "OUTPUT.csv")

data_ed <- data %>%
   count(Phylum)
## same result

## remove following
#Entoprocta
#Cephalorhyncha
#Ctenophora
#Coelenterata
#Gnathostomulida
#Hemichordata
#Platyhelminthes
#Priapulida
#Nemertea
#Rotifera
## PLUs NO PHYLA ASSIGNED- IE 'METAZOA'

## read in edited file
data <- read.csv("temp_4_depth_CCZ_ALL_TAXA_v2_2022-11-06.csv")

ggplot(data = data, mapping = aes(x = Phylum, y = Depth, 
                                  fill = Phylum)) + geom_boxplot() + coord_flip()

## plot HORIZONTAL
ggplot(data = data, mapping = aes(x = Phylum, y = Depth, 
                                  fill = Phylum)) + geom_violin() + coord_flip() + theme_bw() +
   theme(text=element_text(size=12), #change font size of all text
         axis.text=element_text(size=14), #change font size of axis text
         axis.title=element_text(size=14), #change font size of axis titles
         plot.title=element_text(size=14), #change font size of plot title
         legend.text=element_text(size=14), #change font size of legend text
         legend.title=element_text(size=14)) + scale_fill_viridis_d() + 
   theme(legend.position = "none")


## vertical VERSION
ggplot(data = data, mapping = aes(x = Phylum, y = Depth, 
                                  fill = Phylum)) + geom_violin() + theme_bw() +
   scale_fill_viridis_d() +
   theme(text=element_text(size=14), #change font size of all text
         axis.text=element_text(size=14), #change font size of axis text
         axis.title=element_text(size=14), #change font size of axis titles
         plot.title=element_text(size=14), #change font size of plot title
         legend.text=element_text(size=14), #change font size of legend text
         legend.title=element_text(size=14)) + theme(legend.position = "none") +
   theme(axis.title.x=element_blank(),
         axis.text.x = element_text(angle = 55, hjust = 1))


#################################

## 8TH NOV
## CLEANED UP FILES

## INCL-ALL-TAXA	INCL-ALL SPP

## MERGE IN TAXONOMY FOR FAMILY CURVE

## REDO ALL ESTIMATES AND CURVES - FIG 3 + TABLE 1

## SUBTABLE OF NAMES BY HIGHER TAXONOMY- IE FIG 2 BUT ALL LEVELS

data <- read.csv("CCZ_CHECKLIST_2022-11-06.csv")
data2 <- data %>% 
   group_by(Phylum, Class, Order, Family, Genus) %>%
   count(short_spp_cl)


## merge taxonomy data to all taxa data
data1 <- read.csv("temp_checklist_2022-11-06.csv")
data2 <- read.csv("CCZ_ALL_TAXA_DATA_FIN_2022-11-08.csv")
merga <- merge(data1, data2, by="scientificName", all.y = T)
dim(merga) ## 511   7
dim(data1) ## 268   6
dim(data2) ## 511   2
write.csv(merga, "OUTPUT.csv")
