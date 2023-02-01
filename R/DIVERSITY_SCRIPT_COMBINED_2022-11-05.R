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
library(tidyverse)

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
 # theme(legend.position = "none") +
  scale_fill_manual('Position', values=c('steelblue', 'coral2')) +
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
data <- read.csv("temp_fam_v1_abundance_2022-11-08.csv") ## 4410 SPP 


## by size class
## by region

## 8TH NOV- LIT DATA ONLY- NO DUPLIC ACROSS SAMPLING UNITS ACROSS DATASETS- all zero rows incl
#data <- (read.csv("TEMP_ALL_SPP_LIT_ONLY_NO_ZERO_2022-11-08.csv",header=T,sep=",",fileEncoding="latin1")[-2,])

head(data)
dim(data) ## 39116 ## no nas- 36 842
sum(data$Abundance)

## MAKE SPP MATRIX

data.matrix <- sample2matrix(data)
data.matrix[1:5, 1:5]

length(unique(data$Species[data$Abundance > 0])) ## 4372


## CLEANUP

## length(names(apply(data.matrix, is.numeric))) == length(colnames(data.matrix))
## data.matrix$Abundance <- as.numeric(data.matrix$Abundance)

## SPECIES ACCUMULATION CURVE, SPEC POOOL- VEGAN

data3 <- specpool(data.matrix)
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
write.csv(data.curve.df, "OUTPUT.csv")

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
sum(data$Adundance, na.rm = T)
data.matrix <- sample2matrix(data)
data.matrix[1:5, 1:5]

species_no <- specnumber(data.matrix)
raremax <- min(rowSums(data.matrix))
Srare <- rarefy(data.matrix, raremax)
plot(species_no, Srare, xlab = "Observed No. of Species", ylab = "Species")
abline(0, 1)

rarecurve(data.matrix, step = 20, sample = raremax, cex = 1, col = viridis(n=4)[-4], 
          lwd = 2, lty = c("solid", rep("dashed", 2)), label = F,xlab = "No. of individuals", 
          ylab = "Species Richness", cex.axis = 1)
# legend("bottomright", legend = rownames(com_mat), col = viridis(n=4)[-4], pch=19, cex = 0.5, y.intersp = 0.8, pt.cex = 1)
rarecurve(data.matrix, step = 20)

## estimateR for CCZ only data

test <- estimateR(data.matrix) 


############################

## CHAO1 AND CHAO2 IN INEXT- 

#### CHAO1 
data <- read.csv("temp_fam_v1_abundance_2022-11-08.csv") %>% mutate(Site = "CCZ")
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


Hills_q <- iNEXT(com_matT, q=c(0,1,2,3), datatype = "abundance", nboot = 100)

#####################################

#### CHAO2 in INEXT - INCIDENCE

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

m_dat <- read_csv("temp_fam_v1_presence_2022-11-08.csv")%>%
  filter(!is.na(Site))
dim(m_dat)

com_dat <- dcast(m_dat, Site~Species, fun.aggregate = sum, value.var = "Abundance", na.rm=T)
com_mat <- as.matrix(com_dat[,-1])
rownames(com_mat) <- com_dat[,1]

com_mat_inc <- com_mat
com_mat_inc[which(com_mat_inc>1)] <- 1
inc_freq <- as.incfreq(t(com_mat_inc))

Hills_q_0_inc <- iNEXT(inc_freq, q=0, datatype = "incidence_freq", nboot = 2 )
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
data <- read.csv("temp_log_2022-11-06.csv")

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

background on deep-sea biodiversity

knowledge of deep-sea biodiversity

#################################

## 8TH NOV
## CLEANED UP FILES

## INCL-ALL-TAXA	INCL-ALL SPP

## MERGE IN TAXONOMY FOR FAMILY CURVE

## REDO ALL ESTIMATES AND CURVES - FIGure 3 + TABLE 1

## SUBTABLE OF NAMES BY HIGHER TAXONOMY- IE FIG 2 BUT ALL LEVELS

data <- read.csv("CCZ_CHECKLIST_2022-11-06.csv")
names(data)
data2 <- data %>%
  filter(short_spp_cl == "yes")
dim(data2)

data3 <- data2 %>% 
  group_by(Phylum, Class, Order, Family, Genus) %>%
  count(taxonRank)
write.csv(data3, "OUTPUT.csv")

data3 <- data2 %>% 
  group_by(Phylum, Class, Order) %>%
  count(taxonRank)

data3 <- data2 %>% 
  group_by(Phylum, Class) %>%
  count(taxonRank)

## merge taxonomy data to all taxa data
data1 <- read.csv("temp_checklist_2022-11-06.csv")
data2 <- read.csv("CCZ_ALL_TAXA_DATA_FIN_2022-11-08.csv")
merga <- merge(data1, data2, by="scientificName", all.y = T)
dim(merga) ## 511   7
dim(data1) ## 268   6
dim(data2) ## 511   2
write.csv(merga, "OUTPUT.csv")

#####################################

## UPSET- other script

#upset <- read.csv("temp_DD_ed_4_diversity_all_spp_region_2022-04-14.csv")
#data.matrix <- as.matrix(upset)
#test <- make_comb_mat(upset)

#upset <- read.csv("DD_pres-abs_table_v1_2021-10-26.csv")

#data.matrix <- as.matrix(upset)
#head(data.matrix)
#write.csv(data.matrix, "OUTPUT.csv")

## transpose- species in rows, contractors in columns
## t - transpose

#test <- t(data.matrix)
#test <- t(data.presabs)
#head(test)


### PLOTTING

basin_COLS <- viridis(n=8)[-8]
names(basin_COLS) <- colnames(test)
morph_comb <- make_comb_mat(test, mode = "intersect")

# pdf("Shared_wCCZ_species-Upset_plot.pdf", width = 10, height = 3.5)
morph_plot <- UpSet(morph_comb, set_order = colnames(morph_comb), 
                    comb_order = rev(order(comb_size(morph_comb))), 
                    pt_size = unit(2, "mm"), lwd = 1, top_annotation = HeatmapAnnotation(
                      
                      "Shared species (all)" = anno_barplot(comb_size(morph_comb),
                                                            ylim = c(0, max(comb_size(morph_comb))*1.1),
                                                            border = FALSE,
                                                            gp = gpar(fill = "black"),
                                                            height = unit(10, "cm")  ## adjust bar to dot matrix ratio: 10 REG   
                      ),
                      annotation_name_side = "left",
                      annotation_name_rot = 90), ## angle of legend
                    left_annotation = rowAnnotation(
                      "No. species per region (all)" = anno_barplot(-set_size(morph_comb),
                                                                    baseline = 0,
                                                                    axis_param = list(
                                                                      at = c(0, -250, -500, -1000, -1500),
                                                                      labels = c(0, 250, 500, 1000, 1500),
                                                                      labels_rot = 0),
                                                                    border = FALSE,
                                                                    gp = gpar(fill = basin_COLS),
                                                                    width = unit(7, "cm"),  ## width of lhs of y axis: 7: CA; 5: SA; 10 REG
                                                                    height = unit(0.5, "cm")),
                      set_name = anno_text(set_name(morph_comb),
                                           location = 0.5, ## labels for side bars
                                           just = "center",
                                           width = max_text_width(set_name(morph_comb)) + 
                                             unit(1, "mm"))), right_annotation = NULL, show_row_names=FALSE)

morph_plot = draw(morph_plot)


###########################
## reset ggINEXT

plot_q <- ggiNEXT(Hills_q_0_inc, type=1)+
  theme_clean() +
  scale_fill_manual(values = viridis(n=4)[-4])+
  scale_colour_manual(values = viridis(n=4)[-4])+
  facet_wrap(~order, scales = "free", ncol = 1)


'save(Hills_q_CCZ, file = 'data-processed/Hills_q_CCZ.RData’)'

###################################################

## 19th Jan 2023
## do plot of sampling completeness- add to supplementary
##


