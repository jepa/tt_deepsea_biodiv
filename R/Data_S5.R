library(tidyverse)
library(patchwork)
library(vegan)
library(ggthemes)
library(ggstance)
library(picante)
library(iNEXT)
library(reshape)
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
library(reshape2)
library(fishualize)
library(patchwork)
theme_set(theme_bw(base_size = 14)) ## was 12

# -----------------------------------------------------------------------------------------------------------------
# Data import
# -----------------------------------------------------------------------------------------------------------------
# All raw data
full_data <- read.csv("CCZ_ALL_TAXA_DATA_FIN_2023-02-19.csv") %>%
  rename_all(tolower) %>% 
  mutate(abundance = as.numeric(abundance))
dim(full_data) #91996

# Raw checklist data
checklist_raw <- read.csv("CCZ_CHECKLIST_2023-02-17.csv") %>% 
  dplyr::select(Family, ScientificName_accepted) %>% 
  unique() %>% 
  rename_all(tolower) 

# Add family column to full_data
#full_data <- merge(full_data, checklist_raw, all.x = T)

# Species abundance by site
species_abundance <- full_data %>% 
  filter(!is.na(species), species != "") %>% 
  group_by(site, species) %>% 
  summarise(abundance = sum(abundance))%>% 
  filter(!is.na(abundance))

# Family abundance by site
family_abundance <- full_data %>% 
  filter(!is.na(family), family != "") %>% 
  group_by(site, family) %>% 
  summarise(abundance = sum(abundance)) %>% 
  filter(!is.na(abundance))
## COMMENT OUT LAST LINE

# genera abundance by site
genera_abundance <- full_data %>% 
  filter(!is.na(genus), genus != "") %>% 
  group_by(site, genus) %>% 
  summarise(abundance = sum(abundance)) %>% 
  filter(!is.na(abundance))

# Species presence-absence by site
species_presAbs <- full_data %>% 
  filter(!is.na(species), species != "") %>% 
  group_by(site, species) %>% 
  summarise(presence = ifelse(abundance > 0, 1, 0)) %>% 
  filter(!is.na(presence))

#Warning message:
#   Returning more (or less) than 1 row per `summarise()` group was deprecated in dplyr 1.1.0.
# Please use `reframe()` instead.
# When switching from `summarise()` to `reframe()`, remember that `reframe()` always returns an
#ungrouped data frame and adjust accordingly.
#Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated. 

# family presence-absence by site ## RENAME AS FAMILY
family_presAbs <- full_data %>% 
  filter(!is.na(family), family != "") %>% 
  group_by(site, family) %>% 
  summarise(presence = ifelse(abundance > 0, 1, 0)) %>% 
  filter(!is.na(presence))

# Species description data
species_descriptions <- read.csv("TEMP_papers_table_1980_on_2023-01-19.csv") %>% 
  pivot_longer(cols = 6:8, names_to = "var", values_to = "number")

# Phyla data
phyla_overview <- read.csv("TEMP_SUMMARY_FIG2_ALL_PHYLA_2023-01-20.csv") %>% 
  group_by(Phylum) %>% 
  mutate(total_phylum = sum(Total)) %>% 
  group_by(Data) %>% 
  mutate(perc = paste0(round((Total/total_phylum)*100), "%"))
phyla_overview$ypos <- ifelse(phyla_overview$Data == "named species", phyla_overview$total_phylum, phyla_overview$total_phylum + 70)

# -----------------------------------------------------------------------------------------------------------------
# Figure 3: Taxonomic work
# -----------------------------------------------------------------------------------------------------------------
p1 <- ggplot(data = species_descriptions, aes(x = Year)) + theme_bw() +
  geom_line(aes(y = All_CCZ_Descriptions, colour = "All_CCZ_Descriptions"), size = 1.5) +
  geom_line(aes(y = cumul_desc, colour = "cumul_desc"), size = 1.5) +
  geom_line(aes(y = cumul_spp, colour = "cumul_spp"), size = 1.5) +
  geom_line(aes(y = cumul_pubs, colour = "cumul_pubs"), size = 1.5) +
  labs(x = "Year", y = "Total") +
  #  ggtitle("Cumulative totals of publications and descriptions (year 2000- 2021)") +
  theme(text=element_text(size=14), #change font size of all text
        axis.text=element_text(size=14), #change font size of axis text
        axis.title=element_text(size=14), #change font size of axis titles
        plot.title=element_text(size=14), #change font size of plot title
        legend.text=element_text(size=14), #change font size of legend text
        legend.title=element_text(size=14)) + 
  theme(legend.justification = c(-0.35, 1), 
        legend.position = c(0, 1),
        legend.box.margin=margin(c(20, 20, 20, 20)),
        legend.title = element_blank(),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
  scale_colour_manual("", 
                      breaks = c("cumul_desc","cumul_spp","cumul_pubs","All_CCZ_Descriptions"),
                      values = c("black", "coral", "steelblue","pink"),
                      labels = c("Cumulative new descriptions", "Cumulative new species", "Cumulative publications","Yearly descriptions"))

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
                    labels = c('Unnamed species', 'Named species')) +
  scale_color_manual(values = c('steelblue', 'coral2'), guide = "none") +
  ylim(0, 2000); phyla_figure

figure_3 <- descriptions_figure / phyla_figure + plot_annotation(tag_levels = 'A')

ggsave(figure_3,
       filename = 'output-figures/figure_2.tiff', 
       width = 8.5, height = 10, units = 'in', dpi = 150)


# -----------------------------------------------------------------------------------------------------------------
# Figure 4: Relative abundance of Checklists
# -----------------------------------------------------------------------------------------------------------------

library(fishualize)
data <- read.csv("Checklist_summaries_2023-02-08.csv")

rel_abundance_bplot <-  ggplot(data, aes(x = Checklist, y = Proportion, fill = Phylum))+
  geom_bar(stat = "identity", width = 0.8, aes(group = as.factor(Checklist)))+
  scale_fill_fish_d(option = "Coris_gaimard") +
  theme(text=element_text(size=12), #change font size of all text
        axis.text=element_text(size=12), #change font size of axis text
        axis.title=element_text(size=12), #change font size of axis titles
        plot.title=element_text(size=12), #change font size of plot title
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) +
  #scale_fill_manual(values = cols, name = "Species")+
  # geom_text(aes(y = ypos, label = paste0(round(prop, digits = 0), "%")), color = "black", size=2)+
  theme(legend.key.size = unit(0.5, "cm"), legend.key.height =  unit(0.2, "cm"))+
  xlab("Checklist")+
  ylab("Relative abundance (%)") 

# -----------------------------------------------------------------------------------------------------------------
# Figure 5: Species accumulation curves for CCZ
# -----------------------------------------------------------------------------------------------------------------
# Community matrix for species abundance by single site: CCZ
data_CCZ_only <- species_abundance %>% 
  mutate(site = "CCZ")

community_matrix_CCZ <- data_CCZ_only %>%
  dplyr::group_by(site, species) %>%
  dplyr::summarise(cover = sum(abundance)) %>%
  reshape::cast(.,  site ~ species, value = "cover")
community_matrix_CCZ[is.na(community_matrix_CCZ)] <-  0

com_matT <- as.matrix(t(community_matrix_CCZ))
Hills_q_CCZ <- iNEXT::iNEXT(com_matT, q = 0, datatype = "abundance", nboot = 1)
Hills_q_CCZ_df <- as.data.frame(Hills_q_CCZ$iNextEst$size_based)
write.csv(Hills_q_CCZ_df, "OUTPUT.csv")   
length(unique(Hills_q_CCZ$species))

ggiNEXT(Hills_q_CCZ, type=2)

A_Chao1 <- ggplot(Hills_q_CCZ_df %>% filter(Method != 'Observed'), aes(x = m)) +
  geom_ribbon(aes(x = m, ymin = qD.LCL, ymax = qD.UCL), fill = "coral2", alpha = 0.2) +
  geom_line(aes(y = qD, lty = Method),
            col = "coral2", cex = 2) +
  geom_point(data = Hills_q_CCZ_df %>% filter(Method == 'Observed'), aes(x = m, y = qD), col = "coral2", cex = 2) +
  theme(legend.justification = c(0, 1), 
        legend.position = c(.5, .4),
        legend.box.margin=margin(c(10, 10, 10, 10)),
        legend.margin = margin(2, 2, 2, 2),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
  theme(text=element_text(size=12), #change font size of all text
        axis.text=element_text(size=20), #change font size of axis text
        axis.title=element_text(size=20), #change font size of axis titles
        plot.title=element_text(size=20), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20)) + 
  xlab("Number of Individuals") +
  ylab("Species Diversity") +
  scale_y_continuous(breaks = seq(0, 5000, 1000)) +
  scale_linetype_manual(values = c("dashed", "solid")); A_Chao1

############################
# Community matrix for family abundance by single site: CCZ
family_CCZ_only <- family_abundance %>% 
  mutate(site = "CCZ")

family_matrix_CCZ <- family_CCZ_only %>%
  dplyr::group_by(site, family) %>%
  dplyr::summarise(cover = sum(abundance)) %>%
  reshape::cast(.,  site ~ family, value = "cover")
family_matrix_CCZ[is.na(family_matrix_CCZ)] <-  0

family_matT <- as.matrix(t(family_matrix_CCZ))
Hills_q_CCZ_family <- iNEXT::iNEXT(family_matT, q = 0, datatype = "abundance", nboot = 1)
Hills_q_CCZ_family_df <- as.data.frame(Hills_q_CCZ_family$iNextEst$size_based)
write.csv(Hills_q_CCZ_family_df, "OUTPUT.csv")  

ggiNEXT(Hills_q_CCZ_family, type=2)

C_Chao1_family <- ggplot(Hills_q_CCZ_family_df %>% filter(Method != 'Observed'), aes(x = m)) +
  geom_ribbon(aes(x = m, ymin = qD.LCL, ymax = qD.UCL), fill = "coral2", alpha = 0.2) +
  geom_line(aes(y = qD, lty = Method),
            col = "coral2", cex = 2) +
  geom_point(data = Hills_q_CCZ_family_df %>% filter(Method == 'Observed'), aes(x = m, y = qD), col = "coral2", cex = 2) +
  theme(legend.justification = c(0, 1), 
        legend.position = c(.5, .4),
        legend.box.margin=margin(c(20, 20, 20, 20)),
        legend.margin = margin(1, 1, 1, 1),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
  theme(text=element_text(size=20), #change font size of all text
        axis.text=element_text(size=20), #change font size of axis text
        axis.title=element_text(size=20), #change font size of axis titles
        plot.title=element_text(size=20), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20)) + 
  xlab("Number of Individuals") +
  ylab("Family Diversity") +
  # scale_y_continuous(breaks = seq(0, 5000, 1000)) +
  scale_linetype_manual(values = c("dashed", "solid")); C_Chao1_family

####################################

# Community matrix for species presence-absence by all sites (the sampling units)
community_matrix_pres <- species_presAbs %>%
  dplyr::group_by(site, species) %>%
  dplyr::summarise(presence = ifelse(sum(presence) > 0, 1, 0)) %>%
  reshape::cast(.,  site ~ species, value = "presence") ## this causing issue
community_matrix_pres[is.na(community_matrix_pres)] <-  0
inc_freq <- as.incfreq(as.matrix(t(community_matrix_pres)))
Hills_q_0_inc_CCZ <- iNEXT(inc_freq, q = 0, datatype = "incidence_freq", nboot = 1)
Hills_q0_inc_df <- as.data.frame(Hills_q_0_inc_CCZ$iNextEst$size_based)

community_matrix_pres <- species_presAbs %>%
  dplyr::group_by(site, species) %>%
  dplyr::summarise(presence = ifelse(sum(presence) > 0, 1, 0)) 

community_matrix_pres[is.na(community_matrix_pres)] <-  0
inc_freq <- as.incfreq(as.matrix(t(community_matrix_pres)))
Hills_q_0_inc_CCZ <- iNEXT(inc_freq, q = 0, datatype = "incidence_freq", nboot = 1)
Hills_q0_inc_df <- as.data.frame(Hills_q_0_inc_CCZ$iNextEst$size_based)

B_Chao2 <- ggplot(Hills_q0_inc_df %>% filter(Method != 'Observed')) +
  geom_ribbon(aes(x = t, ymin = qD.LCL, ymax = qD.UCL), fill = "coral2", alpha = 0.2, show.legend = FALSE) +      
  geom_line(aes(x = t, y = qD, lty = Method), col = "coral2", cex = 1) +
  geom_point(data = Hills_q0_inc_df %>% filter(Method == 'Observed'), aes(x = t, y = qD), col = "coral2", cex = 2) +
  theme(legend.justification = c(0, 1), 
        legend.position = c(.5, .4),
        legend.box.margin=margin(c(20, 20, 20, 20)),
        legend.margin = margin(1, 1, 1, 1),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
  theme(text=element_text(size=20), #change font size of all text
        axis.text=element_text(size=20), #change font size of axis text
        axis.title=element_text(size=20), #change font size of axis titles
        plot.title=element_text(size=20), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20)) + 
  xlab("Number of sampling units") +
  ylab("Species Diversity") +
  scale_y_continuous(breaks = seq(0, 6000, 1000)) +
  scale_linetype_manual(values = c("dashed", "solid")); B_Chao2

#########################################
# Community matrix for family presence-absence by all sites (the sampling units)
family_matrix_pres <- family_presAbs %>%
  dplyr::group_by(site, family) %>%
  dplyr::summarise(presence = ifelse(sum(presence) > 0, 1, 0)) %>%
  reshape::cast(.,  site ~ family, value = "presence")
family_matrix_pres[is.na(family_matrix_pres)] <-  0
inc_freq_family <- as.incfreq(as.matrix(t(family_matrix_pres)))
Hills_q_0_inc_CCZ_family <- iNEXT(inc_freq_family, q = 0, datatype = "incidence_freq", nboot = 1)
Hills_q_0_inc_CCZ_family_df <- as.data.frame(Hills_q_0_inc_CCZ_family$iNextEst$size_based)
write.csv(Hills_q_0_inc_CCZ_family_df, "OUTPUT.csv")   

D_Chao2_family <- ggplot(Hills_q_0_inc_CCZ_family_df %>% filter(Method != 'Observed')) +
  geom_ribbon(aes(x = t, ymin = qD.LCL, ymax = qD.UCL), fill = "coral2", alpha = 0.2, show.legend = FALSE) +      
  geom_line(aes(x = t, y = qD, lty = Method), col = "coral2", cex = 1) +
  geom_point(data = Hills_q_0_inc_CCZ_family_df %>% filter(Method == 'Observed'), aes(x = t, y = qD), col = "coral2", cex = 2) +
  theme(legend.justification = c(0, 1), 
        legend.position = c(.5, .4),
        legend.box.margin=margin(c(20, 20, 20, 20)),
        legend.margin = margin(1, 1, 1, 1),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
  theme(text=element_text(size=20), #change font size of all text
        axis.text=element_text(size=20), #change font size of axis text
        axis.title=element_text(size=20), #change font size of axis titles
        plot.title=element_text(size=20), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20)) +
  xlab("Number of sampling units") +
  ylab("Family Diversity") +
  # scale_y_continuous(breaks = seq(0, 6000, 1000)) +
  scale_linetype_manual(values = c("dashed", "solid")); D_Chao2_family

figure_3 <- (A_Chao1_species | B_Chao2_species) / (C_Chao1_family | D_Chao2_family) + plot_annotation(tag_levels = 'A')

ggsave(figure_3,
       filename = 'output-figures/figure_3.tiff', 
       width = 8.5, height = 7, units = 'in', dpi = 150)


# -----------------------------------------------------------------------------------------------------------------
# Figure 6: UpSet plots
# -----------------------------------------------------------------------------------------------------------------

upset <- (read.csv("temp_4_upset_APEI-CA_2023-02-18.csv",header=T,sep=",",fileEncoding="latin1")[-2,])
upset <- (read.csv("TEMP_4_UPSET_REGION_2023-02-18.csv",header=T,sep=",",fileEncoding="latin1")[-2,])

## w presence absense from fuzzy sim
data.presabs <- splist2presabs(upset, sites.col = "Site",
                               sp.col = "Species", keep.n = FALSE)

test <- t(data.presabs)
head(test)
test <- t(as.matrix(data.presabs[,-1]))
colnames(test) <- data.presabs$Site

### PLOTTING

basin_COLS <- viridis(n=3)[-1] ## ADJUST BY TOT SITES - 3 - 1 for CA- APEI, 4 - 1 for region
names(basin_COLS) <- colnames(test)
morph_comb <- make_comb_mat(test, mode = "intersect")

# pdf("Shared_wCCZ_species-Upset_plot.pdf", width = 10, height = 3.5)
morph_plot <- UpSet(morph_comb, set_order = colnames(morph_comb), 
                    comb_order = rev(order(comb_size(morph_comb))), 
                    pt_size = unit(2, "mm"), lwd = 1, top_annotation = HeatmapAnnotation(
                      
                      "Shared species" = anno_barplot(comb_size(morph_comb),
                                                      ylim = c(0, max(comb_size(morph_comb))*1.1),
                                                      border = FALSE,
                                                      gp = gpar(fill = "black"),
                                                      height = unit(11, "cm")  ## adjust bar to dot matrix ratio: 10-12 REG   
                      ),
                      annotation_name_side = "left",
                      annotation_name_rot = 90), ## angle of legend
                    left_annotation = rowAnnotation(
                      "Species per region" = anno_barplot(-set_size(morph_comb),
                                                          baseline = 0,
                                                          axis_param = list(
                                                            at = c(0, -1000, -2000),## del -500, -1500
                                                            labels = c(0, 1000, 2000),## del 500,1500
                                                            labels_rot = 0),
                                                          border = FALSE,
                                                          gp = gpar(fill = basin_COLS),
                                                          width = unit(10, "cm"),  ## width of lhs of y axis: 7: CA; 5: SA; 10 REG
                                                          height = unit(0.9, "cm")),
                      set_name = anno_text(set_name(morph_comb),
                                           location = 0.5, ## labels for side bars
                                           just = "center",
                                           width = max_text_width(set_name(morph_comb)) + 
                                             unit(1, "mm"))), right_annotation = NULL, show_row_names=FALSE)

morph_plot = draw(morph_plot)

# -----------------------------------------------------------------------------------------------------------------
# Supplementary figures
# -----------------------------------------------------------------------------------------------------------------

## Figure S1: species accumulation curves/rarefaction curves in vegan

data.matrix <- sample2matrix(data)
data.matrix[1:5, 1:5]
length(unique(data$Species[data$Abundance > 0])) 
length(unique(data$Species)) 
df <- as.data.frame(data.matrix)
data3 <- specpool(data.matrix)

data.curve <- specaccum(data.matrix, method = "random", permutations = 1000)

data.curve.df <- data.frame(Sites = data.curve$sites,
                            Richness = data.curve$richness,
                            sd = data.curve$sd)
head(data.curve.df)

ggplot(data.curve.df, aes(x = Sites, y = Richness)) +
  # Add line
  geom_line() +
  # Add confidence intervals
  geom_ribbon(aes(ymin = Richness - 1.96*sd, ymax = Richness + 1.96*sd), 
              alpha = 0.2, colour = "lightgrey") +
  theme(text=element_text(size=20), #change font size of all text
        axis.text=element_text(size=20), #change font size of axis text
        axis.title=element_text(size=20), #change font size of axis titles
        plot.title=element_text(size=20), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20))

## rarefaction curve

sum(data$Adundance, na.rm = T)
data.matrix <- sample2matrix(data)
data.matrix[1:5, 1:5]

species_no <- specnumber(data.matrix) #4324
raremax <- min(rowSums(data.matrix))
Srare <- rarefy(data.matrix, raremax)
plot(species_no, Srare, xlab = "Observed No. of Species", ylab = "Species")
abline(0, 1)

rarecurve(data.matrix, step = 20, sample = raremax, cex = 1, col = viridis(n=4)[-4], 
          lwd = 2, lty = c("solid", rep("dashed", 2)), label = F,xlab = "No. of individuals", 
          ylab = "Species Richness", cex.axis = 1)
# legend("bottomright", legend = rownames(com_mat), col = viridis(n=4)[-4], pch=19, cex = 0.5, y.intersp = 0.8, pt.cex = 1)
rarecurve(data.matrix, step = 20)

estimateR(data.matrix) 

# -----------------------------------------------------------------------------------------------------------------
# Figure S2: Diversity of Genera
# -----------------------------------------------------------------------------------------------------------------

data_CCZ_only <- genera_abundance %>% 
  mutate(site = "CCZ")

community_matrix_CCZ <- data_CCZ_only %>%
  dplyr::group_by(site, genus) %>%
  dplyr::summarise(cover = sum(abundance)) %>%
  reshape::cast(.,  site ~ genus, value = "cover")
community_matrix_CCZ[is.na(community_matrix_CCZ)] <-  0

com_matT <- as.matrix(t(community_matrix_CCZ))
Hills_q_CCZ <- iNEXT::iNEXT(com_matT, q = 0, datatype = "abundance", nboot = 1)

ggiNEXT(Hills_q_CCZ, type=2)

Hills_q_CCZ_df <- as.data.frame(Hills_q_CCZ$iNextEst$size_based)
write.csv(Hills_q_CCZ_df, "OUTPUT.csv")   

AB_Chao1 <- ggplot(Hills_q_CCZ_df %>% filter(Method != 'Observed'), aes(x = m)) +
  geom_ribbon(aes(x = m, ymin = qD.LCL, ymax = qD.UCL), fill = "coral2", alpha = 0.2) +
  geom_line(aes(y = qD, lty = Method),
            col = "coral2", cex = 2) +
  geom_point(data = Hills_q_CCZ_df %>% filter(Method == 'Observed'), aes(x = m, y = qD), col = "coral2", cex = 2) +
  theme(legend.justification = c(0, 1), 
        legend.position = c(.5, .4),
        legend.box.margin=margin(c(10, 10, 10, 10)),
        legend.margin = margin(2, 2, 2, 2),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
  theme(text=element_text(size=20), #change font size of all text
        axis.text=element_text(size=20), #change font size of axis text
        axis.title=element_text(size=20), #change font size of axis titles
        plot.title=element_text(size=20), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20)) + 
  xlab("Number of Individuals") +
  ylab("Diversity of Genera") +
  # scale_y_continuous(breaks = seq(0, 5000, 1000)) +
  scale_linetype_manual(values = c("dashed", "solid")); AB_Chao1

gen_matrix_pres <- gen_presAbs %>%
  dplyr::group_by(site, genus) %>%
  dplyr::summarise(presence = ifelse(sum(presence) > 0, 1, 0)) %>%
  reshape::cast(.,  site ~ family, value = "presence")
family_matrix_pres[is.na(family_matrix_pres)] <-  0
inc_freq_gen <- as.incfreq(as.matrix(t(gen_matrix_pres)))
Hills_q_0_inc_CCZ_gen <- iNEXT(inc_freq_gen, q = 0, datatype = "incidence_freq", nboot = 1)
Hills_q_0_inc_CCZ_gen_df <- as.data.frame(Hills_q_0_inc_CCZ_gen$iNextEst$size_based)
write.csv(Hills_q_0_inc_CCZ_gen_df, "OUTPUT.csv")   

D_Chao2_gen <- ggplot(Hills_q_0_inc_CCZ_gen_df %>% filter(Method != 'Observed')) +
  geom_ribbon(aes(x = t, ymin = qD.LCL, ymax = qD.UCL), fill = "coral2", alpha = 0.2, show.legend = FALSE) +      
  geom_line(aes(x = t, y = qD, lty = Method), col = "coral2", cex = 1) +
  geom_point(data = Hills_q_0_inc_CCZ_gen_df %>% filter(Method == 'Observed'), aes(x = t, y = qD), col = "coral2", cex = 2) +
  theme(legend.justification = c(0, 1), 
        legend.position = c(.5, .4),
        legend.box.margin=margin(c(20, 20, 20, 20)),
        legend.margin = margin(1, 1, 1, 1),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
  theme(text=element_text(size=20), #change font size of all text
        axis.text=element_text(size=20), #change font size of axis text
        axis.title=element_text(size=20), #change font size of axis titles
        plot.title=element_text(size=20), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20)) +
  xlab("Number of sampling units") +
  ylab("Diversity of Genera") +
  # scale_y_continuous(breaks = seq(0, 6000, 1000)) +
  scale_linetype_manual(values = c("dashed", "solid")); D_Chao2_family


# -----------------------------------------------------------------------------------------------------------------
# Figure S5: Density of sampling by depth
# -----------------------------------------------------------------------------------------------------------------

full_data <- read.csv('CCZ_ALL_TAXA_DATA_FIN_2023-02-19.csv') %>% 
  filter(!is.na(LONGITUDE), !is.na(LATITUDE)) %>% 
  rename(site = SITE, 
         lat = LATITUDE, 
         long = LONGITUDE,
         depth = DEPTH) %>% 
  mutate(year = sub('.*(\\d{4}).*', '\\1', DATE),
         depth = as.numeric(depth),
         year = as.integer(year),
         phylum = Phylum,
         depth_group_width = as.numeric(cut_width(depth, width = 100, center = 0)),
         sample_unit = paste0(site, "_", year)) %>% 
  drop_na()

full_data <- full_data[nzchar(full_data$phylum), ]

effort_data <- full_data  %>% 
  dplyr::select(site, depth, year, phylum) %>% 
  unique() %>% 
  # How should we subdivide this?
  mutate(depth_group_number = as.numeric(cut_number(depth, 10)),
         depth_group_interval = as.numeric(cut_interval(depth, n = 10)),
         depth_group_width = as.numeric(cut_width(depth, width = 100, center = 0)),
         sample_unit = paste0(site, "_", year)) %>% 
  group_by(depth_group_number) %>% 
  mutate(depth_group_number_min = min(depth),
         depth_group_number_max = max(depth)) %>% 
  ungroup() %>% 
  group_by(depth_group_interval) %>% 
  mutate(depth_group_interval_min = min(depth),
         depth_group_interval_max = max(depth)) %>% 
  ungroup() %>% 
  group_by(depth_group_width) %>% 
  mutate(depth_group_width_min = min(depth),
         depth_group_width_max = max(depth)) %>% 
  ungroup()

depth_width_effort <- effort_data %>% 
  dplyr::select(depth_group_width, sample_unit) %>% 
  unique() %>% 
  group_by(depth_group_width) %>% 
  summarise(width_effort = n())

samples_by_depth <- ggplot(full_data) +
  geom_density(aes(y = depth), fill = "lightblue") +
  # geom_hline(aes(yintercept = depth_group_min, colour = depth_group)) +
  scale_y_reverse() +
  #    ggtitle("Number of samples by depth, all phyla combined") +
  scale_colour_gradient(low = "lightblue", high = "dodgerblue4") +
  theme_bw(); samples_by_depth

depth_group_number <- ggplot(effort_data) +
  geom_rect(data = effort_data %>% select(depth_group_number_min, depth_group_number_max, depth_group_number) %>% unique(),
            aes(ymin = depth_group_number_min, ymax = depth_group_number_max,
                xmin = 0, xmax = 0.003, fill = depth_group_number), col = "black", alpha = 0.6) +
  geom_density(aes(y = depth), fill = "coral1")+
  scale_y_reverse() +
  theme(legend.position = "none") +
  theme(text=element_text(size=18), #change font size of all text
        axis.text=element_text(size=18), #change font size of axis text
        axis.title=element_text(size=18), #change font size of axis titles
        plot.title=element_text(size=18), #change font size of plot title
        legend.text=element_text(size=18), #change font size of legend text
        legend.title=element_text(size=18)) 
#      ggtitle("Number of samples by depth, grouped by 10 sample quantiles") +
scale_fill_gradient(low = "lightblue", high = "dodgerblue4") +
  theme_bw(); depth_group_number

depth_group_interval <- ggplot(effort_data) +
  geom_rect(data = effort_data %>% select(depth_group_interval_min, depth_group_interval_max, depth_group_interval) %>% unique(),
            aes(ymin = depth_group_interval_min, ymax = depth_group_interval_max,
                xmin = 0, xmax = 0.003, fill = depth_group_interval), col = "black", alpha = 0.6) +
  geom_density(aes(y = depth), fill = "coral1")+
  scale_y_reverse() +
  ggtitle("Number of samples by depth, grouped by 10 depth quantiles") +
  scale_fill_gradient(low = "lightblue", high = "dodgerblue4") +
  theme_bw(); depth_group_interval

#depth_group_width <- ggplot(effort_data) +
#  geom_rect(data = effort_data %>% select(depth_group_width_min, depth_group_width_max, depth_group_width) %>% unique(),
#            aes(ymin = depth_group_width_min, ymax = depth_group_width_max,
#                xmin = 0, xmax = 0.003, fill = depth_group_width), col = "black", alpha = 0.6) +
#  geom_density(aes(y = depth), fill = "coral1")+
#  scale_y_reverse() +
#  ggtitle("Number of samples by depth, grouped by 100m intervals") +
#  scale_fill_gradient(low = "lightblue", high = "dodgerblue4") +
#  theme_bw(); depth_group_width


depth_groups <- samples_by_depth / depth_group_number / depth_group_interval / depth_group_width

samples_by_depth_phyla <- ggplot(full_data) +
  geom_density(aes(y = depth), fill = "lightblue") +
  facet_wrap(facets = vars(phylum)) +
  scale_y_reverse() +
  ggtitle("Number of samples by depth and phyla") +
  theme_bw(); samples_by_depth_phyla

ggsave(depth_groups,
       filename = 'output-figures/samples_by_depth.jpg', 
       width = 9, height = 14, units = 'in', dpi = 150)

ggsave(samples_by_depth_phyla,
       filename = 'output-figures/samples_by_depth_phyla.jpg', 
       width = 8, height = 6, units = 'in', dpi = 150)

full_with_effort <- merge(full_data, depth_width_effort, all.x = T) 
stand_abundance <- full_with_effort %>% 
  group_by(phylum, depth_group_width) %>% 
  summarise(stand_abundance = round(sum(ABUNDANCE)/width_effort)) %>% 
  unique()

ggplot(stand_abundance) +
  geom_point(aes(x = depth_group_width, y = stand_abundance, col = phylum)) 

# -----------------------------------------------------------------------------------------------------------------
# Figure S6: Contract area by depth
# -----------------------------------------------------------------------------------------------------------------

data <- read.csv("test_boxplot_depth_2023-02-18.csv")
data2 <- data %>%
  filter(Area != "")
dim(data2)
dim(data)
ggplot(data, mapping = aes(x = Area, y = Depth, 
                                   fill = Area)) + geom_boxplot() + coord_flip() + theme_bw() +
  scale_fill_fish_d(option = "Cirrhilabrus_solorensis") +
  theme(text=element_text(size=14), #change font size of all text
        axis.text=element_text(size=14), #change font size of axis text
        axis.title=element_text(size=14), #change font size of axis titles
        plot.title=element_text(size=14), #change font size of plot title
        legend.text=element_text(size=14), #change font size of legend text
        legend.title=element_text(size=14)) + 
  theme(legend.position = "none")

###########################
