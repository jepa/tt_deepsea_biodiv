library(tidyverse)
library(patchwork)
library(vegan)
library(ggthemes)
library(ggstance)
library(picante)
library(iNEXT)
options(scipen=999)
theme_set(theme_bw(base_size = 12))

# -----------------------------------------------------------------------------------------------------------------
# Data import
# -----------------------------------------------------------------------------------------------------------------
# All raw data
full_data <- read.csv('data-raw/CCZ_ALL_TAXA_DATA_FIN_2023-02-01.csv') %>%
      rename_all(tolower) %>% 
      mutate(abundance = as.numeric(abundance))

# Species abundance by site
species_abundance <- full_data %>% 
      filter(!is.na(species), species != "") %>% 
      group_by(site, species) %>% 
      summarise(abundance = sum(abundance)) %>% 
      filter(!is.na(abundance))

# Family abundance by site
family_abundance <- full_data %>% 
      filter(!is.na(family), family != "") %>% 
      group_by(site, family) %>% 
      summarise(abundance = sum(abundance)) %>% 
      filter(!is.na(abundance))

# Genus abundance by site
genus_abundance <- full_data %>% 
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

# Family presence-absence by site
family_presAbs <- full_data %>% 
      filter(!is.na(family), family != "") %>% 
      group_by(site, family) %>% 
      summarise(presence = ifelse(abundance > 0, 1, 0)) %>% 
      filter(!is.na(presence))

# Genus presence-absence by site
genus_presAbs <- full_data %>% 
      filter(!is.na(genus), genus != "") %>% 
      group_by(site, genus) %>% 
      summarise(presence = ifelse(abundance > 0, 1, 0)) %>% 
      filter(!is.na(presence))

# Species description data
species_descriptions <- read.csv("data-raw/ARCHIVED_DATA/TEMP_papers_table_1980_on_2022-11-05.csv") %>% 
      pivot_longer(cols = 6:8, names_to = "var", values_to = "number")

# Phyla data
phyla_overview <- read.csv("data-raw/ARCHIVED_DATA/TEMP_SUMMARY_FIG2_ALL_PHYLA_2022-11-05.csv") %>% 
      group_by(Phylum) %>% 
      mutate(total_phylum = sum(Total)) %>% 
      group_by(Data) %>% 
      mutate(perc = paste0(round((Total/total_phylum)*100), "%"))
phyla_overview$ypos <- ifelse(phyla_overview$Data == "named species", phyla_overview$total_phylum, phyla_overview$total_phylum + 70)

# -----------------------------------------------------------------------------------------------------------------
# Figure 2: Taxonomic work
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

# -----------------------------------------------------------------------------------------------------------------
# Figure 3: Species accumulation curves for CCZ
# -----------------------------------------------------------------------------------------------------------------
# Community matrix for species abundance by single site: CCZ
species_CCZ_only <- species_abundance %>% 
      mutate(site = "CCZ")

species_matrix_CCZ <- species_CCZ_only %>%
dplyr::group_by(site, species) %>%
      dplyr::summarise(cover = sum(abundance)) %>%
      reshape::cast(.,  site ~ species, value = "cover")
species_matrix_CCZ[is.na(species_matrix_CCZ)] <-  0

species_matT <- as.matrix(t(species_matrix_CCZ))
Hills_q_CCZ_species <- iNEXT::iNEXT(species_matT, q = 0, datatype = "abundance", nboot = 10)
Hills_q_CCZ_species_df <- as.data.frame(Hills_q_CCZ_species$iNextEst$size_based)
   
A_Chao1_species <- ggplot(Hills_q_CCZ_species_df %>% filter(Method != 'Observed'), aes(x = m)) +
      geom_ribbon(aes(x = m, ymin = qD.LCL, ymax = qD.UCL), fill = "coral2", alpha = 0.2) +
      geom_line(aes(y = qD, lty = Method),
                col = "coral2", cex = 1) +
      geom_point(data = Hills_q_CCZ_species_df %>% filter(Method == 'Observed'), aes(x = m, y = qD), col = "coral2", cex = 2) +
      theme(legend.justification = c(0, 1), 
            legend.position = c(.5, .4),
            legend.box.margin=margin(c(20, 20, 20, 20)),
            legend.margin = margin(1, 1, 1, 1),
            legend.spacing.x = unit(0, "mm"),
            legend.spacing.y = unit(0, "mm"),
            legend.title = element_blank(),
            legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
      xlab("Number of Individuals") +
      ylab("Species Diversity") +
      # scale_y_continuous(breaks = seq(0, 5000, 1000)) +
      scale_linetype_manual(values = c("dashed", "solid")); A_Chao1_species

# Community matrix for species presence-absence by all sites (the sampling units)
species_matrix_pres <- species_presAbs %>%
      dplyr::group_by(site, species) %>%
      dplyr::summarise(presence = ifelse(sum(presence) > 0, 1, 0)) %>%
      reshape::cast(.,  site ~ species, value = "presence")
species_matrix_pres[is.na(species_matrix_pres)] <-  0
inc_freq_spec <- as.incfreq(as.matrix(t(species_matrix_pres)))
Hills_q_0_inc_CCZ_spec <- iNEXT(inc_freq_spec, q = 0, datatype = "incidence_freq", nboot = 2)
Hills_q_0_inc_CCZ_spec_df <- as.data.frame(Hills_q_0_inc_CCZ_spec$iNextEst$size_based)

B_Chao2_species <- ggplot(Hills_q_0_inc_CCZ_spec_df %>% filter(Method != 'Observed')) +
      geom_ribbon(aes(x = t, ymin = qD.LCL, ymax = qD.UCL), fill = "coral2", alpha = 0.2, show.legend = FALSE) +      
      geom_line(aes(x = t, y = qD, lty = Method), col = "coral2", cex = 1) +
      geom_point(data = Hills_q_0_inc_CCZ_spec_df %>% filter(Method == 'Observed'), aes(x = t, y = qD), col = "coral2", cex = 2) +
      theme(legend.justification = c(0, 1), 
            legend.position = c(.5, .4),
            legend.box.margin=margin(c(20, 20, 20, 20)),
            legend.margin = margin(1, 1, 1, 1),
            legend.spacing.x = unit(0, "mm"),
            legend.spacing.y = unit(0, "mm"),
            legend.title = element_blank(),
            legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
      xlab("Number of sampling units") +
      ylab("Species Diversity") +
      # scale_y_continuous(breaks = seq(0, 6000, 1000)) +
      scale_linetype_manual(values = c("dashed", "solid")); B_Chao2_species

# Community matrix for family abundance by single site: CCZ
family_CCZ_only <- family_abundance %>% 
      mutate(site = "CCZ")

family_matrix_CCZ <- family_CCZ_only %>%
      dplyr::group_by(site, family) %>%
      dplyr::summarise(cover = sum(abundance)) %>%
      reshape::cast(.,  site ~ family, value = "cover")
family_matrix_CCZ[is.na(family_matrix_CCZ)] <-  0

family_matT <- as.matrix(t(family_matrix_CCZ))
Hills_q_CCZ_family <- iNEXT::iNEXT(family_matT, q = 0, datatype = "abundance", nboot = 10)
Hills_q_CCZ_family_df <- as.data.frame(Hills_q_CCZ_family$iNextEst$size_based)

C_Chao1_family <- ggplot(Hills_q_CCZ_family_df %>% filter(Method != 'Observed'), aes(x = m)) +
      geom_ribbon(aes(x = m, ymin = qD.LCL, ymax = qD.UCL), fill = "coral2", alpha = 0.2) +
      geom_line(aes(y = qD, lty = Method),
                col = "coral2", cex = 1) +
      geom_point(data = Hills_q_CCZ_family_df %>% filter(Method == 'Observed'), aes(x = m, y = qD), col = "coral2", cex = 2) +
      theme(legend.justification = c(0, 1), 
            legend.position = c(.5, .4),
            legend.box.margin=margin(c(20, 20, 20, 20)),
            legend.margin = margin(1, 1, 1, 1),
            legend.spacing.x = unit(0, "mm"),
            legend.spacing.y = unit(0, "mm"),
            legend.title = element_blank(),
            legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
      xlab("Number of Individuals") +
      ylab("Family Diversity") +
      # scale_y_continuous(breaks = seq(0, 5000, 1000)) +
      scale_linetype_manual(values = c("dashed", "solid")); C_Chao1_family

# Community matrix for family presence-absence by all sites (the sampling units)
family_matrix_pres <- family_presAbs %>%
      dplyr::group_by(site, family) %>%
      dplyr::summarise(presence = ifelse(sum(presence) > 0, 1, 0)) %>%
      reshape::cast(.,  site ~ family, value = "presence")
family_matrix_pres[is.na(family_matrix_pres)] <-  0
inc_freq_family <- as.incfreq(as.matrix(t(family_matrix_pres)))
Hills_q_0_inc_CCZ_family <- iNEXT(inc_freq_family, q = 0, datatype = "incidence_freq", nboot = 2)
Hills_q_0_inc_CCZ_family_df <- as.data.frame(Hills_q_0_inc_CCZ_family$iNextEst$size_based)

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
      xlab("Number of sampling units") +
      ylab("Family Diversity") +
      # scale_y_continuous(breaks = seq(0, 6000, 1000)) +
      scale_linetype_manual(values = c("dashed", "solid")); D_Chao2_family

# Community matrix for genus abundance by single site: CCZ
genus_CCZ_only <- genus_abundance %>% 
      mutate(site = "CCZ")

genus_matrix_CCZ <- genus_CCZ_only %>%
      dplyr::group_by(site, genus) %>%
      dplyr::summarise(cover = sum(abundance)) %>%
      reshape::cast(.,  site ~ genus, value = "cover")
genus_matrix_CCZ[is.na(genus_matrix_CCZ)] <-  0

genus_matT <- as.matrix(t(genus_matrix_CCZ))
Hills_q_CCZ_genus <- iNEXT::iNEXT(genus_matT, q = 0, datatype = "abundance", nboot = 10)
Hills_q_CCZ_genus_df <- as.data.frame(Hills_q_CCZ_genus$iNextEst$size_based)

E_Chao1_genus <- ggplot(Hills_q_CCZ_genus_df %>% filter(Method != 'Observed'), aes(x = m)) +
      geom_ribbon(aes(x = m, ymin = qD.LCL, ymax = qD.UCL), fill = "coral2", alpha = 0.2) +
      geom_line(aes(y = qD, lty = Method),
                col = "coral2", cex = 1) +
      geom_point(data = Hills_q_CCZ_genus_df %>% filter(Method == 'Observed'), aes(x = m, y = qD), col = "coral2", cex = 2) +
      theme(legend.justification = c(0, 1), 
            legend.position = c(.5, .4),
            legend.box.margin=margin(c(20, 20, 20, 20)),
            legend.margin = margin(1, 1, 1, 1),
            legend.spacing.x = unit(0, "mm"),
            legend.spacing.y = unit(0, "mm"),
            legend.title = element_blank(),
            legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
      xlab("Number of Individuals") +
      ylab("Genus Diversity") +
      # scale_y_continuous(breaks = seq(0, 5000, 1000)) +
      scale_linetype_manual(values = c("dashed", "solid")); E_Chao1_genus

# Community matrix for genus presence-absence by all sites (the sampling units)
genus_matrix_pres <- genus_presAbs %>%
      dplyr::group_by(site, genus) %>%
      dplyr::summarise(presence = ifelse(sum(presence) > 0, 1, 0)) %>%
      reshape::cast(.,  site ~ genus, value = "presence")
genus_matrix_pres[is.na(genus_matrix_pres)] <-  0
inc_freq_genus <- as.incfreq(as.matrix(t(genus_matrix_pres)))
Hills_q_0_inc_CCZ_genus <- iNEXT(inc_freq_genus, q = 0, datatype = "incidence_freq", nboot = 2)
Hills_q_0_inc_CCZ_genus_df <- as.data.frame(Hills_q_0_inc_CCZ_genus$iNextEst$size_based)

F_Chao2_genus <- ggplot(Hills_q_0_inc_CCZ_genus_df %>% filter(Method != 'Observed')) +
      geom_ribbon(aes(x = t, ymin = qD.LCL, ymax = qD.UCL), fill = "coral2", alpha = 0.2, show.legend = FALSE) +      
      geom_line(aes(x = t, y = qD, lty = Method), col = "coral2", cex = 1) +
      geom_point(data = Hills_q_0_inc_CCZ_genus_df %>% filter(Method == 'Observed'), aes(x = t, y = qD), col = "coral2", cex = 2) +
      theme(legend.justification = c(0, 1), 
            legend.position = c(.5, .4),
            legend.box.margin=margin(c(20, 20, 20, 20)),
            legend.margin = margin(1, 1, 1, 1),
            legend.spacing.x = unit(0, "mm"),
            legend.spacing.y = unit(0, "mm"),
            legend.title = element_blank(),
            legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
      xlab("Number of sampling units") +
      ylab("Genus Diversity") +
      # scale_y_continuous(breaks = seq(0, 6000, 1000)) +
      scale_linetype_manual(values = c("dashed", "solid")); F_Chao2_genus


figure_3 <- (A_Chao1_species | B_Chao2_species) / (C_Chao1_family | D_Chao2_family)  / (E_Chao1_genus | F_Chao2_genus) + plot_annotation(tag_levels = 'A')

ggsave(figure_3,
       filename = 'output-figures/figure_3.png', 
       width = 10, height = 10, units = 'in', dpi = 150)


# Can you output the curves and also the results e.g. the chao2 richness estimates?
#       I’m getting >8000 for spp, >550 for fam and >1000 for genera, 
# NB had wrong input data for family- reran and is 485
Hills_q_0_inc_CCZ_spec$AsyEst
Hills_q_0_inc_CCZ_family$AsyEst
Hills_q_0_inc_CCZ_genus$AsyEst
