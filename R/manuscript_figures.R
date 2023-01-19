library(tidyverse)
library(patchwork)
library(vegan)
library(ggthemes)
library(ggstance)
library(picante)
library(iNEXT)

theme_set(theme_bw(base_size = 12))

# -----------------------------------------------------------------------------------------------------------------
# Data import
# -----------------------------------------------------------------------------------------------------------------
# All raw data
full_data <- read.csv('data-raw/CCZ_ALL_SPP_DATA_FIN_2022-11-08.csv') %>% 
      filter(!is.na(LONGITUDE), !is.na(LATITUDE))
names(full_data) <- tolower(names(full_data))

# Species abundance by site
species_abundance <- full_data %>% 
      group_by(site, species) %>% 
      summarise(abundance = sum(abundance)) %>% 
      drop_na()

# Species presence-absence by site
species_presAbs <- full_data %>% 
      group_by(site, species) %>% 
      summarise(presence = ifelse(abundance > 0, 1, 0)) %>% 
      drop_na()

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
data_CCZ_only <- species_abundance %>% 
      mutate(site = "CCZ")

community_matrix_CCZ <- data_CCZ_only %>%
dplyr::group_by(site, species) %>%
      dplyr::summarise(cover = sum(abundance)) %>%
      reshape::cast(.,  site ~ species, value = "cover")
community_matrix_CCZ[is.na(community_matrix_CCZ)] <-  0

com_matT <- as.matrix(t(community_matrix_CCZ))
Hills_q_CCZ <- iNEXT::iNEXT(com_matT, q = 0, datatype = "abundance", nboot = 10)
Hills_q_CCZ_df <- as.data.frame(Hills_q_CCZ$iNextEst$size_based)
   

A_Chao1 <- ggplot(Hills_q_CCZ_df %>% filter(Method != 'Observed'), aes(x = m)) +
      geom_ribbon(aes(x = m, ymin = qD.LCL, ymax = qD.UCL), fill = "coral2", alpha = 0.2) +
      geom_line(aes(y = qD, lty = Method),
                col = "coral2", cex = 1) +
      geom_point(data = Hills_q_CCZ_df %>% filter(Method == 'Observed'), aes(x = m, y = qD), col = "coral2", cex = 2) +
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
      scale_y_continuous(breaks = seq(0, 5000, 1000)) +
      scale_linetype_manual(values = c("dashed", "solid")); A_Chao1

# Community matrix for species presence-absence by all sites (the sampling units)
community_matrix_pres <- species_presAbs %>%
      dplyr::group_by(site, species) %>%
      dplyr::summarise(presence = ifelse(sum(presence) > 0, 1, 0)) %>%
      reshape::cast(.,  site ~ species, value = "presence")
community_matrix_pres[is.na(community_matrix_pres)] <-  0
inc_freq <- as.incfreq(as.matrix(t(community_matrix_pres)))
Hills_q_0_inc_CCZ <- iNEXT(inc_freq, q = 0, datatype = "incidence_freq", nboot = 2)
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
      xlab("Number of sampling units") +
      ylab("Species Diversity") +
      scale_y_continuous(breaks = seq(0, 6000, 1000)) +
      scale_linetype_manual(values = c("dashed", "solid")); B_Chao2






