library(vegan)
library(tidyverse)
library(iNEXT)
library(ggthemes)
theme_set(ggthemes::theme_few(base_size = 12))

# -----------------------------------------------------------------------------------------------------------------
# Data import
# -----------------------------------------------------------------------------------------------------------------
species_descriptions <- read.csv("data-raw/TEMP_papers_table_1980_on_2022-09-21.csv") %>% 
      pivot_longer(cols = 6:8, names_to = "var", values_to = "number")
community_matrix <- read.table('data-processed/CCZ_community_matrix.txt')
load('data-processed/CCZ_com_matrix_standardised.RData') 
specaccum_df <- read.csv('data-processed/CCZ_specaccum.csv')

# -----------------------------------------------------------------------------------------------------------------
# Figure: Raw data
# -----------------------------------------------------------------------------------------------------------------
descriptions_figure <- ggplot(species_descriptions, aes(x = Year, y = number, colour = var, fill = var)) +
      geom_line(size = 1) +
      labs(x = "Year", y = "Cumulative totals") +
      theme(legend.justification = c(0, 1), 
            legend.position = c(0, 1),
            legend.box.margin=margin(c(50, 50, 50, 50)))+
      scale_colour_colorblind("", 
                          breaks = c("cumul_desc", "cumul_spp", "cumul_pubs"),
                          labels = c("all descriptions", "new species", " publications")); descriptions_figure

# ggsave(descriptions_figure, filename = 'output-figures/descriptions_figure.jpg', 
#        width = 15, height = 15, units = 'cm', dpi = 150)

# -----------------------------------------------------------------------------------------------------------------
# Figure: Family/species accumulation by sampling effort
# -----------------------------------------------------------------------------------------------------------------
specaccum_stand <- specaccum(com_matrix_standardised, method = 'random', permutations = 100)  
specaccum_stand_df <-  data.frame(sites = specaccum_stand$sites, richness = specaccum_stand$richness, sd = specaccum_stand$sd)

sample_based_accum <- ggplot(specaccum_df) +
      geom_ribbon(aes(sites, ymin = richness-sd, ymax = richness + sd), alpha = .3, fill = 'blue') +
      geom_line(aes(sites, richness)) +
      xlab("Number of samples") +
      ylab("Species richness"); sample_based_accum

sample_based_accum_stand <- ggplot(specaccum_stand_df) +
      geom_ribbon(aes(sites, ymin = richness-sd, ymax = richness + sd), alpha = .3, fill = 'blue') +
      geom_line(aes(sites, richness)) +
      xlab("Number of samples") +
      ylab("Species richness"); sample_based_accum_stand

# -----------------------------------------------------------------------------------------------------------------
# Figure: Rarefaction curves
# -----------------------------------------------------------------------------------------------------------------
# Overall estimate for CZZ



