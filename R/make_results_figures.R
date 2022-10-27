library(vegan)
library(tidyverse)
library(iNEXT)
library(ggthemes)
theme_set(ggthemes::theme_few(base_size = 12))

# -----------------------------------------------------------------------------------------------------------------
# Data import
# -----------------------------------------------------------------------------------------------------------------
species_descriptions <- read.csv("TEMP_papers_table_1980_on_2022-09-21.csv")  # Missing
community_matrix <- read.table('data-processed/CCZ_community_matrix.txt')
specaccum_df <- read.csv('data-processed/CCZ_specaccum.csv')
load('data-processed/iNEXT_abundance.RData')

# -----------------------------------------------------------------------------------------------------------------
# Figure: Raw data
# -----------------------------------------------------------------------------------------------------------------
ggplot(species_descriptions, aes(x = Year)) +
      geom_line(aes(y = cumul_desc, colour = "cumul_desc"), size = 1) +
      geom_line(aes(y = cumul_spp, colour = "cumul_spp"), size = 1) +
      geom_line(aes(y = cumul_pubs, colour = "cumul_pubs"), size = 1) +
      labs(x = "Year", y = "Cumulative totals") +
      theme(legend.position = "none") +
      scale_colour_manual("", 
                          breaks = c("cumul_desc", "cumul_spp", "cumul_pubs"),
                          values = c("black", "coral2", "steelblue"),
                          labels = c("all descriptions", "new species", " publications"))




# -----------------------------------------------------------------------------------------------------------------
# Figure: Family/species accumulation by sampling effort
# -----------------------------------------------------------------------------------------------------------------
sample_based_accum <- ggplot(specaccum_df) +
      geom_ribbon(aes(sites, ymin = richness-sd, ymax = richness + sd), alpha = .3, fill = 'blue') +
      geom_line(aes(sites, richness)) +
      xlab("Number of samples") +
      ylab("Species richness"); sample_based_accum


# -----------------------------------------------------------------------------------------------------------------
# Figure: Rarefaction curves
# -----------------------------------------------------------------------------------------------------------------
# Overall estimate for CZZ

ggiNEXT(i.out, type = 3)


