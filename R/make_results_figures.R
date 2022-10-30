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
specaccum_df <- read.csv('data-processed/CCZ_specaccum.csv')
load('data-processed/iNEXT_abundance.RData')

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
                          labels = c("all descriptions", "new species", " publications"))

ggsave(descriptions_figure, filename = 'output-figures/descriptions_figure.jpg', 
       width = 15, height = 15, units = 'cm', dpi = 150)



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

# rarefy from vegan- test

S <- specnumber(community_matrix)
raremax <- min(rowSums(community_matrix))
Srare <- rarefy(community_matrix, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Species")
abline(0, 1)
rarecurve(community_matrix, step = 20, sample = raremax)


