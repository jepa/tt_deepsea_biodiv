library(tidyverse)
library(patchwork)
library(vegan)
library(ggthemes)
library(ggstance)
library(picante)

theme_set(theme_bw(base_size = 12))

# -----------------------------------------------------------------------------------------------------------------
# Data import
# -----------------------------------------------------------------------------------------------------------------
full_data <- read.csv('data-raw/CCZ_ALL_SPP_DATA_FIN_2022-11-08.csv') %>% 
      filter(!is.na(LONGITUDE), !is.na(LATITUDE)) %>% 
      rename(site = SITE, 
             lat = LATITUDE, 
             long = LONGITUDE) 

species_descriptions <- read.csv("data-raw/ARCHIVED_DATA/TEMP_papers_table_1980_on_2022-11-05.csv") %>% 
      pivot_longer(cols = 6:8, names_to = "var", values_to = "number")

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
data_CCZ_only <- full_data %>% mutate(Site = "CCZ")
community_matrix_CCZ <- picante::sample2matrix(data_CCZ_only) 









