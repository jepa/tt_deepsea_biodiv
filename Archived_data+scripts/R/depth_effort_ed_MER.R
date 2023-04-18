library(tidyverse)
library(patchwork)

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

depth_group_width <- ggplot(effort_data) +
      geom_rect(data = effort_data %>% select(depth_group_width_min, depth_group_width_max, depth_group_width) %>% unique(),
                aes(ymin = depth_group_width_min, ymax = depth_group_width_max,
                    xmin = 0, xmax = 0.003, fill = depth_group_width), col = "black", alpha = 0.6) +
      geom_density(aes(y = depth), fill = "coral1")+
      scale_y_reverse() +
      ggtitle("Number of samples by depth, grouped by 100m intervals") +
      scale_fill_gradient(low = "lightblue", high = "dodgerblue4") +
      theme_bw(); depth_group_width


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

