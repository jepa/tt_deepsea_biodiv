library(tidyverse)
library(patchwork)

data_coords <- read.csv('data-raw/CCZ_ALL_SPP_DATA_FIN_2022-11-08.csv') %>% 
      filter(!is.na(LONGITUDE), !is.na(LATITUDE)) %>% 
      rename(site = SITE, 
             lat = LATITUDE, 
             long = LONGITUDE,
             depth = DEPTH) %>% 
      mutate(year = sub('.*(\\d{4}).*', '\\1', DATE),
             depth = as.numeric(depth),
             year = as.integer(year),
             phylum = Phylum.) %>% 
      drop_na() %>% 
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

data_coords <- data_coords[nzchar(data_coords$phylum), ]


samples_by_depth <- ggplot(data_coords) +
      geom_density(aes(y = depth), fill = "lightblue") +
      # geom_hline(aes(yintercept = depth_group_min, colour = depth_group)) +
      scale_y_reverse() +
      ggtitle("Number of samples by depth, all phyla combined") +
      scale_colour_gradient(low = "lightblue", high = "dodgerblue4") +
      theme_bw(); samples_by_depth

depth_group_number <- ggplot(data_coords) +
      geom_rect(data = data_coords %>% select(depth_group_number_min, depth_group_number_max, depth_group_number) %>% unique(),
                aes(ymin = depth_group_number_min, ymax = depth_group_number_max,
                    xmin = 0, xmax = 0.003, fill = depth_group_number), col = "black", alpha = 0.6) +
      geom_density(aes(y = depth), fill = "coral1")+
      scale_y_reverse() +
      ggtitle("Number of samples by depth, grouped by 10 sample quantiles") +
      scale_fill_gradient(low = "lightblue", high = "dodgerblue4") +
      theme_bw(); depth_group_number

depth_group_interval <- ggplot(data_coords) +
      geom_rect(data = data_coords %>% select(depth_group_interval_min, depth_group_interval_max, depth_group_interval) %>% unique(),
                aes(ymin = depth_group_interval_min, ymax = depth_group_interval_max,
                    xmin = 0, xmax = 0.003, fill = depth_group_interval), col = "black", alpha = 0.6) +
      geom_density(aes(y = depth), fill = "coral1")+
      scale_y_reverse() +
      ggtitle("Number of samples by depth, grouped by 10 depth quantiles") +
      scale_fill_gradient(low = "lightblue", high = "dodgerblue4") +
      theme_bw(); depth_group_interval

depth_group_width <- ggplot(data_coords) +
      geom_rect(data = data_coords %>% select(depth_group_width_min, depth_group_width_max, depth_group_width) %>% unique(),
                aes(ymin = depth_group_width_min, ymax = depth_group_width_max,
                    xmin = 0, xmax = 0.003, fill = depth_group_width), col = "black", alpha = 0.6) +
      geom_density(aes(y = depth), fill = "coral1")+
      scale_y_reverse() +
      ggtitle("Number of samples by depth, grouped by 100m intervals") +
      scale_fill_gradient(low = "lightblue", high = "dodgerblue4") +
      theme_bw(); depth_group_width


depth_groups <- samples_by_depth / depth_group_number / depth_group_interval / depth_group_width



samples_by_depth_phyla <- ggplot(data_coords) +
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
