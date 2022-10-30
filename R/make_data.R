library(picante)
library(tidyverse)
library(rgeos)
library(rgdal)
library(raster)
library(spatialEco)

# Community matrix
data <- read.csv("data-raw/TEST_LIT+DD_ALL_SPP_2022-10-10.csv", header = T, sep = ",", fileEncoding = "latin1")[-2, ] %>% 
      drop_na(Abundance)
data <- data[nzchar(data$Site), ]
community_matrix <- picante::sample2matrix(data) 

# Remove sites without any species records
empty_sites <- rownames(community_matrix[rowSums(community_matrix) == 0, ])
community_matrix <- community_matrix[!rownames(community_matrix) %in% empty_sites, ]

# Check for non-numeric columns
length(names(sapply(community_matrix, is.numeric))) == length(colnames(community_matrix))

# # iNEXT data
# i.out <- iNEXT(community_matrix, q = 0, datatype = "abundance")

# Species data with coords etc.
data_coords <- read.csv('data-raw/TEMP_SPECIES_ALL_v3_2022-09-22.csv') %>% 
      filter(!is.na(decimalLongitude), !is.na(decimalLatitude), !is.na(individualCount)) %>% 
      rename(lat = decimalLatitude, 
             long = decimalLongitude) %>% 
      group_by(long, lat) %>% 
      mutate(site_ID = cur_group_id())

sites_df <- data.frame(long = data_coords$long, lat = data_coords$lat, site_ID = data_coords$site_ID) %>% 
      unique()

sites_spdf <- SpatialPointsDataFrame(coords = cbind(sites_df$long, sites_df$lat), data = data.frame(sites_df$site_ID))
proj4string(sites_spdf) <- CRS("+init=epsg:4326")

# CCZ extent shapefile
CCZ_outline <- readWKT(text = "POLYGON ((-160.3477 12.5906, -154.0898 -5.4997, 
                        -107.5430 8.6544, -121.8867 23.3343, -160.3477 12.5906))")
proj4string(CCZ_outline) <- CRS("+init=epsg:4326")
grid <- raster(extent(CCZ_outline), resolution = c(5, 5), crs = proj4string(CCZ_outline))
gridPolygon <- rasterToPolygons(grid)
gridPolygon$grid_ID <- 1:nrow(gridPolygon)
grid_clipped <- raster::intersect(gridPolygon, CCZ_outline)
intersectGrid <- gridPolygon[gridPolygon$grid_ID %in% grid_clipped$grid_ID, ]

plot(CCZ_outline)
plot(intersectGrid, add = T)
plot(sites_spdf, add = T)


points_intersected <- data.frame(point.in.poly(sites_spdf, intersectGrid)) %>% 
      rename(site_ID = sites_df.site_ID, long = coords.x1, lat = coords.x2) %>% 
      dplyr::select(site_ID, long, lat, grid_ID) %>% 
      group_by(grid_ID) %>% 
      mutate(pseudo_effort = length(site_ID))

data_merged <- merge(points_intersected, data_coords, by = c("site_ID", "long", "lat"), all.y = T)

# Community matrix by grid, standardised abundance
com_matrix_standardised <- data_merged %>% 
      dplyr::group_by(grid_ID, NAME) %>% 
      dplyr::summarise(cover = sum(individualCount)/pseudo_effort) %>%  
      unique() %>% 
      reshape::cast(.,  grid_ID ~ NAME, value = "cover")
com_matrix_standardised[is.na(com_matrix_standardised)] <-  0
rownames(com_matrix_standardised) <- com_matrix_standardised$grid_ID

# Species curve data
specaccum_result <- specaccum(community_matrix, method = 'random', permutations = 100)  
specaccum_df <-  data.frame(sites = specaccum_result$sites, richness = specaccum_result$richness, sd = specaccum_result$sd)


# Export
write.table(community_matrix, file = 'data-processed/CCZ_community_matrix.txt')
write.table(com_matrix_standardised, file = 'data-processed/CCZ_com_matrix_standardised.txt')  # Some strange characters causing error, need to clea data

write.csv(specaccum_df, file = 'data-processed/CCZ_specaccum.csv')
writeOGR(intersectGrid, dsn = 'data-processed', layer = 'CCZ_grid_5degree', driver = "ESRI Shapefile")
# save(i.out, file = 'data-processed/iNEXT_abundance.RData')
write.csv(data_merged, file = 'data-processed/CCZ_specdata_pseudoeffort.csv')


