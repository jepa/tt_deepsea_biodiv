library(picante)
library(tidyverse)
library(rgeos)
library(rgdal)
library(raster)
library(spatialEco)
library(iNEXT)

# -----------------------------------------------------------------------------------------------------------------
# Raw data
# -----------------------------------------------------------------------------------------------------------------
data_all_species <- read.csv("data-raw/CCZ_ALL_SPP_DATA_FIN_2022-11-08.csv", header = T, sep = ",", fileEncoding = "latin1") %>% 
      drop_na(ABUNDANCE)
data_all_species <- data_all_species[nzchar(data_all_species$SITE), ]
data_all_species <- data_all_species[nzchar(data_all_species$SPECIES), ]

data_CCZ_only <- data_all_species %>% mutate(site = "CCZ")

# Species data with coords etc.
data_coords <- read.csv('data-raw/CCZ_ALL_SPP_DATA_FIN_2022-11-08.csv') %>% 
      filter(!is.na(LONGITUDE), !is.na(LATITUDE)) %>% 
      rename(site = SITE, 
             lat = LATITUDE, 
             long = LONGITUDE) 

# -----------------------------------------------------------------------------------------------------------------
# Make community matrices
# -----------------------------------------------------------------------------------------------------------------
# By site
# community_matrix <- picante::sample2matrix(data_all_species) 
community_matrix <- data_all_species %>% 
      dplyr::group_by(SITE, SPECIES) %>% 
      dplyr::summarise(cover = sum(ABUNDANCE)) %>%  
      reshape::cast(.,  SITE ~ SPECIES, value = "cover")
community_matrix[is.na(community_matrix)] <-  0
rownames(community_matrix) <- community_matrix$SITE
community_matrix <- community_matrix[, !(names(community_matrix) %in% "SITE")]

empty_sites <- rownames(community_matrix[rowSums(community_matrix) == 0, ])
community_matrix <- community_matrix[!rownames(community_matrix) %in% empty_sites, ]

length(names(sapply(community_matrix, is.numeric))) == length(colnames(community_matrix))

# community_matrix_CCZ <- picante::sample2matrix(data_CCZ_only) 
community_matrix_CCZ <- data_CCZ_only %>% 
dplyr::group_by(site, SPECIES) %>% 
      dplyr::summarise(cover = sum(ABUNDANCE)) %>%  
      reshape::cast(.,  site ~ SPECIES, value = "cover")
community_matrix_CCZ[is.na(community_matrix_CCZ)] <-  0
rownames(community_matrix_CCZ) <- community_matrix_CCZ$site
community_matrix_CCZ <- community_matrix_CCZ[, !(names(community_matrix_CCZ) %in% "site")]

# By grid, abundance standardised by pseudo-effort (number of sampling events in grid)
sites_df <- data.frame(long = data_coords$long, lat = data_coords$lat, site = data_coords$site) %>% 
      unique()

sites_spdf <- SpatialPointsDataFrame(coords = cbind(sites_df$long, sites_df$lat), data = data.frame(sites_df$site))
proj4string(sites_spdf) <- CRS("+init=epsg:4326")

CCZ_outline <- readWKT(text = "POLYGON ((-160.3477 12.5906, -154.0898 -5.4997, 
                        -107.5430 8.6544, -121.8867 23.3343, -160.3477 12.5906))")
proj4string(CCZ_outline) <- CRS("+init=epsg:4326")
grid <- raster(extent(CCZ_outline), resolution = c(5, 5), crs = proj4string(CCZ_outline))
gridPolygon <- rasterToPolygons(grid)
gridPolygon$grid_ID <- 1:nrow(gridPolygon)
grid_clipped <- raster::intersect(gridPolygon, CCZ_outline)
intersectGrid <- gridPolygon[gridPolygon$grid_ID %in% grid_clipped$grid_ID, ]
test <- point.in.poly(sites_spdf, intersectGrid)
plot(CCZ_outline)
plot(grid_clipped, add = T)
plot(sites_spdf, add = T)

points_intersected <- data.frame(point.in.poly(sites_spdf, intersectGrid)) %>% 
      rename(site = sites_df.site, long = coords.x1, lat = coords.x2) %>% 
      dplyr::select(site, long, lat, grid_ID) %>% 
      group_by(grid_ID) %>% 
      mutate(pseudo_effort = length(site))

data_merged <- merge(points_intersected, data_coords, by = c("site", "long", "lat"), all.y = T)

# Community matrix by grid, standardised abundance
com_matrix_standardised <- data_merged %>% 
      dplyr::group_by(grid_ID, SPECIES) %>% 
      dplyr::summarise(cover = sum(ABUNDANCE)/pseudo_effort) %>%  
      unique() %>% 
      reshape::cast(.,  grid_ID ~ SPECIES, value = "cover")
com_matrix_standardised[is.na(com_matrix_standardised)] <-  0
rownames(com_matrix_standardised) <- com_matrix_standardised$grid_ID

# -----------------------------------------------------------------------------------------------------------------
# Species curve data - some take long to compile
# -----------------------------------------------------------------------------------------------------------------
specaccum_sites <- specaccum(community_matrix, method = 'random', permutations = 100)  
specaccum_sites_df <-  data.frame(sites = specaccum_sites$sites, richness = specaccum_sites$richness, sd = specaccum_sites$sd)

CCZ_rarecurve <- as.data.frame(rarecurve(community_matrix_CCZ, step = 20))
CCZ_rarecurve$Individuals <- as.integer(gsub('N', '', rownames(CCZ_rarecurve)))
colnames(CCZ_rarecurve)[1] <- "Species"

com_matT <- as.matrix(t(community_matrix_CCZ))
Hills_q_CCZ <- iNEXT(com_matT, q = 0, datatype = "abundance", nboot = 2)
ChaoRichness((com_matT))
com_mat_inc <- community_matrix
com_mat_inc <- ifelse(com_mat_inc > 0, 1, 0)
inc_freq <- as.incfreq(t(com_mat_inc))
ChaoRichness((inc_freq))
Hills_q_0_inc_CCZ <- iNEXT(inc_freq, q=0, datatype = "incidence_freq", nboot = 2)

# -----------------------------------------------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------------------------------------------
write.table(community_matrix, file = 'data-processed/CCZ_community_matrix.txt')
write.table(community_matrix_CCZ, file = 'data-processed/community_matrix_CCZ.txt')

save(com_matrix_standardised, file = 'data-processed/CCZ_com_matrix_standardised.RData')  
write.csv(specaccum_sites_df, file = 'data-processed/CCZ_specaccum_sites.csv')
write.csv(CCZ_rarecurve, file = 'data-processed/CCZ_rarecurve.csv')
save(Hills_q_CCZ, file = 'data-processed/Hills_q_CCZ.RData')
save(Hills_q_0_inc_CCZ, file = 'data-processed/Hills_q_0_inc.RData')

# writeOGR(intersectGrid, dsn = 'data-processed', layer = 'CCZ_grid_5degree', driver = "ESRI Shapefile")
# write.csv(data_merged, file = 'data-processed/CCZ_specdata_pseudoeffort.csv')


