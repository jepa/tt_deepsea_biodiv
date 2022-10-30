library(picante)
library(tidyverse)
library(rgeos)
library(rgdal)
library(raster)

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

# Species curve data
specaccum_result <- specaccum(community_matrix, method = 'random', permutations = 100)  
specaccum_df <-  data.frame(sites = specaccum_result$sites, richness = specaccum_result$richness, sd = specaccum_result$sd)

# iNEXT data
i.out <- iNEXT(community_matrix, q = 0, datatype = "abundance")

# CCZ extent shapefile
CCZ_outline <- readWKT(text = "POLYGON ((-160.3477 12.5906, -154.0898 -5.4997, 
                        -107.5430 8.6544, -121.8867 23.3343, -160.3477 12.5906))")
proj4string(CCZ_outline) <- CRS("+init=epsg:4326")
grid <- raster(extent(CCZ_outline), resolution = c(5, 5), crs = proj4string(CCZ_outline))
gridPolygon <- rasterToPolygons(grid)
gridPolygon$id <- 1:nrow(gridPolygon)
grid_clipped <- raster::intersect(gridPolygon, CCZ_outline)
intersectGrid <- gridPolygon[gridPolygon$id %in% grid_clipped$id, ]

plot(CCZ_outline)
plot(intersectGrid, add = T)

# Export
write.table(community_matrix, file = 'data-processed/CCZ_community_matrix.txt')
write.csv(specaccum_df, file = 'data-processed/CCZ_specaccum.csv')
writeOGR(intersectGrid, dsn = 'data-processed', layer = 'CCZ_grid_5degree', driver = "ESRI Shapefile")
save(i.out, file = 'data-processed/iNEXT_abundance.RData')



test <- iNEXT::ChaoRichness(community_matrix, datatype = "abundance", conf = 0.95)
test <- rowSums(community_matrix)
