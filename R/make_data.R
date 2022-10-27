library(picante)
library(tidyverse)

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

# Export
write.table(community_matrix, file = 'data-processed/CCZ_community_matrix.txt')
write.csv(specaccum_df, file = 'data-processed/CCZ_specaccum.csv')
save(i.out, file = 'data-processed/iNEXT_abundance.RData')


test <- iNEXT::ChaoRichness(community_matrix, datatype = "abundance", conf = 0.95)
test <- rowSums(community_matrix)
