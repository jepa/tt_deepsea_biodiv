library(picante)
library(vegan)
library(tidyverse)

# -----------------------------------------------------------------------------------------------------------------
# Data import & preparation
# -----------------------------------------------------------------------------------------------------------------
data <- read.csv("data-raw/TEST_LIT+DD_ALL_SPP_2022-10-10.csv",header=T,sep=",",fileEncoding="latin1")[-2,]
community_matrix <- picante::sample2matrix(data)


test <- data %>% 
      dplyr::group_by(Site, Species) %>% 
      dplyr::summarise(cover = Abundance) %>%  
      reshape::cast(.,  Site ~ Species, value = "cover")
test[is.na(test)] <-  0
rownames(test) <- test$Site
test <- test[-1, ]
# -----------------------------------------------------------------------------------------------------------------
# Figure: Raw diversity
# -----------------------------------------------------------------------------------------------------------------







# -----------------------------------------------------------------------------------------------------------------
# Figure: Family/species accumulation by sampling effort
# -----------------------------------------------------------------------------------------------------------------
curve <-  data.frame(sites = specaccum(test[, -1])$sites, richness = specaccum(test[, -1])$richness, sd = specaccum(test[, -1])$sd)
ggplot(curve) +
      geom_ribbon(aes(sites, ymin = richness-sd, ymax = richness + sd), alpha = .5) +
      geom_line(aes(sites, richness)) +
      xlab("Site") +
      ylab("Species richness") +
      theme(legend.title = element_blank(),
            legend.position = c(.9, .5))


# -----------------------------------------------------------------------------------------------------------------
# Figure: Rarefaction curves
# -----------------------------------------------------------------------------------------------------------------


rarecurve_var <- vegan::rarecurve(data.matrix) ## abundance based- whole CCZ
rarefy_var <- vegan::rarefy(data.matrix, 90611) ## arg 'sample missing
rrarefy_var <- vegan::rrarefy(data.matrix, 93456) ## arg 'sample missing
test <- vegan::estimateR(data.matrix) ## abund of single sample - FOR CCZ ONLY DATASET

## gives ACE estimator - what about chao2, ES2000

## sample based
test <- specpool(data.matrix)## sample based
test <- poolaccum(data.matrix) ## try pool accum- this is sample based
write.csv(as.data.frame(pool.aurora, "OUTPUT.csv"))

rarefy(data.matrix, 93456, se= FALSE, MARGIN = 1)
