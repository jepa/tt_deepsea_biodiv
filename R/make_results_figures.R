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
load('data-processed/CCZ_com_matrix_standardised.RData') 
specaccum_df <- read.csv('data-processed/CCZ_specaccum.csv')

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
                          labels = c("all descriptions", "new species", " publications")); descriptions_figure

# ggsave(descriptions_figure, filename = 'output-figures/descriptions_figure.jpg', 
#        width = 15, height = 15, units = 'cm', dpi = 150)

# -----------------------------------------------------------------------------------------------------------------
# Figure: Family/species accumulation by sampling effort
# -----------------------------------------------------------------------------------------------------------------
specaccum_stand <- specaccum(com_matrix_standardised, method = 'random', permutations = 100)  
specaccum_stand_df <-  data.frame(sites = specaccum_stand$sites, richness = specaccum_stand$richness, sd = specaccum_stand$sd)

sample_based_accum <- ggplot(specaccum_df) +
      geom_ribbon(aes(sites, ymin = richness-sd, ymax = richness + sd), alpha = .3, fill = 'blue') +
      geom_line(aes(sites, richness)) +
      xlab("Number of samples") +
      ylab("Species richness"); sample_based_accum

sample_based_accum_stand <- ggplot(specaccum_stand_df) +
      geom_ribbon(aes(sites, ymin = richness-sd, ymax = richness + sd), alpha = .3, fill = 'blue') +
      geom_line(aes(sites, richness)) +
      xlab("Number of samples") +
      ylab("Species richness"); sample_based_accum_stand

# -----------------------------------------------------------------------------------------------------------------
# Figure: Rarefaction curves
# -----------------------------------------------------------------------------------------------------------------
# Overall estimate for CZZ



# -----------------------------------------------------------------------------------------------------------------
# Mu's hodge podge code corner of doom
# -----------------------------------------------------------------------------------------------------------------
# 5TH NOVEMBER = LONG LIVE GUY FAWKES

# log of total number of taxa by taxonomic level- ie no of phlya- to genera - log this and use slope to get a prdicted species no


data <- read.csv("CCZ_CHECKLIST_2022-11-04.CSV")

## get total per taxon level- ie total phyla- order-calss-fam-genera recorded
data2 <- tapply(data[, "scientificName"], data[, "taxonRank"],
                function(x) { length(unique(x))})

## output summary as csv
write.csv(data2, "OUTPUT.csv")

## log of total per taxa level
test <- log(data$Total)
## append this as a column in file (could do within script- here I did in excel)
write.csv(test, "OUTPUT.csv")

## read in file with column of log value and taxon order as number- 1 - phlya, 2- class, 3 - order, 4 fam, 5 genus
data <- read.csv("temp_summary_checklist_taxonrank.csv")

## model
model1 <- lm(log ~ order, data = data)

autoplot(model1, smooth.colour = NA)
anova(model1)
summary(model1)

## PLOT
plot(data$olog, data$order, pch = 16)
## add line to basic plot -  line best fit with linear model
abline(lm(log ~ order, data = data)) ## 

## how to get log value for order = 6 (i.e. next taxon level- sp will be 6)
fitted <- predict(lm(data$log ~ data$order))

data$order[7];data$log[7]
order[7];log[7]; data = data


# -----------------------------------------------------------------------------------------------------------------
# scrappy notes for iNEXT

## transposed input file

com_matT <- t(com_mat)

## richness estimates
test <- ChaoRichness((com_matT))

## full version- q nos 0-5, bootstrap at 100 (will take ages to run- set up for overnight)
Hills_q <- iNEXT(com_matT, q=0, datatype = "abundance", nboot = 2)
# Hills_q <- iNEXT(com_matT, q=c(0,1,2,3,4,5), datatype = "abundance", nboot = 100)
# saveRDS(Hills_q, "Hills_q-5_orders.rds")
# Hills_q <- readRDS("Hills_q-5_orders.rds")
# Hills_df <- bind_rows(Hills_q$iNextEst, .id="site")

# plot
ggiNEXT(Hills_q, type=1)

plot_q <- ggiNEXT(Hills_q, type=1)+
   theme_clean()