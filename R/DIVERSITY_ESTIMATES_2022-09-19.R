###########################

library(vegan)
library(ape)
library(iNEXT)
library(picante)
library(MuMIn)
library(UpSetR)
library(lmodel2)
library(betapart)
library(fuzzySim)
library(viridis)
library(RColorBrewer)
library(kableExtra)
library(VennDiagram)
library(ComplexHeatmap)
library(estimateR)
library(knitr)
library(ade4)
library(ggfortify)
library(ggplot2)
library(spadeR)

###########################

## combined data file for analysis- depth, lat, long, sciname/taxonConceptID etc select
## merge with checklist to get higher taxonomy
##DD_ed filter out duplicates

## COOPER SCRIPT

# kable(data) 
head(data)
dim(data)

# Create a matrix we can use with vegan
data.matrix <- sample2matrix(data)
# Error in FUN(X[[i]], ...) : invalid 'type' (character) of argument
# replaced NA values with 1- GSR and UKSR

list <- lapply(data, class)
as.data.frame(list)
#Abundance      Site   Species
#1   numeric character character
# cooper data- abundance = integer

head(data.matrix)

#The simplest measure of  Î± diversity is just the number of species in each site
data2 <- specnumber(data.matrix)
write.csv(data2, "test.csv")

#Simpsonâs and Shannonâs diversity indices can be estimated using the function diversity.
diversity(data.matrix, index = "simpson")
diversity(data.matrix, index = "shannon")

# Î² diversity. # Jaccard's index

test <- betadiver(data.matrix, method = "j")
#test2 <- as.data.frame(test)
#write.csv(test, "beta-J_DD_ed_4_diversity_analyses_family.csv")
# can't coerce class 'dist' to a data frame

betadiver(data.matrix, method = "sor")

# Note that the outputs here are pairwise matrices, as these indices measure the similarity among each pair of sites. You can estimate Whittakerâs original version using method = "w" (this is a dissimilarity method so completely different communities get a score of 1).

# Whittaker's betadiversity index
test <-betadiver(data.matrix, method = "w")
write.csv(test, "beta-J_DD_ed_4_diversity_analyses_v3__incl_all_TCIDS_BC_only_2021-09-02.CSV")

# 2.5.3   Î³ diversity. In this example, Î³ diversity is the total number of species found across all sites. We can very simply calculate this in R using the following code:

# How many unique species are there?

length(unique(data$Species[data$Abundance > 0]))
## boxcore only- 46
## 467
## incl morphospecies 976
# To view unique species
unique(data$Species)

#Note that the [data$Abundance > 0] bit of the code ensures we donât count species where we have them 
## in the species list, but their abundance at all sites is zero.

#2.6 Species accumulation curves (Colwell & Coddington 1994)

data.curve <- specaccum(data.matrix, method = "random", permutations = 1000)

# Look at the results
data.curve

# Plot the curve
plot(data.curve, ci.type = "poly", col = "blue", ci.col = "lightblue", 
     lwd = 2, ci.lty = 0, xlab = "number of sites", 
     ylab = "cumulative number of species")

# "ci.type = "poly" tells R that you want a shaded area showing the confidence intervals 
# from your randomisations. 
# ggplot version- Note that because ggplot works with dataframes we first need to create a dataframe 
## from the data.curve object. Also we use standard deviation * 1.96 to get the confidence intervals.

## family level
plot(data.curve, ci.type = "poly", col = "blue", ci.col = "lightblue", 
     lwd = 2, ci.lty = 0, xlab = "number of sites", 
     ylab = "cumulative number of families") +
  
  
  # Load ggplot
  library(ggplot2)
# Make a new dataframe
data.curve.df <- data.frame(sites = data.curve$sites,
                            richness = data.curve$richness,
                            sd = data.curve$sd)
head(data.curve.df)
tail(data.curve.df)

# Plot
ggplot(data.curve.df, aes(x = sites, y = richness)) +
  # Add line
  geom_line() +
  # Add confidence intervals
  geom_ribbon(aes(ymin = richness - 1.96*sd, ymax = richness + 1.96*sd), 
              alpha = 0.5, colour = "lightblue") +
  labs(title="Species accumulation curve- DeepData all named species records") 
# Remove grey backgrouund
#  theme_bw(base_size = 14)

#To demonstrate why we need the randomisations, look at two curves for just one permutation each.

# Fit one curve with just one permutation
data.curve1 <- specaccum(data.matrix, method = "random", permutations = 1)
# Fit another curve with just one permutation
data.curve2 <- specaccum(data.matrix, method = "random", permutations = 1)

# Set the plotting window so we can plot two plots
par(mfrow = c(1,2))

# Plot the first curve
plot(data.curve1,  
     xlab = "number of sites", ylab = "cumulative number of data species")

# Plot the second curve
plot(data.curve2, 
     xlab = "number of sites", ylab = "cumulative number of data species")

# Reset the plotting window so we see just one plot again
par(mfrow = c(1,1))

## not working

# Finally to estimate total species richness across all sites we can (again) use many different metrics. 
# Some common ones include Chao 2 (Chao 1987), Jackknife and Bootstrapping approaches and these are easy 
# to estimate using the vegan function specpool.

# Estimate diversity
data3 <- specpool(data.matrix)
write.csv(data3, "OUTPUT.csv")

#Species    chao  chao.se    jack1 jack1.se   jack2     boot  boot.se   n
# All     477 700.697 40.72924 689.1446 74.55021 800.649 572.7805 56.10595 249

#Estimates range from 700-800

## fam
#Species	chao	chao.se	jack1	jack1.se	jack2	boot	boot.se	n
#All	432	547.8193548	27.916189	551.68	48.71687784	609.5355579	487.026591	26.24960531	375
## checklist- 535 unique fam names

###############################################

## short version

data <- read.csv("temp_all_spp_macro_2022-04-24.csv")
data <- read.csv("temp_all_spp_mega_2022-04-24.csv")
data <- read.csv("temp_all_spp_meio_2022-04-24.csv")
data <- read.csv("temp_all_spp_SITE_2022-04-24.csv")
data <- read.csv("temp_all_fam_SITE_2022-04-24.csv")
data <- read.csv("temp_named_spp_SITE_2022-04-24.csv")

data <- read.csv("temp_DeepData_spp_v_abundance_2022-09-23.csv")
data <- read.csv("temp_lit+DD_spp_v_abundance_2022-09-23.csv")
## not working- file with spp and abundance only
data <- read.csv("TEMP_lit_all_spp_2022-09-23.csv")
data <- read.csv("TEMP_lit+DD_all_spp_site_by_coords_2022-09-23.csv")

## 27th sept
data <- read.csv("TEMP_lit+DD_all_spp_2022-09-27.csv")
data <- read.csv("TEMP_lit+DD_families_2022-09-27.csv")

head(data)
dim(data)

data.matrix <- sample2matrix(data)
data.matrix[1:5, 1:5]

length(unique(data$Species[data$Abundance > 0]))

data3 <- specpool(data.matrix)
write.csv(data3, "OUTPUT.csv")

##################

## do curve

data.curve <- specaccum(data.matrix, method = "random", permutations = 1000)

# Look at the results
data.curve

# Make a new dataframe
data.curve.df <- data.frame(Sites = data.curve$sites,
                            Richness = data.curve$richness,
                            sd = data.curve$sd)
head(data.curve.df)
tail(data.curve.df)

ggplot(data.curve.df, aes(x = Sites, y = Richness)) +
  # Add line
  geom_line() +
  # Add confidence intervals
  geom_ribbon(aes(ymin = Richness - 1.96*sd, ymax = Richness + 1.96*sd), 
              alpha = 0.2, colour = "lightgrey") +
  theme(text=element_text(size=14), #change font size of all text
        axis.text=element_text(size=20), #change font size of axis text
        axis.title=element_text(size=20), #change font size of axis titles
        plot.title=element_text(size=20), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20)) + theme_bw() #
labs(title="Species accumulation curve- DeepData (named species and morphospecies)")


##########################################################################################

##########################################################################################

## rarefaction

library(iNEXT)

ChaoShannon(data.matrix, datatype = "incidence", transform = FALSE,
            conf = 0.95, B = 200)
#Error in if (Q0.hat == 0) { : missing value where TRUE/FALSE needed
#In addition: Warning message:
#  In digamma(t) : NaNs produced

iNEXT(data2, q = 0, datatype = "incidence_freq", size = NULL,
      endpoint = NULL, knots = 40, se = TRUE, conf = 0.95,
      nboot = 50)

#Error in iNEXT.Sam(Spec = x, q = q, t = size, endpoint = ifelse(is.null(endpoint),  : 
#invalid data structure!, first element should be number of sampling units

iNEXT(data2, q = 0, datatype = "incidence_freq", size = NULL,
      endpoint = NULL, knots = 40, se = TRUE, conf = 0.95,
      nboot = 50)
#Error in if (t > sum(y)) { : missing value where TRUE/FALSE needed
#In addition: Warning message:
#In Fun(x, q) : NAs introduced by coercion
###########################

## iNEXT R documentation
data(ciliates)
lapply(ciliates, as.abucount)
head(ciliates)

data(ciliates)
lapply(ciliates, as.incfreq)

ChaoShannon(ciliates, datatype = "incidence", transform = FALSE,
            +             conf = 0.95, B = 200)
#Error in FUN(X[[i]], ...) : invalid data structure
#> iNEXT(ciliates, q = 0, datatype = "incidence_freq", size = NULL,
#+       endpoint = NULL, knots = 40, se = TRUE, conf = 0.95,
#+       nboot = 50)
#Error in iNEXT.Sam(Spec = x, q = q, t = size, endpoint = ifelse(is.null(endpoint),  : 
#invalid data structure!, first element should be number of sampling units

data(bird)
library(ggplot2)
data(bird)
head(bird)
out <- iNEXT(bird, datatype="abundance")
ggiNEXT(out)

out <- iNEXT(data2, datatype="incidence_freq")
ggiNEXT(out)
## no joy

data(spider)
ChaoRichness(spider$Girdled, datatype="abundance")

############################
head(data.matrix)
names(data.matrix)
write.csv(data.matrix, "OUTPUT.csv")
## need to transpose
data2 <- t(data.matrix)
data2[1:5, 1:5]
write.csv(data2, "OUTPUT.csv")

## still not working

data <- read.csv("temp_all_spp_MEIO_no_dupl_2022-06-13.csv")
data <- read.csv("temp_named_spp_REG-CENTRAL-EAST_2022-04-24.csv")
data <- read.csv("temp_named_spp_REG-CENTRAL-EAST_2022-04-24_V2.csv")
data <- read.csv("temp_DD+lit_data_4_spp_curve_SPECIES_2022-08-27.csv")
data <- read.csv("temp_inext_dd+lit_abund_2022-09-23.csv")

head(data)
data.matrix <- sample2matrix(data)
data.matrix[1:5, 1:5]
## transpose
data2 <- t(data.matrix)
data2[1:5, 1:5]
head(data2)
write.csv(data2, "OUTPUT.csv")

out <- iNEXT(data2, datatype="abundance")
ggiNEXT(out)

out <- iNEXT(data2, datatype="incidence_freq")
ggiNEXT(out)

## more luck- try more straightforward data- named only and region

ChaoRichness(data.matrix$Central, datatype="abundance")
ChaoRichness(data2$Central, datatype="abundance") 
#Error in data2$Central : $ operator is invalid for atomic vectors
ChaoRichness(bird$North.site, datatype="abundance")
ChaoRichness(data2$Central, datatype="abundance")
#Error in data2$Central : $ operator is invalid for atomic vectors

#a vector of species abundances or incidence frequencies. If datatype = "incidence",
#then the first entry of the input data must be total number of sampling units, followed
#by species incidence frequencies.

data(spider)
ChaoShannon(spider$Girdled, datatype="abundance")

## incidence data
data <- read.csv("temp_named_spp_REG-CENTRAL-EAST_PA_2022-04-24.csv")

##############

data(bird)
out2 <- iNEXT(bird, q=0, datatype="abundance")
ggiNEXT(out2)

#######################################################################

## rarefaction from vegan

#rarefy
#rarecurve
#estimateR
#specaccum (method = “random”)
#specaccum(method = ”rarefaction”)

## spp accumulation
#xvar Variable used for the horizontal axis: "individuals" can be used only with
#method = "rarefaction".

###########################################

## RAREFACTION

## 27th sept

data <- read.csv("TEMP_lit+DD_all_spp_2022-09-23.csv")

test2 <- rarefy(data.matrix, sample, se= FALSE, MARGIN = 1)
#Error in sample > minsample : 
#  comparison (6) is possible only for atomic and list types

data <- (read.csv("TEMP_lit+DD_all_spp_2022-09-27.csv",header=T,sep=",",fileEncoding="latin1")[-2,])
data <- (read.csv("TEMP_lit+DD_all_spp_site-as-CCZ-only_2022-09-27.csv",header=T,sep=",",fileEncoding="latin1")[-2,])
data <- (read.csv("TEMP_all_spp_DD+lit_site-as-area_2022-10-02.csv",header=T,sep=",",fileEncoding="latin1")[-2,])
data <- (read.csv("TEMP_lit+DD_all_spp_site-as-CCZ-only_2022-09-30-no-site.csv",header=T,sep=",",fileEncoding="latin1")[-2,])

head(data)
dim(data)

data.matrix <- sample2matrix(data)
data.matrix[1:5, 1:5]

## try rare
rarecurve(data.matrix, col = "blue") 
## works but looks awful- each site is a new curve- 
## all sites renamed as 'CCZ- works- but is this valid= sample size x axis- up to 80,000- treating
## each row as a separate site? dimensions- only 23,929

data.presabs <- splist2presabs(data, sites.col = "Site",
                               sp.col = "Species", keep.n = FALSE)
#Error in make.names(vnames, unique = TRUE) : invalid multibyte string 2

head(data.presabs)

test <- t(data.presabs)
head(test)
test <- t(as.matrix(data.presabs[,-1]))
colnames(test) <- data.presabs$Site
head(test)

rarecurve(test, col = "pink") ## didn't work- long straight linke

rarefy(test, sample, se= FALSE, MARGIN = 1)
#Error in sample > minsample : 
#  comparison (6) is possible only for atomic and list types

## orig test
data <- read.csv("DD_no_pelagic_all_spp_4matrix_2021-11-05-SITE.csv")

data.presabs <- splist2presabs(data, sites.col = "Site",
                               sp.col = "Species", keep.n = FALSE)
#Error in make.names(vnames, unique = TRUE) : invalid multibyte string 2

head(data.presabs)

test <- t(data.presabs)
head(test)
test <- t(as.matrix(data.presabs[,-1]))
colnames(test) <- data.presabs$Site
head(test)

rarefy(test, sample, se= FALSE, MARGIN = 1)
#Error in sample > minsample : 
#  comparison (6) is possible only for atomic and list types

RARE_TEST <- rarefy(test, 1, se= FALSE, MARGIN = 1)

plot(RARE_tEST, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(test, step = 20, sample = raremax, col = "blue", cex = 0.6)

########################################
## TRY ABUNDANCE ONLY

data <- (read.csv("TEMP_lit+DD_all_spp_site-as-CCZ-only_2022-09-30-no-site.csv",header=T,sep=",",fileEncoding="latin1")[-2,])

data.presabs <- splist2presabs(data, sites.col = "Abundance",
                               sp.col = "Species", keep.n = FALSE)
#Error in make.names(vnames, unique = TRUE) : invalid multibyte string 2

head(data.presabs)
test <- t(data.presabs)
head(test)
test <- t(as.matrix(data.presabs[,-1]))
colnames(test) <- data.presabs$Site
head(test)

rarecurve(test, col = "pink") ## didn't work- long straight linke


###########################################

data(package = "vegan") ## names of data sets in the package
data(dune) # Vegetation and Environment in Dutch Dune Meadows
str(dune) #a data frame of observations of 30 species at 20 sites
head(dune)
class(dune) ## data frame
diversity(dune,index = "simpson") 
# calculate Simpson's 1-D Index of Diversity for each site. # closer to 1 = greater divers
simpson <- diversity(dune, "simpson") # or assign to var.
simpson
shannon <- diversity(dune) # note that Shannon's is default
shannon #Typically ranges from 1.5 - 3.4, higher = more diverse 
# lets compare the two
par(mfrow = c(1, 2))  # use par to generate panels with 1 row of 2 graphs
hist(simpson)
hist(shannon)
par(mfrow = c(1, 2))
bray = vegdist(dune, "bray") 
gower = vegdist(dune, "gower")
hist(bray, xlim = range(0.0,1.0))
hist(gower, xlim = range(0.0,1.0))
# r allows for multiple iterations of each dissimilarity index to examine #freqeuncy of differences
spAbund <- rowSums(dune)  #gives the number of individuals found in each plot
spAbund # view observations per plot 
raremin <- min(rowSums(dune))  #rarefaction uses the smallest number of observations per sample to extrapolate the expected number if all other samples only had that number of observations
raremin # view smallest # of obs (site 17)
sRare <- rarefy(dune, raremin) # now use function rarefy
sRare #gives an "expected"rarefied" number of species (not obs) if only 15 individuals were present
rarecurve(dune, col = "blue") 
rarefy(dune, 10, se= FALSE, MARGIN = 1)
## number- add in total samples

# produces rarefaction curves # squares are site numbers positioned at observed space. To "rarefy" a larger site, follow the rarefaction curve until the curve corresponds with the lesser site obs. This gives you rarefied species richness

###########################################

## scaling of higher taxon
# library(ggplot2)

## checklist summaries all taxon levels
data <- read.csv("CCZ_CHECKLIST_2022-09-26.CSV")

data <- read.csv("CCZ_MSPP_CHECKLIST_2022-10-06.CSV")


names(data)
data2 <- tapply(data[, "scientificName"], data[, "taxonRank"],
                function(x) { length(unique(x))})

write.csv(data2, "OUTPUT.csv")

## log of total
test <- log(data$Total)
write.csv(test, "OUTPUT.csv")

data <- read.csv("temp_summary_checklist_taxonrank_2022-09-26.csv")
data <- read.csv("temp_summary_checklist_taxonrank_no_spp_2022-10-01.csv")

names(data)

## model
model1 <- lm(order ~ log, data = data)
#Coefficients:
#  (Intercept)          log  
#-2.267        1.029  
model2 <- lm(order ~ Total, data = data)
#Coefficients:
#  (Intercept)        Total  
#1.833771     0.003135
autoplot(model1, smooth.colour = NA)
autoplot(model2, smooth.colour = NA)
anova(model1)
summary(model1)

## PLOT
plot(data$order, data$log, pch = 16)
## add line to basic plot -  line best fit with linear model
abline(lm(data$order, data$log)) ## invalid
abline(lm(data$log, data$order)) ## invalid
abline(lm(order ~ log, data = data)) ## invalid
abline(lm(log ~ order, data = data)) ## worked
fitted <- predict(lm(data$order ~ data$log))
#1        2        3        4        5 
#1.196720 1.756986 2.966092 4.127307 4.952894 

## predicting order values - wrong

## reverse - log- order
model1 <- lm(Total ~ order, data = data)
model1 <- lm(log ~ order, data = data)
#Call:
#  lm(formula = log ~ order, data = data)
#
#Coefficients:
#  (Intercept)        order  
# 2.2380       0.9609  
summary(model1)
#Call:
#  lm(formula = log ~ order, data = data)
#Residuals:
#  1        2        3        4        5 
#0.16845 -0.24768 -0.03297  0.13518 -0.02298 
## gap between value and fitted line

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  2.23799    0.20048   11.16 0.001541 ** 
#  order        0.96086    0.06045   15.90 0.000541 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.1911 on 3 degrees of freedom
#Multiple R-squared:  0.9883,	Adjusted R-squared:  0.9844 
#F-statistic: 252.7 on 1 and 3 DF,  p-value: 0.0005413

summary.aov(model1)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#order        1  9.232   9.232   252.7 0.000541 ***
#  Residuals    3  0.110   0.037                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary.lm(model1)

plot(model1)
autoplot(model1, smooth.colour = NA)

data$order[7];data$log[7]
order[7];log[7]; data = data

## looking at influential point- tannin at 6%, growth at 2
model2 <- update(model1,subset=(order !=6))
summary(model2)

## PREDICT VALUES
plot(data$order, data$log, pch = 16)
abline(lm(log ~ order, data = data)) ## worked
summary(model1)
fitted <- predict(lm(data$log ~ data$order))
#1        2        3        4        5 
#3.198847 4.159705 5.120564 6.081422 7.042281 
## predicting log values

lines(c(0,0),c(3.36,2.23))
for (i in 1:5) lines (c(order[i],order[i],c(log[i],fitted[i])))
for (i in 1:5) lines (c(order[i],order[i],c(Total[i],fitted[i])))
#Error in order[i] : object of type 'closure' is not subsettable
sum(data$order*data$Total)
sum(data$order*data$log)

fitted <- predict(lm(data$order ==6))


data <- read.csv("temp_summary_checklist_taxonrank_test_2022-10-06.csv")
test <- log(data$Total)
write.csv(test, "OUTPUT.csv")
model1 <- lm(log ~ order, data = data)


########################

## plot non- log as exponential
plot(data$order, data$Total, pch = 16)
expon <- lm(log(data$Total) ~ data$order) ## or other way round?
yv2 <- exp(predict(expon,list(x=data$Total)))
yv2 <- exp(predict(expon,list()))## same result

data <- read.csv("temp_summary_checklist_taxonrank_no_spp_2022-10-01.csv")
names(data)
plot(data$order, data$Total, pch = 16)
expon <- lm(log(data$Total) ~ data$order) ## or other way round?
yv2 <- exp(predict(expon,list()))

newX <- expand.grid(data$order == 6)
prediction <- predict(model1, newdata = newX, interval = "confidence")


ggplot(data, aes(x = order, y = log)) + 
  geom_point(col = "cornflowerblue", size = 3) +
  labs(x = "taxon rank", y = "total") +
  theme_bw() +
  geom_smooth(method = "lm", se = TRUE)

ggplot(data, aes(x = order, y = Total)) + 
  geom_point(col = "cornflowerblue", size = 3) +
  labs(x = "taxon rank", y = "total") +
  theme_bw() +
  geom_smooth(method = "lm", se = TRUE)
labs(x = "Soil Moisture", y = "Growth Rate") 
newX <- expand.grid(soil.moisture.content = seq(from = 0.25, to = 2, length = 100))

newY <- predict(model1, newdata = newX, interval = "confidence")

# Look at newY
head(newY)


###########

data(BCI)
sp1 <- specaccum(BCI)
sp2 <- specaccum(BCI, "random")
sp2
summary(sp2)
plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(sp2, col="yellow", add=TRUE, pch="+")
## Fit Lomolino model to the exact accumulation
mod1 <- fitspecaccum(sp1, "lomolino")


#########################################

## 2ND OCTOBER- 

install.packages("estimateR") ## not avail for this version

## RARE CURVE- SITES AS AREA
## redo spp accumulation curve- ggtheme- remove grey background

## redo DD data extraction- abundance= zero- remove, NA abundance- as zero  

## estimates from log- linear regresssion
## spp estimate from rarefaction- ES(2000)


rm(x,y)
par(mfrow=c(1,1))
curve(read.table("c"))

x2<-x^2
quadratic<- lm(y~x + x2)
summary(quadratic)
XV <- seq(0,30,0.1)
yv <- predict(quadratic,list(x=xv,x2=xv^2))
exponential <- lm(log(y)~x)
yv2 <- exp(predict(exponential,list(x=xv)))

#################################################

## 6th OCTOBER 2022

undisturbed <- read.csv(choose.files(),row.names=1,sep=";")
disturbed <- read.csv(choose.files(),row.names=1,sep=";")

rarecurve(undisturbed)
rarecurve(disturbed)

# We can also do rarefaction curves comparing all undisturbed (considering them as one) with all 
#disturbed forest fragments (also as one). So, let us sum all the individuals of each species for 
#both undisturbed...

sumundisturbed<-colSums(undisturbed)

# and disturbed forest fragments.

sumdisturbed<-colSums(disturbed)

# Now, let us combine both objects and do rarefaction curves comparing undisturbed and disturbed 
#forest fragments.

data<-rbind(sumundisturbed, sumdisturbed)
rarecurve(data)

plot(specaccum(undisturbed))
plot(specaccum(disturbed))
estimateR(undisturbed)

specpool(undisturbed)

##################################

## disturbed- matrix

#sp1 sp2 sp3 sp4 sp5 sp6 sp7 sp8 sp9 sp10 sp11 sp12 sp13 sp14 sp15
#sample1   1   0   0   2   0   1   0   0   0    2    0    0    3    3    0



#####################################
## 2 columns- species and abundance- convert to matrix
## sample2matrix- 3 columns, presabs- pres abs data

## TRY pool.aurora <- poolaccum(auroraT)

## this file- species and abundance only
data <- (read.csv("TEMP_lit+DD_all_spp_site-as-CCZ-only_2022-09-30-no-site.csv",header=T,sep=",",fileEncoding="latin1")[-2,])
## species abundance and site
data <- (read.csv("TEMP_lit+DD_all_spp_2022-09-27.csv",header=T,sep=",",fileEncoding="latin1")[-2,])
data <- (read.csv("TEMP_lit+DD_all_spp_site_as_CCZ_2022-10-06.csv",header=T,sep=",",fileEncoding="latin1")[-2,])

####################
## 10TH OCT
data <- (read.csv("TEST_LIT+DD_ALL_SPP_2022-10-10.csv",header=T,sep=",",fileEncoding="latin1")[-2,])
data <- (read.csv("TEST_LIT+DD_ALL_SPP_SITE-AS-CCZ-ONLY_2022-10-10.csv",header=T,sep=",",fileEncoding="latin1")[-2,])

head(data)

data.matrix <- sample2matrix(data)
data.matrix[1:5, 1:5]


rarecurve(data.matrix) ## abundance based- whole CCZ
rarefy(data.matrix, 90611) ## arg 'sample missing
rrarefy(data.matrix, 93456) ## arg 'sample missing
test <- estimateR(data.matrix) ## abund of single sample - FOR CCZ ONLY DATASET
## gives ACE estimator - what about chao2, ES2000
write.csv(test, "OUTPUT.csv")

## sample based
test <- specpool(data.matrix)## sample based
test <- poolaccum(data.matrix) ## try pool accum- this is sample based
write.csv(as.data.frame(pool.aurora, "OUTPUT.csv"))

rarefy(data.matrix, 93456, se= FALSE, MARGIN = 1)
## larger value- 
#CCZ 
#3813 
#attr(,"Subsample")
#[1] 934560
#Warning message:
#  In rarefy(data.matrix, 934560, se = FALSE, MARGIN = 1) :
#  requested 'sample' was larger than smallest site maximum (93456)

## incidence

data.presabs <- splist2presabs(data, sites.col = "Site",
                               sp.col = "Species", keep.n = FALSE)
#Error in make.names(vnames, unique = TRUE) : invalid multibyte string 2
data.presabs <- splist2presabs(data, sites.col = "Abundance",
                               sp.col = "Species", keep.n = FALSE)
head(data.presabs)

test <- t(data.presabs)
head(test)
test <- t(as.matrix(data.presabs[,-1]))
colnames(test) <- data.presabs$Site
head(test)

rarecurve(test)
## straight lin

## 2 COLUMN DATA- TURN INTO MATRIX..
test <- as.matrix(data) ## not working
test <- data.matrix(data) ## turns species into numbers

###################

## 9TH OCTOBER- RERUN ANALYSIS WITH NEW SITE DATA

install.packages("SpadeR")
library(SpadeR)
