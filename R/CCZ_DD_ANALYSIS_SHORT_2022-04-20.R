###########################################################
########### CCZ AND DEEPDATA ANALYSIS
setwd("C:\\Users\\murr\\Desktop\\EVERYTHING\\CURRENT_MANUSCRIPTS\\CCZ_and_DEEPDATA\\Data_analysis_CCZ_DeepData\\")

######################################

## load general data analysis packages
install.packages("obistools")
#library(tidyverse)
library(dplyr)
library(tidyr)
library(robis)
library(curl) 
library(obistools)
library(devtools)
library(lubridate)
library(glmmTMB)
library(ggthemes)
library(DataExplorer)
library(gridExtra)
library(lme4)
library(car)
library(reshape2)
library(ggfortify)
library(boot)
library(aod)
library(mgcv)
library(pscl)
library(ggplot2)
library(RCurl)
library(stringr)
library(forcats)

library(readr)
library(glue)
library(xml2)
library(jsonlite)
library(purrr)
library(knitr)
library(uuid)
library(stringr)
library(worms)

#diversity estimate packages + plotting
library(vegan)
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

install.packages("remotes")
remotes::install_github("jokergoo/ComplexHeatmap")

## load mapping packages + associated
library(GADMTools)
library(sp)
library(spatialEco)
library(spData)
library(maptools)
library(rgdal)
library(rgeos)
library(sf)
library(mapr)
library(leaflet)
library(rgbif)
library(rfishbase)# library(taxizesoap) deprecated?
library(taxize)
library(scrubr)
library(mapr)
library(mregions)

###############################

## OBIS data collections
library(robis)

data.ccz <- occurrence(geometry = "POLYGON ((-160.3477 12.5906, -154.0898 -5.4997, 
                       -107.5430 8.6544, -121.8867 23.3343, -160.3477 12.5906))")

data.newer <-
  data.ccz %>%
  filter(depth >= 3000)

write.csv(data.newer, "OBIS_CCZ_3000m+_at_2021-06-28.csv")


## read in files

#g_all <- read.csv("GBIF_CCZ_all_depths_2021-07-12_v10.csv")
#o_all <- read.csv("OBIS_CCZ_all_depths_2021-07-12_v12.csv")
#DD <- read.csv("DD_PUBLISHED_main_file_2021-07-12_v9.csv")
DD_ed <- read.csv("DD_PUBLISHED_4analysis_ed_2022-05-24.csv")
gbif <- read.csv("GBIF_CCZ_4analysis_2022-09-11.csv")
#obis_all <- read.csv("OBIS_CCZ_3000m_only_4analysis_2021-08-10.csv")
obis <- read.csv("OBIS_only_4_analysis_2022-09-11.csv") ## DD records removed
obis_DD <- read.csv("OBIS_DD_4_analysis_2022-04-20.csv")
lit_papers <- read.csv("CCZ_LITERATURE_PAPERS_2022-08-26.csv") ## Lit minus 'include'= N- unpublished datasets
lit <- read.csv("CCZ_LITERATURE_RECORDS_2022-08-26.csv")
checklist <- read.csv("CCZ_CHECKLIST_2022-08-26.csv")
tcid <- read.csv("CCZ_MSPP_CHECKLIST_2022-08-23.csv")## morphospecies checklist
genbank <- read.csv("CCZ_GENBANK_2022-06-24.csv")

###############################

dim(DD_ed)           # 40518  65 PREV 40521
dim(gbif)            # 2405   49 PREV 2403   51
#dim(obis_all)       # 50695  59
dim(obis_DD)         # 48554  60 prev 48526  60
dim(obis)            # 2185   60 PREV 2169 
dim(lit)             # 8260   69 PREV 4235
dim(lit_papers)      # 163    48
dim(checklist)       # 2525   53 PREV 2354   29
dim(check_short)     # 2371   60
dim(tcid)            # 5277   17  PREV 4803  11
dim(genbank)         # 4719   30

########################################

## SELECT COLUMNS

data1 <- read.csv("DD_PUBLISHED_4analysis_ed_2021-08-18.csv") ## 3676
data2 <- read.csv("DD_published_main_file_2021-07-12_v8.csv") ## 606
merga <- merge(data1, data2, by="rec", all = T)
dim(merga) ## 511   7
dim(data2) ## 268   6
dim(data1) ## 511   2
write.csv(merga, "OUTPUT.csv")

DD <- read.csv("DD_published_main_file_2021-07-12_v10.csv")
names(DD)

DD_ed <- read.csv("DD_published_main_file_2021-07-12_v9.csv") %>% 
  dplyr::select(rec, comp_key, comp_key_temp, occurrenceID, occurrenceID_ed, Match.type, 
                Taxon.status, ScientificName_accepted, Phylum_ed,	Class_ed,	Order_ed,	Family_ed,	Genus_ed,
                Subgenus, Species_ed,	put_sp_ed,	put_cat,	put_valid,	scientificName,	taxonRank,
                TR_w_TCID,	qualifier,	rank_qual,	column_notes,	non_specimen,	incl_in_taxon_analysis,
                incl_strict_taxon_records, INCL,	CL_short_incl, sample_time,	SampleDate_ed,	sample_month,	
                sample_year, SampleCollectionMethod_ed, WaterDepth_ed, StorageLocation_ed, 
                HabitatType_ed, HabitatDescription_ed, ResearchVessel_ed,	Nominal.size.Category_ed, 
                Number.of.individuals_ed, ContractorID,	SubArea,	StationID, StationDesc,	SampleID,
                ActualLatitude,	ActualLongitude, AreaSampled, Phylum, Class, Order, Family, Subfamily,
                Species, Putative.species.name.or.number, Notes.on.taxanomic.identification,
                Taxonomic.Status, Description.of.DNA.Sequence,
                Voucher.code, Voucher.institution.code, Voucher.status,
                Relative.abundance...., Relative.abundance....) %>%
  dplyr::filter(INCL != "NO") %>% 
  mutate(SampleDate_ed = lubridate::dmy(SampleDate_ed), 
         #         WaterDepth = replace(WaterDepth, WaterDepth == "-9", "NA"),
         decimalLongitude = as.numeric(ActualLongitude),
         decimalLatitude = as.numeric(ActualLatitude)
  )

names(DD_ed)
dim(DD_ed) # 40521    65
head(DD_ed)

##################################################

## FINAL DATASET FOR ANALYSIS- remove rows to exclude, select relevant fields

obis.all <- read.csv("OBIS_CCZ_all_depths_2021-07-12_v12.csv")
dim(obis.all) ##    126329    216
names(obis.all)
obis <- obis.all%>% 
  dplyr::select(rec, layer, inclusion_criteria, inc_notes_mar2022, include, taxon_note, scientificName_ed,
                DD_record, date_year, year, ScientificName_accepted, qualif_ed, 
                TR_w_TCID, CCZ_endemic, other_rank, taxonRank_ed,
                Authority, AphiaID_accepted, occurrenceID, occurrenceID_ed,
                scientificName, decimalLatitude, originalScientificName, basisOfRecord, 
                minimumDepthInMeters, maximumDepthInMeters, depth, id, 
                class, order, family, decimalLongitude, decimalLatitude, occurrenceID, kingdom, phylum, rightsHolder,
                institutionID, individualCount, type, catalogNumber, samplingProtocol, eventID,
                taxonRank, eventTime, taxonomicStatus, verbatimDepth, locationID, accessRights, 
                habitat,recordedBy, identifiedBy, taxonRemarks, taxonomicStatus, genus, species, 
                subclass, taxonConceptID, identificationQualifier, previousIdentifications,
                institutionCode, recordNumber, datasetName, associatedSequences)%>% 
  dplyr::filter(include != "no") 
names(obis)
dim(obis) ## 50851    56

obis_DD <- obis %>% 
  filter(DD_record =="DD ISA records")
dim(obis_DD) ## 48554    56
write.csv(obis_DD, "OBIS_DD_4_analysis_2022-04-03.csv")

obis_only <- obis %>% 
  filter(DD_record !="DD ISA records")
dim(obis_only) ## 2297    56
write.csv(obis_only, "OBIS_only_4_analysis_2022-04-21.csv")

###########

gbif.all <- read.csv("GBIF_CCZ_all_depths_2021-07-12_v9.csv")
names(gbif.all)
dim(gbif.all) # 138414    260
gbif <- gbif.all%>% 
  dplyr::select(rec, layer, inclusion_criteria, incl_notes_march22, include, m_taxon_note, acceptedScientificName_ed,
                gbifID, rightsHolder, institutionID,institutionCode, datasetID, datasetName, 
                basisOfRecord, occurrenceID, catalogNumber, qual_ed,
                recordNumber, recordedBy, individualCount, associatedSequences, previousIdentifications,
                eventID, eventDate, eventTime, year, habitat, samplingProtocol, locality, verbatimDepth,
                decimalLongitude, decimalLatitude, typeStatus, identificationQualifier, identifiedBy,
                taxonConceptID, scientificName, kingdom, phylum, class, order, family, genus, species,
                taxonRank, taxonomicStatus, depth, depthAccuracy, species, 
                verbatimScientificName, acceptedScientificName)%>% 
  dplyr::filter(include != "no") 
dim(gbif) # 2416 49

write.csv(gbif, "GBIF_CCZ_4analysis_2022-04-06.csv")

##################################################

## MAPPING -OBIS- GBIF spatial join
## IN QGIS- create shpfile join of all the contract areas, reserved areas and APEIs
## convert to RDS file? not necessary- do readOGR
## could also do join directly in QGIS

zones <- readOGR(dsn = ".", layer = "CCZ_ISA_only_merged_file")
#points <- read.csv("OBIS_CCZ_records_coords_2021-04-17.csv")
points <- read.csv("GBIF_CCZ_records_coords_2021-04-17.csv") ## for GBIF data
coordinates(points) <- ~longitude + latitude ### converted into SPDF, renamed decimalLatitude to latitude
class(points) ## SpatialPointsDataFrame # sp
head(points@coords)
head(points) 
points@bbox 
## gives coords of bounding box
head(zones)
names(zones)
proj4string(points)
## NA 
proj4string(zones) 
## Warning message: In proj4string(zones) : CRS object has comment, which is lost in output
proj4string(points) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
#Warning message:#In showSRID(uprojargs, format = "PROJ", multiline = "NO") :
#  Discarded datum WGS_1984 in CRS definition, but +towgs84= values preserved
points <- spTransform(points, CRSobj = CRS(proj4string(zones)))
#Warning message: In proj4string(zones) : CRS object has comment, which is lost in output
points_with_zones1 <- over(points, zones[,"layer"]) 
## produce a data frame with the same number of rows as your layer of points. 
## connect those two:

points$zonelayer <- points_with_zones1$layer
## replace "name" with the name of column containing zone names. 
## can also use column index
head(points_with_zones1)
write.csv(points_with_zones1, "GBIF_CCZ_SHP_MERGE_2021-04-17_v2.csv")
## worked

points <- read.csv("GBIF_CCZ_SHP_MERGE_2021-04-17_v2.csv")
## dataset with if in or out of TT, composite key, rec, lat long
points <- read.csv("OBIS_CCZ_SHP_MERGE_2021-04-17.csv") ## ie out CCZ
head(zones@data)  # polygons
head(points@data) # points
#Error in h(simpleError(msg, call)) : 
#error in evaluating the argument 'x' in selecting a method for function 'head': trying to get slot "data" from an object (class "data.frame") that is not an S4 object 
plot(zones)
points(points, pch=20)

##############################################

## EXPLORATORY DATA ANALYSIS
## Explore of continuous variables 
## coord_cartesian, xlim, ylim= optimise display
## ggplot(data = survey)+ geom_histogram(mapping = aes(x = Bulinus_tot), binwidth = 0.6) +
## coord_cartesian((ylim = c(0, 50)))

#' Overview over missing values
plot_missing(DD) + theme_few()

## look at the data structure
plot_str(DD) + theme_few()
introduce(DD)

#   rows columns discrete_columns continuous_columns all_missing_columns total_missing_values complete_rows
#1 40521      67               57                  9                   1               343717             0
#total_observations memory_usage
#1            2714907     36818416

## look at the variables
plot_histogram(DD) 

plot_missing(data) + theme_few()

##############################################

## metadata
list <-lapply(DD_ed, class)
as.data.frame(list)
write.csv(list, "DD_list_for_meta_2021-10-02.csv")

meta <- read.csv("DD_list_for_meta_2021-10-02.csv")  
names(meta)
DD_meta <-
  meta %>%
  gather(names, values, X:Total.Biomass.collected..g.m2.)
write.csv(DD_meta, "DD_meta_2021-10-02.csv")

##############################################

## SUMMARIES

## group by- contractorID, depth, researchvessel, sample_year, samplingMethod_ed.
## WaterDepth ; Nominal.size.Category; decimalLatitude; HabitatType" ## waterColumn- look at
## Number.ofIndividuals_ed, habitat_type
## TR_w_TCID; taxonRank; class, order, 

LIT_1 <- lit_recs %>% 
  group_by(Phylum) %>%
  count(category)
write.csv(LIT_1, "OUTPUT.csv")

##############################################

## PLOTS

## violin plots https://www.r-graph-gallery.com/piechart-ggplot2.html

p <- ggplot(data, aes(x=name, y=value, fill=name)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin()

## PIE

ggplot(pie_data, aes(x="", y=value, fill=Phylum)) +
  geom_bar(stat="identity", width= 0.1, color="white") + ## colour adds line between segments
  coord_polar("y", start=0) + theme_void() ## theme- removes numeric labels

## BARPLOT

ggplot(data = DD_art) + 
  geom_bar(mapping = aes(x = Order_ed, fill = Order_ed)) + coord_flip()

## BOXPLOT

ggplot(DD_ed, mapping = aes(x = Nominal.size.Category_ed, y = Number.of.individuals_ed, 
                            fill = Nominal.size.Category_ed)) + geom_boxplot() + coord_flip()

## LINE CHART

names(lit)
ggplot(lit, aes(x=year_published, y=cumul_tot)) +
  geom_line( color="#69b3a2", size=1, alpha=0.9) +
  ggtitle("test publications")

plot(lit$year_published, lit$desc_tot, type = "l") 
## still stepping stone effect

## add mulitple lines
plot(x, y1, type = "l")                                 # Draw first line
lines(x, y2, type = "l", col = "red")                   # Add second line
lines(x, y3, type = "l", col = "green")                 # Add third line

## UPSET

upset(fromExpression(input), 
      nintersects = 40, 
      nsets = 6, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1
)

## lolipop

# Create data
data <- data.frame(x=seq(1,30), y=abs(rnorm(30)))

# Plot
ggplot(data, aes(x=x, y=y)) +
  geom_point() + 
  geom_segment( aes(x=x, xend=x, y=0, yend=y))


####################

## Table 1 publications and descriptions

data <- read.csv("CCZ_MSPP_CHECKLIST_2022-08-23.csv")

## total pubs by year
data <- read.csv("CCZ_LITERATURE_PAPERS_2022-08-17.csv")
names(data)
data2 <- tapply(data[, "reference"], data[, "year_published"],
                function(x) { length(unique(x))})

## total descriptions by year
data <- read.csv("CCZ_LITERATURE_PAPERS_2022-04-06.csv")
data2 <- data %>%
  filter(has_description == "y")
dim(data2)

data2 <- tapply(data2[, "reference"], data2[, "year_published"],
                function(x) { length(unique(x))})

## new spp by year
lit2 <- lit_recs%>%
  filter(category == "new species")
dim(lit2)
data2 <- tapply(lit2[, "scientificName"], lit2[, "year_published"],
                function(x) { length(unique(x))})

## all descriptions by year
lit2 <- lit_recs%>%
  filter(CAT == "CCZ description")
dim(lit2)
data2 <- tapply(lit2[, "scientificName"], lit2[, "year_published"],
                function(x) { length(unique(x))})
write.csv(data2, "OUTPUT.csv")

## add cumulative totals in excel- =SUM($X$2:X10) # X being column in question e.g. total descr

###########################################################

## Publication plots (fig 3- rates of taxonomic work in the CCZ)

data <- read.csv("TEMP_papers_table_1980_on_2022-09-21.csv")

## cumulative totals
p1 <- ggplot(data = data, aes(x = Year)) +
  geom_line(aes(y = cumul_desc, colour = "cumul_desc"), size = 1) +
  geom_line(aes(y = cumul_spp, colour = "cumul_spp"), size = 1) +
  geom_line(aes(y = cumul_pubs, colour = "cumul_pubs"), size = 1) +
  labs(x = "Year", y = "Cumulative totals") +
  #  ggtitle("Cumulative totals of publications and descriptions (year 2000- 2021)") +
  theme(text=element_text(size=12), #change font size of all text
            axis.text=element_text(size=12), #change font size of axis text
            axis.title=element_text(size=12), #change font size of axis titles
            plot.title=element_text(size=12), #change font size of plot title
            legend.text=element_text(size=12), #change font size of legend text
            legend.title=element_text(size=12)) + 
  theme(legend.position = "none") +
  scale_colour_manual("", 
                      breaks = c("cumul_desc", "cumul_spp", "cumul_pubs"),
                      values = c("black", "coral2", "steelblue"),
                      labels = c("all descriptions", "new species", " publications"))

## actual totals
p2 <- ggplot(data = data, aes(x = Year)) +
  geom_line(aes(y = all_desc, colour = "cumul_desc"), size = 0.7) +
  geom_line(aes(y = new_spp, colour = "cumul_species"), size = 0.7) +
  geom_line(aes(y = Publications, colour = "cumul_pubs"), size = 0.7) +
  labs(x = "Year", y = "Total") +
  #  ggtitle("Totals of publications and descriptions (year 2000- 2021)") +
  theme(text=element_text(size=12), #change font size of all text
        axis.text=element_text(size=12), #change font size of axis text
        axis.title=element_text(size=12), #change font size of axis titles
        plot.title=element_text(size=12), #change font size of plot title
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) + 
  scale_colour_manual("", 
                      breaks = c("cumul_desc", "cumul_species", "cumul_pubs"),
                      values = c("black", "red", "blue"),
                      labels = c("all descriptions", "new species", " publications"))

#grid.arrange(p2, p1, ncol = 2)
grid.arrange(p2, p1, nrow = 2)

#########################################

## Fig 4 CCZ descriptions by phyla

data  <- checklist %>%
  filter(description != "")
dim(data)

ggplot(data = data) + 
  geom_bar(mapping = aes(x = description, fill = Phylum)) + 
  labs(x = "Descriptions", y = "Total")+ scale_fill_viridis_d() +
  theme(text=element_text(size=12), #change font size of all text
        axis.text=element_text(size=12), #change font size of axis text
        axis.title=element_text(size=12), #change font size of axis titles
        plot.title=element_text(size=12), #change font size of plot title
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) # +
#  labs(title="CCZ descriptions (2000- present) by Phyla") 

################

## Table 2: New species of the major macrofaunal groups

data <- read.csv("CCZ_CHECKLIST_2022-04-07.csv")

data <- checklist %>%
  filter(SPP_CHECKLIST == "yes")
dim(data) ## 654

data2 <- data %>%
  filter(Order == "Tanaidacea") %>%
  count(Family)
dim(data2)
data2

data2 <- data %>%
  filter(Order == "Tanaidacea") %>%
  filter(description == "new species") %>%
  count(Family)
data2
dim(data2) ## 654

##polychaetes
data2 <- data %>%
  filter(Class == "Polychaeta") %>%
 # filter(description == "new species") %>% #
  group_by(Order) %>%
  count(Family)
data2
dim(data2) ## 654

write.csv(data2, "OUTPUT.csv")

####################################

## Table 3 descriptions by size class and id method

data <- read.csv("CCZ_CHECKLIST_2022-04-07.csv")

data2 <- data %>%
  filter(description == "new species") %>% #
  count(size_cat)

data2 <- data %>%
  filter(description == "new genus") %>% #
  count(size_cat)

data <- read.csv("CCZ_LITERATURE_RECORDS_2022-04-06.csv")
data2 <- data %>%
  filter(CAT == "CCZ description") %>% 
  group_by(descr_id_with) %>%
  count(size_cat)
data2

data <- read.csv("CCZ_LITERATURE_PAPERS_2022-04-06.csv")
data2 <- data %>%
  filter(has_description == "y") %>% 
#  group_by(descr_id_with) %>% #
  count(size_cat)
data2

data2 <- data %>%
  #  group_by(descr_id_with) %>% #
  count(size_cat)
data2
write.csv(data2, "OUTPUT.csv")

#######################################

## FIG 5

data <- read.csv("CCZ_CHECKLIST_2022-04-07.csv")

data2 <- checklist%>%
  filter(description != "")
dim(data2)

ggplot(data = data2) + 
  geom_bar(mapping = aes(x = size_cat, fill = description)) +
  scale_fill_viridis_d() + 
  labs(x = "Size class", y = "Total") +
  theme(text=element_text(size=12), #change font size of all text
        axis.text=element_text(size=12), #change font size of axis text
        axis.title=element_text(size=12), #change font size of axis titles
        plot.title=element_text(size=12), #change font size of plot title
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) # +
  #    labs(title="Checklist: descriptions by size class and Phlyum")

########################################

## Fig 6: records by year of sampling

## SUMMARY LEVEL DATA

data <- lit_recs %>%
  count(year_published)

data <- DD_ed %>%
  count(sample_year)

data <- gbif %>%
  count(year)

data <- obis %>%
  count(date_year)

write.csv(data, "OUTPUT.csv")

data <- read.csv("temp_all_recs_by_year_v4_2000_on_2022-04-10.csv")
#data <- read.csv("temp_all_recs_by_year_v2_2022-04-10.csv") ## all years

head(data)
ggplot(data, aes(year, n, fill = Source)) +
  geom_bar(stat="identity", position = "dodge") + 
  labs(x = "Year", y = "Total") +
  scale_fill_viridis_d() + 
  theme(text=element_text(size=12), #change font size of all text
        axis.text=element_text(size=12), #change font size of axis text
        axis.title=element_text(size=12), #change font size of axis titles
        plot.title=element_text(size=12), #change font size of plot title
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) # +
#  labs(title="Records by year of sampling: DeepData and the literature")
  
data2 <- tapply(data[, "reference"], data[, "size_Cat"],
                function(x) { length(unique(x))})

#####################################

## Fig 9: taxonomic resolution of records BY DATA SOURCE

data <- lit_recs %>%
  filter(other_rank != "y")

data <- DD_ed%>%
  filter(other_rank != "y")

data <- obis%>%
  filter(other_rank != "y")

data <- gbif%>%
  filter(other_rank != "y")

p <- ggplot(data = data) + 
  geom_bar(mapping = aes(x = TR_w_TCID, fill =TR_w_TCID)) +
  labs(x = "Taxonomic rank", y = "Records")+ 
  #  labs(title="Taxonomic resolution of records: CCZ literature") +
  scale_fill_viridis_d() +
  theme(text=element_text(size=12), #change font size of all text
        axis.text=element_text(size=12), #change font size of axis text
        axis.title=element_text(size=12), #change font size of axis titles
        plot.title=element_text(size=12), #change font size of plot title
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12))
p_ed <- p + theme(legend.position = "none")

#####################################

## Fig 10: pie chart- phyla by data source

data <- lit_recs %>% 
  count(Phylum)
write.csv(data, "OUTPUT.csv")

data <- DD_ed %>% 
  count(Phylum_ed)

data <- read.csv("temp_fig10_lit_phylum_pie_data_2022-04-12.csv")
data <- read.csv("temp_fig10_DD_phylum_pie_data_2022-04-12.csv")
data <- read.csv("temp_fig10_obis_phylum_pie_data_2022-04-12.csv")
data <- read.csv("temp_fig10_gbif_phylum_pie_data_2022-04-12.csv")

ggplot(data, aes(x="", y=value, fill=Phylum)) +
  scale_fill_viridis_d() +
  geom_bar(stat="identity", width= 0.1, color="white") + ## colour adds line between segments
  coord_polar("y", start=0) + theme_void() ## theme- removes numeric labels

## Sipuncula- edit taxon rank to order- now updated in WoRMS
## DD - ScientificName_accepted ## phylum_ed
## gbif acceptedScientificName_ed (no 'Sipuncula' records in sci name)
## obis scientificName_ed (no 'Sipuncula' records in sci name)

## Table 7

data <- DD_ed %>% 
  group_by(ContractorID) %>%
  count(taxonRank)

data <- DD_ed %>% 
  group_by(ContractorID) %>%
  count(TR_w_TCID)
write.csv(data, "OUTPUT.csv")

## table 8 records by contractor and taxon level 

DD_ed <- read.csv("DD_PUBLISHED_4analysis_ed_2022-04-23.csv")

data2 <- tapply(data[, "ScientificName_accepted"], data[, "size_cat"],
                function(x) { length(unique(x))})

DD2 <- DD_ed %>% 
  filter(pelagic != "pelagic")
## also filter pelagic for those totals

data <- DD_ed %>% 
  filter(taxonRank == "phylum")
dim(data)

data <- DD_ed %>% 
  filter(taxonRank == "class")

data <- DD_ed %>% 
  filter(taxonRank == "order")

data <- DD_ed %>% 
  filter(taxonRank == "family")

data <- DD2 %>% 
  filter(taxonRank == "genus")

data <- DD2 %>% 
  filter(taxonRank == "species")

data2 <- tapply(data[, "ScientificName_accepted"], data[, "ContractorID"],
                function(x) { length(unique(x))})

length(unique(data$ScientificName_accepted))

length(unique(data$put_sp_ed))

## do these tables by contract area - with and without pelagics
## new section in report- after 'trends by contractor'

## morphospp
data <- DD_ed %>% 
  filter(TR_w_TCID == "temporary name")

data2 <- tapply(data[, "ScientificName_accepted"], data[, "contract_area"],
                function(x) { length(unique(x))})
## this counting sci names tho- want to count put_sp_ed

DD_ed <- read.csv("DD_PUBLISHED_4analysis_ed_2022-04-23.csv")

data <- DD_ed %>% 
  filter(pelagic != "pelagic") %>% 
  filter (TR_w_TCID == "temporary name") %>%
  filter(put_valid != "no")

dim(data)
length(unique(data$put_sp_ed))

data2 <- tapply(data[, "put_sp_ed"], data[, "contract_area"],
                function(x) { length(unique(x))})
data2 <- tapply(data[, "put_sp_ed"], data[, "ContractorID"],
               function(x) { length(unique(x))})
write.csv(data2, "OUTPUT.csv")

##########################

## table 12

data <- read.csv("CCZ_MSPP_CHECKLIST_2022-04-24.csv")

length(unique(data$taxonConceptID)) ## 4436

dim(data)
data2 <- tapply(data[, "taxonConceptID"], data[, "source"],
                function(x) { length(unique(x))})

data2 <- tapply(data[, "taxonConceptID"], data[, "Phylum"],
                function(x) { length(unique(x))})

data2 <- tapply(data[, "taxonConceptID"], data[, "size_cat"],
                function(x) { length(unique(x))})
write.csv(data2, "OUTPUT.csv")

data2 <- data %>% 
  filter(source == "DeepData")

data2 <- data %>% 
  filter(source == "literature")

data2 <- data %>% 
  filter(source == "OBIS")

data2 <- data %>% 
  filter(source == "GBIF")

write.csv(data3, "OUTPUT.csv")

## Table 11 morphospp by size cat and phylum 

data <- read.csv("CCZ_MSPP_CHECKLIST_2022-04-24.csv")
dim(data)

data2 <- data %>% 
  filter(size_cat == "macrofauna")
dim(data2)

data2 <- data %>% 
  filter(size_cat == "megafauna")

data2 <- data %>% 
  filter(size_cat == "meiofauna")

data3 <- tapply(data2[, "taxonConceptID"], data2[, "Phylum"],
                function(x) { length(unique(x))})
write.csv(data3, "OUTPUT.csv")

## look at undescribed spp ## tot by class - order etc

data <- read.csv("CCZ_MSPP_CHECKLIST_2022-04-24.csv")
dim(data)
data2 <- data %>% 
  filter(new_spp_gen != "")
length(unique(data2$taxonConceptID)) ## 1221

data3 <- tapply(data2[, "taxonConceptID"], data2[, "new_spp_gen"],
                function(x) { length(unique(x))})

data3 <- tapply(data2[, "taxonConceptID"], data2[, "source"],
                function(x) { length(unique(x))})

data3 <- tapply(data2[, "taxonConceptID"], data2[, "Phylum"],
                function(x) { length(unique(x))})

data3 <- tapply(data2[, "taxonConceptID"], data2[, "Class"],
                function(x) { length(unique(x))})

data3 <- tapply(data2[, "taxonConceptID"], data2[, "Order"],
                function(x) { length(unique(x))})

data3 <- tapply(data2[, "taxonConceptID"], data2[, "size_cat"],
                function(x) { length(unique(x))})


write.csv(data3, "OUTPUT.csv")

## class - nest within phylum
data3 <- data2 %>% 
  group_by(Phylum, Class, Order) %>%
  count(new_spp_gen)
write.csv(sum, "OUTPUT.csv")


#################################

## Fig 12 morphospp names

length(unique(data$taxonConceptID))

ggplot(data = data) + 
  geom_bar(mapping = aes(x = Phylum, fill = source)) +
  labs(x = "taxonomic rank", y = "Records")+ scale_fill_viridis_d() +
  #  labs(title="New species and all recorded species names by Phylum") +
  theme(text=element_text(size=12), #change font size of all text
        axis.text=element_text(size=12), #change font size of axis text
        axis.title=element_text(size=12), #change font size of axis titles
        plot.title=element_text(size=12), #change font size of plot title
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 55, hjust = 1))

length(unique(data$scientificName))


## Fig 11 records by size class 

data <- read.csv("CCZ_LITERATURE_RECORDS_2022-04-10.csv")

ggplot(data) + 
  geom_bar(mapping = aes(x = size_cat, fill = Record)) +
  labs(x = "", y = "")+ 
#  labs(title="Literature: records by size class") + 
  scale_fill_viridis_d() +
  theme(text=element_text(size=14), #change font size of all text
        axis.text=element_text(size=14), #change font size of axis text
        axis.title=element_text(size=14), #change font size of axis titles
        plot.title=element_text(size=14), #change font size of plot title
        legend.text=element_text(size=14), #change font size of legend text
        legend.title=element_text(size=14)) 

data <-read.csv("DD_PUBLISHED_4analysis_ed_2022-04-14.csv")

ggplot(data) + 
  geom_bar(mapping = aes(x = size_cat, fill = Record)) +
  labs(x = "", y = "")+ 
  #  labs(title="DeepData: records by size class") + 
  scale_fill_viridis_d() +
  theme(text=element_text(size=14), #change font size of all text
        axis.text=element_text(size=14), #change font size of axis text
        axis.title=element_text(size=14), #change font size of axis titles
        plot.title=element_text(size=14), #change font size of plot title
        legend.text=element_text(size=14), #change font size of legend text
        legend.title=element_text(size=14)) 

###############################################

## Table 10

data <- read.csv("CCZ_LITERATURE_RECORDS_2022-04-10.csv")

length(unique(data$scientificName))

## record count by id method
data2 <- tapply(data[, "rec"], data[, "descr_id_with"],
                function(x) { length(unique(x))})

## total by size class
data2 <- tapply(data[, "rec"], data[, "size_cat"],
                function(x) { length(unique(x))})

## sep totals by size class
data2 <- data %>% 
  filter(size_cat == "macrofauna") ## do for mega and meio too
data3 <- tapply(data2[, "rec"], data2[, "descr_id_with"],
                function(x) { length(unique(x))})

## do totals by record type ## already filtered by size class
data <- read.csv("CCZ_LITERATURE_RECORDS_2022-04-10.csv")

data2 <- data %>% 
  filter(size_cat == "macrofauna")

data3 <- data2 %>% 
  filter(Record == "CCZ description")

data4 <- tapply(data3[, "rec"], data3[, "descr_id_with"],
                function(x) { length(unique(x))})
length(unique(data3$scientificName)) 
## 23 mega
## 41 meio
## 150

## or in tidyverse
sum <- data %>% 
  group_by(size_cat, Record) %>%
  count(descr_id_with)
write.csv(sum, "OUTPUT.csv")

## name count by size cat + id method?

data <- read.csv("CCZ_LITERATURE_RECORDS_2022-04-10.csv")

data2 <- data %>% 
  filter(size_cat == "macrofauna")

data2 <- data %>% 
  filter(size_cat == "megafauna")

data2 <- data %>% 
  filter(size_cat == "meiofauna")

## name count by id method
data3 <- tapply(data2[, "scientificName"], data2[, "descr_id_with"],
                function(x) { length(unique(x))})
write.csv(data3, "OUTPUT.csv")

data3 <- data2 %>% 
  filter(Record == "CCZ description")

## name count by size cat
data2 <- tapply(data[, "scientificName"], data[, "size_cat"],
                function(x) { length(unique(x))})


#####################################

## Fig 12 data records by size class and phylum
## filter out records where phylum not recorded

data <- lit_recs %>% 
  filter(Phylum != "")
ggplot(data) + 
  geom_bar(mapping = aes(x = size_cat, fill = Phylum)) +
  labs(x = "", y = "")+ 
  #  labs(title="Literature: records by size class") + 
  scale_fill_viridis_d() +
  theme(text=element_text(size=14), #change font size of all text
        axis.text=element_text(size=14), #change font size of axis text
        axis.title=element_text(size=14), #change font size of axis titles
        plot.title=element_text(size=14), #change font size of plot title
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) 

DD_ed <-read.csv("DD_PUBLISHED_4analysis_ed_2022-04-11.csv")

data <- DD_ed %>% 
  filter(Phylum. != "")
ggplot(data) + 
  geom_bar(mapping = aes(x = size_cat, fill = Phylum.)) +
  labs(x = "", y = "")+ 
  #  labs(title="DeepData: records by size class") + 
  scale_fill_viridis_d() +
  theme(text=element_text(size=14), #change font size of all text
        axis.text=element_text(size=14), #change font size of axis text
        axis.title=element_text(size=14), #change font size of axis titles
        plot.title=element_text(size=14), #change font size of plot title
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) 


## fig 13 and 14- do by name counts

## or in tidyverse
data <- DD_ed %>% 
  group_by(size_cat, Record) %>%
  count(ScientificName_accepted)
write.csv(sum, "OUTPUT.csv")

#############################

## table 12 spp by phylum

check <- read.csv("CCZ_CHECKLIST_2022-08-26.csv")
data <- check %>% 
  filter(SPP_CHECKLIST == "yes")

data2 <- tapply(data[, "scientificName"], data[, "Phylum"],
                function(x) { length(unique(x))})
length(unique(data$scientificName)) 
write.csv(data2, "OUTPUT.csv")

## FIG 15 species names in checklist by phylum

ggplot(data) + 
  geom_bar(mapping = aes(x = Phylum, fill = Phylum)) +
  labs(x = "", y = "Total")+ 
  #  labs(title="Literature: records by size class") + 
  scale_fill_viridis_d() +
  theme(text=element_text(size=12), #change font size of all text
        axis.text=element_text(size=12), #change font size of axis text
        axis.title=element_text(size=12), #change font size of axis titles
        plot.title=element_text(size=12), #change font size of plot title
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) +
        theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 55, hjust = 1))

###########################

## table 12 ## total species names by phylum and data source

data2 <- data %>% 
  filter(source_lit == "literature")

data2 <- data %>% 
  filter(source_gbif == "GBIF")

data2 <- data %>% 
  filter(source_DD == "DD")

data2 <- data %>% 
  filter(source_OBIS_only == "OBIS only")

data3 <- tapply(data2[, "scientificName"], data2[, "Phylum"],
                function(x) { length(unique(x))})

write.csv(data3, "OUTPUT.csv")

## FIG 16 tot species by data soruce- use summary table from above

data <- read.csv("Temp_checklist_summary_table_by_phylum+source_2022-04-13.csv")
names(data)

ggplot(data, aes(Phylum, Total, fill = Source)) +
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_viridis_d() +
  theme(axis.title.x=element_blank(),
        text=element_text(size=12),
        legend.text=element_text(size=12),
        axis.text=element_text(size=12),
        axis.text.x = element_text(angle = 55, hjust = 1))

#######################
### table 13 taxon res by phylum

data <- read.csv("CCZ_CHECKLIST_2022-04-18.csv")

data2 <- tapply(data[, "taxonRank"], data[, "Phylum"],
                function(x) { length(unique(x))})

data2 <- data %>% 
  group_by(taxonRank) %>%
  count(Phylum)
write.csv(data2, "OUTPUT.csv")

data2 <- data %>% 
  count(taxonRank)

#######################
### fig 17 total names by rank and phylum

## filter out other ranks- subgenus etc
data2 <- data %>% 
  filter(other_rank != "y") %>% 
  filter(Phylum != "")

ggplot(data = data2) + 
  geom_bar(mapping = aes(x = Phylum, fill = taxonRank)) +
  labs(x = "", y = "Total") + scale_fill_viridis_d() +
  #  labs(title="Checklist: names by Phylum") + coord_flip()
  theme(axis.title.x=element_blank(),
        text=element_text(size=12),
        legend.text=element_text(size=12),
        axis.text=element_text(size=12),
        axis.text.x = element_text(angle = 55, hjust = 1))

####################################

data <- read.csv("CCZ_CHECKLIST_2022-04-24.csv")
data2 <- data %>% 
  filter(SPP_CHECKLIST == "yes")
dim(data2)

data3 <- tapply(data2[, "scientificName"], data2[, "size_cat"],
                function(x) { length(unique(x))})
length(unique(data2$scientificName)) 

data3 <- data2 %>% 
  group_by(size_cat) %>%
  count(Phylum)
write.csv(data3, "OUTPUT.csv")

data3 <- tapply(data2[, "scientificName"], data2[, "Phylum"],
                function(x) { length(unique(x))})

#########################################

## fig 18 total species names in checklist by size class and phylum

ggplot(data2) + 
  geom_bar(mapping = aes(x = size_cat, fill = Phylum)) +
  labs(x = "", y = "Total") +
 # labs(title="Checklist: species names by size class") 
scale_fill_viridis_d() +
  theme(text=element_text(size=12), #change font size of all text
        axis.text=element_text(size=12), #change font size of axis text
        axis.title=element_text(size=12), #change font size of axis titles
        plot.title=element_text(size=12), #change font size of plot title
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) 

#######################################

## new fig: total names in DD - lit by size cat - w record type / phylum
data <- read.csv("CCZ_CHECKLIST_2022-04-18.csv")

## lit_recs
## DD_ed
## HERE

data <- lit_recs <- read.csv("CCZ_LITERATURE_RECORDS_2022-04-10.csv")
data <- read.csv("DD_PUBLISHED_4analysis_ed_2022-04-14.csv")

data2 <- data %>% 
  group_by(size_cat) %>%
  count(Phylum)
write.csv(data2, "OUTPUT.csv")

data <- read.csv("OUTPUT.csv")
data2 <- data %>% 
  filter(Phylum != "")
dim(data2)

ggplot(data2) + 
  geom_bar(mapping = aes(x = size_cat, fill = Phylum)) +
  labs(x = "", y = "Total") +
  # labs(title="Checklist: species names by size class") 
  scale_fill_viridis_d() +
  theme(text=element_text(size=12), #change font size of all text
        axis.text=element_text(size=12), #change font size of axis text
        axis.title=element_text(size=12), #change font size of axis titles
        plot.title=element_text(size=12), #change font size of plot title
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) 


data2 <- data %>% 
  group_by(size_cat) %>%
  count(Record)

data2<- tapply(data[, "scientificName"], data[, "size_cat"],
                function(x) { length(unique(x))})

data2<- tapply(data[, "scientificName"], data[, "Record"],
               function(x) { length(unique(x))})

data2<- tapply(data[, "scientificName"], data[, "Phylum"],
               function(x) { length(unique(x))})


#######################################

## fig 19 plot no spp by phylum versus CCZ descriptions/newly described
## summary table- total species, and ccz

data <- read.csv("Temp_CCZ_spp_v_named_2022-04-13.csv")
ggplot(data, aes(Phylum, Total, fill = Record)) +
  geom_bar(stat="identity", position = "dodge") + 
  #  labs(title="New species and all recorded species names by Phylum") +
  # scale_fill_viridis_d() +
  theme(axis.title.x=element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        axis.title=element_text(size=12), 
        axis.text=element_text(size=12),
        axis.text.x = element_text(angle = 55, hjust = 1))

## do summary for genera
data <- read.csv("CCZ_CHECKLIST_2022-04-07.csv")
data2 <- data %>% 
  filter(taxonRank == "genus")
dim(data2)

data3 <- data2 %>% 
  count(Phylum)
write.csv(data3, "OUTPUT.csv")

data3 <- data2 %>% 
  filter(description == "new genus")

data4 <- data3 %>% 
  count(Phylum)
write.csv(data4, "OUTPUT.csv")

## S fig 4

data2 <- read.csv("Temp_summary_table_CCZ_gen_v_named_2022-04-13.csv")
ggplot(data2, aes(Phylum, Total, fill = Record), show.legend = FALSE) +
  geom_bar(stat="identity", position = "dodge") + 
  #  labs(title="New genera and all recorded genera by Phylum") +
  theme(axis.title.x=element_blank(),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        axis.title=element_text(size=14), 
        axis.text=element_text(size=14),
        axis.text.x = element_text(angle = 55, hjust = 1))

#########################

data <- tcid %>% 
  count(size_cat)

## table 14

## summary morphospp by phylum
## prop morpho pf all spp- sum all sp- div morpho by tot

data <- read.csv("CCZ_MSPP_CHECKLIST_2022-04-07.csv")

data2 <- tapply(data[, "taxonConceptID"], data[, "Phylum"],
                function(x) { length(unique(x))})
write.csv(data2, "OUTPUT.csv")

###############################################

## diversity estimates

## SPP ACCUM CURVES (figs 20 and 21)

## DD data select: Species4analysis (pelagics stripped, spp level data only- tcid + named combined)
## ScientificName_accepted put_sp_ed	put_valid
##	pelagic	size_cat	TR_w_TCID	Number.of.individuals_test	site
## Species_4analysis now in 1 column- pelagics removed- temp names + named spp
## remove sub and superfam from family dataset?

data <- read.csv("temp_all_spp_macro_2022-04-24.csv")
data <- read.csv("temp_all_spp_mega_2022-04-24.csv")
data <- read.csv("temp_all_spp_meio_2022-04-24.csv")
data <- read.csv("temp_all_spp_SITE_2022-04-24.csv")
data <- read.csv("temp_all_fam_SITE_2022-04-24.csv")
data <- read.csv("temp_named_spp_SITE_2022-04-24.csv")


head(data)
dim(data)

data.matrix <- sample2matrix(data)
data.matrix[1:5, 1:5]

## transpose
#data2 <- t(data.matrix)
#data2[1:5, 1:5]
#write.csv(data2, "OUTPUT.csv")
## write to file

##############

alpha <- specnumber(data.matrix)
write.csv(alpha, "OUTPUT.csv")

#Simpson’s and Shannon’s diversity indices 
data2 <- diversity(data.matrix, index = "simpson")
data2 <- diversity(data.matrix, index = "shannon")

#test <- as.data.frame(data2)
#write.csv(test, "OUTPUT.csv")

# β diversity. # Jaccard's index

data2 <- betadiver(data.matrix, method = "j")
data2[1:5]

data2 <- betadiver(data.matrix, method = "sor")

# Note that the outputs here are pairwise matrices, as these indices measure the similarity 
#among each pair of sites. You can estimate Whittaker’s original version using method = "w" 

#Whittakers dissimilarity method (completely different communities=1).

# Whittaker's betadiversity index
data2 <-betadiver(data.matrix, method = "w")

## convert matrices to table

test <- data.frame(t(combn(rownames(data.matrix),2)), as.numeric(data2))
names(test) <- c("c1", "c2", "distance")
#test2 <- melt(test)[melt(upper.tri(test))$value,]
#names(test2) <- c("c1", "c2", "distance")

write.csv(test, "OUTPUT.csv")

##########################

## total number of species found across all sites.- 
## How many unique species by site

length(unique(data$Species[data$Abundance > 0]))

unique(data$Species)

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
  theme(text=element_text(size=12), #change font size of all text
        axis.text=element_text(size=14), #change font size of axis text
        axis.title=element_text(size=14), #change font size of axis titles
        plot.title=element_text(size=14), #change font size of plot title
        legend.text=element_text(size=14), #change font size of legend text
        legend.title=element_text(size=14)) #+
labs(title="Species accumulation curve- DeepData (named species and morphospecies)")


## Species richness estimates

data3 <- specpool(data.matrix)
write.csv(data3, "OUTPUT.csv")

############################

### UPSET PLOTS (figs 25- 28)

## all species
## region plots- region name = 'site', same for subarea etc
upset <- read.csv("temp_all_spp_REG_2022-04-24.csv")
upset <- read.csv("temp_all_spp_CA_2022-04-24.csv") 
## remove APEI 9- too few records- 
upset <- read.csv("temp_all_spp_CA_NO_APEI-9_2022-04-24.csv")
upset <- read.csv("temp_all_spp_CA_NO_APEIS_2022-04-24.csv")
upset <- read.csv("temp_all_spp_SA_2022-04-24.csv")
upset <- read.csv("temp_named_spp_CA_2022-04-24.csv")

upset <- read.csv("TEMP_LIT_NAMED_SPP_UPSET_TEST_2022-05-01.csv")
upset <- read.csv("TEMP_LIT_ALL_SPP_UPSET_TEST_2022-05-01.csv")

## w presence absense from fuzzy sim
data.presabs <- splist2presabs(upset, sites.col = "Site",
                               sp.col = "Species", keep.n = FALSE)

test <- t(data.presabs)
head(test)
test <- t(as.matrix(data.presabs[,-1]))
colnames(test) <- data.presabs$Site
head(test)

write.csv(test, "OUTPUT.csv")

#####################################

## UPSET- other script

#upset <- read.csv("temp_DD_ed_4_diversity_all_spp_region_2022-04-14.csv")
#data.matrix <- as.matrix(upset)
#test <- make_comb_mat(upset)

#upset <- read.csv("DD_pres-abs_table_v1_2021-10-26.csv")

#data.matrix <- as.matrix(upset)
#head(data.matrix)
#write.csv(data.matrix, "OUTPUT.csv")

## transpose- species in rows, contractors in columns
## t - transpose

#test <- t(data.matrix)
#test <- t(data.presabs)
#head(test)

#test <- t(as.matrix(data.presabs[,-1]))
#colnames(test) <- data.presabs$Site
#head(test)
#write.csv(test, "OUTPUT.csv")

### PLOTTING

basin_COLS <- viridis(n=8)[-8]
names(basin_COLS) <- colnames(test)
morph_comb <- make_comb_mat(test, mode = "intersect")

# pdf("Shared_wCCZ_species-Upset_plot.pdf", width = 10, height = 3.5)
morph_plot <- UpSet(morph_comb, set_order = colnames(morph_comb), 
                    comb_order = rev(order(comb_size(morph_comb))), 
                    pt_size = unit(2, "mm"), lwd = 1, top_annotation = HeatmapAnnotation(
  
  "Shared species (all)" = anno_barplot(comb_size(morph_comb),
                                        ylim = c(0, max(comb_size(morph_comb))*1.1),
                                        border = FALSE,
                                        gp = gpar(fill = "black"),
                                        height = unit(7, "cm")  ## adjust bar to dot matrix ratio: 10 REG   
  ),
  annotation_name_side = "left",
  annotation_name_rot = 90), ## angle of legend
  left_annotation = rowAnnotation(
    "No. species per region (all)" = anno_barplot(-set_size(morph_comb),
                                            baseline = 0,
                                            axis_param = list(
                                              at = c(0, -250, -500, -1000, -1500),
                                              labels = c(0, 250, 500, 1000, 1500),
                                              labels_rot = 0),
                                            border = FALSE,
                                            gp = gpar(fill = basin_COLS),
                                            width = unit(7, "cm"),  ## width of lhs of y axis: 7: CA; 5: SA; 10 REG
                                            height = unit(0.5, "cm")),
    set_name = anno_text(set_name(morph_comb),
                         location = 0.5, ## labels for side bars
                         just = "center",
                         width = max_text_width(set_name(morph_comb)) + 
                           unit(1, "mm"))), right_annotation = NULL, show_row_names=FALSE)

morph_plot = draw(morph_plot)

######################################

## PLOT - NAMED ONLY 

## do for named species only
upset <- read.csv("temp_named_spp_REG_2022-04-24.csv")
upset <- read.csv("temp_named_spp_CA_2022-04-24.csv") ## removed APEI9 - 1 record only
upset <- read.csv("temp_named_spp_CA_NO_APEIS_2022-04-24.csv")
upset <- read.csv("temp_named_spp_SA_2022-04-24.csv")

data.presabs <- splist2presabs(upset, sites.col = "Site",
                               sp.col = "Species", keep.n = FALSE)

test <- t(data.presabs)
head(test)
test <- t(as.matrix(data.presabs[,-1]))
colnames(test) <- data.presabs$Site
head(test)


basin_COLS <- viridis(n=8)[-8]
names(basin_COLS) <- colnames(test)
morph_comb <- make_comb_mat(test, mode = "intersect")

# pdf("Shared_wCCZ_species-Upset_plot.pdf", width = 10, height = 3.5)
morph_plot <- UpSet(morph_comb, set_order = colnames(morph_comb), 
                    comb_order = rev(order(comb_size(morph_comb))), 
                    pt_size = unit(2, "mm"), lwd = 1, top_annotation = HeatmapAnnotation(
                      
                      "Shared named species" = anno_barplot(comb_size(morph_comb),
                                                            ylim = c(0, max(comb_size(morph_comb))*1.1),
                                                            border = FALSE,
                                                            gp = gpar(fill = "black"),
                                                            height = unit(9, "cm")  ## adjust bar to dot matrix ratio    
                      ),
                      annotation_name_side = "left",
                      annotation_name_rot = 90), ## angle of legend
                    left_annotation = rowAnnotation(
                      "No. named species per region" = anno_barplot(-set_size(morph_comb),
                                                                    baseline = 0,
                                                                    axis_param = list(
                                                                      at = c(0, -20, -40, -60, -80, -100),
                                                                      labels = c(0, 20, 40, 60, 80, 100),
                                                                      labels_rot = 0),
                                                                    border = FALSE,
                                                                    gp = gpar(fill = basin_COLS),
                                                                    width = unit(7, "cm"),  ## width of lhs of y axis: 7: CA; 5: SA; 10 REG
                                                                    height = unit(0.5, "cm")),
                      set_name = anno_text(set_name(morph_comb),
                                           location = 0.5, ## labels for side bars
                                           just = "center",
                                           width = max_text_width(set_name(morph_comb)) + 
                                             unit(1, "mm"))), right_annotation = NULL, show_row_names=FALSE)

morph_plot = draw(morph_plot)

######################

## adjusted

baseline = 0,
axis_param = list(
  at = c(0, -5, -10, -15),
  labels = c(0, 5, 10, 15),
  labels_rot = 0),
border = FALSE,
gp = gpar(fill = basin_COLS),
width = unit(10, "cm"),  ## width of lhs of y axis: 7: CA; 5: SA; 10 REG
height = unit(0.5, "cm")),

  at = c(0, -20, -40, -60, -80, -100),
  labels = c(0, 20, 40, 60, 80, 100),
  

## rarefaction

data(package = "vegan") ## names of data sets in the package
data(dune) # Vegetation and Environment in Dutch Dune Meadows
str(dune) 

spAbund <- rowSums(dune)
raremin <- min(rowSums(dune))  #rarefaction uses the smallest number of observations per sample to e
## xtrapolate the expected number if all other samples only had that number of observations
raremin # view smallest # of obs (site 17)

sRare <- rarefy(dune, raremin) # now use function rarefy
sRare
rarecurve(dune, col = "blue") # produces rarefaction curves # squares are site numbers positioned at 
## observed space. To "rarefy" a larger site, follow the rarefaction curve until the curve corresponds with the lesser site obs. This gives you rarefied species richness


###################################
## fig 20

## merge in TR_w_TCID from old version - missing capsula galatea as had come up as chromista
data <- read.csv("temp_OBIS_dd_2021_2022-04-20.csv")
merga <- merge(data, obis_DD, by="rec", all = T)
dim(merga) ## 511   7
dim(obis_DD) ## 268   6
dim(data) ## 511   2
write.csv(merga, "OUTPUT.csv")

data <- read.csv("OBIS_DD_4_analysis_2022-04-20.csv")

names(obis_DD)
p <- ggplot(obis_DD) + 
  geom_bar(mapping = aes(x = TR_w_TCID, fill = TR_w_TCID)) +
  labs(x = "", y = "Total")+ 
  scale_fill_viridis_d() +
  theme(text=element_text(size=12), #change font size of all text
        axis.text=element_text(size=12), #change font size of axis text
        axis.title=element_text(size=12), #change font size of axis titles
        plot.title=element_text(size=12), #change font size of plot title
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) 
p_ed <- p + theme(legend.position = "none")

data1 <- data                                                 # Replicate original data
data1$x <- factor(data1$x,                                    # Change ordering manually
                  levels = c("phylum", "class", "order", "family", "genus", "species", "temp name"))

ggplot(data1, aes(x, y)) +                                    # Manually ordered barchart
  geom_bar(stat = "identity")

reorder()

install.packages("worms")
library(worms)    
help(worms)    

data(northseamacrozoobenthos)

#########################################

## SPP MATRIX TABLES REDO

#data <- read.csv("temp_all_spp_SITE_2022-04-24.csv")
#data <- read.csv("temp_all_spp_SA_2022-04-24.csv")
#data <- read.csv("temp_all_spp_CA_2022-04-24.csv")
#data <- read.csv("temp_all_spp_REG_2022-04-24.csv")

#data <- read.csv("temp_named_spp_SITE_2022-04-24.csv")
#data <- read.csv("temp_named_spp_SA_2022-04-24.csv")
data <- read.csv("temp_named_spp_CA_2022-04-24.csv")
data <- read.csv("temp_named_spp_REG_2022-04-24.csv")

head(data)
dim(data)

## abundance matrix
data.matrix <- sample2matrix(data)
data.matrix[1:5, 1:5]

## transpose
data2 <- t(data.matrix)
data2[1:5, 1:5]
write.csv(data2, "OUTPUT.csv")

##########################
## w presence absense from fuzzy sim

data.presabs <- splist2presabs(data, sites.col = "Site",
                               sp.col = "Species", keep.n = FALSE)
head(data.presabs)

## transpose
test <- t(data.presabs)
head(test)
test <- t(as.matrix(data.presabs[,-1]))
colnames(test) <- data.presabs$Site
head(test)
write.csv(test, "OUTPUT.csv")
## worked

## as.matrix not working
data.matrix <- as.matrix(data)
test <- make_comb_mat(data)
## error Can not find columns which are logical or only contain 0 or 1.
write.csv(data.matrix, "OUTPUT.csv")
## not worked

##########################

## relative abundance polychaetes fig 23

DD_ed <- read.csv("DD_PUBLISHED_4analysis_ed_2022-04-23.csv")
DD_poly <- DD_ed %>% 
  filter(pelagic =="") %>%
  filter(Class_ed == "Polychaeta")
dim(DD_poly) ##7899

DD_poly2 <- DD_poly %>%
  group_by(Family_ed, contract_area) %>%
  summarise(tot = sum(Number.of.individuals_ed, na.rm = T))

DD_poly3 <- DD_poly2 %>%
  group_by(contract_area) %>%
  summarise(tot.region = sum(tot))

## do proportions and arrange so consistency in display
merga <- merge(DD_poly2, DD_poly3, by="ContractorID", all= T) %>%
  mutate(prop = tot/tot.region) %>%
  arrange(desc(tot.region))

## use this version- families are same colour
merga <- merge(DD_poly2, DD_poly3, by="ContractorID", all= T) %>%
  mutate(prop = tot/tot.region) %>%
  arrange(desc(Family_ed))

table(merga$ContractorID, merga$prop)

## plotting

rel_abundance_bplot <- ggplot(merga, aes(x = ContractorID, y = prop, fill = Family_ed))+
  geom_bar(stat = "identity", width = 0.8, aes(group = as.factor(ContractorID)))+
  theme_minimal()+
  #scale_fill_manual(values = cols, name = "Species")+
  # geom_text(aes(y = ypos, label = paste0(round(prop, digits = 0), "%")), color = "black", size=2)+
  theme(legend.key.size = unit(0.4, "cm"), legend.key.height =  unit(0.1, "cm"), legend.title = element_text(face = "bold"))+
  xlab("Contractor ID")+
  ylab("Relative abundance (%)") +
  labs(title="DeepData: Polychaete family records: relative abundance by contractor area")  +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 55, hjust = 1))


