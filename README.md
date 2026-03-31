# A Baseline Assessment of Deep-Sea Biodiversity, and Guide for Future Actions to Inform Stewardship, in Trinidad and Tobago

This repository contains the code for a project on the Baseline Assessment of Deep-Sea Biodiversity in Trinidad and Tobago lead by Jamie Edghill. It is partially based on [howlerMoonkey/CCZ_BIODIVERSITY:main](https://github.com/howlerMoonkey/CCZ_BIODIVERSITY).

# Authors

Jaime-Leigh Lue Chin^1^, Muriel Rabone^2^, Juliano Palacios-Abrantes^3^, Judith Gobin^1^, Alana Jute^4^, La Daana Kanhai^1^, Diva J. Amon^1,5,6^

## Affiliations

1.  Department of Life Sciences, The University of the West Indies, St. Augustine, Trinidad and Tobago
2.  Natural History Museum, London, UK
3.  Institute for the Oceans and Fisheries, The University of British Columbia, Vancouver, Canada
4.  Institute of Marine Affairs, Chaguaramas, Trinidad and Tobago
5.  Marine Science Institute, University of California, Santa Barbara, Santa Barbara, CA, USA
6.  SpeSeas, D’Abadie, Trinidad and Tobago

*Corresponding author: Jaime-Leigh Lue Chin, jaimeleigh.edghill[at]my.uwi.edu

# Reference 

Chin, J.-L. L., Rabone, M., Palacios-Abrantes, J., Gobin, J., Jute, A., Kanhai, L. D. K. and Amon, D. J., 2026. A baseline assessment of deep-sea biodiversity, and guide for future actions to inform stewardship, in Trinidad and Tobago. Discover Oceans, 3 (1), 18.  

[![Paper -- DOI](https://img.shields.io/badge/Paper-Article-FF69B4)](https://link.springer.com/article/10.1007/s44289-026-00119-3) [![Code -- DOI](https://zenodo.org/badge/908994126.svg)](https://doi.org/10.5281/zenodo.19354802) [![Data -- DOI](https://img.shields.io/badge/Data--DOI-10.5683%2FSP3%2FJMBHER-orange)](https://obis.org/dataset/7eca348a-9a68-4440-9333-c9b685527f96) 

# Files and folders organization

```         
tt_deepsea_biodiv/                    
├── R/                              # Custom R functions used throughout the project
│   ├── *Data_S5.R*                 # original script from `howlerMoonkey/CCZ_BIODIVERSITY:main` only for reference
│   ├── *biodiv_script.Rmd*         # script with the current Chao2 analysis
│   └──  *tt_spatial_analysis.Rmd*   # script for producing maps figures
├── .gitignore                      # Specifies which files/folders Git should ignore
├── LICENSE                         # License file (e.g., MIT, CC BY-NC 4.0)
├── tt_deepsea_biodiv.Rproj         # RStudio project file
└── README.md
```

# Data Availability 

## Primary data
The primary dataset for this article will be provided in the Supplementary Information Table 1 (*Pending*). Additionally, new records currently not included on OBIS will eventually be uploaded.

[![Data -- DOI](https://img.shields.io/badge/Data--DOI-10.5683%2FSP3%2FJMBHER-orange)](https://obis.org/dataset/7eca348a-9a68-4440-9333-c9b685527f96) 

## Support data
- *SAUEEZ_July2015.shp*, This is the Sea Around Us (SAU) shapefile for the world EEZs. Contact the [SAU](http://www.seaaroundus.org/) for a version of it. 
- *bathy_map*, Data avialble from NOAA National Centers for Environmental Information. 2022: ETOPO 2022 15 Arc-Second Global Relief Model. NOAA National Centers for Environmental Information. tools:::Rd_expr_doi("doi.org/10.25921/fd45-gt74")
- *202106_TT_MEEI_dWGS84_offshore_blocks_Trinidad_and_Tobago.shp*, shapefile of oil leases in the T&T EEZ. Contact main author for data.


# External References

-   Rabone, M. and others. 2023. How many metazoan species live in the world’s largest mineral exploration region? Curr. Biol. 33: 2383-2396.e5. <doi:10.1016/j.cub.2023.04.052>

-   [Code for Rabone et al., 2023](https://github.com/howlerMoonkey/CCZ_BIODIVERSITY)


_**As of March of 2026 the paper has been published and this repository has been archived**_

