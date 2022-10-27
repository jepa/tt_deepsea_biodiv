library(picante)

# Make and export community matrix
data <- read.csv("data-raw/TEST_LIT+DD_ALL_SPP_2022-10-10.csv",header=T,sep=",",fileEncoding="latin1")[-2,]
community_matrix <- picante::sample2matrix(data)