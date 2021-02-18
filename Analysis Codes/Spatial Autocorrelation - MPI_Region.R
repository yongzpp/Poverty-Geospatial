# Spatial Autocorrelation (on sub-national regions)

# Import packages
library(sf)
library(sp)
library(tmap)
library(maptools)
library(GISTools)
library(tidyverse)
library(raster)
library(spData)
library(countrycode)
library(dplyr)
library(spdep)
library(spatialEco)

# From the spatial autocorrelation on MPI, we can tell that there is a cluster at Niger
# Niger MPI value: 0.605
# Niger local Moran statistic: 1.042341
# Niger local Moran p-value: 0.001868499
# The Niger local Moran statistic is statistically significant.

# Neighbours of Niger: Chad, Mali, Benin, Nigeria, Cameroon, Burkina Faso

# Load dataset
MPI <- read.csv("MPIData_augmented.csv")

# Filtering points in Sub-Saharan Africa only
MPI_africa <- MPI[MPI$World.region=='Sub-Saharan Africa',]
# Create a new column of 'ISO Country Code 2 digit'
MPI_africa$'iso_a2' = countrycode(MPI_africa$ISO.country.code, origin = "iso3c", destination = "iso2c")
# Getting the DataFrame containing every country MPI level
region_mpi_df <- data.frame(unique(MPI_africa[,c("Country", "Sub.national.region", "MPI_Region", "iso_a2")]))

# Loading the various sp objects
# Found in the "sp objects - Regional MPI" folder
niger <- readRDS(file = "sp objects - Regional MPI/gadm36_NER_1_sp.rds")
chad <- readRDS(file = "sp objects - Regional MPI/gadm36_TCD_1_sp.rds")
mali <- readRDS(file = "sp objects - Regional MPI/gadm36_MLI_1_sp.rds")
benin <- readRDS(file = "sp objects - Regional MPI/gadm36_BEN_1_sp.rds")
nigeria <- readRDS(file = "sp objects - Regional MPI/gadm36_NGA_1_sp.rds")
cameroon <- readRDS(file = "sp objects - Regional MPI/gadm36_CMR_1_sp.rds")
burkina_faso <- readRDS(file = "sp objects - Regional MPI/gadm36_BFA_1_sp.rds")

niger0 <- readRDS(file = "sp objects - Regional MPI/gadm36_NER_0_sp.rds")
chad0 <- readRDS(file = "sp objects - Regional MPI/gadm36_TCD_0_sp.rds")
mali0 <- readRDS(file = "sp objects - Regional MPI/gadm36_MLI_0_sp.rds")
benin0 <- readRDS(file = "sp objects - Regional MPI/gadm36_BEN_0_sp.rds")
nigeria0 <- readRDS(file = "sp objects - Regional MPI/gadm36_NGA_0_sp.rds")
cameroon0 <- readRDS(file = "sp objects - Regional MPI/gadm36_CMR_0_sp.rds")
burkina_faso0 <- readRDS(file = "sp objects - Regional MPI/gadm36_BFA_0_sp.rds")

cluster <- rbind(niger, chad, mali, benin, nigeria, cameroon, burkina_faso)
country_border <- rbind(niger0, chad0, mali0, benin0, nigeria0, cameroon0, burkina_faso0)

# Plotting the map
tmap_mode('plot')

# Map of the regions in the cluster identified
tm_shape(cluster) + tm_polygons()

# Map of the country borders
tm_shape(country_border) + tm_polygons("NAME_0", legend.show = FALSE) + tm_text("NAME_0", size = 1.0)

# Data cleaning
cluster$NAME_1[7] <- "Tillaberi"
cluster$NAME_1[9] <- "Barh El Gazal"
cluster$NAME_1[11] <- "Borkou/Tibesti"
cluster$NAME_1[12] <- "Chari Baguirmi"
cluster$NAME_1[13] <- "Ennedi"
cluster$NAME_1[14] <- "Ennedi"
cluster$NAME_1[15] <- "Guéra"
cluster$NAME_1[22] <- "Mayo Kebbi Est"
cluster$NAME_1[23] <- "Mayo Kebbi Ouest"
cluster$NAME_1[24] <- "Moyen Chari"
cluster$NAME_1[25] <- "Ouaddaï"
cluster$NAME_1[28] <- "Tandjilé"
cluster$NAME_1[29] <- "Borkou/Tibesti"
cluster$NAME_1[30] <- "N'Djamena"
#cluster$NAME_1[33] <- "S�gou"
#cluster$NAME_1[35]
cluster$NAME_1[38] <- "S�gou"
#cluster$NAME_1[40] 
cluster$NAME_1[42] <- "Atacora"
cluster$NAME_1[47] <- "Couffo"
cluster$NAME_1[50] <- "Ouémé"
#cluster$NAME_1[68] 
#cluster$NAME_1[79]
#cluster$NAME_1[97]
#cluster$NAME_1[112]
#cluster$NAME_1[120]
#cluster$NAME_1[123]
#cluster$NAME_1[128]

# Merge cluster and region_mpi_df
cluster_mpi <- sp::merge(x = cluster, y = region_mpi_df, by.x = "NAME_1", by.y = "Sub.national.region", duplicateGeoms = T)

# Map of the regions in the cluster identified
tm_shape(cluster_mpi) + tm_polygons("MPI_Region", legend.show = TRUE, title = "MPI_Region in Cluster")

# Plot of the African countries with available MPI
# Removing African countries without any MPI values
sub_cluster_mpi <- sp.na.omit(cluster_mpi, col.name = "MPI_Region")
tm_shape(sub_cluster_mpi) + tm_polygons("MPI_Region", legend.show = TRUE, title = "MPI_Region in Cluster")

# Neighbour Selection Method: Queen's case
# Any countries touching the borders will be considered as neighbours
# Neighbour list
sub_cluster_mpi.nb <- poly2nb(sub_cluster_mpi) 
sub_cluster_mpi.nb # Average number of links: 5.666667

# Convert the neighbour list to a listw object 
sub_cluster_mpi.lw <- nb2listw(sub_cluster_mpi.nb) 
sub_cluster_mpi.lw

# Plotting the lagged means
sub_cluster_mpi$MPI_Region.lagged.means <- lag.listw(sub_cluster_mpi.lw, sub_cluster_mpi$MPI_Region) 
tm_shape(sub_cluster_mpi) + tm_polygons(col = 'MPI_Region.lagged.means', title= 'MPI_Region Lagged Means') +
  tm_layout(legend.bg.color = "white")

# Computing Global Moran's I statistic
# 1. Assume MPI_Country is normally distributed
moran.test(sub_cluster_mpi$MPI_Region, sub_cluster_mpi.lw, randomisation=FALSE)

# 2. Simulation-based approach (Monte Carlo)
moran.mc(sub_cluster_mpi$MPI_Region, sub_cluster_mpi.lw, 10000)

# Computing local Moran's I statistic
# Compute the local Moran's I
sub_cluster_mpi$localI <- localmoran(sub_cluster_mpi$MPI_Region, sub_cluster_mpi.lw)[, 1]

# Plot on local Moran's I values
tm_shape(sub_cluster_mpi) + tm_polygons(col= 'localI',title= "Local Moran's
I",legend.format=list(flag= "+")) + tm_style('col_blind')

# Computing the p-value of local Moran's I
sub_cluster_mpi$pval <- localmoran(sub_cluster_mpi$MPI_Region, sub_cluster_mpi.lw)[, 5]

# Plot on local Moran's I p-values
tm_shape(sub_cluster_mpi) + tm_polygons(col= 'pval',title= "p-value", border.col = "black", breaks= c(0, 0.01, 0.05, 0.10, 1), palette = "-Greens")

# Spatial Autoregression
# 1. Conditional autoregressive (CAR) model, First-order dependency
# Only depends on the neighbours
sar.res <- spautolm(MPI_Region ~ 1, listw=sub_cluster_mpi.lw, data=sub_cluster_mpi, family = "CAR") 
summary(sar.res)

sar.res$lambda.se
sar.res$lambda + c(-2,2)*sar.res$lambda.se

# 2. Conditional autoregressive (CAR) model, First-order dependency
# Only depends on the neighbours
sar.res2 <- spautolm(MPI_Region ~ 1, listw=sub_cluster_mpi.lw, data=sub_cluster_mpi, family = "SAR") 
summary(sar.res2)

sar.res2$lambda.se
sar.res2$lambda + c(-2,2)*sar.res2$lambda.se

# Identifying the clusters
list_of_MPI_region <- cluster_mpi$MPI_Region
# identifying countries with MPI_region > 0.60
index <- which(list_of_MPI_region > 0.60)
# removing the index 96 as it represents Cameroon
index <- index[index != 96]
# Chad represents index 9 to 31
central_cluster_index <- index[!index %in% 9:31]
east_cluster_index <- index[index %in% 9:31]

central_cluster <- cluster_mpi[central_cluster_index,]
east_cluster <- cluster_mpi[east_cluster_index,]

# Map of the central cluster
tm_shape(central_cluster) + tm_polygons("NAME_1", legend.show = FALSE) + tm_text("NAME_1", size = 1.0) +
  tm_layout(title = "Central cluster")

# Map of the eastern cluster
tm_shape(east_cluster) + tm_polygons("NAME_1", legend.show = FALSE) + tm_text("NAME_1", size = 1.0) +
  tm_layout(title = "Eastern cluster")
