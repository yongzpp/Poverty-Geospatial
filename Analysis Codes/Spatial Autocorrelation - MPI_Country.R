# BT4015 Project
# Spatial Autocorrelation

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

# Load dataset
MPI <- read.csv("MPIData_augmented.csv")

# Filtering points in Sub-Saharan Africa only
MPI_africa <- MPI[MPI$World.region=='Sub-Saharan Africa',]
# Create a new column of 'ISO Country Code 2 digit'
MPI_africa$'iso_a2' = countrycode(MPI_africa$ISO.country.code, origin = "iso3c", destination = "iso2c")
# Getting the DataFrame containing every country MPI level
country_mpi_df <- data.frame(unique(MPI_africa[,c("Country", "MPI_Country", "iso_a2")]))

# Plotting the map
tmap_mode('plot')

# Plot on MPI_Country
africa = world %>% 
  filter(continent == "Africa", !is.na(iso_a2)) %>% 
  left_join(country_mpi_df, by = "iso_a2") %>% 
  dplyr::select(Country, MPI_Country)
# Plot of the entire Africa
tm_shape(africa) + tm_polygons("MPI_Country")

# Plot of the African countries with available MPI
# Removing African countries without any MPI values
sub_africa = na.omit(africa)
tm_shape(sub_africa) + tm_polygons("MPI_Country", alpha = 0.7)

# Drop Lesotho and Madagascar as they have no neighbours
sub_africa <- sub_africa[-c(5,32),]
# Converting to sp object
sub_africa.sp <- as_Spatial(sub_africa)
tm_shape(sub_africa.sp) + tm_polygons("MPI_Country", alpha = 0.7)

# Neighbour Selection Method: Queen's case
# Any countries touching the borders will be considered as neighbours
# Neighbour list
sub_africa.sp.nb <- poly2nb(sub_africa.sp) 
sub_africa.sp.nb # Average number of links: 4.057143

# Convert the neighbour list to a listw object 
sub_africa.sp.lw <- nb2listw(sub_africa.sp.nb) 
sub_africa.sp.lw 

# Plotting the lagged means
sub_africa.sp$MPI_Country.lagged.means <- lag.listw(sub_africa.sp.lw, sub_africa.sp$MPI_Country) 
tm_shape(sub_africa.sp) + tm_polygons(col = 'MPI_Country.lagged.means', title= 'MPI Lagged Means') +
  tm_layout(legend.bg.color = "white")

# Computing Global Moran's I statistic
# 1. Assume MPI_Country is normally distributed
moran.test(sub_africa.sp$MPI_Country, sub_africa.sp.lw, randomisation=FALSE)

# 2. Simulation-based approach (Monte Carlo)
moran.mc(sub_africa.sp$MPI_Country, sub_africa.sp.lw, 10000)

# Computing local Moran's I statistic
# Compute the local Moran's I
sub_africa.sp$localI <- localmoran(sub_africa.sp$MPI_Country, sub_africa.sp.lw)[, 1]

# Plot on local Moran's I values
tm_shape(sub_africa.sp) + tm_polygons(col= 'localI',title= "Local Moran's
I",legend.format=list(flag= "+")) + tm_style('col_blind')

# Computing the p-value of local Moran's I
sub_africa.sp$pval <- localmoran(sub_africa.sp$MPI_Country, sub_africa.sp.lw)[, 5]

# Plot on local Moran's I p-values
tm_shape(sub_africa.sp) + tm_polygons(col= 'pval',title= "p-value", border.col = "black", breaks= c(0, 0.01, 0.05, 0.10, 1), palette = "-Greens") 

# Spatial Autoregression
# 1. Conditional autoregressive (CAR) model, First-order dependency
# Only depends on the neighbours
sar.res <- spautolm(MPI_Country ~ 1, listw=sub_africa.sp.lw, data=sub_africa.sp, family = "CAR") 
sar.res

sar.res$lambda.se
sar.res$lambda + c(-2,2)*sar.res$lambda.se

# 2. Simultaneous autoregressive (SAR) model, Second-order dependency
# Depends on neighbours' neighbours
sar.res2 <- spautolm(MPI_Country ~ 1, listw=sub_africa.sp.lw, data=sub_africa.sp, family = "SAR") 
sar.res2

sar.res2$lambda.se
sar.res2$lambda + c(-2,2)*sar.res2$lambda.se







