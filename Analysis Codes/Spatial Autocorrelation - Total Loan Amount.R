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

# Load datasets
kiva <- read.csv("kivaData_augmented.csv")

# Extracting the world region from country
# Creating a new column called "world_region"
kiva$world_region <- countrycode(sourcevar = kiva[,"country"], origin = "country.name", 
                                 destination = "continent")

# Filtering points in Africa continent
kiva_africa <- kiva[kiva$world_region=='Africa',]
# Removing data points from Egypt (not Sub-Saharan Africa)
kiva_africa <- kiva_africa[!(kiva_africa$country=='Egypt'),]

# Total loan amount over the years
loan_amount <- data.frame(aggregate(kiva_africa$loan_amount, by=list(country=kiva_africa$country), FUN=sum))
colnames(loan_amount) <- c("country", "total_loan_amount")
loan_amount$total_loan_amount <- loan_amount$total_loan_amount/1000000 # divide by 1 million
loan_amount$iso_a2 = countrycode(sourcevar = loan_amount[,"country"], origin = "country.name", 
                                        destination = "iso2c")

# Plot on total loan amounts
africa = world %>% 
  filter(continent == "Africa", !is.na(iso_a2)) %>% 
  left_join(loan_amount, by = "iso_a2") %>% 
  dplyr::select(country, total_loan_amount) 
tm_shape(africa) + tm_polygons("total_loan_amount", title = "Total Loan Amount")

# Plot of the African countries with available loan amounts
# Removing African countries without any loan amount
sub_africa = na.omit(africa)
tm_shape(sub_africa) + tm_polygons("total_loan_amount", title = "Total Loan Amount")

# Drop Madagascar as she has no neighbours
sub_africa <- sub_africa[-c(26),]
# Converting to sp object
sub_africa.sp <- as_Spatial(sub_africa)
tm_shape(sub_africa.sp) + tm_polygons("total_loan_amount", title = "Total Loan Amount")

# Neighbour Selection Method: Queen's case
# Any countries touching the borders will be considered as neighbours
# Neighbour list
sub_africa.sp.nb <- poly2nb(sub_africa.sp) 
sub_africa.sp.nb # Average number of links: 3.357143 

# Convert the neighbour list to a listw object 
sub_africa.sp.lw <- nb2listw(sub_africa.sp.nb) 
sub_africa.sp.lw

# Plotting the lagged means
sub_africa.sp$total_loan_amount.lagged.means <- lag.listw(sub_africa.sp.lw, sub_africa.sp$total_loan_amount) 
tm_shape(sub_africa.sp) + tm_polygons(col = 'total_loan_amount.lagged.means', title= 'Loan Amount Lagged Means') +
  tm_layout(legend.bg.color = "white")

# Computing Global Moran's I statistic
# 1. Assume MPI_Country is normally distributed
moran.test(sub_africa.sp$total_loan_amount, sub_africa.sp.lw, randomisation=FALSE)

# 2. Simulation-based approach (Monte Carlo)
moran.mc(sub_africa.sp$total_loan_amount, sub_africa.sp.lw, 10000)

# Computing local Moran's I statistic
# Compute the local Moran's I
sub_africa.sp$localI <- localmoran(sub_africa.sp$total_loan_amount, sub_africa.sp.lw)[, 1]

# Plot on local Moran's I values
tm_shape(sub_africa.sp) + tm_polygons(col= 'localI',title= "Local Moran's
I",legend.format=list(flag= "+")) + tm_style('col_blind')

# Computing the p-value of local Moran's I
sub_africa.sp$pval <- localmoran(sub_africa.sp$total_loan_amount, sub_africa.sp.lw)[, 5]

# Plot on local Moran's I p-values
tm_shape(sub_africa.sp) + tm_polygons(col= 'pval',title= "p-value", border.col = "black", breaks= c(0, 0.01, 0.05, 0.10, 1), palette = "-Greens") 

# Spatial Autoregression
# 1. Conditional autoregressive (CAR) model, First-order dependency
# Only depends on the neighbours
sar.res <- spautolm(total_loan_amount ~ 1, listw=sub_africa.sp.lw, data=sub_africa.sp, family = 'CAR') 
sar.res

sar.res$lambda.se
sar.res$lambda + c(-2,2)*sar.res$lambda.se

# 2. Simultaneous autoregressive (SAR) model, Second-order dependency
# Depends on neighbours' neighbours
sar.res2 <- spautolm(total_loan_amount ~ 1, listw=sub_africa.sp.lw, data=sub_africa.sp, family = 'SAR') 
sar.res2

sar.res2$lambda.se
sar.res2$lambda + c(-2,2)*sar.res2$lambda.se