# BT4015 Project
# Import libraries
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

# Loading the datasets
MPI <- read.csv("MPIData_augmented.csv")

# Descriptive analysis
# Removing null values
# summary(MPI)
sum(is.na(MPI$latitude))
sum(is.na(MPI$longitude))
MPI = MPI[!is.na(MPI$latitude),]
MPI = MPI[!is.na(MPI$longitude),]

# Filtering points in Sub-Saharan Africa only
MPI_africa <- MPI[MPI$World.region=='Sub-Saharan Africa',]

# Cleaning up the data points which clearly do not lie in Africa
# latitude (positive - north, negative - south)
MPI_africa <- MPI_africa[(MPI_africa$latitude < 37.21 & MPI_africa$latitude > -34.41),]
# longitude (positive - east, negative north)
MPI_africa <- MPI_africa[(MPI_africa$longitude < 57.27 & MPI_africa$longitude > -17.31),]
# Wiki source - geography of africa

# Create a new column of 'ISO Country Code 2 digit'
MPI_africa$'iso_a2' = countrycode(MPI_africa$ISO.country.code, origin = "iso3c", destination = "iso2c")

# Creating sp dataframe
MPI_africa.sp = SpatialPointsDataFrame(coords = MPI_africa[c("longitude", "latitude")], # longitude and latitude is reversed
                                       data = MPI_africa, proj4string = CRS("+init=epsg:4326"))
# ==========================================================================================
# Plotting the map
tmap_mode('plot')

# Interpolation
library(rgdal)
library(tmap)
# Remove rows with Elevation that have NA values for MPI_africa
MPI_africa = MPI_africa[!is.na(MPI_africa$Elevation),]

# Remove rows with TimeToCity that have NA values for MPI_africa
#MPI_africa = MPI_africa[!is.na(MPI_africa$TimeToCity),]

# Remove rows with precipitation that have NA values for MPI_africa
MPI_africa = MPI_africa[!is.na(MPI_africa$precipitation),]

# Remove rows with temperature that have NA values for MPI_africa
MPI_africa = MPI_africa[!is.na(MPI_africa$Temperature),]

# Creating sp dataframe for elev
africa_elev.sp = SpatialPointsDataFrame(coords = MPI_africa[c("longitude", "latitude")], # longitude and latitude is reversed
                                        data = data.frame(MPI_africa$Elevation))

# Creating sp dataframe for TimeToCity
#africa_ttc.sp = SpatialPointsDataFrame(coords = MPI_africa[c("longitude", "latitude")], # longitude and latitude is reversed
                                       #data = data.frame(MPI_africa$TimeToCity))

# Creating sp dataframe for precipitation
africa_prep.sp = SpatialPointsDataFrame(coords = MPI_africa[c("longitude", "latitude")], # longitude and latitude is reversed
                                        data = data.frame(MPI_africa$precipitation))

# Creating sp dataframe for temp
africa_temp.sp = SpatialPointsDataFrame(coords = MPI_africa[c("longitude", "latitude")], # longitude and latitude is reversed
                                        data = data.frame(MPI_africa$Temperature))

africa = world %>% filter(continent == "Africa") 
africa_bound <- st_union(africa)
africa_bound <- st_transform(africa_bound, 4326)
africa_bound.sp <- as(africa_bound, "Spatial")

africa_elev.sp@bbox <- africa_bound.sp@bbox
#africa_ttc.sp@bbox <- africa_bound.sp@bbox
africa_prep.sp@bbox <- africa_bound.sp@bbox
africa_temp.sp@bbox <- africa_bound.sp@bbox
#IDW
library(gstat) # Use gstat's idw routine
library(sp)    # Used for the spsample function

#=========== Elevation raster==================
# Create an empty grid where n is the total number of cells
grd              <- as.data.frame(spsample(africa_elev.sp, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object

# Add P's projection information to the empty grid
proj4string(africa_elev.sp) <- proj4string(africa_elev.sp) # Temp fix until new proj env is adopted
proj4string(grd) <- proj4string(africa_elev.sp)

# Interpolate the grid cells using a power value of 2 (idp=2.0)
africa_elev.sp.idw <- gstat::idw(MPI_africa.Elevation  ~ 1, africa_elev.sp, newdata=grd, idp=2.0)

# Convert to raster object then clip to Africa
r       <- raster(africa_elev.sp.idw)
r.m.elev     <- mask(r, africa_bound.sp)

# #=========== Time to city raster==================
# # Create an empty grid where n is the total number of cells
# grd              <- as.data.frame(spsample(africa_ttc.sp, "regular", n=50000))
# names(grd)       <- c("X", "Y")
# coordinates(grd) <- c("X", "Y")
# gridded(grd)     <- TRUE  # Create SpatialPixel object
# fullgrid(grd)    <- TRUE  # Create SpatialGrid object
# 
# # Add P's projection information to the empty grid
# proj4string(africa_ttc.sp) <- proj4string(africa_ttc.sp) # TimeToCity fix until new proj env is adopted
# proj4string(grd) <- proj4string(africa_ttc.sp)
# 
# 
# # Interpolate the grid cells using a power value of 2 (idp=2.0)
# africa_ttc.sp.idw <- gstat::idw(MPI_africa.TimeToCity ~ 1, africa_ttc.sp, newdata=grd, idp=2.0)
# 
# # Convert to raster object then clip to Africa
# r       <- raster(africa_ttc.sp.idw)
# r.m.ttc     <- mask(r, africa_bound.sp)


#=========== Precipitation raster==================
# Create an empty grid where n is the total number of cells
grd              <- as.data.frame(spsample(africa_prep.sp, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object

# Add P's projection information to the empty grid
proj4string(africa_prep.sp) <- proj4string(africa_prep.sp) # precipitation fix until new proj env is adopted
proj4string(grd) <- proj4string(africa_prep.sp)

#africa_prep.sp <- as(africa_precipitation, "Spatial")
# Interpolate the grid cells using a power value of 2 (idp=2.0)
africa_prep.sp.idw <- gstat::idw(MPI_africa.precipitation ~ 1, africa_prep.sp, newdata=grd, idp=2.0)

# Convert to raster object then clip to Africa
r       <- raster(africa_prep.sp.idw)
r.m.precip     <- mask(r, africa_bound.sp)


#=========== Temperature raster==================
# Create an empty grid where n is the total number of cells
grd              <- as.data.frame(spsample(africa_temp.sp, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object

# Add P's projection information to the empty grid
proj4string(africa_temp.sp) <- proj4string(africa_temp.sp) # Temp fix until new proj env is adopted
proj4string(grd) <- proj4string(africa_temp.sp)

#africa_temp.sp <- as(africa_temp, "Spatial")
# Interpolate the grid cells using a power value of 2 (idp=2.0)
africa_temp.sp.idw <- gstat::idw(MPI_africa.Temperature ~ 1, africa_temp.sp, newdata=grd, idp=2.0)

# Convert to raster object then clip to Texas
r       <- raster(africa_temp.sp.idw)
r.m.temp     <- mask(r, africa_bound.sp)

#===============Scaling==========================
r.m.elev.scaled <- (r.m.elev-minValue(r.m.elev))/(maxValue(r.m.elev)-minValue(r.m.elev))
#r.m.ttc.scaled <- (r.m.ttc-minValue(r.m.ttc))/(maxValue(r.m.ttc)-minValue(r.m.ttc))
r.m.precip.scaled <- (r.m.precip-minValue(r.m.precip))/(maxValue(r.m.precip)-minValue(r.m.precip))
r.m.temp.scaled <- (r.m.temp-minValue(r.m.temp))/(maxValue(r.m.temp)-minValue(r.m.temp))

#================= unifying dimensions and resolutions of rasters============
raster_elev <- resample(r.m.elev.scaled,r.m.temp.scaled)
raster_precip <- resample(r.m.precip.scaled,r.m.temp.scaled)
#raster_ttc <- resample(r.m.ttc.scaled,r.m.temp.scaled)

#================ Combine Rasters=================
raster_combine <- raster_elev - raster_precip - r.m.temp.scaled

# Plot
tm_shape(raster_combine) + 
  tm_raster(n=10,palette = "-Spectral", auto.palette.mapping = FALSE,
            title="Raster combine") + 
  tm_legend(legend.outside=TRUE)


