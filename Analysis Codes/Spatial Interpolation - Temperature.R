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
# Remove rows with temperature that have NA values for MPI_africa
MPI_africa = MPI_africa[!is.na(MPI_africa$Temperature),]

# Creating sp dataframe for temp
africa_temp.sp = SpatialPointsDataFrame(coords = MPI_africa[c("longitude", "latitude")], # longitude and latitude is reversed
                                       data = data.frame(MPI_africa$Temperature))

# Create Africa boundary
africa = world %>% filter(continent == "Africa") 
africa_bound <- st_union(africa)
africa_bound <- st_transform(africa_bound, 4326)
africa_bound.sp <- as(africa_bound, "Spatial")

africa_temp.sp@bbox <- africa_bound.sp@bbox
# Plot points Temperature 
tm_shape(africa_bound.sp) + tm_polygons()+
  tm_shape(africa_temp.sp)+
  tm_dots()

tm_shape(africa_bound.sp) + tm_polygons(col='white') +
  tm_shape(africa_temp.sp) +
  tm_dots(col="MPI_africa.Temperature", palette = "YlOrRd", auto.palette.mapping = FALSE,
          title="Sampled temperature") +
  tm_legend(legend.outside=FALSE)

# Thiessen polygons for temperature
library(spatstat)  # Used for the dirichlet tessellation function
library(maptools)  # Used for conversion from SPDF to ppp
library(raster)    # Used to clip out thiessen polygons

# Create a tessellated surface
th  <-  as(dirichlet(as.ppp(africa_temp.sp)), "SpatialPolygons")
proj4string(africa_temp.sp) <- proj4string(MPI_africa.sp)
proj4string(th) <- proj4string(africa_temp.sp)

crs(africa_temp.sp) <- crs(th)
th.z     <- over(th, africa_temp.sp, fn=mean)
th.spdf  <-  SpatialPolygonsDataFrame(th, th.z)
th.clp   <- raster::intersect(as_Spatial(africa),th.spdf)

# Predicted temperature in celcius
tm_shape(th.clp) + 
  tm_polygons(col="MPI_africa.Temperature", palette="YlOrRd", auto.palette.mapping=FALSE,
              title="Predicted temperature") +
  tm_legend(legend.outside=TRUE)


#IDW
library(gstat) # Use gstat's idw routine
library(sp)    # Used for the spsample function

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
r.m     <- mask(r, africa_bound.sp)

# Plot
tm_shape(r.m) + 
  tm_raster(n=10,palette = "YlOrRd", auto.palette.mapping = FALSE,
            title="Predicted temperature") + 
  tm_shape(P) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

# Leave-one-out validation routine
IDW.out <- vector(length = length(africa_temp.sp))
for (i in 1:length(africa_temp.sp)) {
  IDW.out[i] <- idw(MPI_africa.Temperature ~ 1, africa_temp.sp[-i,], africa_temp.sp[i,], idp=2.0)$var1.pred
}

# Plot the differences
OP <- par(pty="s", mar=c(4,3,0,0))
plot(IDW.out ~ africa_temp.sp$MPI_africa.Temperature, asp=1, xlab="Observed", ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5))
abline(lm(IDW.out ~ africa_temp.sp$MPI_africa.Temperature), col="red", lw=2,lty=2)
abline(0,1)
par(OP)

# Compute RMSE
sqrt( sum((IDW.out - africa_temp.sp$MPI_africa.Temperature)^2) / length(africa_temp.sp))

# Define the 1st order polynomial equation
f.1 <- as.formula(MPI_africa.Temperature ~ X + Y) 

# Add X and Y to P
africa_temp.sp$X <- coordinates(africa_temp.sp)[,1]
africa_temp.sp$Y <- coordinates(africa_temp.sp)[,2]

# Run the regression model
lm.1 <- lm( f.1, data=africa_temp.sp)

# Use the regression model output to interpolate the surface
dat.1st <- SpatialGridDataFrame(grd, data.frame(var1.pred = predict(lm.1, newdata=grd))) 

# Clip the interpolated raster to Texas
r   <- raster(dat.1st)
r.m <- mask(r, africa_bound.sp)

# Plot the map
tm_shape(r.m) + 
  tm_raster(n=10, palette="YlOrRd", auto.palette.mapping=FALSE, 
            title="Predicted temperature") +
  tm_shape(africa_temp.sp) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)


# Define the 2nd order polynomial equation
f.2 <- as.formula(MPI_africa.Temperature ~ X + Y + I(X*X)+I(Y*Y) + I(X*Y))

# Add X and Y to P
africa_temp.sp$X <- coordinates(africa_temp.sp)[,1]
africa_temp.sp$Y <- coordinates(africa_temp.sp)[,2]

# Run the regression model
lm.2 <- lm( f.2, data=africa_temp.sp)

# Use the regression model output to interpolate the surface
dat.2nd <- SpatialGridDataFrame(grd, data.frame(var1.pred = predict(lm.2, newdata=grd))) 

# Clip the interpolated raster to Texas
r   <- raster(dat.2nd)
r.m <- mask(r, africa_bound.sp)

# Plot the map
tm_shape(r.m) + 
  tm_raster(n=10, palette="YlOrRd", auto.palette.mapping=FALSE,
            title="Predicted temperature") +
  tm_shape(africa_temp.sp) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)
