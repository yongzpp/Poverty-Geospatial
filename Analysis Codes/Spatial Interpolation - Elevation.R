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

# Creating sp dataframe for temp
africa_elev.sp = SpatialPointsDataFrame(coords = MPI_africa[c("longitude", "latitude")], # longitude and latitude is reversed
                                        data = data.frame(MPI_africa$Elevation))

# Create boundary Africa
africa = world %>% filter(continent == "Africa") 
africa_bound <- st_union(africa)
africa_bound <- st_transform(africa_bound, 4326)
africa_bound.sp <- as(africa_bound, "Spatial")

africa_elev.sp@bbox <- africa_bound.sp@bbox
# Plot points Elevation 
# Plot points Elevation 
tm_shape(africa_bound.sp) + tm_polygons()+
  tm_shape(africa_elev.sp)+
  tm_dots()

tm_shape(africa_bound.sp) + tm_polygons(col='grey') +
  tm_shape(africa_elev.sp) +
  tm_dots(col="MPI_africa.Elevation", palette = "YlGn", auto.palette.mapping = FALSE,size=0.1,
          title="Sampled Elevation\n(in metres)") +
  tm_legend(legend.outside=FALSE)
# Thiessen polygons for Elevation
library(spatstat)  # Used for the dirichlet tessellation function
library(maptools)  # Used for conversion from SPDF to ppp
library(raster)    # Used to clip out thiessen polygons

# Create a tessellated surface
th  <-  as(dirichlet(as.ppp(africa_elev.sp)), "SpatialPolygons")
proj4string(africa_elev.sp) <- proj4string(MPI_africa.sp)
proj4string(th) <- proj4string(africa_elev.sp)

crs(africa_elev.sp) <- crs(th)
th.z     <- over(th, africa_elev.sp, fn=mean)
th.spdf  <-  SpatialPolygonsDataFrame(th, th.z)
th.clp   <- raster::intersect(as_Spatial(africa),th.spdf)

# Predicted Elevation in metres
tm_shape(th.clp) + 
  tm_polygons(col="MPI_africa.Elevation", palette="YlGn", auto.palette.mapping=FALSE,
              title="Predicted Elevation\n(in metres)") +
  tm_legend(legend.outside=TRUE)


#IDW
library(gstat) # Use gstat's idw routine
library(sp)    # Used for the spsample function

# Create an empty grid where n is the total number of cells
grd              <- as.data.frame(spsample(africa_elev.sp, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object

# Add P's projection information to the empty grid
proj4string(africa_elev.sp) <- proj4string(africa_elev.sp) # Temp fix until new proj env is adopted
proj4string(grd) <- proj4string(africa_elev.sp)

#africa_elev.sp <- as(africa_temp, "Spatial")
# Interpolate the grid cells using a power value of 2 (idp=2.0)
africa_elev.sp.idw <- gstat::idw(MPI_africa.Elevation  ~ 1, africa_elev.sp, newdata=grd, idp=2.0)

# Convert to raster object then clip to Africa
r       <- raster(africa_elev.sp.idw)
r.m     <- mask(r, africa_bound.sp)

# Plot
tm_shape(r.m) + 
  tm_raster(n=10,palette = "YlGn", auto.palette.mapping = FALSE,
            title="Predicted elevation \n(in metres)") + 
  tm_shape(africa_elev.sp) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

# Leave-one-out validation routine
IDW.out <- vector(length = length(africa_elev.sp))
for (i in 1:length(africa_elev.sp)) {
  IDW.out[i] <- idw(MPI_africa.Elevation ~ 1, africa_elev.sp[-i,], africa_elev.sp[i,], idp=2.0)$var1.pred
}

# Plot the differences
OP <- par(pty="s", mar=c(4,3,0,0))
plot(IDW.out ~ africa_elev.sp$MPI_africa.Elevation, asp=1, xlab="Observed", ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5))
abline(lm(IDW.out ~ africa_elev.sp$MPI_africa.Elevation), col="red", lw=2,lty=2)
abline(0,1)
par(OP)

# Compute RMSE
sqrt( sum((IDW.out - africa_temp.sp$MPI_africa.Temperature)^2) / length(africa_temp.sp))
#========================================================
# Define the 1st order polynomial equation
f.1 <- as.formula(MPI_africa.Elevation  ~ X + Y) 

# Add X and Y to P
africa_elev.sp$X <- coordinates(africa_elev.sp)[,1]
africa_elev.sp$Y <- coordinates(africa_elev.sp)[,2]

# Run the regression model
lm.1 <- lm( f.1, data=africa_elev.sp)

# Use the regression model output to interpolate the surface
dat.1st <- SpatialGridDataFrame(grd, data.frame(var1.pred = predict(lm.1, newdata=grd))) 

# Clip the interpolated raster to Texas
r   <- raster(dat.1st)
r.m <- mask(r, africa_bound.sp)
plot(r)
# Plot the map
tm_shape(r.m) + 
  tm_raster(n=10, palette="YlGn", auto.palette.mapping=FALSE, 
            title="Predicted Elevation \n(in metres)") +
  tm_shape(africa_elev.sp) + tm_dots(size=0.1) +
  tm_legend(legend.outside=TRUE)
#==========================================================
# Define the 2nd order polynomial equation
f.2 <- as.formula(MPI_africa.Elevation  ~ X + Y + I(X*X)+I(Y*Y) + I(X*Y))

# Add X and Y to P
africa_elev.sp$X <- coordinates(africa_elev.sp)[,1]
africa_elev.sp$Y <- coordinates(africa_elev.sp)[,2]

# Run the regression model
lm.2 <- lm( f.2, data=africa_elev.sp)

# Use the regression model output to interpolate the surface
dat.2nd <- SpatialGridDataFrame(grd, data.frame(var1.pred = predict(lm.2, newdata=grd))) 

# Clip the interpolated raster to Texas
r   <- raster(dat.2nd)
r.m <- mask(r, africa_bound.sp)

# Plot the map
tm_shape(r.m) + 
  tm_raster(n=10, palette="YlGn", auto.palette.mapping=FALSE,
            title="Predicted Elevation \n(in metres)") +
  tm_shape(africa_elev.sp) + tm_dots(size=0.1) +
  tm_legend(legend.outside=TRUE)
