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
library(spatstat)
library(rgdal)

#preparing the data
df <- read.csv("all_loan_theme_merged_with_geo_mpi_regions.csv")
df$continent <- countrycode(sourcevar = df[,"Country"], origin = "country.name", destination = "continent") #assigning continent
df_africa <- df[df$continent=='Africa',] #filtering african countries

#converting geo coordinates
tmp1 <- gsub("\\(","",as.character(df_africa$geo))  
df_africa$geo <- gsub("\\)","",as.character(tmp1))
df_africa <- df_africa %>% separate(geo, c("Y", "X"), sep=',')
df_africa$X <- as.numeric(df_africa$X)
df_africa$Y <- as.numeric(df_africa$Y)

#cleaning & filtering by coordinates
df_africa <- drop_na(df_africa)
min(df_africa$X)
max(df_africa$X)
min(df_africa$Y)
max(df_africa$Y)

# tessalation of african countries
africa = world %>% filter(continent == "Africa") 
x <- as(africa,"Spatial")
regions <- slot(x, 'polygons')
regions <- lapply(regions, function(x) {SpatialPolygons(list(x))})
windows <- lapply(regions, as.owin)
te <- tess(tiles=windows)

#creating window boundary for ppp
africa_bound <- st_union(africa)
y <- as(africa_bound, "Spatial")
regions <- slot(y, 'polygons')
regions <- lapply(regions, function(y) {SpatialPolygons(list(y))})
windows <- lapply(regions, as.owin)

#create ppp spatstat obj
window <- owin(xrange=c(-17.62, 51.13), yrange=c(-34.13,37.35))
pp1 <- ppp(df_africa$X, df_africa$Y, window=windows[[1]])
pp1 <- as.ppp(pp1)
plot(pp1, main=NULL, cols=rgb(0,0,0,.2), pch=20)

#Density based analysis
#quadratic density
Q <- quadratcount(pp1, nx= 10, ny=10)
plot(pp1, pch=20, cols="grey70", main=NULL)  # Plot points
plot(Q, add=TRUE)  # Add quadrat grid

# Compute the density for each quadrat
Q.d <- intensity(Q)

# Plot the density
plot(intensity(Q, image=TRUE), main=NULL, las=1)  # Plot density raster
plot(pp1, pch=20, cex=0.6, col=rgb(0,0,0,.5), add=TRUE)  # Add points

#quadrat density on a tessallated surface
plot(te, main="", las=1)
Q   <- quadratcount(pp1, tess = te)  # Tally counts
Q.d <- intensity(Q)  # Compute density
Q.d
plot(intensity(Q, image=TRUE), las=1, main=NULL)
plot(pp1, pch=20, cex=0.6, col=rgb(1,1,1,.5), add=TRUE)
cl <-  interp.colours(c("lightblue", "orange" ,"red"), te$n)
plot( intensity(Q, image=TRUE), las=1, col=cl, main=NULL)
plot(pp1, pch=20, cex=0.6, col=rgb(0,0,0,.5), add=TRUE)

#TEMPERATURE AS COVARIATE
#plot rasterlayer histogram
temp <- r.m
hist(temp, main=NULL, las=1)
temp.lg <- log(temp)
hist(temp.lg, main=NULL, las=1)

#quadrat density on a tessallated surface for temperature raster
brk  <- c( -Inf, 4.5, 5.5, 6.5, Inf)  # Define the breaks
Zcut <- cut(temp.lg, breaks=brk)  # Classify the raster
E    <- tess(image=Zcut)  # Create a tesselated surface

plot(E, main="", las=1)
Q   <- quadratcount(pp1, tess = E)  # Tally counts
Q.d <- intensity(Q)  # Compute density
Q.d
plot(intensity(Q, image=TRUE), las=1, main=NULL)
plot(pp1, pch=20, cex=0.6, col=rgb(1,1,1,.5), add=TRUE)
cl <-  interp.colours(c("lightblue", "orange" ,"red"), te$n)
plot( intensity(Q, image=TRUE), las=1, col=cl, main=NULL)
plot(pp1, pch=20, cex=0.6, col=rgb(0,0,0,.5), add=TRUE)

#Kernal density raster
K1 <- density(pp1) # Using the default bandwidth
plot(K1, main=NULL, las=1)
contour(K1, add=TRUE)

K2 <- density(pp1, sigma=1) 
plot(K2, main=NULL, las=1)
contour(K2, add=TRUE)

K3 <- density(pp1, kernel = "disc", sigma=1) 
plot(K3, main=NULL, las=1)
contour(K3, add=TRUE)

#Kernal density adjusted for covariate
# Compute rho using the ratio method
rho <- rhohat(pp1, as.im(temp.lg),  method="ratio")
# Generate rho vs covariate plot
plot(rho, las=1, main=NULL, xlab='Logged Temperature')
pred <- predict(rho)
cl   <- interp.colours(c("lightyellow", "orange" ,"red"), 100) # Create color scheme
plot(pred, col=cl, las=1, main=NULL, gamma = 0.25)
K1_vs_pred <- pairs(K1, pred, plot = FALSE)
plot(K1_vs_pred$pred ~ K1_vs_pred$K1, pch=20,
     xlab = "Observed intensity", 
     ylab = "Predicted intensity", 
     col = rgb(0,0,0,0.1))
abline(a=0, b = 1, col = "red")
K1_vs_K1 <- pairs(K1, K1, labels = c("K1a", "K1b"), plot = FALSE)
plot(K1_vs_K1$K1a ~ K1_vs_K1$K1b, pch=20,
     xlab = "Observed intensity", 
     ylab = "Observed intensity")
summary(as.data.frame(K1_vs_pred))
plot(K1_vs_pred$pred ~ K1_vs_pred$K1, pch=20,
     xlab = "Observed intensity", 
     ylab = "Predicted intensity", 
     col = rgb(0,0,0,0.1),
     ylim = c(0, 5))
abline(a=0, b = 1, col = "red")
K1_vs_K1 <- pairs(K1, K1, labels = c("K1a", "K1b"), plot = FALSE)

#Distance based analysis
#average neighbour analysis
mean(nndist(pp1, k=1)) #To compute the average first nearest neighbor distance
mean(nndist(pp1, k=2)) #To compute the average second nearest neighbor distance
ANN <- apply(nndist(pp1, k=1:1000),2,FUN=mean)
plot(ANN ~ eval(1:1000), type="b", main=NULL, las=1)

#K and L functions
K <- Kest(pp1, correction="border")
plot(K, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))
plot(K)
kf.env <- envelope(pp1,Kest,correction="border")
plot(kf.env)

L <- Lest(pp1, main=NULL)
plot(L, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))
plot(L, . -r ~ r, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))
lf.env <- envelope(pp1,Lest,correction="border")
plot(lf.env)

#pair correlation function g
g  <- pcf(pp1)
plot(g, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))
gf.env <- envelope(pp1,Gest,correction="border")
plot(gf.env)


#Hypothesis tests
#test for clustering/dispersion
ann.p <- mean(nndist(pp1, k=1))
ann.p

n     <- 299L               # Number of simulations
ann.r <- vector(length = n) # Create an empty object to be used to store simulated ANN values
for (i in 1:n){
  rand.p   <- rpoint(n=pp1$n, win=windows[[1]])  # Generate random point locations
  ann.r[i] <- mean(nndist(rand.p, k=1))  # Tally the ANN values
}
plot(rand.p, pch=16, main=NULL, cols=rgb(0,0,0,0.5))
hist(ann.r, main=NULL, las=1, breaks=40, col="bisque", xlim=range(ann.p, ann.r))
abline(v=ann.p, col="blue")

n     <- 299L
ann.r <- vector(length=n)
for (i in 1:n){
  Z<- setcov(as.im(temp), W=windows[[1]])
  rand.p   <- rpoint(n=pp1$n, f=Z)  
  ann.r[i] <- mean(nndist(rand.p, k=1))
}

Window(rand.p) <- windows[[1]]  
plot(rand.p, pch=16, main=NULL, cols=rgb(0,0,0,0.5))
hist(ann.r, main=NULL, las=1, breaks=40, col="bisque", xlim=range(ann.p, ann.r))
abline(v=ann.p, col="blue")

N.greater <- sum(ann.r > ann.p)
p <- min(N.greater + 1, n + 1 - N.greater) / (n +1)
p


#Modelling intensity as a function of a covariate
# Create the Poisson point process model
temp.lg.im <- as.im(temp.lg)
PPM1 <- ppm(pp1 ~ temp.lg.im)
plot(effectfun(PPM1, "temp.lg.im", se.fit=TRUE), main=NULL, las=1)
PPM1
PPM0 <- ppm(pp1 ~ 1)
PPM0
pp1$n / area(windows[[1]]) 
anova(PPM0, PPM1, test="LRT")


#PRECIPIRATION AS COVARIATE
#plot rasterlayer histogram
precip <- r.m  #get from interpolation script
hist(precip, main=NULL, las=1)
precip.lg <- log(precip)
hist(precip.lg, main=NULL, las=1)

#quadrat density on a tessallated surface for precipitation raster
brk  <- c( -Inf, 4.5, 5.5, 6.5, Inf)  # Define the breaks
Zcut <- cut(precip.lg, breaks=brk)  # Classify the raster
E    <- tess(image=Zcut)  # Create a tesselated surface

plot(E, main="", las=1)
Q   <- quadratcount(pp1, tess = E)  # Tally counts
Q.d <- intensity(Q)  # Compute density
Q.d
plot(intensity(Q, image=TRUE), las=1, main=NULL)
plot(pp1, pch=20, cex=0.6, col=rgb(1,1,1,.5), add=TRUE)
cl <-  interp.colours(c("lightblue", "orange" ,"red"), te$n)
plot( intensity(Q, image=TRUE), las=1, col=cl, main=NULL)
plot(pp1, pch=20, cex=0.6, col=rgb(0,0,0,.5), add=TRUE)

#Kernal density raster
K1 <- density(pp1) # Using the default bandwidth
plot(K1, main=NULL, las=1)
contour(K1, add=TRUE)

K2 <- density(pp1, sigma=1) 
plot(K2, main=NULL, las=1)
contour(K2, add=TRUE)

K3 <- density(pp1, kernel = "disc", sigma=1) 
plot(K3, main=NULL, las=1)
contour(K3, add=TRUE)

#Kernal density adjusted for covariate
# Compute rho using the ratio method
rho <- rhohat(pp1, as.im(precip.lg),  method="ratio")
# Generate rho vs covariate plot
plot(rho, las=1, main=NULL, xlab='Logged Precipitation')
pred <- predict(rho)
cl   <- interp.colours(c("lightyellow", "orange" ,"red"), 100) # Create color scheme
plot(pred, col=cl, las=1, main=NULL, gamma = 0.25)
K1_vs_pred <- pairs(K1, pred, plot = FALSE)
plot(K1_vs_pred$pred ~ K1_vs_pred$K1, pch=20,
     xlab = "Observed intensity", 
     ylab = "Predicted intensity", 
     col = rgb(0,0,0,0.1))
abline(a=0, b = 1, col = "red")
K1_vs_K1 <- pairs(K1, K1, labels = c("K1a", "K1b"), plot = FALSE)
plot(K1_vs_K1$K1a ~ K1_vs_K1$K1b, pch=20,
     xlab = "Observed intensity", 
     ylab = "Observed intensity")
summary(as.data.frame(K1_vs_pred))
plot(K1_vs_pred$pred ~ K1_vs_pred$K1, pch=20,
     xlab = "Observed intensity", 
     ylab = "Predicted intensity", 
     col = rgb(0,0,0,0.1),
     ylim = c(0, 5))
abline(a=0, b = 1, col = "red")
K1_vs_K1 <- pairs(K1, K1, labels = c("K1a", "K1b"), plot = FALSE)

#Distance based analysis
#average neighbour analysis
mean(nndist(pp1, k=1)) #To compute the average first nearest neighbor distance
mean(nndist(pp1, k=2)) #To compute the average second nearest neighbor distance
ANN <- apply(nndist(pp1, k=1:1000),2,FUN=mean)
plot(ANN ~ eval(1:1000), type="b", main=NULL, las=1)

#K and L functions
K <- Kest(pp1, correction="border")
plot(K, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))
plot(K)
kf.env <- envelope(pp1,Kest,correction="border")
plot(kf.env)

L <- Lest(pp1, main=NULL)
plot(L, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))
plot(L, . -r ~ r, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))
lf.env <- envelope(pp1,Lest,correction="border")
plot(lf.env)

#pair correlation function g
g  <- pcf(pp1)
plot(g, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))
gf.env <- envelope(pp1,Gest,correction="border")
plot(gf.env)


#Hypothesis tests
#test for clustering/dispersion
ann.p <- mean(nndist(pp1, k=1))
ann.p

n     <- 299L               # Number of simulations
ann.r <- vector(length = n) # Create an empty object to be used to store simulated ANN values
for (i in 1:n){
  rand.p   <- rpoint(n=pp1$n, win=windows[[1]])  # Generate random point locations
  ann.r[i] <- mean(nndist(rand.p, k=1))  # Tally the ANN values
}
plot(rand.p, pch=16, main=NULL, cols=rgb(0,0,0,0.5))
hist(ann.r, main=NULL, las=1, breaks=40, col="bisque", xlim=range(ann.p, ann.r))
abline(v=ann.p, col="blue")

n     <- 299L
ann.r <- vector(length=n)
for (i in 1:n){
  Z<- setcov(as.im(precip), W=windows[[1]])
  rand.p   <- rpoint(n=pp1$n, f=Z)  
  ann.r[i] <- mean(nndist(rand.p, k=1))
}

Window(rand.p) <- windows[[1]]  
plot(rand.p, pch=16, main=NULL, cols=rgb(0,0,0,0.5))
hist(ann.r, main=NULL, las=1, breaks=40, col="bisque", xlim=range(ann.p, ann.r))
abline(v=ann.p, col="blue")

N.greater <- sum(ann.r > ann.p)
p <- min(N.greater + 1, n + 1 - N.greater) / (n +1)
p


#Modelling intensity as a function of a covariate
# Create the Poisson point process model
precip.lg.im <- as.im(precip.lg)
PPM1 <- ppm(pp1 ~ precip.lg.im)
plot(effectfun(PPM1, "precip.lg.im", se.fit=TRUE), main=NULL, las=1)
PPM1
PPM0 <- ppm(pp1 ~ 1)
PPM0
pp1$n / area(windows[[1]]) 
anova(PPM0, PPM1, test="LRT")



#ELEVATION AS COVARIATE
#plot rasterlayer histogram
elev <- r.m
hist(elev, main=NULL, las=1)
elev.lg <- log(elev)
hist(elev.lg, main=NULL, las=1)

#quadrat density on a tessallated surface for elevation raster
brk  <- c( -Inf, 4.5, 5.5, 6.5, Inf)  # Define the breaks
Zcut <- cut(elev.lg, breaks=brk)  # Classify the raster
E    <- tess(image=Zcut)  # Create a tesselated surface

plot(E, main="", las=1)
Q   <- quadratcount(pp1, tess = E)  # Tally counts
Q.d <- intensity(Q)  # Compute density
Q.d
plot(intensity(Q, image=TRUE), las=1, main=NULL)
plot(pp1, pch=20, cex=0.6, col=rgb(1,1,1,.5), add=TRUE)
cl <-  interp.colours(c("lightblue", "orange" ,"red"), te$n)
plot( intensity(Q, image=TRUE), las=1, col=cl, main=NULL)
plot(pp1, pch=20, cex=0.6, col=rgb(0,0,0,.5), add=TRUE)

#Kernal density raster
K1 <- density(pp1) # Using the default bandwidth
plot(K1, main=NULL, las=1)
contour(K1, add=TRUE)

K2 <- density(pp1, sigma=1) 
plot(K2, main=NULL, las=1)
contour(K2, add=TRUE)

K3 <- density(pp1, kernel = "disc", sigma=1) 
plot(K3, main=NULL, las=1)
contour(K3, add=TRUE)

#Kernal density adjusted for covariate
# Compute rho using the ratio method
rho <- rhohat(pp1, as.im(elev.lg),  method="ratio")
# Generate rho vs covariate plot
plot(rho, las=1, main=NULL, xlab='Logged Elevation')
pred <- predict(rho)
cl   <- interp.colours(c("lightyellow", "orange" ,"red"), 100) # Create color scheme
plot(pred, col=cl, las=1, main=NULL, gamma = 0.25)
K1_vs_pred <- pairs(K1, pred, plot = FALSE)
plot(K1_vs_pred$pred ~ K1_vs_pred$K1, pch=20,
     xlab = "Observed intensity", 
     ylab = "Predicted intensity", 
     col = rgb(0,0,0,0.1))
abline(a=0, b = 1, col = "red")
K1_vs_K1 <- pairs(K1, K1, labels = c("K1a", "K1b"), plot = FALSE)
plot(K1_vs_K1$K1a ~ K1_vs_K1$K1b, pch=20,
     xlab = "Observed intensity", 
     ylab = "Observed intensity")
summary(as.data.frame(K1_vs_pred))
plot(K1_vs_pred$pred ~ K1_vs_pred$K1, pch=20,
     xlab = "Observed intensity", 
     ylab = "Predicted intensity", 
     col = rgb(0,0,0,0.1),
     ylim = c(0, 5))
abline(a=0, b = 1, col = "red")
K1_vs_K1 <- pairs(K1, K1, labels = c("K1a", "K1b"), plot = FALSE)

#Distance based analysis
#average neighbour analysis
mean(nndist(pp1, k=1)) #To compute the average first nearest neighbor distance
mean(nndist(pp1, k=2)) #To compute the average second nearest neighbor distance
ANN <- apply(nndist(pp1, k=1:1000),2,FUN=mean)
plot(ANN ~ eval(1:1000), type="b", main=NULL, las=1)

#K and L functions
K <- Kest(pp1, correction="border")
plot(K, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))
plot(K)
kf.env <- envelope(pp1,Kest,correction="border")
plot(kf.env)

L <- Lest(pp1, main=NULL)
plot(L, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))
plot(L, . -r ~ r, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))
lf.env <- envelope(pp1,Lest,correction="border")
plot(lf.env)

#pair correlation function g
g  <- pcf(pp1)
plot(g, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))
gf.env <- envelope(pp1,Gest,correction="border")
plot(gf.env)


#Hypothesis tests
#test for clustering/dispersion
ann.p <- mean(nndist(pp1, k=1))
ann.p

n     <- 299L               # Number of simulations
ann.r <- vector(length = n) # Create an empty object to be used to store simulated ANN values
for (i in 1:n){
  rand.p   <- rpoint(n=pp1$n, win=windows[[1]])  # Generate random point locations
  ann.r[i] <- mean(nndist(rand.p, k=1))  # Tally the ANN values
}
plot(rand.p, pch=16, main=NULL, cols=rgb(0,0,0,0.5))
hist(ann.r, main=NULL, las=1, breaks=40, col="bisque", xlim=range(ann.p, ann.r))
abline(v=ann.p, col="blue")

n     <- 299L
ann.r <- vector(length=n)
for (i in 1:n){
  Z<- setcov(as.im(elev), W=windows[[1]])
  rand.p   <- rpoint(n=pp1$n, f=Z)  
  ann.r[i] <- mean(nndist(rand.p, k=1))
}

Window(rand.p) <- windows[[1]]  
plot(rand.p, pch=16, main=NULL, cols=rgb(0,0,0,0.5))
hist(ann.r, main=NULL, las=1, breaks=40, col="bisque", xlim=range(ann.p, ann.r))
abline(v=ann.p, col="blue")

N.greater <- sum(ann.r > ann.p)
p <- min(N.greater + 1, n + 1 - N.greater) / (n +1)
p


#Modelling intensity as a function of a covariate
# Create the Poisson point process model
elev.lg.im <- as.im(elev.lg)
PPM1 <- ppm(pp1 ~ elev.lg.im)
plot(effectfun(PPM1, "elev.lg.im", se.fit=TRUE), main=NULL, las=1)
PPM1
PPM0 <- ppm(pp1 ~ 1)
PPM0
pp1$n / area(windows[[1]]) 
anova(PPM0, PPM1, test="LRT")

