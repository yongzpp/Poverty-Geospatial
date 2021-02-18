# BT4015 Project 
# Temporal Analysis

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
library(readxl)

# Load datasets
kiva <- read.csv("kivaData_augmented.csv")
MPI <- read.csv("MPIData_augmented.csv")
MPI_old <- read.csv("MPI_old.csv")

# Descriptive analysis
# summary(kiva)

# 1st dataset: kivaData_augmented
# Extracting the world region from country
# Creating a new column called "world_region"
kiva$world_region <- countrycode(sourcevar = kiva[,"country"], origin = "country.name", 
                                 destination = "continent")

# Filtering points in Africa continent
kiva_africa <- kiva[kiva$world_region=='Africa',]
# Removing data points from Egypt (not Sub-Saharan Africa)
kiva_africa <- kiva_africa[!(kiva_africa$country=='Egypt'),]

# Extracting the year from the 'funded time'
# Timeframe: 2014-2017
kiva_africa$year_funded <- substr(kiva_africa$funded_time, start = 1, stop = 4)

# Conduct before-after analysis
# Before: 2014-2015 (sum)
# After: 2016-2017 (sum)
#kiva_africa_before <- kiva_africa[(kiva_africa$year_funded == '2014')|(kiva_africa$year_funded == '2015'),]
kiva_africa_before <- kiva_africa[(kiva_africa$year_funded == '2014'),]
#kiva_africa_after <- kiva_africa[(kiva_africa$year_funded == '2016')|(kiva_africa$year_funded == '2017'),]
kiva_africa_after <- kiva_africa[(kiva_africa$year_funded == '2017'),]

# Total loan amount over the years
before_loan_amount <- data.frame(aggregate(kiva_africa_before$loan_amount, by=list(country=kiva_africa_before$country), FUN=sum))
colnames(before_loan_amount) <- c("country", "total_loan_amount_before")
before_loan_amount$total_loan_amount_before <- before_loan_amount$total_loan_amount_before/1000000 # divide by 1 million
before_loan_amount$iso_a2 = countrycode(sourcevar = before_loan_amount[,"country"], origin = "country.name", 
                                             destination = "iso2c")
#before_loan_amount$hover <- paste(before_loan_amount$country, before_loan_amount$total_loan_amount_before, sep=": ")

# Plot on total loan amounts
africa = world %>% 
  filter(continent == "Africa", !is.na(iso_a2)) %>% 
  left_join(before_loan_amount, by = "iso_a2") %>% 
  dplyr::select(country, total_loan_amount_before) 
tm_shape(africa) + tm_polygons("total_loan_amount_before", title = "Total Loan Amount (Before)")

# Total loan amount over the years
after_loan_amount <- data.frame(aggregate(kiva_africa_after$loan_amount, by=list(country=kiva_africa_after$country), FUN=sum))
colnames(after_loan_amount) <- c("country", "total_loan_amount_after")
after_loan_amount$total_loan_amount_after <- after_loan_amount$total_loan_amount_after/1000000 # divide by 1 million
after_loan_amount$iso_a2 = countrycode(sourcevar = after_loan_amount[,"country"], origin = "country.name", 
                                        destination = "iso2c")
#after_loan_amount$hover <- paste(after_loan_amount$country, after_loan_amount$total_loan_amount_after, sep=": ")

# Plot on total loan amounts
africa = world %>% 
  filter(continent == "Africa", !is.na(iso_a2)) %>% 
  left_join(after_loan_amount, by = "iso_a2") %>% 
  dplyr::select(country, total_loan_amount_after) 
tm_shape(africa) + tm_polygons("total_loan_amount_after", title = "Total Loan Amount (After)", breaks = c(0,2,4,6,8))

# Table of difference between before and after
# Using an outer join
df_compare <- merge(x = before_loan_amount, y = after_loan_amount, by = "country", all = TRUE)
# Fill all NAs with zeros
df_compare$total_loan_amount_before[is.na(df_compare$total_loan_amount_before)] <- 0
df_compare$total_loan_amount_after[is.na(df_compare$total_loan_amount_after)] <- 0

# Creating a new column called 'difference'
# After - Before
df_compare$loan_difference <- df_compare$total_loan_amount_after - df_compare$total_loan_amount_before
# Getting the table of difference
df_compare[,c("country", "loan_difference")]

# Plot on difference in total loan amounts
df_compare <- df_compare %>% rename(iso_a2 = iso_a2.x)

# Filling the missing values
df_compare[df_compare$country == "Madagascar", "iso_a2"] <- "MG"
df_compare[df_compare$country == "Lesotho", "iso_a2"] <- "LS"
df_compare[df_compare$country == "Cote D'Ivoire", "iso_a2"] <- "CI"

africa = world %>% 
  filter(continent == "Africa", !is.na(iso_a2)) %>% 
  left_join(df_compare, by = "iso_a2") %>% 
  dplyr::select(country, loan_difference) 
tm_shape(africa) + tm_polygons("loan_difference", title = "Difference in Loan Amount", legend.format=list(flag= "+")) + 
  tm_layout(legend.title.size = 0.8)

# Retrieving the list of countries with loans
# Extracting the ISO country code from the country column
df_compare$iso_a3 = countrycode(sourcevar = df_compare[,"country"], origin = "country.name", 
                                        destination = "iso3c")
list_of_countries <- df_compare$iso_a3

# 2nd dataset: MPIData_augmented
# Filtering points in Sub-Saharan Africa only
MPI_africa <- MPI[MPI$World.region=='Sub-Saharan Africa',]
#MPI_africa <- MPI[MPI$ISO.country.code %in% list_of_countries,]

# Create a new column of 'ISO Country Code 2 digit'
MPI_africa$'iso_a2' = countrycode(MPI_africa$ISO.country.code, origin = "iso3c", destination = "iso2c")

# Keeping the Country and MPI_Country columns
sub_MPI_africa <- MPI_africa[,c("Country", "MPI_Country", "ISO.country.code", "iso_a2")]
sub_MPI_africa <- unique(sub_MPI_africa)
sub_MPI_africa <- sub_MPI_africa %>% rename(MPI_Country_after = MPI_Country)

# Merging the datasets
merged_MPI <- merge(x = MPI_old, y = sub_MPI_africa, by = "ISO.country.code", all.x = TRUE)
merged_MPI <- na.omit(merged_MPI)

# Since the MPI_Country is not available for every year, we need to pick the most suitable year.
# 2013, 2014, 2011, then 2015
MPI_2013 <- merged_MPI[merged_MPI$Year == '2013',]
countries <- unique(MPI_2013$ISO.country.code)

MPI_2014 <- merged_MPI[merged_MPI$Year == '2014',]
MPI_2014 <- MPI_2014[!(MPI_2014$ISO.country.code %in% countries),]
all_MPI <- rbind(MPI_2013, MPI_2014)
countries <- unique(all_MPI$ISO.country.code)

MPI_2011 <- merged_MPI[merged_MPI$Year == '2011',]
MPI_2011 <- MPI_2011[!(MPI_2011$ISO.country.code %in% countries),]
all_MPI <- rbind(all_MPI, MPI_2011)
countries <- unique(all_MPI$ISO.country.code)

MPI_2015 <- merged_MPI[merged_MPI$Year == '2015',]
MPI_2015 <- MPI_2015[!(MPI_2015$ISO.country.code %in% countries),]
all_MPI <- rbind(all_MPI, MPI_2015)

# Plot the Country MPI (before)
all_MPI <- all_MPI %>% rename(Country = Country.x)
africa = world %>% 
  filter(continent == "Africa", !is.na(iso_a2)) %>% 
  left_join(all_MPI, by = "iso_a2") %>% 
  dplyr::select(Country, MPI_Country_before) 
tm_shape(africa) + tm_polygons("MPI_Country_before", title = "Country MPI (before)")

# Plot the Country MPI (after)
africa = world %>% 
  filter(continent == "Africa", !is.na(iso_a2)) %>% 
  left_join(all_MPI, by = "iso_a2") %>% 
  dplyr::select(Country, MPI_Country_after) 
tm_shape(africa) + tm_polygons("MPI_Country_after", title = "Country MPI (after)")

# Plot the difference in MPI
all_MPI$mpi_difference <- all_MPI$MPI_Country_after - all_MPI$MPI_Country_before
africa = world %>% 
  filter(continent == "Africa", !is.na(iso_a2)) %>% 
  left_join(all_MPI, by = "iso_a2") %>% 
  dplyr::select(Country, mpi_difference) 
tm_shape(africa) + tm_polygons("mpi_difference", title = "Difference in Country MPI", legend.format=list(flag= "+"), palette = "seq") +
  tm_layout(aes.palette = list(seq = "-RdYlGn"))

# Computing a binary variable to indicate if kiva has distributed the loans effectively
total_loan_amt_diff <- df_compare[, c("iso_a3", "country", "loan_difference")]
MPI_diff <- all_MPI[, c("ISO.country.code", "mpi_difference", "iso_a2")]
MPI_diff <- MPI_diff %>% rename(iso_a3 = ISO.country.code)
df_effective <- merge(x = total_loan_amt_diff, y = MPI_diff, by = "iso_a3")

# The following are considered to be effective distribution
# 1. positive loan_difference, negative mpi_difference (becomes better)
# 2. negative loan_difference, negative mpi_difference (becomes better)
df_effective$effective <- with(df_effective, ifelse(mpi_difference < 0, "Improved", 
                                                    ifelse(loan_difference < 0, "Worsened (reduced loan)",
                                                           "Worsened (increased loan)")))

africa = world %>% 
  filter(continent == "Africa", !is.na(iso_a2)) %>% 
  left_join(df_effective, by = "iso_a2") %>% 
  dplyr::select(country, effective) 
tm_shape(africa) + tm_polygons("effective", title = "Effectiveness of loans", style = "cat", palette = c("green", "red", "blue")) +
  tm_layout(legend.width = 1.1)











