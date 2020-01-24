#### Introduction to raster data ----
## Goals: 
# Understand the difference between raster and vector data
# Introduce raster parameters - dim(), crs(), res(), extent(), names()
# Use basic raster operations - subset(), stack(), calc(), extract(), crop(), mask(), cellStats()
# Visualization of rasters - plot(), `col = `, `xlim = `, `ylim = `
# Access and download data - getData()

## If time permits:
# Data access: discussion of other sources? 
# - EarthEngine, NASA, CA Geospatial portal


#### Load packages, set working directory ----
install.packages(c("raster","rgdal", "here"))

library(raster) # For working with raster datasets
library(rgdal)  # For working with vector datsets


#### Acquire data ----
# Let's start off with two sets of data. First a vector object with information about California's shape, and second, some historical temperature data in raster format.

## Vector:
# https://data.ca.gov/dataset/ca-geographic-boundaries/resource/3db1e426-fb51-44f5-82d5-a54d7c6e188b
# Download the shapefile from the above source.
CA = readOGR(dsn = "/PATH/TO/DOWNLOADS/ca-state-boundary/",
             layer = "CA_State_TIGER2016")
# Always look at your data!
str(CA)
plot(CA)
# We'll come back to the California shapefile later.


## Raster:
?getData

# Now, let's work with some raster data First, download historical (also called "current" = 1960-1990) maximum monthly temperatures. Each layer of the raster stack is the maximum temperature for a month (nLayers = 12). Note that temperature is stored in units of 0.1 degrees Celcius. The resolution argument (`res`) tells `getData()` that we want coarse data (resolution = 10 minutes means that a pixel in Sacramento is ~11.5 miles wide). 
maxTempHistorical = getData(name = 'worldclim', 
                         var='tmax', 
                         res=10)
maxTempHistorical



#### Elements of a raster* object ----

# Raster* objects include rasters, rasterStacks, and rasterBricks. A raster has one layer, while rasterStacks and rasterBricks can have more than one layer.
maxTempHistorical
# How many rows and columns (and optionally layers) are in the object?
dim(maxTempHistorical)
# How big are the pixels?
res(maxTempHistorical) 
# Where are the geographic limits of this dataset?
extent(maxTempHistorical) 
# How are these data being projected from a 3D surface onto a 2D map?
crs(maxTempHistorical)
# What are the names of the layers in the object?
names(maxTempHistorical) 

## An important note on geospatial data:
# Calling `res(maxTempHistorical)` returns 0.16666... This is because we downloaded worldclim data at a 10-minute scale. There are 60 minutes in a degree; 10/60 = 1/6 = 0.16666... The resolution that is reported for a raster depends on its crs!!! So, if you have an object with crs() = "...longlat..." you're measuring your data in degrees. Another common crs is UTM (California is in UTM zones 10 and 11), which is measured in meters! Comparing the resolution of data in longlat and UTM units is like comparing a space shuttle's speed in miles per hour and kilometers per hour.




#### Plotting ----
# Using typical graphics functions in R we can visualize raster datasets.
# Note that when we call `plot()` on a raster dataset with many layers, each layer is plotted.
plot(maxTempHistorical)



#### Subsetting layers ----
# We can subset a layer of a raster stack in two ways:
# 1, using the `subset()` function
# 2, using indexing

# Method 1: Can use layer name or numerical index.
?raster::subset
names(maxTempHistorical)
maxJanHisto = subset(maxTempHistorical, 
                     subset = "tmax1")

# Method 2: Also works for layer name or numerical index
sample2 = maxTempHistorical[[1]]
sample3 = maxTempHistorical[["tmax1"]]

# Check that both methods do the same thing:
all.equal(maxJanHisto, sample2) # TRUE
all.equal(sample2, sample3) # TRUE

# Let's look at what we've got now:
plot(maxJanHisto, 
     main = "Historical max Jan temps")



#### Plotting (reprise) ----
# We can specify a color ramp using the `col` argument:
plot(maxJanHisto, 
     col = topo.colors(100), 
     main = "Historical max Jan temps")

# We can zoom in on a region using the xlim and ylim arguments:
plot(maxJanHisto,
     xlim = c(-180,-50),
     ylim = c(10,90), 
     main = "Historical max Jan temps, North America")



#### Subsetting in space ----

# Often, we only want to look at a specific region within a raster dataset. This can be accomplished by cropping and/or masking, and is akin to spatially subsetting our data. The raster functions for this are `crop()` and `mask()`. To crop or mask a raster, we will need a raster and a secondary object, typically a Spatial* object. The raster and secondary objects must be projected using the same crs arguments. Now we'll return to the California shapefile we downloaded a while ago.
str(CA) # CA is a SpatialPolygonsDataFrame
proj4string(CA)
proj4string(maxJanHisto)

# Since the projections aren't the same, we need to transform a dataset. Reprojecting rasters can degrade the data quality (and is computationally expensive), so let's transform the California shapefile
CA_proj = spTransform(CA, CRSobj = crs(maxJanHisto))

# Check the transformation using plot, appending the last plot
lines(CA_proj)
# Does that look right?

# Now, we can finally do some subsetting! Let's compare cropping and masking.
# First, crop:
maxJanCaliCrop = crop(maxJanHisto,
                      CA_proj)
# Next, mask:
maxJanCaliMask = mask(maxJanHisto,
                      CA_proj)

# Visualize the output of these functions
plot(maxJanCaliCrop, main = "cropped dataset")
plot(maxJanCaliMask, main = "masked dataset")

# So, what's the difference?
extent(maxJanCaliCrop)
extent(maxJanCaliMask)
# Remember we can always check the documentation
?mask
?crop

# If we want to crop to the extent, *and* mask to the bounds of the region, we need to use both `crop()` and `mask()`:
maxJanCaliClip = crop(maxJanHisto,
                      CA_proj)
maxJanCaliClip = mask(maxJanCaliClip,
                      CA_proj)
plot(maxJanCaliClip)
lines(CA_proj)



#### Raster calculations ----
# Now, let's use raster calculations to compute values within and across rasters.

## A "within" example:

# Suppose we want to visualize Jan temps with meaningful values. Remember, temperature is reported in 0.1-degree Celcius units. So, let's divide by 10 and visualize:
janRealUnits = calc(maxJanHisto,
                    fun = function(x){x/10})
plot(maxJanHisto,
     main = "Raw data (units = 0.1 deg C)")
plot(janRealUnits,
     main = "Scaled data (units = 1 deg C)")
# Values in the sampleRasCalc raster are double that of the maxJanHistorical.
# This was the same as the following vectorized approach:
plot(maxJanHisto/10,
     main = "Vectorized method (units = 1 deg C)")


## An "across" example:

# Suppose we want to estimate the mean January temperature across space. To do that, we'll also need a tmin raster layer. The `var` argument can be set to `tmin` for minimum monthly temperatures. While it isn't perfect, a rough estimate for mean temperature can be computed by taking the mean of the tmax and tmin values.
minTempHisto = getData(name = 'worldclim', 
                            var='tmin', 
                            res=10)
# Select the layer with January
minJanHisto = minTempHisto[[1]]
# Visualize - do these data make sense?
plot(minJanHisto)

# Next, stack the raster datsets we want to use into one object using the `stack()` function:
janHistoStack = stack(maxJanHisto,
                      minJanHisto)

# Finally, we can calculate a mean value using the `calc()` function, and specifying that we want to calculate a mean:
meanJanHisto = calc(janHistoStack,
                    fun = mean)
# Look at the data
plot(meanJanHisto)



#### A practical application ----

## Let's make things a bit more complex. At this point, we can answer the following question: how much warming can we expect by the year 2070 under the current Carbon emissions regime?
# To answer: Download projected maximum and minimum monthly temperatures under a climate warming scenario. We'll use the RCP 8.5 ("business as usual") in the year 2070. Predicted climate data came from the Coupled Model Intercomparisin Project (CMIP5) - here we're just using one model from that project (AC). We'll use the same resolution as the other data we downloaded, 10-minutes.
tmax85 = getData(name = 'CMIP5', 
                 var='tmax', 
                 res=10, 
                 rcp=85, 
                 model='AC', 
                 year=70)
tmin85 = getData(name = 'CMIP5', 
                 var='tmin', 
                 res=10, 
                 rcp=85, 
                 model='AC', 
                 year=70)
tmax85
tmin85

## Let's see how much warmer we can expect mean January temperatures to be in 2070, under the RCP 8.5 scenario. 
# First, stack the layers of interest into a single object
jan85 = stack(tmax85[[1]],
              tmin85[[1]])
# Use `calc()` to find the mean
meanJan85 = calc(jan85,
                 fun = mean)

# An alternative to `calc()` is to use basic arithmetic directly on two raster objects. Let's find the difference between predicted and historical January temperatures using an arithmetical approach.
diffJan85 = meanJan85 - meanJanHisto

# Let's look at the difference:
plot(diffJan85, main = "Change in Jan temp, historical to 2070")

# Where do we see the greatest warming? How much warmer do we expect those regions to be compared to the 1960-1990 baseline?



#### Extracting values ----

## Suppose we want to know what the situation is for a particular location - say, Sacramento. To do that, we'll use the `extract()` function:
?extract

# First, identify where Sacramento is 
# NOTE: be sure to use appropriate units given the CRS arguments!!
sacLon = -121.4944  # Longitude = east/west
sacLat = 38.5816    # Latitude = north/south

# We need to convert the coordinates for Sacramento into matrix form to use the `extract()` function
sacMatrix =  matrix(c(sacLon, sacLat),
                    ncol = 2)

# Use the `extract()` function to identify the value at a specific location
sacWarming = extract(diffJan85, y = sacMatrix)
sacWarming

# How much is the mean January temperature in Sacramento expected to increase by 2070?


#### Raster statistics ----

## We can get broad-scale aggregated summary statistics from rasters using the `cellStats()` function and specifying a summarizing function. Available functions include mean, standard deviation, min, etc. See documentation details for more.
?cellStats
meanTempJan = cellStats(meanJanHisto,
                        stat = "mean")
meanTempJan

# We can also use cellStats for raster bricks, returning a value for each layer. Let's use California temperature as an example. First, crop and mask the full historical max temp raster stack to CA.
maxHistoCaliClip = crop(maxTempHistorical,
                        CA_proj)
maxHistoCaliClip = mask(maxHistoCaliClip,
                        CA_proj)

# Find the mean value of each layer
avgMaxtempCA = cellStats(maxHistoCaliClip,
                         stat = "mean")
# Visualize
plot(avgMaxtempCA/10, 
     type = "b",
     main = "Average maximum temperature by month\nin California, 1960-1990",
     ylab = "Temperature (deg C)",
     xlab = "Month")



#### Challenge ---- 

## Compare predicted temperature in a month of your choice between the two different RCP scenarios, RCP 4.5 and RCP 8.5. How does Carbon mitigation affect expected warming? Does mitigation affect warming uniformly over space?

# Hint: First, develop a workflow with comments. What data do you need to answer the question? What operations need to be performed on the data? Are there any unit corrections to consider?



#### Additional resources ----

## Coordinate reference system and spatial projection
# https://www.earthdatascience.org/courses/earth-analytics/spatial-data-r/intro-to-coordinate-reference-systems/

## What is raster data?
# https://desktop.arcgis.com/en/arcmap/10.3/manage-data/raster-and-images/what-is-raster-data.htm

## Raster package documentation 
# https://cran.r-project.org/web/packages/raster/raster.pdf

