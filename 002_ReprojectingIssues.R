#### Background info ----
# Managing raster data: https://desktop.arcgis.com/en/arcmap/10.3/manage-data/raster-and-images/raster-coordinate-systems.htm
# Datums: https://gisgeography.com/geodetic-datums-nad27-nad83-wgs84/
# Datums vs Reference Systems: https://gis.stackexchange.com/questions/664/difference-between-projection-and-datum
# Resampling methods: https://gis.stackexchange.com/questions/2587/deciding-what-interpolation-method-to-use-for-resampling-raster-data

#### Load packages, set seed ----
## Mandatory package
library(raster)

## If you're not plotting the context polygons you don't need to import the following
library(USAboundaries)
library(sf)
library(ggplot2)

## Set seed for replicability
set.seed(-356) # The start of the conquests of Alexander the Great

#### Generate some data ----
# A made-up raster dataset located in northeast Arizona
madeUpData <- matrix(runif(10000,0,100),
                     nrow=100,
                     ncol=100)

# Set up extent information of bounding box
xmn = 300000
xmx = 400000
ymn = 4000000
ymx = 4100000

# Start out with a common coordinate reference system: 
# The Universal Transverse Mercator (units are in meters)
# Zone 13 runs up through Arizona and New Mexico
AZ_UTM = "+proj=utm +zone=13 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# Compile the raster using the data, projection information, and extent information
baseRaster <- raster(madeUpData,
                     crs = AZ_UTM,
                     xmn = xmn, xmx = xmx,
                     ymn = ymn, ymx = ymx)
str(baseRaster)
plot(baseRaster)


#### Context plots ----

# You won't be able to run this without the miscellaneous libraries above
# Access the 4 corners states using the USAboundaries library
myStates <- us_states(resolution = "low", 
                      states = c("Arizona",
                                 "New Mexico",
                                 "Utah",
                                 "Colorado"))
str(myStates)

# Spatially transform the states lines to the baseRaster crs
myStates <- st_transform(myStates,crs(baseRaster))

# Look at where the made up data lives
ggplot() +
  geom_sf(data = myStates) +
  geom_rect(aes(xmin=extent(baseRaster)@xmin,
                xmax=extent(baseRaster)@xmax,
                ymin=extent(baseRaster)@ymin,
                ymax=extent(baseRaster)@ymax))



#### project Ad Nauseum ----

# List of potential CRS systems to project to
# Each of these projections were found on epsg.io by searching "Arizona"
mySystems <- c("+init=epsg:2222",
               "+init=epsg:2223",
               "+init=epsg:2224",
               "+init=epsg:2761",
               "+init=epsg:2762",
               "+init=epsg:2763",
               "+init=epsg:2867",
               "+init=epsg:2868",
               "+init=epsg:2869")

# Make a new raster based on the one we made above
newRaster <- baseRaster

# Use a `for` loop to iterate over projections
for(i in 1:length(mySystems)){ 
  # For each of the projections listed above...
  myCRS <- mySystems[i]
  # ... project the last raster to the new projection
  newRaster <- projectRaster(newRaster,
                             crs = crs(myCRS),
                             method = "bilinear")
}

# Convert newRaster back to baseRaster projection
newRaster <- projectRaster(newRaster,
                           crs = crs(baseRaster),
                           method = "bilinear")



#### Comparing inputs and outputs ----

# Check if all cells are equal
all.equal(newRaster, baseRaster) # Nope!

# Plot the two raster sets side-by-side:
par(mfrow=c(1,2))
plot(baseRaster,
     xlim = c(xmn,xmx),
     ylim = c(ymn,ymx),
     main = "Original raster")
plot(newRaster,
     xlim = c(xmn,xmx),
     ylim = c(ymn,ymx),
     main = "New raster")
# Do they look the same?
# How much change do we see?

# Compare the spatial resolutions:
res(baseRaster)
res(newRaster)

# Plot a histogram:
hist(baseRaster, 
     breaks = 2500, 
     main = "Base raster")
hist(newRaster, 
     breaks = 2500, 
     main = "New raster")
# What happened?
# Hint: https://en.wikipedia.org/wiki/Regression_toward_the_mean
# Think about what bilinear interpolation is doing

# Can you imagine a real-life scenario where this would affect your work?
