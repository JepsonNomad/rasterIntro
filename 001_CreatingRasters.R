## Let's create a raster object from scratch.

#### Load packages ----
library(raster)




#### Define parameters ----

# Let's make up some data here. A raster at its core is a grid of cells, so let's make a 10x10 matrix containing random data. That means our raster will also ultimately have dimensions c(10,10). Let's call this first dataset forest cover.
forestData = matrix(rnorm(100),
                nrow = 10,
                ncol = 10)

# For geospatial work, we need a projection to work in. Let's pick UTM zone 11N (west half of CA)
UTM11N = CRS("+init=epsg:32611")




#### Create raster object ----
# First, create a raster object using only the data we created.
for0 = raster(forestData)
# Next, create a raster object using the data, as well as projection information
for1 = raster(forestData,
            crs = UTM11N)
# Finally, create a raster object using the data, projection information, and parameters defining the bounding box.
for2 = raster(forestData,
            xmn = 0, xmx = 1000,
            ymn = 0, ymx = 1000,
            crs = UTM11N)

## Compare the various resulting raster objects. What is different between them?
for0 # NA CRS, but resolution and extent are auto-generated
for1 # UTM with auto-generated resolution and extent
for2 # UTM with defined extent and auto-calculated resolution

plot(for0)
plot(for1)
plot(for2)

#### Create raster stack ----
## To stack several raster layers on top of each other, use the `stack()` function. 
# Try stacking each of the rasters we created on top of each other.
s0 = stack(for0,for1)
s1 = stack(for0,for2)
s2 = stack(for1,for2)
# This fails; why?

## Remember, the point of raster stacks is to systematically connect data in a *meaningful* way.
# Stacks should only be generated when the geographic parameters of each component layer match. Let's check if these forest cover rasters are comparable.
compareRaster(for0,for1)
compareRaster(for0,for2)
compareRaster(for1,for2)

# Make a new dataset: elevation.
elevData <- matrix(rnorm(100),
                   nrow = 10,
                   ncol = 10)

# Convert it into a raster object that could be stacked with for2 - the most precisely defined object
dem <- raster(elevData,
              xmn = 0, xmx = 1000,
              ymn = 0, ymx = 1000,
              crs = UTM11N)

