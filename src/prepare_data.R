# R script to identify and visualize pixels with a consistent value of 23 (23 = broad-leaved forest)
# across a raster stack of .tif files.
# The script reads .tif files from a specified folder into a raster stack,
# checks each pixel in each layer for the value 23,
# and then identifies those pixels that consistently have the value 23 across all layers.

# Load the setup script for the NDVI Germany project
source("~/edu/NDVI_germany/src/NDVI_germany_setup.R")

# Define the folder path where .tif files are stored
folder_path <- paste0(envrmt$path_data_level0, "/corine_landcover_2000_2018")

# List all .tif files in the specified folder
tif_files <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)

# Load the .tif files into a raster stack
corine_landcover_stack <- rast(tif_files)

# Apply a function to each layer in the stack to check if the pixel values are equal to 23
layers_with_23 <- lapply(corine_landcover_stack, function(x) x == 23)

# Reduce the list of logical rasters to a single raster, where a pixel is TRUE only if it is TRUE in all layers
consistent_23 <- Reduce("&", layers_with_23)

# Write the raster 'consistent_23' to a .tif file, allowing overwriting of existing file if it exists
terra::writeRaster(consistent_23, paste0(envrmt$path_data_level1, "/consitent_broad_leaved_forest.tif"), overwrite = TRUE)


