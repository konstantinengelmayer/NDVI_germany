library(envimaR)

packagesToLoad <- c(
  "terra",
  "reticulate",
  "sf",
  "osmdata",
  "rsample",
  "tfdatasets",
  "purrr",
  "stars",
  "magick",
  "mapview",
  "lubridate",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2"
)

# define a project root folder
rootDir <- "~/edu/NDVI_germany"

# some new paths
projectDirList <- c(
  "data/",
  "data/raster_data/data_level0",
  "data/raster_data/data_level0/NDVI_data",
  "data/raster_data/data_level0/NDVI_data/ndvi",
  "data/raster_data/data_level0/NDVI_data/quality_layer",
  "data/raster_data/data_level1",
  "data/vector_data",
  "docs/",
  "run/",
  "tmp",
  "src/",
  "src/functions/"
)

# Now set automatically root direcory, folder structure and load libraries
envrmt <- envimaR::createEnvi(
  root_folder = rootDir,
  folders = projectDirList,
  path_prefix = "path_",
  libs = packagesToLoad,
  alt_env_id = "COMPUTERNAME",
  alt_env_value = "PCRZP",
  alt_env_root_folder = "F:/BEN/edu"
)

## set terra temp path
terra::terraOptions(tempdir = envrmt$path_tmp)