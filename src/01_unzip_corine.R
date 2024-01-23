# R script for extracting nested zip files in corine landcover folder
# The script first unzips the outer zip file (EEA.zip) to a temporary folder.
# Then, it iterates through all the inner zip files contained within this folder
# and extracts each to a specified destination folder.

# Load required setup from an external R script
source("~/edu/NDVI_germany/src/NDVI_germany_setup.R")

# Path to the outer zip file
outer_zip_file <- paste0(envrmt$path_data_level0, "/corine_landcover_2000_2018/EEA.zip")

# Temporary folder to store inner zip files
temp_folder <- envrmt$path_tmp

# Destination folder for the final extracted files
dest_folder <- paste0(envrmt$path_data_level0, "/corine_landcover_2000_2018")

# Unzip the outer zip file to the temporary folder
unzip(outer_zip_file, exdir = temp_folder)

# List all inner zip files within the temporary folder
zip_files <- list.files(temp_folder, pattern = "\\.zip$", full.names = TRUE)

# Iterate over each inner zip file and unzip it to the destination folder
for (file in zip_files) {
  unzip(file, exdir = dest_folder)
}
