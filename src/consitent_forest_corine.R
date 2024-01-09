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

                         
######################################################################################################################################################
# Den NDVI-Layern das korrekte Dataum geben
year <- substr(basename(tif_files), 35, 38)
year # Jahr extrahiert
day <- substr(basename(tif_files), 39, 41)
day # Tag extrahiert

day <- as.integer(day) # Tag als integer schreiben
datestring <- paste0(year,"-01-01") # einen Datastring daraus machen
datestring
origin <- as.Date(datestring) # als Datum formatieren
origin
newdate <- origin + (day-1) # Neue Datum mit korrekten Daten erstellen
newdate

# Write it out:
outputdir <- "YOUR_PATH"
outext <- "_correct_Date_NDVI.tif" # Korrekte Namen zuweisen
outputfilename <- substr(basename(files), 1, 26) # Filenamen erstellen
outputfilename
outfile <- paste0(outputdir,outputfilename,newdate,outext) # Filenamen zusammensetzen
outfile

# Write raster
for (i in 1:length(outfile)) {
  print(i)
  print(outfile[i])
#  writeRaster(NDVI_roh[[i]], outfile[i], overwrite=TRUE)
}

# TODO: Reading the files with correct Date
#TODO Geting all month with winter:

input <- "C:/Users/Malte/Desktop/NDVI_2023/Data/NDVI_Date"

files_DJF <- list.files(path=input, pattern = "-12-|-01-|-02-.+.tif$", full.names=T)
files_DJF
year <- as.numeric(substr(basename(files_DJF), 27, 30))
year # die mehrfachen Jahre in einen numerischen Vektor packen
outputdir <- "C:/Users/Malte/Desktop/NDVI_2023/Data/NDVI_Date/"
outputfilename <- "MOD13Q1.061__250m_16_days_NDVI_mean_Winter_"
years <- c(2000:2023) # die Jahre 2000 bis 2023 in einen Vektor erstellen
years

# Loop for: NDVI_Winter
for (i in years) { #
  print(i)
  aggrfiles <- subset(files_DJF,year==i)
  print(aggrfiles)
  rstack_mean <- calc((stack(aggrfiles)), mean, na.rm=TRUE)
  yearextension <- as.character(i)
  print(yearextension)
  outfile <- paste0(outputdir,outputfilename,yearextension) 
  print(outfile)
#  writeRaster(rstack_mean, outfile, overwrite=TRUE)
}

winter <- list.files(path=input, pattern="Winter", full.names=T)
winter <- stack(winter)
winter
plot(winter[[2]])



# TODO Merging all month with spring:

files_MAM <- list.files(path=input, pattern = "-03-|-04-|-05-.+.tif$", full.names=T)
files_MAM
year <- as.numeric(substr(basename(files_MAM), 27, 30))
year # die mehrfachen Jahre in einen numerischen Vektor packen
outputfilename <- "MOD13Q1.061__250m_16_days_NDVI_mean_Spring_"
years <- c(2000:2023) # die Jahre 2000 bis 2023 in einen Vektor erstellen
years

# Loop for: NDVI_spring
for (i in years) { #
  print(i)
  aggrfiles <- subset(files_MAM,year==i)
  print(aggrfiles)
  rstack_mean <- calc((stack(aggrfiles)), mean, na.rm=TRUE)
  yearextension <- as.character(i)
  print(yearextension)
  outfile <- paste0(outputdir,outputfilename,yearextension) 
  print(outfile)
  writeRaster(rstack_mean, outfile, overwrite=TRUE)
}

spring <- list.files(path=input, pattern="Spring", full.names=T)
spring <- stack(spring)
spring
plot(spring[[2]])



# TODO Merging all month for Summer:

files_JJA <- list.files(path=input, pattern = "-06-|-07-|-08-.+.tif$", full.names=T)
files_JJA
year <- as.numeric(substr(basename(files_JJA), 27, 30))
year # die mehrfachen Jahre in einen numerischen Vektor packen
outputfilename <- "MOD13Q1.061__250m_16_days_NDVI_mean_Summer_"
years <- c(2000:2023) # die Jahre 2000 bis 2023 in einen Vektor erstellen
years

# Loop for: NDVI_summer
for (i in years) { #
  print(i)
  aggrfiles <- subset(files_JJA,year==i)
  print(aggrfiles)
  rstack_mean <- calc((stack(aggrfiles)), mean, na.rm=TRUE)
  yearextension <- as.character(i)
  print(yearextension)
  outfile <- paste0(outputdir,outputfilename,yearextension) 
  print(outfile)
  writeRaster(rstack_mean, outfile, overwrite=TRUE)
}

summer <- list.files(path=input, pattern="Summer", full.names=T)
summer <- stack(Summer)
summer
plot(summer[[39]])


                         
# TODO Merging all month for Autumn:

files_SON <- list.files(path=input, pattern = "-09-|-10-|-11-.+.tif$", full.names=T)
files_SON
year <- as.numeric(substr(basename(files_SON), 27, 30))
year # die mehrfachen Jahre in einen numerischen Vektor packen
outputfilename <- "MOD13Q1.061__250m_16_days_NDVI_mean_Autumn_"
years <- c(2000:2023) # die Jahre 2000 bis 2023 in einen Vektor erstellen
years

# Loop for: NDVI_autumn
for (i in years) { #
  print(i)
  aggrfiles <- subset(files_SON,year==i)
  print(aggrfiles)
  rstack_mean <- calc((stack(aggrfiles)), mean, na.rm=TRUE)
  yearextension <- as.character(i)
  print(yearextension)
  outfile <- paste0(outputdir,outputfilename,yearextension) 
  print(outfile)
  writeRaster(rstack_mean, outfile, overwrite=TRUE)
}

autumn <- list.files(path=input, pattern="Autumn", full.names=T)
autumn <- stack(autumn)
autumn
plot(autumn[[39]])
