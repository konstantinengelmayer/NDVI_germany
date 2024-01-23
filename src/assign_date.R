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
  writeRaster(NDVI_roh[[i]], outfile[i], overwrite=TRUE)
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