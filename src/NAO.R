# NAO

# Get the whole dataframe
df <- read.table("~/edu/NDVI_germany/data/NAO/norm.nao.monthly.b5001.current.ascii") 
head(df)




# Source: https://www.ufz.de/index.php?de=37937
install.packages("ncdf4")
library("ncdf4")
install.packages("RColorBrewer")
library(RColorBrewer)


drought <- nc_open("~/edu/NDVI_germany/data/NAO/271359_Dürreintensität_Gesamtboden_1952-2022_Apr-Okt.nc")
print(drought)

# Load the NetCDF file
nc_file_path <- "~/edu/NDVI_germany/data/NAO/271357_Dürremagnitude_Gesamtboden_1952-2022_Apr-Okt.nc"
smi_raster <- rast(nc_file_path, subds="SMI")
print(smi_raster)
# Replace -9999 with NA
values(smi_raster)[values(smi_raster) == -9999] <- NA

# Define a custom color palette (inverted)
color_palette <- (heat.colors(100)) # Using terrain.colors for example, but you can choose any palette
# Create a green-to-red color palette
green_to_red <- colorRampPalette(rev(c("green", "yellow", "red")))

# Generate a palette with a specified number of colors, for example, 100
palette <- green_to_red(100)

# Function to add the year to the plot title
plot_with_year <- function(raster_layer, timestep, start_year=1952) {
  # Calculate the year for the given timestep
  year <- start_year + timestep - 1
  
  # Plot the raster layer with the custom color palette
  plot(raster_layer, main=paste("Soil Moisture Index (SMI) - Year:", year), col=palette)
}

# Example: Plot the first time step with the year in the title
plot_with_year(smi_raster[[52]], 52)

