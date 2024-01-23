# try with dataframe



# Load necessary libraries
library(raster)
library(zoo)

# Assuming 'ndvi_data' is a RasterStack or RasterBrick of NDVI data over time
# And 'dates' is a vector of corresponding dates

# Convert to a time series data frame
ndvi_ts <- as.data.frame(getValues(ndvi_data))
colnames(ndvi_ts) <- dates

# Apply a rolling mean to smooth the data (optional, can adjust window size)
ndvi_smooth <- rollapply(ndvi_ts, width = 3, FUN = mean, by.column = TRUE, fill = NA)

# Define a function to find SOS and EOS based on a threshold
find_season <- function(ndvi_series, threshold) {
  above_thresh <- ndvi_series > threshold
  sos <- which(diff(above_thresh) == 1)[1]
  eos <- which(diff(above_thresh) == -1)[1]
  return(c(sos, eos))
}

# Apply the function to your NDVI time series
season_indices <- apply(ndvi_smooth, 2, find_season, threshold = 0.5)

# Extracting SOS and EOS dates
sos_dates <- dates[season_indices[1,]]
eos_dates <- dates[season_indices[2,]]
