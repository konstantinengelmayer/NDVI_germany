source("~/edu/NDVI_germany/src/NDVI_germany_setup.R")
quality_raster <- rast(file.path(envrmt$path_quality_layer, "MOD13Q1.061__250m_16_days_VI_Quality_doy2013193_aid0001.tif"))

to_binary_string <- function(number, bits) {
  if (is.na(number)) {
    return(NA_character_)  # Return NA for NA inputs
  } else {
    bit_vector <- rev(as.integer(intToBits(number)))
    bit_string <- paste(bit_vector[(32-bits+1):32], collapse = "")
    return(bit_string)
  }
}

# Convert raster values to binary with enough bits (16 bits for MODIS), handling NA values
quality_raster_binary <- sapply(values(quality_raster), to_binary_string, bits = 15)

# Function to check if the first 6 bits are 0, handling NA values
check_first_6_bits <- function(binary_val) {
  if (is.na(binary_val)) {
    return(NA)  # Preserve NA values
  } else {
    return(substr(binary_val, 1, 3) == "000")
  }
}

# Apply the function
selected_pixels_values <- sapply(quality_raster_binary, check_first_6_bits)


# Convert back to raster
selected_pixels <- rast(matrix(selected_pixels_values, nrow=nrow(quality_raster), ncol=ncol(quality_raster), byrow = TRUE))

plot(selected_pixels)
# selected_pixels is now a logical raster where TRUE indicates pixels with the first 6 bits as 0
