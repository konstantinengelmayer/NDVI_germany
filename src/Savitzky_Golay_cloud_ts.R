source("~/edu/NDVI_germany/src/NDVI_germany_setup.R")
df <- readRDS("~/edu/NDVI_germany/data/raster_data/data_level1/long_ndvi_quality_df.rds")
head(df)
nrow(df) # 18600582


# Define the conversion function
to_binary_string <- function(number, bits = 16) {
  if (is.na(number)) {
    return(NA_character_)  # Handle NA values
  } else {
    bit_vector <- rev(as.integer(intToBits(number)))
    bit_string <- paste(bit_vector[(33-bits):32], collapse = "")
    return(bit_string)
  }
}

# Apply the function to the 'Quality' column of your dataframe
df$Quality_Binary = sapply(df$Quality, to_binary_string, bits = 15)
head(df)
tail(df)
nrow(df)
str(df)


# Filter the good values
result_df <- df %>%
  filter(
    str_sub(Quality_Binary, 1, 2) != "01" &
      str_sub(Quality_Binary, 1, 2) != "11" & 
      str_sub(Quality_Binary, 3, 4) != "11" &
      str_sub(Quality_Binary, 3, 4) != "1" &
      str_sub(Quality_Binary, 15, 15) == "0" 
  )


tail(result_df)
nrow(result_df) # Old: 11874448 / New: 14773579
14773579-11874448 # 2899131 cloud pixel integrated
18600582-14773579 # 3827003 Ovservations less
18600582-11874448 # 6726134 Observations less

# result_df$Quality_Binary <- NULL # Deleting the Quality bits
# result_df$Quality <- NULL # Deleting the coloumn with the Quality

print(colnames(result_df)) # Get the names = we should rename the Year coloumn to the Date coloumn.
colnames(result_df)[1] <- "Date" # This makes the next steps more intuitive: Date and Year coloumn
result_df$Year <- year(result_df$Date)
head(result_df) # Have a look
str(result_df)


# TODO: Filter pixel groups with less than 3 observations per season
# Convert Date to Date 
result_df$Date <- as.Date(result_df$Date)
str(result_df)

# Function to determine the season based on the month
getSeason <- function(month) {
  if (month %in% 3:5) {
    return("Spring")
  } else if (month %in% 6:8) {
    return("Summer")
  } else if (month %in% 9:11) {
    return("Autumn")
  } else {
    return("Winter")
  }
}

# Apply the function to create a new Season column
result_df$Season <- sapply(month(result_df$Date), getSeason)

# Group by Year, x, y, and Season, then filter groups with at least 3 observations
result <- result_df %>%
  group_by(Date = year(Date), x, y, Season) %>%
  filter(n() >= 3) %>%
  ungroup()

# View the first few rows of the result
head(result)
nrow(result)
14773579-14471973 # 301606 filtered out because of too less observations!


############# CAUTION: IS 8 REALLY THE RIGHT BIT ???? ##########################

# TODO: Detect clouds and filter pixel with more than 3 cloud observations together. 
# TODO: Detect the clouds first: 
detect_clouds <- function(quality_binary) {
  if (str_sub(quality_binary, 8, 8) == "0" & 
      str_sub(quality_binary, 10, 10) == "0" &
      str_sub(quality_binary, 14, 14) == "0") {
    return(0) # No clouds
  } else {
    return(1) # Clouds present
  }
}

# Apply the function to create a new Clouds column
result$Clouds <- sapply(result$Quality_Binary, detect_clouds)




# TODO: Exclude all pixel with more than X cloudy observations. 
# Function to check for consecutive cloudy observations
check_consecutive_clouds <- function(clouds) {
  rle_clouds <- rle(clouds)
  if(any(rle_clouds$values == 1 & rle_clouds$lengths > 4)) {
    return(TRUE)  # Indicates exclusion due to more than 4 consecutive cloudy observations
  } else {
    return(FALSE) # Indicates keeping the pixel group
  }
}

# Apply the function to each pixel group and filter out those with more than 4 consecutive clouds
df_filtered <- result %>%
  group_by(x, y) %>%
  mutate(Exclude = check_consecutive_clouds(Clouds)) %>%
  ungroup() %>%
  filter(Exclude == FALSE) %>%
  select(-Exclude)  # Remove the helper column if not needed

# View the result
head(df_filtered)



# Safe and read the new file with all good values
saveRDS(df_filtered, file="~/edu/NDVI_germany/data/raster_data/data_level1/long_ndvi_clouds.rds")



################################################################################
######################## Start of calculation ##################################
################################################################################
# TODO: Load the dataframe
# Note: The smoother value affects the result of SOS and EOS!
df <- readRDS("~/edu/NDVI_germany/data/raster_data/data_level1/long_ndvi_clouds.rds")# Read the RDS file
head(df)


