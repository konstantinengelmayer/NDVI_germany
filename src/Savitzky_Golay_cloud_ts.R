source("~/edu/NDVI_germany/src/NDVI_germany_setup.R")
df <- readRDS("~/edu/NDVI_germany/data/raster_data/data_level1/long_ndvi_quality_df.rds")
head(df)
nrow(df) # 18600582
df <- subset(df, Year < 2001)
tail(df)
head(df)
str(df)


df$Year <- as.Date(df$Year) # Converting 'Year' to Date format if it's not already

# Group by pixel coordinates and slice the first few observations per group
# Assuming each pixel group can be uniquely identified by 'x' and 'y' coordinates
grouped_df <- df %>%
  group_by(x, y) %>%
  slice_head(n = 1) %>%
  ungroup()

# Now, if you want to select the first 4 unique pixel groups based on 'x' and 'y'
unique_pixels <- distinct(df, x, y) %>% slice_tail(n = 1)

# Filter based on matching 'x' and 'y' to those in 'unique_pixels'
df <- df %>%
  semi_join(unique_pixels, by = c("x", "y"))

# Show the structure and head of the extracted data
str(df)
head(df)
backup <- df

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

# filter the values
result_df <- df %>%
  dplyr::filter(
    str_sub(Quality_Binary, 1, 2) != "01" & 
      str_sub(Quality_Binary, 1, 2) != "11" & 
      str_sub(Quality_Binary, 3, 4) != "11" & 
      str_sub(Quality_Binary, 3, 4) != "10" & 
      str_sub(Quality_Binary, 7, 8) != "11"
  )

# str_sub(Quality_Binary, 15) == "0" # possible ice/snow

tail(result_df)
nrow(result_df)

# result_df$Quality_Binary <- NULL # Deleting the Quality bits
# result_df$Quality <- NULL # Deleting the coloumn with the Quality

print(colnames(result_df)) # Get the names = we should rename the Year coloumn to the Date coloumn.
colnames(result_df)[1] <- "Date" # This makes the next steps more intuitive: Date and Year coloumn
result_df$Date <- as.Date(result_df$Date)
result_df$Year <- format(result_df$Date, "%Y")
head(result_df) # Have a look
str(result_df)



# TODO: Filter pixel groups with less than 3 observations per season
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
head(result_df)
str(result_df)


##################### Winter Spanning ##########################################
Year <- result_df$Year

result_df$SeasonYear <- ifelse(
  result_df$Season == "Winter" & month(result_df$Date) %in% c(11, 12),
  paste(Year, as.numeric(Year) + 1, sep="-"),
  ifelse(
    result_df$Season == "Winter" & month(result_df$Date) == 1,
    paste(as.numeric(Year) - 1, Year, sep="-"),
    Year
  )
)

# Then group by this new SeasonYear column along with x, y, and Season
# Filter the groups
result <- result_df %>%
  dplyr::group_by(Year, x, y, SeasonYear) %>%  
  dplyr::filter(n() >= 4) %>%
  dplyr::ungroup()


# View the first few rows of the result
head(result)
str(result)
nrow(result)


############# Detect clouds/ice ################################################

# TODO: Detect clouds/ice and filter pixel
detect_clouds <- function(Quality_Binary) {
  if (str_sub(Quality_Binary, 11, 11) == "1" # Mixed clouds
      ) { # possible ice/snow
    return(0) # No clouds/ice/snow
  } else {
    return(1) # Clouds/ice/snow present
  }
}

# str_sub(Quality_Binary, 9, 9) == "0" & # exclude adjacent clouds 

# Apply the function to create a new Clouds column
result$Clouds <- sapply(result$Quality_Binary, detect_clouds)
head(result)
str(result)

result$Clouds <- as.numeric(result$Clouds)  # Ensuring Clouds is numeric



################################################################################

# Exclude all pixel with more than X cloudy observations. 
# Function to check for consecutive cloudy observations

check_consecutive_clouds <- function(Clouds) {
  rle_clouds <- rle(Clouds)
  if(any(rle_clouds$values == 1 & rle_clouds$lengths > 3)) {
    return(TRUE)  # Indicates exclusion due to more than 4 consecutive cloudy observations
  } else {
    return(FALSE) # Indicates keeping the pixel group
  }
}


df_filtered <- result %>%
  group_by(x, y) %>%
  mutate(Exclude = check_consecutive_clouds(Clouds)) %>%
  ungroup() %>%
  filter(Exclude == FALSE) %>%
  dplyr::select(-Exclude)  # Use dplyr's select explicitly


# Now, view the result
head(df_filtered)
tail(df_filtered)
nrow(df_filtered)
################################################################################


# TODO: Safe and read the new file with all good values
saveRDS(df_filtered, file="~/edu/NDVI_germany/data/raster_data/data_level1/long_ndvi_clouds.rds")


# Identify the first 20 unique pixel groups based on x and y coordinates
unique_pixels <- distinct(df_filtered, x, y) %>% head(20)

# Filter the main dataframe for only these pixel groups
filtered_result <- df_filtered %>%
  semi_join(unique_pixels, by = c("x", "y"))

# Plotting
ggplot(filtered_result, aes(x = Date, y = NDVI_Value, group = interaction(x, y), color = interaction(x, y))) +
  geom_line() +
  geom_point(size = 0.5) +
  theme_minimal() +
  labs(title = "NDVI Values Over Time for First 20 Pixel Groups",
       x = "Date",
       y = "NDVI Value") +
  scale_x_date(date_labels = "%b %d", date_breaks = "1 month") +
  theme(legend.position = "none") # Remove legend to avoid clutter

# Note: Adjust `scale_x_date` for different time formats as needed



################################################################################
######################## Start of calculation ##################################
################################################################################
# TODO: Load the dataframe
# Note: The smoother value affects the result of SOS and EOS!
df_prep <- readRDS("~/edu/NDVI_germany/data/raster_data/data_level1/long_ndvi_clouds.rds")# Read the RDS file
head(df_prep)


# Step 1: Interpolate NDVI values for cloudy observations
result_interpolated <- result %>%
  group_by(x, y) %>%
  mutate(NDVI_Value_Interpolated = if_else(Clouds == 1, NA_real_, NDVI_Value)) %>%
  mutate(NDVI_Value_Interpolated = na.approx(NDVI_Value_Interpolated, na.rm = FALSE)) %>%
  ungroup()

head(result_interpolated)

# Step 2: Apply the Savitzky-Golay filter to the interpolated NDVI time series
result_smoothed <- result_interpolated %>%
  group_by(x, y) %>%
  mutate(NDVI_SG_Filtered = sgolayfilt(NDVI_Value_Interpolated, p = 3, n = 5)) %>%
  ungroup()

head(result_smoothed)



ppp <- subset(result_smoothed, Year ==2022)

ppp

# Ensure the Date column is in Date format
ppp$Date <- as.Date(ppp$Date)

# Extract Month from Date
ppp$Month <- month(ppp$Date, label = TRUE, abbr = TRUE)

# Group by pixel coordinates (x, y) and filter for the first 20 unique groups
ppp_grouped <- ppp %>%
  group_by(x, y) %>%
  mutate(GroupID = cur_group_id()) %>%
  ungroup()

# Plot
ggplot(ppp, aes(x = Month, y = NDVI_SG_Filtered, group = interaction(x, y), color = as.factor(GroupID))) +
  geom_line() +
  geom_point() +
  labs(title = "Monthly NDVI Values for last 4 Pixel Groups",
       x = "Month",
       y = "NDVI Value") +
  theme_minimal() +
  theme(legend.position = "none") # Hide legend to avoid clutter

plot(ppp$Month, ppp$NDVI_Value)
plot(ppp$Month, ppp$NDVI_Value_Interpolated)
plot(ppp$Month, ppp$NDVI_SG_Filtered, type="b")







################################################################################
library(zoo)

result_smoothed <- result_interpolated %>%
  group_by(x, y) %>%
  mutate(NDVI_Value_Interpolated = na.locf(NDVI_Value_Interpolated, na.rm = FALSE), # Fill NAs forward
         NDVI_Value_Interpolated = na.locf(NDVI_Value_Interpolated, fromLast = TRUE, na.rm = FALSE)) %>% # Fill NAs backward
  mutate(NDVI_SG_Filtered = ifelse(is.na(NDVI_Value_Interpolated), NA, 
                                   sgolayfilt(NDVI_Value_Interpolated, p = 3, n = 5))) %>%
  ungroup()

# Replace the original NDVI values with the smoothed values, if desired
# result_smoothed$NDVI_Value <- result_smoothed$NDVI_SG_Filtered
################################################################################
# View the result
