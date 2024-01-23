# NDVI Data Processing for Consistent Broad-Leaved Forest Areas in Germany

# This script processes NDVI (Normalized Difference Vegetation Index) data specifically for Hessen, 
# focusing on areas with consistent broad-leaved forests. The script involves several key steps:
# - Setting up the environment and loading necessary libraries and functions.
# - Creating a stack of NDVI rasters from provided TIFF files.
# - Similarly processing quality layer data, which is essential for assessing the NDVI data's reliability.
# - Loading and processing consistent broad-leaved forest data, including resampling to match the NDVI stack's resolution.
# - Masking the NDVI data using the processed broad-leaved forest data to focus the analysis on these forested areas.
# - Finally, converting the masked NDVI data into a DataFrame with spatial coordinates, enabling further analysis 
#   and comparison with the quality layer. This DataFrame is then saved for use in subsequent analytical processes.


# Load necessary functions and libraries
source("~/edu/NDVI_germany/src/NDVI_germany_setup.R")

# List all NDVI files in the specified directory
tif_files <- list.files(path = envrmt$path_ndvi, pattern = "*.tif", full.names = TRUE)

# Create a stack of NDVI rasters
ndvi_stack <- rast(tif_files)
ndvi_stack

# Extract file names for later date conversion
tif_names <- list.files(path = envrmt$path_ndvi, pattern = "*.tif")

# Function to convert file names to dates
convert_date <- function(tif_names) {
  year <- substr(tif_names, 35, 38)
  day_of_year <- as.integer(substr(tif_names, 39, 41))
  date <- as.Date(paste(year, "-01-01", sep = ""), format = "%Y-%m-%d") + 
    days(day_of_year - 1)
  return(as.character(date))
}

# Convert date format for all file names
converted_dates <- sapply(tif_names, convert_date)
names(ndvi_stack) <- converted_dates
ndvi_stack

# List all quality layer files in the specified directory
quality_files <- list.files(path = envrmt$path_quality_layer, pattern = "*.tif", full.names = TRUE)

# Create a stack of quality layer rasters
quality_stack <- rast(quality_files)
quality_stack

# Extract file names for later date conversion
quality_names <- list.files(path = envrmt$path_quality_layer, pattern = "*.tif")

# Function to convert quality layer file names to dates
convert_date_qual <- function(quality_names) {
  year <- substr(quality_names, 41, 44)
  day_of_year <- as.integer(substr(quality_names, 45, 47))
  date <- as.Date(paste(year, "-01-01", sep = ""), format = "%Y-%m-%d") + 
    days(day_of_year - 1)
  return(as.character(date))
}

# Convert date format for all quality layer file names
converted_dates <- sapply(quality_names, convert_date_qual)
names(quality_stack) <- converted_dates
quality_stack

# Load and process consistent broad-leaved forest data
consistent_broad_leaved <- rast(paste0(envrmt$path_data_level1, "/consistent_broad_leaved_forest.tif"))
resampled_broad_leaved <- resample(consistent_broad_leaved, ndvi_stack, method='bilinear')
coreccted_broad_leaved <- resampled_broad_leaved == 1
coreccted_broad_leaved <- crop(coreccted_broad_leaved, ndvi_stack)
writeRaster(coreccted_broad_leaved, paste0(envrmt$path_data_level1, "/consistent_broad_leaved_resampled.tif"))
consistent_broad_leaved <- rast(paste0(envrmt$path_data_level1, "/consistent_broad_leaved_resampled.tif"))

# Mask the NDVI stack with consistent broad-leaved forest data
masked_ndvi <- mask(ndvi_stack, consistent_broad_leaved, maskvalue=0)

# Step 2: Convert the masked NDVI stack to a DataFrame
ndvi_df <- as.data.frame(masked_ndvi, xy=TRUE, na.rm=TRUE)
saveRDS(ndvi_df, paste0(envrmt$path_data_level1, "/ndvi_df.rds"))

ndvi_df <- readRDS(paste0(envrmt$path_data_level1, "/ndvi_df.rds"))

# Reshape the DataFrame to a long format
long_ndvi_df <- ndvi_df %>%
  pivot_longer(
    cols = -c(x, y), # Exclude the position columns from reshaping
    names_to = "Year", # Name of the new column for years
    values_to = "NDVI_Value" # Name of the new column for NDVI values
  )

saveRDS(long_ndvi_df, paste0(envrmt$path_data_level1, "/long_ndvi_df.rds"))
long_ndvi_df <- readRDS(paste0(envrmt$path_data_level1, "/long_ndvi_df.rds"))

quality_data <- data.frame(Year = integer(), x = numeric(), y = numeric(), Quality = integer())

# Loop through each year
for(year in unique(long_ndvi_df$Year)) {
  # Extract the quality layer for the year
  quality_layer <- quality_stack[[year]]
  
  # Subset long_ndvi_df for the current year
  year_data <- subset(long_ndvi_df, Year == year)
  
  # Extract the quality values for the corresponding positions
  extracted_values <- terra::extract(quality_layer, cbind(year_data$x, year_data$y))
  
  # Create a temporary data frame for the current year's data
  temp_df <- data.frame(Year = rep(year, nrow(year_data)), 
                        x = year_data$x, 
                        y = year_data$y, 
                        Quality = extracted_values)
  names(temp_df)[4] <- c("Quality")
  # Combine the temporary data frame with the main quality data frame
  quality_data <- rbind(quality_data, temp_df)
}

# merge resulting df with the long ndvi df
final_df <- merge(long_ndvi_df, quality_data, by = c("Year", "x", "y"))
saveRDS(final_df, paste0(envrmt$path_data_level1, "/long_ndvi_quality_df.rds"))
long_ndvi_df <- readRDS(paste0(envrmt$path_data_level1, "/long_ndvi_df.rds"))
final_df <- readRDS(paste0(envrmt$path_data_level1, "/long_ndvi_quality_df.rds"))

# function to convert integer to binarys
int_to_15bit_binary <- function(value) {
  # Convert the integer to a vector of bits and take the first 15 bits
  bits <- as.integer(intToBits(value))[1:15]
  
  # Convert the bits vector to a binary string
  binary_str <- paste(rev(bits), collapse = "")
  
  return(binary_str)
}
# creatae df of convertet unique Quality values to binary
binarys <- sapply(unique(final_df$Quality), int_to_15bit_binary)
binarys_df <- data.frame(integer = unique(final_df$Quality), binarys = binarys)

# merge over unique quality values
final_df_with_binary <- merge(final_df, binarys_df, by.x = "Quality", by.y = "integer")
final_df_with_binary <- final_df_with_binary[,c(2,3,4,5,1,6)]

saveRDS(final_df_with_binary, paste0(envrmt$path_data_level1, "/master_df.rds"))
test <- readRDS(paste0(envrmt$path_data_level1, "/master_df.rds"))

result_df <- test %>%
  filter(
    str_sub(binarys, 1, 2) == "00" & 
      str_sub(binarys, 3, 4) != "11" &
      str_sub(binarys, 3, 4) != "1" &
      str_sub(binarys, 8, 8) == "0" & 
      str_sub(binarys, 10, 10) == "0" &
      str_sub(binarys, 14, 14) == "0" &
      str_sub(binarys, 15, 15) == "0" 
  )

str(result_df)
result_df$Year <- as.Date(result_df$Year, format = "%Y-%m-%d" )
colnames(result_df)[1] <- "date"

result_df_Year <- result_df %>%
  mutate(Year = year(date)) %>%
  group_by(Year) %>%
  summarize(Mean_NDVI = mean(NDVI_Value, na.rm = TRUE))

result_df_year_pixel <- result_df %>%
  mutate(Year = year(date)) %>%
  group_by(x, y, Year) %>%
  summarize(Mean_NDVI = mean(NDVI_Value, na.rm = TRUE), .groups = 'drop')

ggplot(result_df_year_pixel, aes(x = Year, y = Mean_NDVI))+
  geom_point()

get_season <- function(date) {
  month <- month(date)
  ifelse(month %in% c(12, 1, 2), 'Winter',
         ifelse(month %in% 3:5, 'Spring',
                ifelse(month %in% 6:8, 'Summer', 'Autumn')))
}

result_df_season_pixel <- result_df %>%
  mutate(
    Year = year(date),  # Extract year
    Season = get_season(date)  # Assign season
  ) %>%
  group_by(x, y, Year, Season) %>%  # Group by pixel location, year, and season
  summarize(Mean_NDVI = mean(NDVI_Value, na.rm = TRUE), .groups = 'drop')  # Calculate mean NDVI for each group

ggplot(result_df_season_pixel, aes(x = Year, y = Mean_NDVI, col = Season))+
  geom_point()

spring_data <- result_df_season_pixel %>% 
  filter(Season == "Spring")

ggplot(spring_data, aes(x = Year, y = Mean_NDVI)) +
  geom_point()

summer_data <- result_df_season_pixel %>% 
  filter(Season == "Summer")

ggplot(summer_data) +
  geom_boxplot(aes(y=Mean_NDVI, group = Year))


result_df_season <- result_df %>%
  mutate(
    Year = year(date),  # Extract year
    Season = get_season(date)  # Assign season
  ) %>%
  group_by(Year, Season) %>%  # Group by pixel location, year, and season
  summarize(Mean_NDVI = mean(NDVI_Value, na.rm = TRUE), .groups = 'drop')  # Calculate mean NDVI for each group

ggplot(result_df_season, aes(x = Year, y = Mean_NDVI, col = Season)) +
  geom_line()
result_df_season  


wide_df <- result_df_season %>%
  pivot_wider(names_from = Season, values_from = Mean_NDVI)

head(wide_df)

model <- lm(wide_df$Autumn~wide_df$Spring) 
summary(model)
# Base plot
plot(wide_df$Spring, wide_df$Autumn, col = 'blue', pch = 19,
     xlab = "Spring NDVI", ylab = "Autumn NDVI", main = "NDVI Spring vs Autumn")

# Add regression line
abline(model, col = "black", lwd = 2)

# Annotate each point with the corresponding year
# Adjusting text position by adding a small offset to x and y coordinates
text_offset <- 0.006
text(wide_df$Spring + text_offset, wide_df$Autumn, labels = wide_df$Year, cex = 0.8, col = 'darkgrey')

