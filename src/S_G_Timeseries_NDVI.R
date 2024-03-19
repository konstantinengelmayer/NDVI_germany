# Written by Malte Date: 16.03.2024
## NDVI timeseries with Savitzky-Golay Filter & Cloud interpolation ############

################################################################################
# Get the whole dataframe
df <- readRDS("~/edu/NDVI_germany/data/raster_data/data_level1/long_ndvi_quality_df.rds") 
str(df)
################################################################################
# TODO: Preprocessing
# Step 1: Preprocessing

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
df$Quality <- NULL
str(df)


################################################################################
# TODO: Exclude bad pixel values
# filter the values, but let the clouds in the data
result_df <- df %>%
  dplyr::filter(
    str_sub(Quality_Binary, 1, 2) != "11" & # exclude all not produced observations
      str_sub(Quality_Binary, 1, 2) != "01" & # exclude all problematic observations 
      str_sub(Quality_Binary, 3, 4) != "11" & # delete all really bad pixel
      str_sub(Quality_Binary, 3, 4) != "10" # delete bad quality pixel
  )

str(result_df)
nrow(df)-nrow(result_df) # -591943 Observations
############# Detect clouds ################################################

# TODO: Detect whole clouds
detect_clouds <- function(Quality_Binary) {
  if (str_sub(Quality_Binary, 1, 2) == "10") { # Clouds
    return(1) # Most probably cloudy
  } else {
    return(0) # No clouds
  }
}

# Apply the function to create a new Clouds column
result_df$Clouds <- sapply(result_df$Quality_Binary, detect_clouds)
str(result_df)
result_df$Clouds <- as.numeric(result_df$Clouds)


################################################################################
# Step 2: 
# TODO: Delete all pixel groups with years with too many clouds or too less observations
# Convert Year to Date
result_df$Year <- as.Date(result_df$Year)

# Filter pixel groups with at least 17 observations in a year
df_filtered <- result_df %>%
  group_by(x, y, Year_group = floor_date(Year, "year")) %>%
  filter(n() >= 17) %>%
  ungroup() %>%
  select(-Year_group)

str(df_filtered)

# Remove pixel groups with two or more consecutive cloud observations from March to November
df_final <- df_filtered %>%
  group_by(x, y) %>%
  mutate(Clouds_lag = lag(Clouds), # Create a lag column for Clouds to identify consecutive cloud observations
         Month = month(Year)) %>%
  filter(!(Month >= 3 & Month <= 11 & Clouds == 1 & Clouds_lag == 1)) %>%
  select(-c(Clouds_lag, Month)) %>%
  ungroup()

str(df_final)

################################################################################
# Step 3
# TODO: Replace NDVI_Value with NA for cloudy observations and create lin_int for linear interpolation
df_final <- df_final %>%
  group_by(x, y) %>%
  mutate(lin_int = ifelse(Clouds == 1, NA, NDVI_Value)) %>% # Replace cloudy NDVI values with NA
  mutate(lin_int = na.approx(lin_int, na.rm = FALSE)) %>% # Perform linear interpolation
  ungroup()

str(df_final)

# saveRDS(df_final, "~/edu/NDVI_germany/data/raster_data/data_level1/cloud_interpolation.rds")
# df_final <- readRDS("~/edu/NDVI_germany/data/raster_data/data_level1/cloud_interpolation.rds")


################################################################################
# Step 4
# TODO: Apply the long term Savitzky-Golay Filter
# Function to apply Savitzky-Golay filter to a vector

apply_savgol <- function(vec, m = 5, d = 3) {
  # Calculate window size
  window_size <- 2 * m + 1
  
  # Apply Savitzky-Golay filter
  sg_filtered <- sgolayfilt(vec, n = window_size, p = d, ts = 1)
  
  # Return the filtered vector
  return(sg_filtered)
}

# Applying the Savitzky-Golay filter to the lin_int column for each pixel group
df_final <- df_final %>%
  group_by(x, y) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(Savitzky = ifelse(is.na(lin_int), NA, apply_savgol(lin_int, m = 5, d = 3))) %>%
  ungroup()

# Handle potential NA values at the edges if needed
df_final$Savitzky <- na.approx(df_final$Savitzky, na.rm = FALSE)


str(df_final)


################################################################################
# Step 5
# TODO: Weights
# Create the weights for the NDVI Values that has been changed through the previous Savitzky-Golay filter

df_final <- df_final %>%
  group_by(x, y) %>% # Group data by each unique pixel
  mutate(weight = ifelse(lin_int >= Savitzky, 
                         1, 
                         1 - (Savitzky - lin_int) / max(Savitzky - lin_int, na.rm = TRUE))) %>%
  ungroup() # Ensure that the data is not grouped for subsequent operations


################################################################################
# Step 6
# TODO: 16 day time-series
# Create a new NDVI time-series with the weights and the Savitzky-Golay Filter

df_final <- df_final %>%
  mutate(Adjusted_NDVI = (lin_int * weight + Savitzky * (1 - weight)))

# saveRDS(df_final, "~/edu/NDVI_germany/data/raster_data/data_level1/Savitzky_long_term.rds")
# df_final <- readRDS("~/edu/NDVI_germany/data/raster_data/data_level1/Savitzky_long_term.rds")
str(df_final)
################################################################################
# Step 7
# TODO: Fit the adjusted NDVI values
df_final$Quality_Binary <- NULL # delete the quality coloumn

# Iterative fitting process where k = the number of iterations

# Apply the Savitzky-Golay filter to Adjusted_NDVI for the iterative fitting process
df_final <- df_final %>%
  group_by(x, y) %>%
  # Apply sgolayfilt with parameters for the iterative fitting (m=4, d=6 for the first iteration)
  mutate(Fitted_NDVI = sgolayfilt(Adjusted_NDVI, p = 6, n = 9)) %>%
  ungroup()

str(df_final)

################################################################################
# Step 8
# TODO: Fit the fitted NDVI values iteratively
# Iterative fitting with exit criteria

library(pracma) # For sgolayfilt function

# Define the custom function for iterative fitting
iterative_fitting <- function(data) {
  # Parameters for Savitzky-Golay filter
  p <- 6
  n <- 9
  max_iterations <- 10 # Maximum number of iterations
  
  # Initialize Final_NDVI with Fitted_NDVI
  data$Final_NDVI <- data$Fitted_NDVI
  
  # Loop for iterative fitting
  for (iter in 1:max_iterations) {
    # Apply the Savitzky-Golay filter for the current iteration
    data$Final_NDVI <- sgolayfilt(data$Final_NDVI, p, n)
    
    # Calculate the fitting-effect index (Fk) for the current iteration
    temp_fitting_effect <- abs(data$Final_NDVI - data$lin_int) * data$weight
    Fk <- sum(temp_fitting_effect, na.rm = TRUE) # Summarize over this pixel group
    
    # Check exit condition for this group based on Fk
    if (iter > 1 && Fk >= prev_Fk) {
      break # Exit loop if Fk does not decrease
    }
    prev_Fk <- Fk # Store Fk for comparison in the next iteration
  }
  
  return(data)
}

# Apply the iterative fitting process to each pixel group separately
df_final <- df_final %>%
  group_by(x, y) %>%
  group_modify(~ iterative_fitting(.x)) %>%
  ungroup()

# After this, `df_final` contains the Final_NDVI after iterative fitting for each pixel group

# saveRDS(df_final, "~/edu/NDVI_germany/data/raster_data/data_level1/Savitzky_iterative.rds")
# df_final <- readRDS("~/edu/NDVI_germany/data/raster_data/data_level1/Savitzky_iterative.rds")


################################################################################
df_final$Date <- as.Date(df_final$Year)
df_final$Year <- format(df_final$Year, "%Y")
str(df_final) # Have a look
################################################################################
############################### EXCLUDE GREEN WINTER ###########################

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

df_final$Season <- sapply(month(df_final$Date), getSeason)

# Step 1: Filter for Winter season and NDVI values above 0.7
winter_high_ndvi <- df_final %>%
  filter(Season == "Winter" & Final_NDVI > 0.6)

# Step 2 & 3: Group by pixel coordinates, then count unique years meeting the condition
pixel_groups_high_ndvi <- winter_high_ndvi %>%
  group_by(x, y) %>%
  summarize(Unique_Winter_Years = n_distinct(Year)) %>%
  filter(Unique_Winter_Years >= 5) %>%
  ungroup()

# Step 4: Create green_winter dataframe to INCLUDE these pixel groups
# This dataframe will include only the records from pixel groups that met the exclusion criteria
green_winter <- df_final %>%
  semi_join(pixel_groups_high_ndvi, by = c("x", "y"))

# Step 5: Update df_final to EXCLUDE these pixel groups
# This dataframe will now exclude the pixel groups that met the criteria
df_savitzky_updated <- df_final %>%
  anti_join(pixel_groups_high_ndvi, by = c("x", "y"))

# Output the updated df_final dataframe without the excluded pixel groups
df_final <- df_savitzky_updated

print(paste("Number of excluded observations during green winter filter:", nrow(df_final)-nrow(green_winter)))

# Save file:
saveRDS(df_final, "~/edu/NDVI_germany/data/raster_data/data_level1/Savitzky_without_green_winter.rds")
# df_final <- readRDS("~/edu/NDVI_germany/data/raster_data/data_level1/Savitzky_without_green_winter.rds")
str(df_final)

################################################################################
################################################################################
# Evaluate the NDVI time-series with one pixel:

set.seed(102)
unique_pixel <- df_final %>% distinct(x, y) %>% sample_n(1)
# Filter the main dataframe to include only the selected pixel group
dfp <- df_final %>%
  filter(x == unique_pixel$x & y == unique_pixel$y)

################## PLOT ########################################################
# Filter the dataframe for the year xxxx
dfp_xxxx <- dfp %>%
  filter(as.Date(Date) >= as.Date("2013-01-01") & as.Date(Date) <= as.Date("2013-12-31"))

# Plotting the Original and the interpolated and the smoothed data
ggplot(dfp_xxxx, aes(x = Date)) +
  geom_line(aes(y = NDVI_Value, color = "Original NDVI"), linewidth=1) +
  geom_line(aes(y = lin_int, color = "Interpolated NDVI"), linewidth=1) +
  geom_line(aes(y = Savitzky, color = "Smoothed NDVI"), linewidth=1) +
  scale_color_manual(values = c("Original NDVI" = "blue", 
                                "Interpolated NDVI" = "green", 
                                "Smoothed NDVI" = "red")) +
  labs(title = "NDVI Values through one Year", 
       y = "NDVI Value", 
       color = "Legend") +
  theme_minimal()

# Plotting the original and the different steps of Savitzky smoothig
ggplot(dfp_xxxx, aes(x = Date)) +
  geom_line(aes(y = NDVI_Value, color ="Original NDVI"), linewidth=0.8) +
  geom_line(aes(y = Adjusted_NDVI, color = "Adjusted_NDVI"), linewidth=0.9) +
  geom_line(aes(y = Fitted_NDVI, color = "Fitted_NDVI"), linewidth=1.1) +
  geom_line(aes(y = Final_NDVI, color = "Final_NDVI"), linewidth=1.2) +
  scale_color_manual(values = c("Adjusted_NDVI" = "blue", 
                                "Fitted_NDVI" = "green", 
                                "Final_NDVI" = "red",
                                "Original NDVI" = "grey")) +
  labs(title = "NDVI Values through one Year", 
       y = "NDVI Value", 
       color = "Legend") +
  theme_minimal()

################################################################################
################################################################################

# TODO: Create a sequence of daily dates covering the range in your data
seq_dates <- seq(from = min(df_final$Date), to = max(df_final$Date), by = "day")
head(seq_dates)

# Initialize a counter for excluded groups
excluded_groups_counter <- 0

# Define a wrapper function for the interpolation step to catch errors
safe_interpolate <- function(dates, values, seq_dates) {
  tryCatch({
    # Perform linear interpolation
    interpolated_values <- approx(x = as.numeric(dates), y = values,
                                  xout = as.numeric(seq_dates), method = "linear", rule = 2)$y
    return(interpolated_values)
  }, error = function(e) {
    # Increment the counter if an error occurs
    excluded_groups_counter <<- excluded_groups_counter + 1
    return(rep(NA, length(seq_dates))) # Return a vector of NAs of the same length as seq_dates
  })
}

# Apply the safe interpolation function within group_modify
grouped <- df_final %>%
  group_by(x, y, Year) %>%
  group_modify(~ {
    # Create a sequence of daily dates for the current year
    start_date <- min(.x$Date)
    end_date <- max(.x$Date)
    seq_dates <- seq(from = start_date, to = end_date, by = "day")
    
    # Interpolate NDVI values for these dates using the safe wrapper function
    daily_ndvi <- safe_interpolate(.x$Date, .x$Final_NDVI, seq_dates)
    
    # Construct the result dataframe
    data.frame(Date = seq_dates, NDVI_Value = daily_ndvi)
  }) %>%
  ungroup()

# Print the number of excluded pixel groups
print(paste("Number of excluded pixel groups during daily interpolation:", excluded_groups_counter))

str(grouped)
unique(is.na(grouped$NDVI_Value))
################################################################################

# Step 1: Identify pixel groups with any NA values in NDVI_Value
groups_with_na <- grouped %>%
  group_by(x, y) %>%
  filter(any(is.na(NDVI_Value))) %>%
  summarise() %>%
  ungroup()

# Count the number of excluded pixel groups
number_of_excluded_groups <- nrow(groups_with_na)

# Step 2: Exclude these groups from the original dataframe
grouped <- grouped %>%
  anti_join(groups_with_na, by = c("x", "y"))

# Print the count of excluded pixel groups
print(paste("Number of excluded pixel groups with NA:", number_of_excluded_groups))

# saveRDS(grouped, "~/edu/NDVI_germany/data/raster_data/data_level1/Savitzky_iterative_daily.rds")
grouped <- readRDS("~/edu/NDVI_germany/data/raster_data/data_level1/Savitzky_iterative_daily.rds")


################################################################################
########################## CALCULATE SOS EOS APG ###############################

################################################################################
# Calculate SOS & EOS with !!!!! pixelwise !!!!! local threshold:
calculate_sos_eos_apg_pgt_v2 <- function(df) {
    # Early return for empty groups
   if(nrow(df) == 0) return(data.frame(SOS = NA, EOS = NA, APG = NA, PGT = NA, sos_ndvi_value = NA, eos_ndvi_value = NA))
    
  NDVImax <- max(df$NDVI_Value, na.rm = TRUE)
  NDVImin <- min(df$NDVI_Value, na.rm = TRUE)
    
    # Guard clause for infinite values or identical min and max
   if(is.infinite(NDVImax) || is.infinite(NDVImin) || NDVImax == NDVImin) {
    return(data.frame(SOS = NA, EOS = NA, APG = NA, PGT = NA, sos_ndvi_value = NA, eos_ndvi_value = NA))
  }
    
    # Normalize NDVI values for the pixel group
  df$NDVIratio <- (df$NDVI_Value - NDVImin) / (NDVImax - NDVImin)
    
    # Calculate APG (Annual Peak Greenness) and PGT (Peak Greenness Time)
  APG <- NDVImax
  PGT <- df$Date[which.max(df$NDVI_Value)]
    
    # Determine SOS based on normalized threshold
   sos_index <- which(df$NDVIratio >= 0.5 & lag(df$NDVIratio, default = first(df$NDVIratio)) < 0.5)[1]
   if (!is.na(sos_index)) {
    sos <- df$Date[sos_index]
    sos_ndvi_value <- df$NDVI_Value[sos_index]
  } else {
    sos <- NA
    sos_ndvi_value <- NA
  }
    
    # Filter data after PGT for EOS calculation and apply the same normalization logic
  df_post_PGT <- df %>% filter(Date > PGT)
  if (nrow(df_post_PGT) > 0) {
    eos_index <- which(df_post_PGT$NDVIratio < 0.5 & lag(df_post_PGT$NDVIratio, default = first(df_post_PGT$NDVIratio)) >= 0.5)[1]
    if (!is.na(eos_index)) {
      eos <- df_post_PGT$Date[eos_index]
      eos_ndvi_value <- df_post_PGT$NDVI_Value[eos_index]
    } else {
      eos <- NA
      eos_ndvi_value <- NA
    }
  } else {
    eos <- NA
    eos_ndvi_value <- NA
    }
  
  return(data.frame(SOS = sos, EOS = eos, APG = APG, PGT = PGT, sos_ndvi_value = sos_ndvi_value, eos_ndvi_value = eos_ndvi_value))
}

# Apply the function to each group
sos <- grouped %>%
  group_by(x, y, Year) %>%
  group_modify(~ calculate_sos_eos_apg_pgt_v2(.x)) %>%

############################## visualize results ###############################
str(sos)
mean(sos$sos_ndvi_value, na.rm=T)
mean(sos$eos_ndvi_value, na.rm=T)
################################################################################

# Count the NA in the coloumns
na_counts <- sos %>%
  summarise_all(~sum(is.na(.)))

print(na_counts)

# Exclude all !rows! with NA
clean_sos <- sos %>%
  # Filter out rows with any NA values in any column
  filter(!if_any(everything(), is.na))

print(paste("Number of excluded Years during NA deleting:", nrow(sos)-nrow(clean_sos)))
nrow(clean_sos)
str(clean_sos)

################################################################################

# Extract month from EOS date
clean_sos <- clean_sos %>%
  mutate(EOS_month = format(EOS, "%m"), # Extract month as a string
         EOS_month = as.integer(EOS_month)) # Convert to integer for ordering

# Extract month from SOS date
clean_sos <- clean_sos %>%
  mutate(SOS_month = format(SOS, "%m"), # Extract month as a string
         SOS_month = as.integer(SOS_month)) # Convert to integer for ordering

# Extract month from PGT date
clean_sos <- clean_sos %>%
  mutate(PGT_month = format(PGT, "%m"), # Extract month as a string
         PGT_month = as.integer(PGT_month)) # Convert to integer for ordering


# Count the number of occurrences for each SOS/EOS/PGT month
sos_month_counts <- clean_sos %>%
  count(SOS_month) %>%
  rename(Month = SOS_month, SOS_Count = n)

eos_month_counts <- clean_sos %>%
  count(EOS_month) %>%
  rename(Month = EOS_month, EOS_Count = n)

pgt_month_counts <- clean_sos %>%
  count(PGT_month) %>%
  rename(Month = PGT_month, PGT_Count = n)

# Sequentially joining the month counts
month_counts <- sos_month_counts %>%
  full_join(eos_month_counts, by = "Month") %>%
  full_join(pgt_month_counts, by = "Month")

# Replace NA with 0 for counts that don't appear in some datasets
month_counts[is.na(month_counts)] <- 0

# Print the result to check
print(month_counts)

################################################################################
# TODO: Exclude outliers

filtered_sos <- clean_sos %>%
  filter(!(EOS_month %in% c(1, 2, 3, 4, 5, 6, 7, 12))) %>%
  filter(!(SOS_month %in% c(1, 2, 6, 7, 8, 9, 10, 11, 12))) %>%
  filter(!(PGT_month %in% c(1, 2, 3, 4, 9, 10, 11, 12)))

str(filtered_sos)


# TODO: Exclude pixel groups that have less than xy Years
# CAUTION: The higher the forced n (=year), the higher the slope!!
filtered_excl <- filtered_sos %>%
  # Group by pixel coordinates
  group_by(x, y) %>%
  # Add count of non-NA SOS and EOS values for each group
  mutate(SOS_count = sum(!is.na(SOS)),
         EOS_count = sum(!is.na(EOS))) %>%
  # Filter out groups with less than xy non-NA values for either SOS or EOS
  filter(SOS_count >= 18 & EOS_count >= 18) %>%
  # Remove the count columns as they are no longer needed
  select(-SOS_count, -EOS_count) %>%
  # Ungroup the data
  ungroup()

nrow(filtered_excl)
print(paste("Number of excluded Years during outliers deleting:", nrow(filtered_sos)-nrow(filtered_excl)))
str(filtered_excl)

################################################################################
# TODO: Converting SOS and EOS to day of the year
filtered_excl <- filtered_excl %>%
  mutate(
    Year = as.numeric(Year),  # Convert Year to numeric
    SOS_julian = yday(SOS)  # Convert SOS to the day of the year
  )

# EOS 
filtered_excl <- filtered_excl %>%
  mutate(
    Year = as.numeric(Year),
    EOS_julian =yday(EOS)
  )



################################################################################

# TODO: Linear Regression
# linear regression with sos and year
lm_result <- lm(SOS_julian ~ Year, data = filtered_excl)
# Display the summary of the linear model
summary(lm_result)

# linear regression with eos and year
lm_result <- lm(EOS_julian ~Year, data = filtered_excl)
summary(lm_result)
################################################################################
# TODO: Check the normal distribution!
#shapiro.test(filtered_excl$SOS_julian) # no normal distribution!
#shapiro.test(filtered_excl$EOS_julian)# no normal distribution!

# Perform Shapiro-Wilk test on the residuals
#shapiro.test(residuals(lm_result))
qqnorm(residuals(lm_result))
qqline(residuals(lm_result))
################################################################################

# Perform linear regression with EOS_julian as the dependent variable
lm_result <- lm(EOS_julian ~ SOS_julian, data = filtered_excl)
# Display the summary of the linear model
summary(lm_result)

################################################################################
# TODO: Look at one pixel
# Define the coordinates for the pixel group to be extracted with a small tolerance for floating-point comparison
x_coord <- 8.663542
y_coord <- 50.01146
tolerance <- 1e-5  # Define a tolerance level, e.g., 0.00001

# Filter the dataframe to extract the specified pixel group within the tolerance range
df_2 <- clean_sos %>%
  filter(abs(x - x_coord) < tolerance & abs(y - y_coord) < tolerance)
head(df_2)
plot(df_2$SOS_julian, df_2$EOS_julian)

################################################################################


cor.test(filtered_excl$SOS_julian, filtered_excl$EOS_julian)

glmer
