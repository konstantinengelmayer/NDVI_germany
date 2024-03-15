################################################################################
# get a random pixel for fast test calculation
set.seed(99)
unique_pixel <- df %>% distinct(x, y) %>% sample_n(1)
# Filter the main dataframe to include only the selected pixel group
df_random_group <- df %>%
  filter(x == unique_pixel$x & y == unique_pixel$y)
nrow(df_random_group)
dfp <- df_random_group


## NDVI timeseries with Savitzky-Golay Filter & Cloud interpolation ############

################################################################################
# Get the whole dataframe
df <- readRDS("~/edu/NDVI_germany/data/raster_data/data_level1/long_ndvi_quality_df.rds") 
str(df)
################################################################################
# Step 0: Preprocessing

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
head(df)
str(df)


################################################################################
# filter the values
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
# Delete all pixel groups with years with too many clouds or too less observations
# Convert Year to Date
result_df$Year <- as.Date(result_df$Year)

# Step 1: Filter pixel groups with at least 17 observations in a year
df_filtered <- result_df %>%
  group_by(x, y, Year_group = floor_date(Year, "year")) %>%
  filter(n() >= 17) %>%
  ungroup() %>%
  select(-Year_group)

str(df_filtered)

# Step 2: Remove pixel groups with two or more consecutive cloud observations from March to November
df_final <- df_filtered %>%
  group_by(x, y) %>%
  mutate(Clouds_lag = lag(Clouds), # Create a lag column for Clouds to identify consecutive cloud observations
         Month = month(Year)) %>%
  filter(!(Month >= 3 & Month <= 11 & Clouds == 1 & Clouds_lag == 1)) %>%
  select(-c(Clouds_lag, Month)) %>%
  ungroup()

str(df_final)

################################################################################
# Step 1
# Replace NDVI_Value with NA for cloudy observations and create lin_int for linear interpolation
df_final <- df_final %>%
  group_by(x, y) %>%
  mutate(lin_int = ifelse(Clouds == 1, NA, NDVI_Value)) %>% # Replace cloudy NDVI values with NA
  mutate(lin_int = na.approx(lin_int, na.rm = FALSE)) %>% # Perform linear interpolation
  ungroup()

str(df_final)

# saveRDS(df_final, "~/edu/NDVI_germany/data/raster_data/data_level1/cloud_interpolation.rds")
df_final <- readRDS("~/edu/NDVI_germany/data/raster_data/data_level1/cloud_interpolation.rds")



################################################################################
# Apply the long term Savitzky-Golay Filter
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
# Create the weights for the NDVI Values that has been changed through the previous Savitzky-Golay filter

df_final <- df_final %>%
  group_by(x, y) %>% # Group data by each unique pixel
  mutate(weight = ifelse(lin_int >= Savitzky, 
                         1, 
                         1 - (Savitzky - lin_int) / max(Savitzky - lin_int, na.rm = TRUE))) %>%
  ungroup() # Ensure that the data is not grouped for subsequent operations


################################################################################
# Create a new NDVI time-series with the weights and the Savitzky-Golay Filter

df_final <- df_final %>%
  mutate(Adjusted_NDVI = (lin_int * weight + Savitzky * (1 - weight)))

# saveRDS(df_final, "~/edu/NDVI_germany/data/raster_data/data_level1/Savitzky_long_term.rds")
df_final <- readRDS("~/edu/NDVI_germany/data/raster_data/data_level1/Savitzky_long_term.rds")
str(df_final)
################################################################################
df_final$Quality_Binary <- NULL

# Iterative fitting process where k = the number of iterations

# Apply the Savitzky-Golay filter to Adjusted_NDVI for the iterative fitting process
df_final <- df_final %>%
  group_by(x, y) %>%
  # Apply sgolayfilt with parameters for the iterative fitting (m=4, d=6 for the first iteration)
  mutate(Fitted_NDVI = sgolayfilt(Adjusted_NDVI, p = 6, n = 9)) %>%
  ungroup()

str(df_final)

################################################################################
# Iterative fitting with exit criteria
# Initialize variables for iteration
# Initialize variables for iteration
# Set parameters for Savitzky-Golay filter
p <- 6
n <- 9

# Initialize Final_NDVI with Fitted_NDVI before the loop
df_final <- df_final %>%
  mutate(Final_NDVI = Fitted_NDVI)

k <- 1
max_iterations <- 10 # To prevent infinite loops
fitting_effect_indices <- numeric(max_iterations) # Store fitting-effect indices

for (iter in 1:max_iterations) {
  # Apply the Savitzky-Golay filter for the current iteration to Final_NDVI
  df_final <- df_final %>%
    group_by(x, y) %>%
    mutate(Final_NDVI = sgolayfilt(Final_NDVI, p, n)) %>%
    ungroup()
  
  # Calculate the fitting-effect index (Fk) for the current iteration using Final_NDVI and lin_int
  df_final <- df_final %>%
    mutate(temp_fitting_effect = abs(Final_NDVI - lin_int) * weight)
  
  Fk <- sum(df_final$temp_fitting_effect, na.rm = TRUE) # Summarize over all pixels
  fitting_effect_indices[iter] <- Fk
  
  # Check exit condition
  if (iter > 1 && (fitting_effect_indices[iter] >= fitting_effect_indices[iter-1] || iter >= max_iterations)) {
    message("Exiting at iteration: ", iter)
    break # Exit the loop if Fk does not decrease or max iterations reached
  }
}

# Cleanup after the loop, removing the temporary fitting effect calculation column
df_final <- df_final %>%
  select(-temp_fitting_effect)

str(df_final)

#saveRDS(df_final, "~/edu/NDVI_germany/data/raster_data/data_level1/Savitzky_iterative.rds")
df_Savitzky <- readRDS("~/edu/NDVI_germany/data/raster_data/data_level1/Savitzky_iterative.rds")
str(df_Savitzky)


################################################################################
df_Savitzky$Date <- as.Date(df_Savitzky$Year)
df_Savitzky$Year <- format(df_Savitzky$Year, "%Y")
head(df_Savitzky) # Have a look
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

df_Savitzky$Season <- sapply(month(df_Savitzky$Date), getSeason)


# Step 1: Filter for Winter season and NDVI values above 0.7
winter_high_ndvi <- df_Savitzky %>%
  filter(Season == "Winter" & Final_NDVI > 0.7)

# Step 2 & 3: Group by pixel coordinates, then count unique years meeting the condition
pixel_groups_high_ndvi <- winter_high_ndvi %>%
  group_by(x, y) %>%
  summarize(Unique_Winter_Years = n_distinct(Year)) %>%
  filter(Unique_Winter_Years >= 5) %>%
  ungroup()

# Step 4: Create green_winter dataframe to INCLUDE these pixel groups
# This dataframe will include only the records from pixel groups that met the exclusion criteria
green_winter <- df_Savitzky %>%
  semi_join(pixel_groups_high_ndvi, by = c("x", "y"))

# Step 5: Update df_Savitzky to EXCLUDE these pixel groups
# This dataframe will now exclude the pixel groups that met the criteria
df_savitzky_updated <- df_Savitzky %>%
  anti_join(pixel_groups_high_ndvi, by = c("x", "y"))

# Output the updated df_Savitzky dataframe without the excluded pixel groups
df_Savitzky <- df_savitzky_updated

# Print the first few rows of each dataframe to verify
print(paste("Number of excluded observations during green winter filter:", nrow(df_Savitzky)-nrow(green_winter)))

################################################################################


# TODO: Create a sequence of daily dates covering the range in your data
seq_dates <- seq(from = min(df_Savitzky$Date), to = max(df_Savitzky$Date), by = "day")
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
grouped <- df_Savitzky %>%
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
# Create a test dataset to calculate just with xx%
unique_pixel_groups <- grouped %>%
  distinct(x, y)

# select xx% of these pixel groups
set.seed(12) # Setting seed for reproducibility
selected_pixel_groups <- unique_pixel_groups %>%
  sample_frac(0.10
  )

# Filter the original dataframe to only include data from the selected pixel groups
selected_df <- grouped %>%
  semi_join(selected_pixel_groups, by = c("x", "y"))

grouped <- selected_df # calculate still with df

str(grouped)

################################################################################
########################## CALCULATE SOS EOS APG ###############################
calculate_sos_eos_apg_pgt_v2 <- function(df) {
  # Initialize return values
  sos <- NA
  eos <- NA
  sos_ndvi_value <- NA
  eos_ndvi_value <- NA
  
  # Calculate APG (Annual Peak Greenness) and PGT (Peak Greenness Time)
  APG <- max(df$NDVI_Value, na.rm = TRUE)
  PGT <- df$Date[which.max(df$NDVI_Value)]
  
  # Calculate SOS
  sos_threshold_exceeded <- rollapply(df$NDVI_Value, width = 2, FUN = function(x) all(x > 0.71056), partial = TRUE, align = "left")
  sos_index <- which(sos_threshold_exceeded == TRUE)[1]  # First occurrence
  if (!is.na(sos_index)) {
    sos <- df$Date[sos_index]
    sos_ndvi_value <- df$NDVI_Value[sos_index]
  }
  
  # Filtering data after PGT for EOS calculation
  df_post_PGT <- df %>% filter(Date > PGT)
  
  if (nrow(df_post_PGT) > 0) {
    eos_threshold_exceeded <- rollapply(df_post_PGT$NDVI_Value, width = 2, FUN = function(x) all(x < 0.70309), partial = TRUE, align = "left")
    eos_index <- which(eos_threshold_exceeded == TRUE)[1]  # First occurrence after PGT
    if (!is.na(eos_index)) {
      eos <- df_post_PGT$Date[eos_index]
      eos_ndvi_value <- df_post_PGT$NDVI_Value[eos_index]
    }
  }
  
  return(data.frame(SOS = sos, EOS = eos, APG = APG, PGT = PGT, sos_ndvi_value = sos_ndvi_value, eos_ndvi_value = eos_ndvi_value))
}

# Apply the function to each group
sos_eos_apg_results_v2 <- grouped %>%
  group_by(x, y, Year) %>%
  group_modify(~ calculate_sos_eos_apg_pgt_v2(.x)) %>%
  ungroup()

sos <- sos_eos_apg_results_v2


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
  ungroup()


############################## visualize results ###############################

str(sos)
mean(sos$sos_ndvi_value, na.rm=T)
mean(sos$eos_ndvi_value, na.rm=T)
################################################################################

clean_sos <- sos %>%
  group_by(x, y) %>%
  # Filter out groups with any NA values in any coloumn
  filter(!any(is.na(SOS), is.na(EOS), is.na(APG), is.na(PGT))) %>%
  ungroup()

print(paste("Number of excluded Pixel during EOS NA deleting filter:", nrow(sos)-nrow(clean_sos)))
nrow(clean_sos)
str(clean_sos)


mean(clean_sos$sos_ndvi_value, na.rm=T)
mean(clean_sos$eos_ndvi_value, na.rm=T)
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

# saveRDS(month_counts, "~/edu/NDVI_germany/docs/month")


filtered_sos <- clean_sos %>%
  filter(!(EOS_month %in% c(1, 2, 3, 4, 5, 12))) %>%
  filter(!(SOS_month %in% c(1, 7, 8, 9, 10, 11)))

head(filtered_sos)
# saveRDS(filtered_sos, "~/edu/NDVI_germany/docs/sos_eos")

mean(filtered_sos$sos_ndvi_value, na.rm=T)
mean(filtered_sos$eos_ndvi_value, na.rm=T)



################################################################################

clean_sos <- clean_sos %>%
  mutate(
    Year = as.numeric(Year),  # Convert Year to numeric
    SOS_julian = yday(SOS)  # Convert SOS to the day of the year
  )

# Perform linear regression
lm_result <- lm(SOS_julian ~ Year, data = clean_sos)

# Display the summary of the linear model
summary(lm_result)

str(clean_eos)
# EOS 
clean_eos <- clean_sos %>%
  mutate(
    Year = as.numeric(Year),
    EOS_julian =yday(EOS)
  )
lm_result <- lm(EOS_julian ~Year, data = clean_eos)
summary(lm_result)
################################################################################

clean_sos <- clean_sos %>%
  mutate(
    SOS_julian = yday(SOS),
    EOS_julian = yday(EOS)
  )

# Perform linear regression with EOS_julian as the dependent variable
lm_result <- lm(EOS_julian ~ SOS_julian, data = clean_sos)
# Display the summary of the linear model
summary(lm_result)


# Define the coordinates for the pixel group to be extracted with a small tolerance for floating-point comparison
x_coord <- 8.663542
y_coord <- 50.01146
tolerance <- 1e-5  # Define a tolerance level, e.g., 0.00001

# Filter the dataframe to extract the specified pixel group within the tolerance range
df_2 <- clean_sos %>%
  filter(abs(x - x_coord) < tolerance & abs(y - y_coord) < tolerance)
head(df_2)
plot(df_2$SOS_julian, df_2$EOS_julian)


cor.test(clean_sos$SOS_julian, clean_sos$EOS_julian)


# CALCULATION WITH BOTH ASPECTS: LOCAL AND GLOBAL THRESHOLD. 

