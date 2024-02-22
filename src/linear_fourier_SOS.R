################################################################################
# Get the whole dataframe
df <- readRDS("~/edu/NDVI_germany/data/raster_data/data_level1/long_ndvi_filtered.rds") 
################################################################################

# Create a test dataset to calculate just with xx%
unique_pixel_groups <- df %>%
  distinct(x, y)

# select xx% of these pixel groups
set.seed(12) # Setting seed for reproducibility
selected_pixel_groups <- unique_pixel_groups %>%
  sample_frac(0.05)

# Filter the original dataframe to only include data from the selected pixel groups
selected_df <- df %>%
  semi_join(selected_pixel_groups, by = c("x", "y"))

df <- selected_df # calculate still with df

str(df)
################################################################################
# TODO: Filter pixel groups with less than 4 observations per season
# Function to determine the season based on the month
# This is a good previous filter. The second problematic group filter from the 
# linear daily interpolation kicked 0.975% less groupes out because of this filter here!


# Your existing functions and code
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

df$Season <- sapply(month(df$Date), getSeason)

Year <- as.numeric(format(df$Date, "%Y"))

df$SeasonYear <- ifelse(
  df$Season == "Winter" & month(df$Date) %in% c(11, 12),
  paste(Year, Year + 1, sep="-"),
  ifelse(
    df$Season == "Winter" & month(df$Date) == 1,
    paste(Year - 1, Year, sep="-"),
    as.character(Year)
  )
)

# Adjusted filtering approach
df_filtered <- df %>%
  group_by(x, y, Season, SeasonYear) %>%
  filter(!(Season == "Winter" & (SeasonYear == "2000-2001" | SeasonYear == "2022-2023"))) %>%
  filter(n() >= 4 | (Season == "Winter" & (SeasonYear == "2000-2001" | SeasonYear == "2022-2023"))) %>%
  ungroup()

print(paste("Number of excluded observations during season filter:", nrow(df)-nrow(df_filtered)))


# The df_filtered dataframe now should contain only the pixel groups 
# with more than four observations per season, excluding specific conditions for Winter 2000 and Winter 2023.
df <- df_filtered

############# Detect clouds ################################################
# TODO: Detect mixed clouds and filter pixel
detect_clouds <- function(Quality_Binary) {
  if (str_sub(Quality_Binary, 11, 11) == "1") { # Mixed clouds
    return(0) # No mixed clouds
  } else {
    return(1) # Mixed Clouds present
  }
}

# Apply the function to create a new Clouds column
df$Clouds <- sapply(df$Quality_Binary, detect_clouds)
str(df)
df$Clouds <- as.numeric(df$Clouds)

################################################################################
############################ DELETE OBSERVATIONDROPS ###########################
################################################################################
# Step 1: Adjust NDVI values based on specified conditions

result_adjusted <- df %>%
  group_by(x, y) %>%
  # Step 1: Set NDVI values below 0.4 to NA
  mutate(NDVI_Value_Adjusted = if_else(NDVI_Value < 0.4, NA_real_, NDVI_Value)) %>%
  # Step 2: Apply logic for sudden drop, ignoring NA values
  mutate(
    Prev_Value = lag(na.locf(NDVI_Value_Adjusted, na.rm = FALSE)),
    Next_Value = lead(na.locf(NDVI_Value_Adjusted, na.rm = FALSE)),
    NDVI_Value_Adjusted = if_else(NDVI_Value_Adjusted < Prev_Value - 0.1 & NDVI_Value_Adjusted < Next_Value - 0.1, NA_real_, NDVI_Value_Adjusted)
  ) %>%
  # Identify first 4 and last 4 rows in each group
  mutate(
    Row = row_number(),
    Max_Row = max(row_number()),
    Adjust_NA = if_else((Row <= 4 | Row > (Max_Row - 4)) & is.na(NDVI_Value_Adjusted), 0.5000, NDVI_Value_Adjusted)
  ) %>%
  # Apply the NA adjustment rule for edge cases
  mutate(NDVI_Value_Adjusted = coalesce(Adjust_NA, NDVI_Value_Adjusted)) %>%
  # Remove the temporary columns used for calculations
  select(-Prev_Value, -Next_Value, -Row, -Max_Row, -Adjust_NA) %>%
  ungroup()


################################################################################
# Step 2: Interpolate adjusted NDVI values
result_interpolated <- result_adjusted %>%
  group_by(x, y) %>%
  mutate(NDVI_Value_Interpolated = na.approx(NDVI_Value_Adjusted, na.rm = FALSE)) %>%
  ungroup()

################################################################################
######################## SAVITZKY-GOLAY FILTER #################################
################################################################################
# Step 3: Apply the Savitzky-Golay filter to the interpolated NDVI time series

# Initialize a list to store dataframes of excluded groups
problematic_groups_list <- list()

# Update the wrapper function for sgolayfilt to catch errors and save the entire group's data
safe_sgolayfilt <- function(data, p, n) {
  tryCatch({
    data$NDVI_SG_Filtered <- sgolayfilt(data$NDVI_Value_Interpolated, p, n)
    return(data)
  }, error = function(e) {
    # Save the entire problematic group's data
    problematic_groups_list <<- c(problematic_groups_list, list(data))
    data$NDVI_SG_Filtered <- NA # Assign NA to NDVI_SG_Filtered to keep the same structure
    return(data)
  })
}

# Apply the function to your dataframe, grouping by x and y
result_smoothed <- result_interpolated %>%
  group_by(x, y) %>%
  group_modify(~ safe_sgolayfilt(.x, p = 3, n = 5)) %>%
  ungroup()

# Combine all problematic groups into one dataframe
problematic_groups_df <- bind_rows(problematic_groups_list)

# Print the number of excluded groups
print(paste("Number of excluded pixel groups during sgolayfilter:", length(problematic_groups_list)))

str(result_smoothed)

############################### EXCLUDE GREEN WINTER ###########################

# Step 1: Filter for Winter season and NDVI values above 0.7
winter_high_ndvi <- result_smoothed %>%
  filter(Season == "Winter" & NDVI_Value > 0.7)

# Step 2 & 3: Group by pixel coordinates, then count unique years meeting the condition
pixel_groups_high_ndvi <- winter_high_ndvi %>%
  group_by(x, y) %>%
  summarize(Unique_Winter_Years = n_distinct(Year)) %>%
  filter(Unique_Winter_Years >= 5) %>%
  ungroup()

# Step 4: Create green_winter dataframe to INCLUDE these pixel groups
# This dataframe will include only the records from pixel groups that met the exclusion criteria
green_winter <- result_smoothed %>%
  semi_join(pixel_groups_high_ndvi, by = c("x", "y"))

# Step 5: Update result_smoothed to EXCLUDE these pixel groups
# This dataframe will now exclude the pixel groups that met the criteria
result_smoothed_updated <- result_smoothed %>%
  anti_join(pixel_groups_high_ndvi, by = c("x", "y"))

# Output the updated result_smoothed dataframe without the excluded pixel groups
result_smoothed <- result_smoothed_updated

# Print the first few rows of each dataframe to verify
print(paste("Number of excluded observations during green winter filter:", nrow(result_smoothed)-nrow(green_winter)))


################################################################################
######################## fourier smoothing #####################################
######################## prepare data ##########################################

# TODO: Create a sequence of daily dates covering the range in your data
seq_dates <- seq(from = min(result_smoothed$Date), to = max(result_smoothed$Date), by = "day")
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
grouped <- result_smoothed %>%
  group_by(x, y, Year) %>%
  group_modify(~ {
    # Create a sequence of daily dates for the current year
    start_date <- min(.x$Date)
    end_date <- max(.x$Date)
    seq_dates <- seq(from = start_date, to = end_date, by = "day")
    
    # Interpolate NDVI values for these dates using the safe wrapper function
    daily_ndvi <- safe_interpolate(.x$Date, .x$NDVI_SG_Filtered, seq_dates)
    
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

# Now, 'grouped_cleaned' contains data without any pixel groups that had NA in NDVI_Value

################################################################################

# grouped$NDVI_SG_Filtered <- sgolayfilt(grouped$NDVI_Value, p = 3, n = 31)

# Note: Ensure n is odd and > p; adjust p and n as needed based on experimentation

############################# create Fourier ###################################
smooth_with_fourier <- function(ndvi_values, frequency = 365, components_to_keep = 6) {
  if(length(ndvi_values) < 2) {  # Ensure there are at least two rows
    return(rep(NA, length(ndvi_values)))
  }
  ndvi_ts <- ts(ndvi_values, frequency = frequency)
  fft_result <- fft(ndvi_ts)
  n <- length(fft_result)
  if(components_to_keep < n / 2) {
    fft_result[(components_to_keep + 1):(n - components_to_keep)] <- 0
  }
  smoothed_ndvi <- Re(fft(fft_result, inverse = TRUE)) / n
  return(smoothed_ndvi)
}

# Apply the smoother to each group
 grouped_smoothed <- grouped %>%
  group_by(x, y, Year) %>%
  mutate(NDVI_Value_Smoothed = smooth_with_fourier(NDVI_Value, 365, 6))


################################################################################
########################## CALCULATE SOS EOS APG ###############################
################################################################################

calculate_sos_eos_apg_pgt <- function(df) {
  if(nrow(df) == 0 || all(is.na(df$NDVI_Value_Smoothed))) {
    return(data.frame(SOS = NA, EOS = NA, APG = NA, PGT = NA)) # Early return for empty groups or all NA NDVI
  }
  
  NDVImax <- max(df$NDVI_Value_Smoothed, na.rm = TRUE)
  NDVImin <- min(df$NDVI_Value_Smoothed, na.rm = TRUE)
  
  # Additional guard clause for cases where NDVImax or NDVImin might still be Inf/-Inf due to all NAs
  if(is.infinite(NDVImax) || is.infinite(NDVImin)) {
    return(data.frame(SOS = NA, EOS = NA, APG = NA, PGT = NA))
  }
  
  # Proceed with existing calculations...
  df$NDVIratio <- (df$NDVI_Value_Smoothed - NDVImin) / (NDVImax - NDVImin)
  
  APG <- NDVImax
  PGT <- df$Date[which.max(df$NDVI_Value_Smoothed)]
  
  df$above_threshold <- df$NDVIratio >= 0.5
  run_length <- rle(df$above_threshold)
  valid_starts <- which(run_length$values & run_length$lengths >= 2)
  
  sos <- if(length(valid_starts) > 0) {
    start_index <- sum(run_length$lengths[1:(valid_starts[1]-1)]) + 1
    df$Date[start_index]
  } else {
    NA
  }
  
  # Filtering data after PGT for EOS calculation
  df_post_PGT <- df %>% filter(Date > PGT)
  if(nrow(df_post_PGT) > 0) {
    df_post_PGT$below_threshold <- df_post_PGT$NDVIratio <= 0.5
    run_length_eos <- rle(df_post_PGT$below_threshold)
    valid_ends <- which(run_length_eos$values & run_length_eos$lengths >= 4)
    
    if(length(valid_ends) > 0) {
      # Find the starting index of the first valid end period
      end_index_start <- sum(run_length_eos$lengths[1:(valid_ends[1]-1)]) + 1
      eos <- df_post_PGT$Date[end_index_start]
    } else {
      eos <- NA
    }
  } else {
    eos <- NA
  }
  
  return(data.frame(SOS = sos, EOS = eos, APG = APG, PGT = PGT))
}


# TODO: Apply the function to the filtered groups
sos_eos_apg_results <- grouped_smoothed %>%
  group_by(x, y, Year) %>%
  group_modify(~ calculate_sos_eos_apg_pgt(.)) %>%
  ungroup()

############################## visualize results ###############################

str(sos_eos_apg_results)
print(sos_eos_apg_results)

################################################################################

clean_sos_eos_apg_results <- sos_eos_apg_results %>%
  group_by(x, y) %>%
  # Filter out groups with any NA values in any column
  filter(!any(is.na(SOS), is.na(EOS), is.na(APG), is.na(PGT))) %>%
  ungroup()

print(paste("Number of excluded Pixel during EOS NA deleting filter:", nrow(sos_eos_apg_results)-nrow(clean_sos_eos_apg_results)))
nrow(clean_sos_eos_apg_results)
str(clean_sos_eos_apg_results)

# Extract month from EOS date
clean_sos_eos_apg_results <- clean_sos_eos_apg_results %>%
  mutate(EOS_month = format(EOS, "%m"), # Extract month as a string
         EOS_month = as.integer(EOS_month)) # Convert to integer for ordering

# Extract month from SOS date
clean_sos_eos_apg_results <- clean_sos_eos_apg_results %>%
  mutate(SOS_month = format(SOS, "%m"), # Extract month as a string
         SOS_month = as.integer(SOS_month)) # Convert to integer for ordering

# Extract month from PGT date
clean_sos_eos_apg_results <- clean_sos_eos_apg_results %>%
  mutate(PGT_month = format(PGT, "%m"), # Extract month as a string
         PGT_month = as.integer(PGT_month)) # Convert to integer for ordering


sos_month_counts <- clean_sos_eos_apg_results %>%
  count(SOS_month) %>%
  rename(Month = SOS_month, SOS_Count = n)

# Count the number of occurrences for each EOS month
eos_month_counts <- clean_sos_eos_apg_results %>%
  count(EOS_month) %>%
  rename(Month = EOS_month, EOS_Count = n)

pgt_month_counts <- clean_sos_eos_apg_results %>%
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

saveRDS(month_counts, "~/edu/NDVI_germany/docs/month_count_fourier_smoother_0.05")


df_f <- readRDS("~/edu/NDVI_germany/docs/month_count_fourier_smoother_0.05")
df_f
df_s <- readRDS("~/edu/NDVI_germany/docs/month_count_Savitzky_0.05")
df_s
df_n <- readRDS("~/edu/NDVI_germany/docs/month_count_no_daily_smoother_0.05")
df_n
str(df_n)

total_eos_count_f <- sum(df_f$EOS_Count)
total_eos_count_s <- sum(df_s$EOS_Count)
total_eos_count_n <- sum(df_n$EOS_Count)
# Print the total count of EOS_Count values
print(paste("Total EOS_Count values:", total_eos_count_f, total_eos_count_s, total_eos_count_n))


filtered_sos_eos_apg_results <- clean_sos_eos_apg_results %>%
  filter(!(EOS_month %in% c(1, 2, 3, 4, 5, 12))) %>%
  filter(!(SOS_month %in% c(1, 7, 8, 9, 10, 11)))
  
filtered_sos_eos_apg_results
saveRDS(filtered_sos_eos_apg_results, "~/edu/NDVI_germany/docs/sos_eos_cleaned_filtered_fourier")

sos_month_counts <- filtered_sos_eos_apg_results %>%
  count(SOS_month) %>%
  rename(Month = SOS_month, SOS_Count = n)

# Count the number of occurrences for each EOS month
eos_month_counts <- filtered_sos_eos_apg_results %>%
  count(EOS_month) %>%
  rename(Month = EOS_month, EOS_Count = n)

pgt_month_counts <- filtered_sos_eos_apg_results %>%
  count(PGT_month) %>%
  rename(Month = PGT_month, PGT_Count = n)

month_counts <- sos_month_counts %>%
  full_join(eos_month_counts, by = "Month") %>%
  full_join(pgt_month_counts, by = "Month")

# Replace NA with 0 for counts that don't appear in some datasets
month_counts[is.na(month_counts)] <- 0

# Print the result to check
print(month_counts)
