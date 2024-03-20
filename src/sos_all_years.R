# Load NDVI Data
source("~/NDVI_germany/src/NDVI_germany_setup.R")
#data <- readRDS(file.path(envrmt$path_data_level1, "value_march_to_october.rds"))
data <- readRDS(file.path(envrmt$path_data_level1, "long_ndvi_quality_df.rds"))

# Prepare date Sequences for Each Year
# This generates a sequence of dates for every year present in the dataset
#years <- unique(data$Year)
#date_sequences <- lapply(years, function(year) {
#  seq(ymd(paste(year, "01-01", sep = "-")), ymd(paste(year, "12-31", sep = "-")), by = "day")
#})

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
data$Quality_Binary = sapply(data$Quality, to_binary_string, bits = 15)
data$Quality <- NULL
head(data)


################################################################################

data <- data %>%
  dplyr::filter(
    str_sub(Quality_Binary, 1, 2) != "11" & # exclude all not produced observations
      str_sub(Quality_Binary, 1, 2) != "01" & # exclude all problematic observations 
      str_sub(Quality_Binary, 3, 4) != "11" & # delete all really bad pixel
      str_sub(Quality_Binary, 3, 4) != "10" # delete bad quality pixel
  )
colnames(data)[1] <- "date"
data$date <- as.Date(data$date)

detect_clouds <- function(Quality_Binary) {
  if (str_sub(Quality_Binary, 1, 2) == "10") { # Clouds
    return(1) # Most probably cloudy
  } else {
    return(0) # No clouds
  }
}

data <- data %>%
  mutate(Month = month(date)) %>%
  mutate(Year = year(date))


# Apply the function to create a new Clouds column
data$Clouds <- sapply(data$Quality_Binary, detect_clouds)
data$Clouds <- as.numeric(data$Clouds)


################################################################################
# Delete all pixel groups with years with too many clouds or too less observations
# Convert Year to Date

# Step 1: Filter pixel groups with at least 17 observations in a year
data <- data %>%
  group_by(x, y, Year) %>%
  dplyr::filter(n() >= 17) %>%
  ungroup()


# Step 2: Remove pixel groups with two or more consecutive cloud observations from March to November
data <- data %>%
  group_by(x, y) %>%
  mutate(Clouds_lag = lag(Clouds)) %>%
  dplyr::filter(!(Month >= 3 & Month <= 11 & Clouds == 1 & Clouds_lag == 1)) %>%
  dplyr::select(-c(Clouds_lag)) %>%
  ungroup()



################################################################################
# Step 1
# Replace NDVI_Value with NA for cloudy observations and create lin_int for linear interpolation
data <- data %>%
  group_by(x, y) %>%
  mutate(NDVI_Value = ifelse(Clouds == 1, NA, NDVI_Value)) %>% # Replace cloudy NDVI values with NA
  mutate(NDVI_Value = na.approx(NDVI_Value, na.rm = FALSE)) %>% # Perform linear interpolation
  ungroup()
data <- data %>%
  group_by(x, y) %>%
  dplyr::filter(n() >= 9) %>%
  ungroup()

data <- data %>%
  arrange(x, y, date) %>%
  group_by(x, y) %>%
  mutate(NDVI_Value = sgolayfilt(NDVI_Value, p = 3, n = 11)) %>%
  ungroup()

#data <- data %>%
#  mutate(Distance = abs(NDVI_Value - NDVI_Smoothed)) %>%
#  mutate(Weight = 1 / (1 + Distance)) # Simple weighting scheme; adjust as needed


#max_iterations <- 10 # Define a maximum number of iterations to prevent infinite loops
#convergence_threshold <- 0.01 # Define a threshold for determining convergence

#iterate_group <- function(df_group) {
#  max_iterations <- 10
#  convergence_threshold <- 0.01
#  prev_fitting_effect_index <- Inf
#  error_count <- 0 # Initialize error counter
#  
#  for (iteration in 1:max_iterations) {
#    # Wrap core logic in tryCatch for error handling
#    tryCatch({
#      # Apply Savitzky-Golay filter for refitting
#      df_group <- df_group %>%
#        mutate(NDVI_Refitted = sgolayfilt(NDVI_Value, p = 6, n = 9))
#      
#      # Recalculate weights
#      df_group <- df_group %>%
#        mutate(Distance = abs(NDVI_Value - NDVI_Refitted),
#               Weight = 1 / (1 + Distance))
#      
#      # Calculate fitting-effect index
#      fitting_effect_index <- mean(df_group$Distance)
#      
#      # Check for convergence
#      if (abs(prev_fitting_effect_index - fitting_effect_index) < convergence_threshold) {
#        break # Convergence achieved
#      }
#      prev_fitting_effect_index <- fitting_effect_index
#    }, error = function(e) {
#      # Increment error counter and print error message
#      error_count <<- error_count + 1
#      message("Error in iteration ", iteration, ": ", e$message)
#    })
#  }
#  
#  # Add error count to the dataframe as a new column (or however you wish to track it)
#  df_group$error_count <- error_count
#  return(df_group)
#}

#grouped_data_list <- data %>%
#  group_by(x, y) %>%
#  group_split()#

# Apply the iterative fitting process to each group
#processed_groups <- map(grouped_data_list, iterate_group)

# Combine the processed groups back into a single data frame
##data_processed <- bind_rows(processed_groups)
#data_processed$NDVI_Value <- data_processed$NDVI_Refitted
#data <- data_processed
#pixels_with_enough_data <- data %>%
#  dplyr::filter(Year == 2020) %>%
#  group_by(x, y) %>%
#  summarise(Entries = n()) %>%
#  dplyr::filter(Entries >= 10) %>%
#  ungroup()
#
# Then, randomly select one of these pixels
#random_pixel <- pixels_with_enough_data %>%
#  sample_n(1)
#
# Now, filter the original dataset for the randomly selected pixel and the year 2006
#data_for_random_pixel_2006 <- data %>%
#  dplyr::filter(x == random_pixel$x & y == random_pixel$y & Year == 2020)

# Proceed with your plotting or data manipulation
#ggplot(data_for_random_pixel_2006, aes(x = date)) +
#  geom_line(aes(y = NDVI_Value_golay, colour = "Golay Filtered NDVI")) +
#  geom_line(aes(y = NDVI_Value, colour = "Original NDVI")) +
#  labs(y = "NDVI", colour = "Legend") +
#  theme_minimal()

data <- data %>%
  mutate(NDVI_Value = if_else(NDVI_Value < 0.3, NA, NDVI_Value)) #%>%
  #group_by(x, y, Year) %>% # Group by pixel and year
  #mutate(
  #  Prev_NDVI = lag(NDVI_Value), # Get previous NDVI value within the same year
  #  Next_NDVI = lead(NDVI_Value), # Get next NDVI value within the same year
  #  Remove = (Prev_NDVI - NDVI_Value >= 0.1) & (Next_NDVI - NDVI_Value >= 0.1), # Condition to check within the same year
  #  NDVI_Value = ifelse(Remove, NA, NDVI_Value) # Mark sudden drops as NA within the same year
  #) %>%
  #dplyr::select(-Prev_NDVI, -Next_NDVI, -Remove) %>% # Remove auxiliary columns
  #ungroup() # Remove the grouping
data <- data[-which(is.na(data$NDVI_Value)),]

pixels_to_exclude <- data %>%
  dplyr::filter(Month %in% c(12, 1, 2), NDVI_Value > 0.65) %>% # Filter for high NDVI in winter
  group_by(x, y) %>%
  summarise(Num_Years = n_distinct(Year)) %>% # Count distinct years meeting the condition
  dplyr::filter(Num_Years >= 5) %>%
  dplyr::select(-Num_Years)

data <- data %>%
  anti_join(pixels_to_exclude, by = c("x", "y"))

df_monthly <- data %>%
  group_by(Year, x, y, Month) %>%
  summarise(Count = n(), .groups = 'drop') 

df_filtered <- df_monthly %>%
  group_by(Year, x, y) %>%
  dplyr::filter(all(c(3, 4, 5, 6, 7, 8, 9, 10) %in% Month)) %>%
  summarise(Months_with_data = n_distinct(Month), .groups = 'drop') 

data <- data %>%
  left_join(df_filtered, by = c("Year", "x", "y")) %>%
  dplyr::filter(!is.na(Months_with_data)) %>%
  dplyr::select(-Months_with_data)
# Adjusted pipeline with progress bar for the date expansion part
data_expanded <- data %>%
  dplyr::select(Year, x, y) %>%
  distinct() %>%
  ungroup() %>%
  mutate(date = pblapply(Year, function(year) {
    seq(ymd(paste0(year, "-01-01")), ymd(paste0(year, "-12-31")), by = "day")
  })) %>%
  unnest(date) %>%
  left_join(data, by = c("Year", "x", "y", "date"))
saveRDS(data_expanded, file.path(envrmt$path_data_level1, "data_expanded_save_15million.rds"))
data_expanded <- readRDS(file.path(envrmt$path_data_level1, "data_expanded_save_15million.rds"))
gc()

# Interpolate NDVI Values for Each Pixel
# Fills in missing NDVI values linearly for days without data
data_interpolated <- data_expanded %>%
  group_by(Year, x, y) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(NDVI_Value = approx(date, NDVI_Value, date, method = "linear")$y) %>%
  ungroup()
rm(data_expanded)
gc()
# Apply Fourier Smoothing Across All Years
# Smoothes NDVI values using Fourier terms and linear regression
#data_smoothed <- data_interpolated %>%
#  group_by(Year, x, y) %>%
#  do({
#    ndvi_values <- .$NDVI_Value
#    time_points <- seq_along(ndvi_values)
#    K <- 3  # Number of harmonics for Fourier transform
#    fourier_terms <- forecast::fourier(ts(ndvi_values, frequency = 365), K = K)
#    fit <- lm(ndvi_values ~ fourier_terms)
#    smoothed_values <- predict(fit, newdata = list(fourier_terms = fourier_terms))
#    data.frame(date = .$date, NDVI_Smooth = as.numeric(smoothed_values), stringsAsFactors = FALSE)
#  }) %>%
#  ungroup()
#rm(data_interpolated)
#gc()

pixel_min_max <- data_interpolated %>%
  group_by(x, y) %>%
  summarise(
    Global_Min_NDVI = min(NDVI_Value, na.rm = TRUE),
    Global_Max_NDVI = max(NDVI_Value, na.rm = TRUE)
  ) %>%
  ungroup()

data_normalized <- data_interpolated %>%
  left_join(pixel_min_max, by = c("x", "y")) %>%
  mutate(
    # Normalize NDVI values based on the global min and max for each pixel
    Norm_NDVI = (NDVI_Value - Global_Min_NDVI) / (Global_Max_NDVI - Global_Min_NDVI)
  ) %>%
  dplyr::select(-Global_Min_NDVI, -Global_Max_NDVI)
# Normalize Smoothed NDVI Values Across All Years
gc()
# Identify Start of Season for Each Pixel Across All Years
# The start of the season is defined as the first date where normalized NDVI exceeds 0.5 for 3 consecutive days
SoS_dates <- data_normalized %>%
  group_by(Year, x, y) %>%
  mutate(
    Above_Threshold_SoS = Norm_NDVI > 0.5,
    Lead1_SoS = lead(Above_Threshold_SoS, 1, default = FALSE),
    Lead2_SoS = lead(Above_Threshold_SoS, 2, default = FALSE),
    Lead3_SoS = lead(Above_Threshold_SoS, 3, default = FALSE),
    Start_Season_Flag = Above_Threshold_SoS & Lead1_SoS & Lead2_SoS & Lead3_SoS
  ) %>%
  summarise(
    Start_of_Season = min(if_else(Start_Season_Flag, date, as.Date(NA)), na.rm = TRUE)
  ) %>%
  dplyr::filter(!is.na(Start_of_Season))

peak_NDVI_dates <- data_normalized %>%
  group_by(Year, x, y) %>%
  # Identify the date of peak NDVI
  summarise(Peak_NDVI_Date = date[which.max(Norm_NDVI)]) %>%
  ungroup()

# Join peak NDVI dates back to the original dataset for reference
data_with_peak <- data_normalized %>%
  left_join(peak_NDVI_dates, by = c("Year", "x", "y"))

# Calculate End of Season based on conditions after peak NDVI
EoS_dates <- data_with_peak %>%
  group_by(Year, x, y) %>%
  mutate(
    # Ensure we only consider dates after the peak NDVI date for EoS calculation
    Post_Peak = date > Peak_NDVI_Date,
    Below_Threshold_EoS = Norm_NDVI < 0.5 & Post_Peak,
    Lead1_EoS = lead(Below_Threshold_EoS, 1, default = FALSE),
    Lead2_EoS = lead(Below_Threshold_EoS, 2, default = FALSE),
    Lead3_EoS = lead(Below_Threshold_EoS, 3, default = FALSE),
    End_Season_Flag = Below_Threshold_EoS & Lead1_EoS & Lead2_EoS & Lead3_EoS
  ) %>%
  summarise(
    # Find the earliest date after the peak where the End_Season_Flag is TRUE
    End_of_Season = if_else(any(End_Season_Flag), min(date[End_Season_Flag], na.rm = TRUE), as.Date(NA))
  ) %>%
  dplyr::filter(!is.na(End_of_Season)) %>%
  ungroup()


season_dates <- left_join(SoS_dates, EoS_dates, by = c("Year", "x", "y"))
season_dates <- season_dates %>%
  mutate(
    SoS_Day_of_Year = yday(Start_of_Season),
    EoS_Day_of_Year = yday(End_of_Season)
  )
cor.test(season_dates$SoS_Day_of_Year,season_dates$EoS_Day_of_Year)

season_dates <- season_dates %>%
  mutate(
    SoS_Month = month(Start_of_Season),
    EoS_Month = month(End_of_Season)
  )
# Summarize Start of Season by month
SoS_summary <- season_dates %>%
  group_by(SoS_Month = month(Start_of_Season, label = TRUE, abbr = FALSE)) %>%
  summarise(SoS_Count = n()) %>%
  ungroup()

# Summarize End of Season by month
EoS_summary <- season_dates %>%
  group_by(EoS_Month = month(End_of_Season, label = TRUE, abbr = FALSE)) %>%
  summarise(EoS_Count = n()) %>%
  ungroup()

# Combine SoS and EoS summaries into one dataframe
month_summary <- full_join(SoS_summary, EoS_summary, by = c("SoS_Month" = "EoS_Month")) %>%
  rename(Month = SoS_Month) %>%
  mutate(SoS_Count = replace_na(SoS_Count, 0),
         EoS_Count = replace_na(EoS_Count, 0)) %>%
  slice(c(1:11, 13:n(), 12))

ndvi_peak <- data_interpolated %>%
  group_by(Year, x, y) %>%
  summarise(Peak_NDVI = max(NDVI_Value, na.rm = TRUE),
            Date_of_Peak = date[which.max(NDVI_Value)],
            .groups = 'drop') # Ensure the dataframe is ungrouped after summarise

# Now, join ndvi_peak with season_dates to add the peak NDVI information
season_dates <- season_dates %>%
  left_join(ndvi_peak, by = c("Year", "x", "y"))

# View the first few rows of the updated season_dates dataframe
head(season_dates)

plot_data <- season_dates %>%
  mutate(
    SoS_Day_of_Year = yday(Start_of_Season),
    EoS_Day_of_Year = yday(End_of_Season)
  ) %>%
  # Reshape data for plotting
  pivot_longer(
    cols = c(SoS_Day_of_Year, EoS_Day_of_Year),
    names_to = "Season_Phase",
    values_to = "Day_of_Year"
  ) %>%
  # Remove NAs
  drop_na(Day_of_Year)

# Create the boxplot
ggplot(plot_data, aes(x = factor(Year), y = Day_of_Year, fill = Season_Phase)) +
  geom_boxplot() +
  labs(x = "Year", y = "Day of Year", title = "Start and End of Season by Year") +
  theme_minimal() +
  scale_fill_manual(values = c("SoS_Day_of_Year" = "blue", "EoS_Day_of_Year" = "red")) +
  theme(legend.title = element_blank(), legend.position = "bottom")



# For interactive map with mapview - assuming results_significant_eos is an sf object with a "Slope" column
mapview(results_significant_eos, 
        zcol = "Slope",  
        map.types = "OpenTopoMap")


saveRDS(data_normalized, file.path(envrmt$path_data_level1, "normalized_all_interpol_golay_3_11.rds"))
saveRDS(season_dates, file.path(envrmt$path_data_level1, "season_date_all_interpol_golay_3_11.rds"))
#season_dates <- readRDS("~/season_date.rds")
#season_dates<- readRDS(file.path(envrmt$path_data_level1, "season_date_with_filter.rds"))


