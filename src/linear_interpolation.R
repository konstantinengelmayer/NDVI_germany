################################################################################
df <- readRDS("~/edu/NDVI_germany/data/raster_data/data_level1/long_ndvi_filtered.rds") 
################################################################################
set.seed(99)
unique_pixel <- df %>% distinct(x, y) %>% sample_n(1)
# Filter the main dataframe to include only the selected pixel group
df_random_group <- df %>%
  filter(x == unique_pixel$x & y == unique_pixel$y)
nrow(df_random_group)
df <- df_random_group
################################################################################
# Get a problematic pixel group:

# Define the coordinates for the pixel group to be extracted with a small tolerance for floating-point comparison
x_coord <- 8.842708
y_coord <- 49.39687
tolerance <- 1e-5  # Define a tolerance level, e.g., 0.00001

# Filter the dataframe to extract the specified pixel group within the tolerance range
df_1 <- df %>%
  filter(abs(x - x_coord) < tolerance & abs(y - y_coord) < tolerance)

# Display the first few rows of df_1 to confirm extraction
head(df_1)
df <- df_1
################################################################################
############# Detect clouds/ice ################################################

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
head(df)
str(df)

df$Clouds <- as.numeric(df$Clouds)  # Ensuring Clouds is numeric


################################################################################
# Step 0: Adjust NDVI values based on specified conditions

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
  # Remove the temporary columns used for calculations
  select(-Prev_Value, -Next_Value) %>%
  ungroup()

############# SECOND TRY WITH THE FILLING OF THE FIRST AND LAST FOUR ROWS ######
library(dplyr)
library(zoo)

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

# Step 1: Interpolate adjusted NDVI values
result_interpolated <- result_adjusted %>%
  group_by(x, y) %>%
  mutate(NDVI_Value_Interpolated = na.approx(NDVI_Value_Adjusted, na.rm = FALSE)) %>%
  ungroup()
head(result_interpolated)
str(result_interpolated)
result_interpolated <- na.omit(result_interpolated)



lambda <- 10 # Smoothing parameter; adjust based on data
result_interpolated$NDVI_Whittaker <- with(result_interpolated, {
  zoo::rollapply(NDVI_Value_Interpolated, width = 7, FUN = function(z) smooth.spline(z, lambda = lambda)$y, 
                 by.column = FALSE, fill = NA, align = "center")
})


# Step 2: Apply the Savitzky-Golay filter to the interpolated NDVI time series
result_smoothed <- result_interpolated %>%
  group_by(x, y) %>%
  mutate(NDVI_SG_Filtered = sgolayfilt(NDVI_Value_Interpolated, p = 3, n = 5)) %>%
  ungroup()


str(result_smoothed)
head(result_smoothed)


################## PLOT ########################################################
# Filter the dataframe for the year xxxx
result_smoothed_2xxx <- result_smoothed %>%
  filter(as.Date(Date) >= as.Date("2005-01-01") & as.Date(Date) <= as.Date("2005-12-31"))

# Plotting
ggplot(result_smoothed_2xxx, aes(x = Date)) +
  geom_line(aes(y = NDVI_Value, color = "Original NDVI"), linewidth=1) +
  geom_line(aes(y = NDVI_Value_Interpolated, color = "Interpolated NDVI"), linewidth=1) +
  geom_line(aes(y = NDVI_SG_Filtered, color = "Smoothed NDVI"), linewidth=1) +
  scale_color_manual(values = c("Original NDVI" = "blue", 
                                "Interpolated NDVI" = "green", 
                                "Smoothed NDVI" = "red")) +
  labs(title = "NDVI Values through one Year", 
       y = "NDVI Value", 
       color = "Legend") +
  theme_minimal()


################################################################################
######################## fourier smoothing #####################################
######################## prepare data ##########################################

str(result_smoothed)

# TODO: Create a sequence of daily dates covering the range in your data
seq_dates <- seq(from = min(result_smoothed$Date), to = max(result_smoothed$Date), by = "day")
head(seq_dates)


# TODO: Group the data and create for every day a daily NDVI-value
grouped <- result_smoothed %>%
  group_by(x, y, Year) %>%
  group_modify(~ {
    # Create a sequence of daily dates for the current year
    start_date <- min(.x$Date)
    end_date <- max(.x$Date)
    seq_dates <- seq(from = start_date, to = end_date, by = "day")
    
    # Prepare for interpolation: Ensure dates are in numeric format for `approx`
    numeric_dates <- as.numeric(.x$Date)
    numeric_seq_dates <- as.numeric(seq_dates)
    
    # Interpolate NDVI values for these numeric dates
    daily_ndvi <- approx(x = numeric_dates, y = .x$NDVI_SG_Filtered,
                         xout = numeric_seq_dates, method = "linear", rule = 2)$y
    
    # Construct the result dataframe
    data.frame(Date = seq_dates, NDVI_Value = daily_ndvi)
  })

print(grouped, n=20)


############################# create ###########################################

# TODO: Create the Fourier Smoother
smooth_with_fourier <- function(df) {
  if(nrow(df) < 2) {  # Ensure there are at least two rows to perform ts and fft
    df$NDVI_Value_Smoothed <- NA
    return(df)
  }
  
  # Assuming 'df' contains 'Date' and 'NDVI_Value' columns
  ndvi_ts <- ts(df$NDVI_Value, frequency = 365)  # Adjust frequency as needed
  
  # Apply Fourier transform
  fft_result <- fft(ndvi_ts)
  
  # Determine the number of components to keep
  components_to_keep <- 6  # Adjust based on your data
  
  # Safely zero out the middle coefficients, if possible
  n <- length(fft_result)
  if(components_to_keep < n / 2) {
    fft_result[(components_to_keep + 1):(n - components_to_keep)] <- 0
  }
  
  # Inverse FFT to get the smoothed series
  smoothed_ndvi <- Re(fft(fft_result, inverse = TRUE)) / length(fft_result)
  
  # Add smoothed NDVI values to the dataframe
  df$NDVI_Value_Smoothed <- smoothed_ndvi
  
  return(df)
}


# TODO: Apply the smoother: Fourier
grouped_smoothed <- grouped %>%
  group_by(x, y, Year) %>%
  group_modify(~ smooth_with_fourier(.))

print(grouped_smoothed, n=20)

# Filter for the year 2012
grouped_smoothed_2012 <- grouped_smoothed %>%
  filter(Year == "2012")

ggplot(grouped_smoothed_2012, aes(x = Date)) +
  geom_line(aes(y = NDVI_Value, color = factor(month(Date))), size = 1) +
  geom_line(aes(y = NDVI_Value_Smoothed, color = factor(month(Date))), size = 1, linetype = "dashed") +
  scale_color_viridis_d(name = "Month", guide = guide_legend(title = "Month")) +
  labs(title = "NDVI and Smoothed NDVI Values in 2012", y = "NDVI Value") +
  theme_minimal() +
  theme(legend.position = "bottom")


################################################################################
################################### calculate ##################################
# Without autumn minimum
calculate_sos_eos_apg_pgt <- function(df) {
  if(nrow(df) == 0) return(data.frame(SOS = NA, EOS = NA, APG = NA, PGT = NA)) # Early return for empty groups
  
  NDVImax <- max(df$NDVI_Value_Smoothed, na.rm = TRUE)
  NDVImin <- min(df$NDVI_Value_Smoothed, na.rm = TRUE)
  
  # Guard clause for infinite values
  if(is.infinite(NDVImax) || is.infinite(NDVImin)) {
    return(data.frame(SOS = NA, EOS = NA, APG = NA, PGT = NA))
  }
  
  # Calculate NDVI ratio for the entire dataset
  df$NDVIratio <- (df$NDVI_Value_Smoothed - NDVImin) / (NDVImax - NDVImin)
  
  # Calculate APG (Annual Peak Growth) as the maximum smoothed NDVI value
  APG <- NDVImax
  # Calculate PGT (Peak Growth Time) as the date corresponding to APG
  PGT <- df$Date[which.max(df$NDVI_Value_Smoothed)]
  
  # Logic to calculate SOS considering two consecutive days above threshold
  df$above_threshold <- df$NDVIratio >= 0.5
  run_length <- rle(df$above_threshold)
  valid_starts <- which(run_length$values & run_length$lengths >= 2)
  
  if(length(valid_starts) > 0) {
    start_index <- sum(run_length$lengths[1:(valid_starts[1]-1)]) + 1
    sos <- df$Date[start_index]
  } else {
    sos <- NA
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



