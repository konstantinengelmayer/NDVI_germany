# This script calculates the SOS (Start of season) and the EOS (End of Season)
# For this, we need to find the one day in every year where a chosen threshold
# is above (spring) or below (autumn) the value.

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

# Filter the good values
result_df <- df %>%
  filter(
    str_sub(Quality_Binary, 1, 2) == "00" & 
      str_sub(Quality_Binary, 3, 4) != "11" &
      str_sub(Quality_Binary, 3, 4) != "1" &
      str_sub(Quality_Binary, 8, 8) == "0" & 
      str_sub(Quality_Binary, 10, 10) == "0" &
      str_sub(Quality_Binary, 14, 14) == "0" &
      str_sub(Quality_Binary, 15, 15) == "0" 
  )

tail(result_df)
nrow(result_df) # 11874448

result_df$Quality_Binary <- NULL # Deleting the Quality bits
result_df$Quality <- NULL # Deleting the coloumn with the Quality

print(colnames(result_df)) # Get the names = we should rename the Year coloumn to the Date coloumn.
colnames(result_df)[1] <- "Date" # This makes the next steps more intuitive: Date and Year coloumn
result_df$Year <- year(result_df$Date)
head(result_df) # Have a last look

# Safe and read the new file with all good values
saveRDS(result_df, file="~/edu/NDVI_germany/data/raster_data/data_level1/long_ndvi_quality_bits_good_df.rds")




################################################################################
# TODO: Load the dataframe
# Note: The smoother value affects the result of SOS and EOS!
df <- readRDS("~/edu/NDVI_germany/data/raster_data/data_level1/long_ndvi_quality_bits_good_df.rds")# Read the RDS file
head(df)


test <- subset(df, Year==2020) # create a test dataset to calculate the first steps faster
str(test) # How is the structure?
unique(is.na(test)) # No NA left
tail(test)

#####################
365/16 # 22.8 Observations per year possible
#####################

# Count NDVI observations per pixel per year
observations_count <- test %>%
  group_by(x, y, Year) %>%
  summarise(Observations = n()) %>%
  ungroup()  # Ensure the result is no longer grouped

# View the first few rows of the counts
print(head(observations_count))

# Visualize the number of observations per pixel per year
ggplot(observations_count, aes(x = factor(Year), y = Observations)) +
  geom_boxplot() +
  labs(title = "NDVI Observations per Pixel per Year",
       x = "Year",
       y = "Number of Observations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Improve x-axis labels readability

###############################

# Assuming the first pixel's coordinates you're interested in are from the first row
first_x <- unique(test$x)[1]
first_y <- unique(test$y)[1]
first_x
first_y

# Filter the dataframe for this specific pixel across all years
specific_pixel_data <- test %>%
  filter(x == first_x, y == first_y)

# View the result
specific_pixel_data

specific_pixel_data$Date <- as.Date(specific_pixel_data$Date)
specific_pixel_data$TimeNumeric <- as.numeric(specific_pixel_data$Date - min(specific_pixel_data$Date))
specific_pixel_data

loess_fit <- loess(NDVI_Value ~ TimeNumeric, data = specific_pixel_data, span = 0.4)


# Create a sequence of days for prediction
prediction_days <- seq(from = min(specific_pixel_data$TimeNumeric), 
                       to = max(specific_pixel_data$TimeNumeric), 
                       by = 1)

predicted_ndvi <- predict(loess_fit, newdata = data.frame(TimeNumeric = prediction_days))

# Define the NDVI threshold for SOS (e.g., 0.5 or another value based on your criteria)
ndvi_threshold <- 0.7

# Find the first day when NDVI exceeds the threshold
sos_index <- which(predicted_ndvi > ndvi_threshold)[1]
sos_date <- min(specific_pixel_data$Date) + sos_index - 1  # Adjusting for the start at 0

# Print the SOS date
if (!is.na(sos_index)) {
  print(paste("Start of Season (SOS):", sos_date))
} else {
  print("No Start of Season (SOS) found based on the threshold.")
}



######################## fourier smoothing #############################
######################## prepare data ##################################

head(test)
str(test)

# Convert Date to Date class if it's not already
test$Date <- as.Date(test$Date) # It is not!

# TODO: Create a sequence of daily dates covering the range in your data
seq_dates <- seq(from = min(test$Date), to = max(test$Date), by = "day")
seq_dates


# TODO: Group the data and create for every day a daily NDVI-value
test_grouped <- test %>%
  group_by(x, y, Year) %>%
  group_modify(~ {
    # Ensure the Date column is correctly formatted as Date
    .x$Date <- as.Date(.x$Date)
    
    # Check for sufficient non-NA NDVI values
    if (sum(!is.na(.x$NDVI_Value)) < 2) {
      return(data.frame(Date = .x$Date, NDVI_Value = rep(NA, length(.x$Date))))
    }
    
    # Create a sequence of daily dates for the current year
    start_date <- min(.x$Date)
    end_date <- max(.x$Date)
    seq_dates <- seq(from = start_date, to = end_date, by = "day")
    
    # Prepare for interpolation: Ensure dates are in numeric format for `approx`
    numeric_dates <- as.numeric(.x$Date)
    numeric_seq_dates <- as.numeric(seq_dates)
    
    # Interpolate NDVI values for these numeric dates
    daily_ndvi <- approx(x = numeric_dates, y = .x$NDVI_Value,
                         xout = numeric_seq_dates, method = "linear", rule = 2)$y
    
    # Construct the result dataframe
    data.frame(Date = seq_dates, NDVI_Value = daily_ndvi)
  })

print(test_grouped, n=50)
str(test_grouped)



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
test_grouped_smoothed <- test_grouped %>%
  group_by(x, y, Year) %>%
  group_modify(~ smooth_with_fourier(.))

print(test_grouped_smoothed, n=80)
str(test_grouped_smoothed)



# TODO: Filter pixel groups with observations in the month March, April, Mai, Agust and Septembre

test_grouped_filtered <- test_grouped_smoothed %>%
  group_by(x, y, Year) %>%
  # Ensure at least one observation in each of the specified months
  mutate(month = month(Date)) %>%
  filter(month %in% c(3, 4, 5, 6, 7, 8, 9, 10)) %>%
  # Summarize to check the presence of each month within groups
  group_by(x, y, Year, month) %>%
  summarise(n = n(), .groups = "drop") %>%
  # Ensure all required months are present within the group
  group_by(x, y, Year) %>%
  filter(all(c(3, 4, 5, 6, 7, 8, 9, 10) %in% month)) %>%
  # Join back to the original smoothed data to get the full dataset
  select(-n) %>%
  distinct(x, y, Year) %>%
  left_join(test_grouped_smoothed, by = c("x", "y", "Year")) %>%
  ungroup()


head(test_grouped_filtered)
str(test_grouped_filtered)
nrow(test_grouped_filtered)
boxplot(test_grouped_filtered$NDVI_Value_Smoothed)
print(test_grouped_filtered, n=900)




# First, ensure the Date column is of type Date if not already
test_grouped_filtered$Date <- as.Date(test_grouped_filtered$Date)

# Create a new column for Month
test_grouped_filtered$Month <- format(test_grouped_filtered$Date, "%m")

# Select the first 20 unique pixel groups
unique_pixel_groups <- test_grouped_filtered %>% 
  distinct(x, y) %>% 
  slice(1:20) # Adjust this number as needed

# Filter the main dataset for these first 20 pixel groups
filtered_for_plot <- test_grouped_filtered %>%
  semi_join(unique_pixel_groups, by = c("x", "y"))

# Plotting
ggplot(filtered_for_plot, aes(x = Month, y = NDVI_Value_Smoothed, group = interaction(x, y), color = as.factor(interaction(x, y)))) +
  geom_line() + 
  geom_point(size = 0.5) +
  theme_minimal() +
  labs(title = "Smoothed NDVI Values Across Months for First 20 Pixel Groups",
       x = "Month",
       y = "Smoothed NDVI Value",
       color = "Pixel Group") +
  theme(legend.position = "none") # Hide legend for clarity




################################### calculate ##################################

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
      end_index <- sum(run_length_eos$lengths[1:(valid_ends[length(valid_ends)])])
      eos <- df_post_PGT$Date[end_index]
    } else {
      eos <- NA
    }
  } else {
    eos <- NA
  }
  
  return(data.frame(SOS = sos, EOS = eos, APG = APG, PGT = PGT))
}


# TODO: Apply the function to the filtered groups
sos_eos_apg_results <- test_grouped_filtered %>%
  group_by(x, y, Year) %>%
  group_modify(~ calculate_sos_eos_apg_pgt(.)) %>%
  ungroup()

############################## visualize results ###############################

str(sos_eos_apg_results)
print(sos_eos_apg_results)
boxplot(sos_eos_apg_results$SOS)
boxplot(sos_eos_apg_results$EOS)
boxplot(sos_eos_apg_results$APG)
boxplot(sos_eos_apg_results$PGT)
hist(sos_eos_apg_results$SOS, breaks=11)
hist(sos_eos_apg_results$EOS, breaks=11)
hist(sos_eos_apg_results$APG, breaks=11)
