source("~/edu/NDVI_germany/src/NDVI_germany_setup.R")
data <- readRDS(file.path(envrmt$path_data_level1, "value_for_12month.rds"))
data <- data[which(data$Year ==2005),]
date_sequence <- seq(ymd("2005-01-01"), ymd("2005-12-31"), by = "day")

# Expand your dataset to include all dates in 2005 for each pixel
data_expanded <- data %>%
  dplyr::select(x, y) %>%
  distinct() %>%
  expand_grid(Date = date_sequence) %>%
  left_join(data, by = c("x", "y", "Date"))
data_expanded <- data_expanded[,-c(5,6)]

# Interpolate NDVI values for each pixel
data_interpolated <- data_expanded %>%
  group_by(x, y) %>%
  arrange(Date, .by_group = TRUE) %>%
  mutate(NDVI_Value = approx(Date, NDVI_Value, Date, method = "linear")$y) %>%
  ungroup()

data_interpolated$DayOfYear <- yday(data_interpolated$Date)

# Apply Fourier smoothing for each pixel
data_smoothed <- data_interpolated %>%
  group_by(x, y) %>%
  do({
    # Extract the NDVI values and corresponding DayOfYear as numeric
    ndvi_values <- .$NDVI_Value
    days <- .$DayOfYear
    
    # Generate Fourier terms for the time series
    K <- 3  # Number of harmonics
    fourier_terms <- forecast::fourier(ts(ndvi_values, frequency = 365), K = K)
    
    # Fit a linear model using the Fourier terms as predictors
    fit <- lm(ndvi_values ~ fourier_terms)
    
    # Predict NDVI values using the fitted model (This smooths the NDVI series)
    smoothed_values <- predict(fit, newdata = list(fourier_terms = fourier_terms))
    
    # Return a data frame with Date and smoothed NDVI values
    data.frame(Date = .$Date, NDVI_Smooth = as.numeric(smoothed_values), stringsAsFactors = FALSE)
  }) %>%
  ungroup()

# Now, proceed with normalizing the smoothed NDVI values
data_normalized <- data_smoothed %>%
  group_by(x, y) %>%
  mutate(
    Min_NDVI = min(NDVI_Smooth, na.rm = TRUE),
    Max_NDVI = max(NDVI_Smooth, na.rm = TRUE),
    Norm_NDVI = (NDVI_Smooth - Min_NDVI) / (Max_NDVI - Min_NDVI)
  ) %>%
  ungroup() %>%
  dplyr::select(-Min_NDVI, -Max_NDVI)

# Step 2: Identify Start of Season for each pixel
start_of_season <- data_normalized %>%
  group_by(x, y) %>%
  mutate(
    Above_Threshold = Norm_NDVI > 0.5,
    Lead1 = lead(Above_Threshold, 1, default = FALSE),
    Lead2 = lead(Above_Threshold, 2, default = FALSE),
    Start_Season_Flag = Above_Threshold & Lead1 & Lead2
  ) %>%
  filter(Start_Season_Flag) %>%
  summarise(Start_of_Season = min(Date)) %>%
  ungroup()
start_of_season <- start_of_season %>%
  mutate(Start_of_Season_Day = yday(Start_of_Season))
boxplot(start_of_season$Start_of_Season_Day)


