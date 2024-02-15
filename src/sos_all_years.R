# Load NDVI Data
source("~/edu/NDVI_germany/src/NDVI_germany_setup.R")
data <- readRDS(file.path(envrmt$path_data_level1, "value_for_12month.rds"))

# Prepare Date Sequences for Each Year
# This generates a sequence of dates for every year present in the dataset
years <- unique(data$Year)
date_sequences <- lapply(years, function(year) {
  seq(ymd(paste(year, "01-01", sep = "-")), ymd(paste(year, "12-31", sep = "-")), by = "day")
})

# Expand Dataset to Include All Dates for Each Pixel Each Year
# Ensures that there's an entry for every day of the year for each pixel
data_expanded <- data %>%
  dplyr::select(Year, x, y) %>%
  distinct() %>%
  ungroup() %>%
  mutate(Date = map(Year, ~seq(ymd(paste0(.x, "-01-01")), ymd(paste0(.x, "-12-31")), by = "day"))) %>%
  unnest(Date) %>%
  left_join(data, by = c("Year", "x", "y", "Date"))

# Interpolate NDVI Values for Each Pixel
# Fills in missing NDVI values linearly for days without data
data_interpolated <- data_expanded %>%
  group_by(Year, x, y) %>%
  arrange(Date, .by_group = TRUE) %>%
  mutate(NDVI_Value = approx(Date, NDVI_Value, Date, method = "linear")$y) %>%
  ungroup()

# Apply Fourier Smoothing Across All Years
# Smoothes NDVI values using Fourier terms and linear regression
data_smoothed <- data_interpolated %>%
  group_by(Year, x, y) %>%
  do({
    ndvi_values <- .$NDVI_Value
    time_points <- seq_along(ndvi_values)
    K <- 3  # Number of harmonics for Fourier transform
    fourier_terms <- forecast::fourier(ts(ndvi_values, frequency = 365), K = K)
    fit <- lm(ndvi_values ~ fourier_terms)
    smoothed_values <- predict(fit, newdata = list(fourier_terms = fourier_terms))
    data.frame(Date = .$Date, NDVI_Smooth = as.numeric(smoothed_values), stringsAsFactors = FALSE)
  }) %>%
  ungroup()

# Normalize Smoothed NDVI Values Across All Years
data_normalized <- data_smoothed %>%
  group_by(Year, x, y) %>%
  mutate(
    Min_NDVI = min(NDVI_Smooth, na.rm = TRUE),
    Max_NDVI = max(NDVI_Smooth, na.rm = TRUE),
    Norm_NDVI = (NDVI_Smooth - Min_NDVI) / (Max_NDVI - Min_NDVI)
  ) %>%
  ungroup() %>%
  dplyr::select(-Min_NDVI, -Max_NDVI)

# Identify Start of Season for Each Pixel Across All Years
# The start of the season is defined as the first date where normalized NDVI exceeds 0.5 for 3 consecutive days
start_of_season <- data_normalized %>%
  group_by(Year, x, y) %>%
  mutate(
    Above_Threshold = Norm_NDVI > 0.5,
    Lead1 = lead(Above_Threshold, 1, default = FALSE),
    Lead2 = lead(Above_Threshold, 2, default = FALSE),
    Start_Season_Flag = Above_Threshold & Lead1 & Lead2
  ) %>%
  filter(Start_Season_Flag) %>%
  summarise(Start_of_Season = min(Date)) %>%
  ungroup()

# Convert Start of Season Date to Day of Year
start_of_season <- start_of_season %>%
  mutate(Day_of_Year = yday(Start_of_Season))

# Visualize Start of Season Variation Across Years
ggplot(start_of_season, aes(x = factor(Year), y = Day_of_Year)) +
  geom_boxplot() +
  labs(title = "Start of Season by Year", x = "Year", y = "Day of Year") +
  theme_minimal()
