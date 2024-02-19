# Load NDVI Data
source("~/edu/NDVI_germany/src/NDVI_germany_setup.R")
data <- readRDS("~/value_for_at_least_10_months.rds")

# Prepare Date Sequences for Each Year
# This generates a sequence of dates for every year present in the dataset
years <- unique(data$Year)
date_sequences <- lapply(years, function(year) {
  seq(ymd(paste(year, "01-01", sep = "-")), ymd(paste(year, "12-31", sep = "-")), by = "day")
})

# Adjusted pipeline with progress bar for the date expansion part
data_expanded <- data %>%
  dplyr::select(Year, x, y) %>%
  distinct() %>%
  ungroup() %>%
  mutate(Date = pblapply(Year, function(year) {
    seq(ymd(paste0(year, "-01-01")), ymd(paste0(year, "-12-31")), by = "day")
  })) %>%
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
season_dates <- data_normalized %>%
  group_by(Year, x, y) %>%
  mutate(
    # Start of Season (SoS) calculation
    Above_Threshold_SoS = Norm_NDVI > 0.5,
    Lead1_SoS = lead(Above_Threshold_SoS, 1, default = FALSE),
    Lead2_SoS = lead(Above_Threshold_SoS, 2, default = FALSE),
    Lead3_SoS = lead(Above_Threshold_SoS, 3, default = FALSE),
    Start_Season_Flag = Above_Threshold_SoS & Lead1_SoS & Lead2_SoS & Lead3_SoS,
    Start_of_Season = if_else(Start_Season_Flag, Date, as.Date(NA)),
    Start_of_Season = min(Start_of_Season, na.rm = TRUE),
    
    # Preparing for End of Season (EoS) calculation
    # Filter within mutate to avoid affecting Start_of_Season calculation
    Below_Threshold_EoS = if_else(month(Date) %in% c(8, 9, 10, 11) & !is.na(Start_of_Season), Norm_NDVI < 0.5, FALSE),
    Lead1_EoS = lead(Below_Threshold_EoS, 1, default = FALSE),
    Lead2_EoS = lead(Below_Threshold_EoS, 2, default = FALSE),
    Lead3_EoS = lead(Below_Threshold_EoS, 3, default = FALSE),
    End_Season_Flag = Below_Threshold_EoS & Lead1_EoS & Lead2_EoS & Lead3_EoS
  ) %>%
  group_by(Year, x, y) %>%
  summarise(
    Start_of_Season = first(Start_of_Season),
    End_of_Season = if_else(any(End_Season_Flag), min(Date[End_Season_Flag]), as.Date(NA))
  ) %>%
  ungroup()
season_dates <- season_dates[-which(is.na(season_dates$Start_of_Season)|is.na(season_dates$End_of_Season)),]
season_dates <- season_dates %>%
  mutate(
    SoS_Day_of_Year = yday(Start_of_Season),
    EoS_Day_of_Year = yday(End_of_Season)
  )
cor.test(season_dates$SoS_Day_of_Year,season_dates$EoS_Day_of_Year)

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

sos_trends <- season_dates %>%
  group_by(x, y) %>%
  do({
    tryCatch({
      mod_sos <- lm(SoS_Day_of_Year ~ Year, data = .)
      data.frame(Slope = coef(mod_sos)["Year"], PValue = summary(mod_sos)$coefficients["Year", "Pr(>|t|)"])
    }, error = function(e) {
      data.frame(Slope = NA, PValue = NA)  # Handle potential errors gracefully
    })
  }) %>%
  ungroup()

# Convert to sf object for spatial plotting
results_sf_sos <- st_as_sf(sos_trends, coords = c("x", "y"), crs = 4326)

# Filter for significant trends
results_significant_sos <- results_sf_sos[which(results_sf_sos$PValue <= 0.05),]

# Plotting significant trends for Start of Season
ggplot(results_significant_sos) +
  geom_sf(aes(color = Slope)) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
  labs(title = "Trend Start of Season over Years for each pixel (p <= 0.05)", color = "Trend Coefficient") +
  theme_minimal()

# For interactive map with mapview - assuming results_significant_sos is an sf object with a "Slope" column
mapview(results_significant_sos, 
        zcol = "Slope",  
        map.types = "OpenTopoMap")

eos_trends <- season_dates %>%
  group_by(x, y) %>%
  do({
    tryCatch({
      mod_eos <- lm(EoS_Day_of_Year ~ Year, data = .)
      data.frame(Slope = coef(mod_eos)["Year"], PValue = summary(mod_eos)$coefficients["Year", "Pr(>|t|)"])
    }, error = function(e) {
      data.frame(Slope = NA, PValue = NA)  # Handle potential errors gracefully
    })
  }) %>%
  ungroup()

# Convert to sf object for spatial plotting
results_sf_eos <- st_as_sf(eos_trends, coords = c("x", "y"), crs = 4326)

# Filter for significant trends
results_significant_eos <- results_sf_eos[which(results_sf_eos$PValue <= 0.05),]

# Plotting significant trends for End of Season
ggplot(results_significant_eos) +
  geom_sf(aes(color = Slope)) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
  labs(title = "Trend End of Season over Years for each pixel (p <= 0.05)", color = "Trend Coefficient") +
  theme_minimal()

# For interactive map with mapview - assuming results_significant_eos is an sf object with a "Slope" column
mapview(results_significant_eos, 
        zcol = "Slope",  
        map.types = "OpenTopoMap")


saveRDS(data_normalized, "~/normalized_NDVI_10_month.rds")
saveRDS(season_dates, "~/season_date.rds")
