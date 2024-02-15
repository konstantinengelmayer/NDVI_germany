source("~/edu/NDVI_germany/src/NDVI_germany_setup.R")
data <- readRDS(file.path(envrmt$path_data_level1, "long_ndvi_quality_bits_good_df.rds"))

data$Date <- as.Date(data$Date)
data$Month <- month(data$Date)

# Step 2: Group by Year, x, y, and Month to ensure unique entries for each month
# Then, count the number of records for each combination
df_monthly <- data %>%
  group_by(Year, x, y, Month) %>%
  summarise(Count = n(), .groups = 'drop') 
df_yearly_count <- df_monthly %>%
  group_by(Year, x, y) %>%
  summarise(Months_with_data = n_distinct(Month), .groups = 'drop')
pixels_with_complete_data <- df_yearly_count %>%
  filter(Months_with_data == 12)
df_filtered <- data %>%
  left_join(pixels_with_complete_data, by = c("Year", "x", "y")) %>%
  filter(Months_with_data == 12) %>%
  select(-Months_with_data)
# Count rows for every year in df_filtered
yearly_counts <- df_filtered %>%
  group_by(Year) %>%
  summarise(Count = n())
pixel_counts_per_year <- df_filtered %>%
  group_by(Year) %>%
  summarise(Pixel_Count = n_distinct(paste(x, y)))
saveRDS(df_filtered, file.path(envrmt$path_data_level1, "value_for_12month.rds"))
