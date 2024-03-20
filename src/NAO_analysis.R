source("~/edu/NDVI_germany/src/NDVI_germany_setup.R")
NAO <- read.table(file.path(envrmt$path_data, "NAO/norm.nao.monthly.b5001.current.ascii"), col.names = c("Year", "Month", "NAO"))
season_dates <- readRDS(paste0(envrmt$path_data_level1, "/season_date_with_filter.rds"))
head(season_dates)
season_dates <- season_dates %>%
  group_by(x, y) %>%
  dplyr::filter(n_distinct(Year) >= 20) %>%
  ungroup()
cor.test(filtered_season_dates$SoS_Day_of_Year, filtered_season_dates$EoS_Day_of_Year)

NAO <- NAO %>%
  mutate(FullDate = as.Date(paste(Year, Month, "01", sep = "-")))
plot(y= NAO$NAO, x= NAO$FullDate)


sos_trends <- season_dates %>%
  group_by(x, y) %>%
  # Filter groups with at least 10 distinct years
  dplyr::filter(n_distinct(Year) >= 20) %>%
  do({
    tryCatch({
      mod_sos <- lm(SoS_Day_of_Year ~ Year, data = .)
      data.frame(Slope = coef(mod_sos)["Year"], PValue = summary(mod_sos)$coefficients["Year", "Pr(>|t|)"])
    }, error = function(e) {
      data.frame(Slope = NA, PValue = NA)  # Handle potential errors gracefully
    })
  }) %>%
  ungroup()
results_sf_sos <- st_as_sf(sos_trends, coords = c("x", "y"), crs = 4326)
results_significant_sos <- results_sf_sos[which(results_sf_sos$PValue <= 0.05),]
ggplot(results_significant_sos) +
  geom_sf(aes(color = Slope)) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
  labs(title = "Trend End of Season over Years for each pixel (p <= 0.05)", color = "Trend Coefficient") +
  theme_minimal()

results_significant_sos <- results_significant_sos %>%
  mutate(x = st_coordinates(geometry)[, 1], # Extracting longitude
         y = st_coordinates(geometry)[, 2]) # Extracting latitude
selected_rows <- season_dates %>%
  inner_join(results_significant_sos, by = c("x", "y"))

#selected_rows <- season_dates %>%
#  rowwise() %>%
#  mutate(
#    Avg_NAO_Before_Season = {
#      # Calculate the start month and year for NAO analysis (3 months before the season start)
#      analysis_start_date <- Start_of_Season %m-% months(3)
#      analysis_start_month <- month(analysis_start_date)
#      analysis_start_year <- year(analysis_start_date)
#      
#      # Adjust for when the analysis spans the new year
#      analysis_end_month <- month(Start_of_Season) - 1
#      analysis_end_year <- if_else(analysis_end_month == 0, Year - 1, Year)
#      analysis_end_month <- if_else(analysis_end_month == 0, 12, analysis_end_month)
#      
#      # Handling case where analysis period spans over two years
#      mean_NAO <- if (analysis_start_year < analysis_end_year) {
#        rbind(
#          NAO %>% dplyr::filter(Year == analysis_start_year & Month >= analysis_start_month),
#          NAO %>% dplyr::filter(Year == analysis_end_year & Month <= analysis_end_month)
#        ) %>% summarise(Avg_NAO = mean(NAO, na.rm = TRUE)) %>% pull(Avg_NAO)
#      } else {
#        NAO %>% 
#          dplyr::filter(Year == analysis_start_year & Month >= analysis_start_month & Month <= analysis_end_month) %>%
#          summarise(Avg_NAO = mean(NAO, na.rm = TRUE)) %>%
#          pull(Avg_NAO)
#      }
#      
#      mean_NAO
#    }
#  ) %>%
#  ungroup() # Ensure further operations are not in rowwise mode

#cor.test(selected_rows$SoS_Day_of_Year, selected_rows$Avg_NAO_Before_Season)

#ggplot(data = selected_rows)+
#  geom_point(aes(x = selected_rows$Avg_NAO_Before_Season, y = SoS_Day_of_Year))


# Reading the data from the file, assuming your data is saved in 'nao_data.txt'
nao_data <- read.table(file.path(envrmt$path_data, "NAO/NAO_daily.txt"), sep = ",", col.names = c("NAO", "Year", "Month", "Day"))

# Adding a full date column
nao_data$FullDate <- as.Date(with(nao_data, paste(Year, Month, Day, sep = "-")), "%Y-%m-%d")

# Display the first few rows to verify
head(nao_data)
selected_rows$Start_of_Season <- as.Date(selected_rows$Start_of_Season)

# Function to calculate mean NAO for 3 months before SoS
calculate_mean_NAO <- function(sos_date) {
  # Calculate the start date 3 months before SoS
  start_date <- sos_date - months(3)
  
  # Filter NAO records between start_date and sos_date
  relevant_nao <- nao_data[nao_data$FullDate >= start_date & nao_data$FullDate < sos_date, ]
  
  # Calculate mean NAO
  mean_nao <- mean(relevant_nao$NAO, na.rm = TRUE)
  
  return(mean_nao)
}

# Apply the function to each Start_of_Season date in selected_rows
season_dates$Mean_NAO_Before_SoS <- sapply(season_dates$Start_of_Season, calculate_mean_NAO)

# Display the updated selected_rows
head(selected_rows)
cor.test(season_dates$SoS_Day_of_Year, season_dates$Mean_NAO_Before_SoS)
ggplot(data = season_dates) +
  geom_point(aes(x = Mean_NAO_Before_SoS, y = SoS_Day_of_Year, col = Year)) +
  scale_color_viridis_c(option = "D") + # option can be "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"
  theme_minimal() +
  labs(colour = "Year")


eos_trends <- season_dates %>%
  group_by(x, y) %>%
  dplyr::filter(n_distinct(Year) >= 20) %>%
  do({
    tryCatch({
      mod_eos <- lm(EoS_Day_of_Year ~ Year, data = .)
      data.frame(Slope = coef(mod_eos)["Year"], PValue = summary(mod_eos)$coefficients["Year", "Pr(>|t|)"])
    }, error = function(e) {
      data.frame(Slope = NA, PValue = NA)  # Handle potential errors gracefully
    })
  }) %>%
  ungroup()
results_sf_eos <- st_as_sf(eos_trends, coords = c("x", "y"), crs = 4326)
results_significant_eos <- results_sf_eos[which(results_sf_eos$PValue <= 0.05),]

results_significant_sos <- results_significant_sos %>%
  mutate(x = st_coordinates(geometry)[, 1], # Extracting longitude
         y = st_coordinates(geometry)[, 2]) # Extracting latitude
selected_rows <- season_dates %>%
  inner_join(results_significant_sos, by = c("x", "y"))




selected_rows$End_of_Season <- as.Date(selected_rows$End_of_Season)

# Function to calculate mean NAO for 3 months before SoS
calculate_mean_NAO <- function(eos_date) {
  # Calculate the start date 3 months before SoS
  start_date <- eos_date - months(3)
  
  # Filter NAO records between start_date and sos_date
  relevant_nao <- nao_data[nao_data$FullDate >= start_date & nao_data$FullDate < eos_date, ]
  
  # Calculate mean NAO
  mean_nao <- mean(relevant_nao$NAO, na.rm = TRUE)
  
  return(mean_nao)
}

# Apply the function to each Start_of_Season date in selected_rows
season_dates$Mean_NAO_Before_EoS <- sapply(season_dates$End_of_Season, calculate_mean_NAO)
cor.test(season_dates$EoS_Day_of_Year, season_dates$Mean_NAO_Before_EoS)
ggplot(data = selected_rows) +
  geom_point(aes(x = Mean_NAO_Before_EoS, y = EoS_Day_of_Year, col = Year)) +
  scale_color_viridis_c(option = "D") + # option can be "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"
  theme_minimal() +
  labs(colour = "Year")

season_dates <- season_dates %>%
  mutate(Vegetation_Period_Length = as.integer(End_of_Season - Start_of_Season))
saveRDS(season_dates, paste0(envrmt$path_data_level1, "/saveSeasonDatesAnalysis.RDS"))
