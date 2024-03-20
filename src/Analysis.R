source("~/NDVI_germany/src/NDVI_germany_setup.R")
season_dates <- readRDS(paste0(envrmt$path_data_level1, "/season_date_all_interpol_golay_3_11.rds"))
season_dates1 <- readRDS(paste0(envrmt$path_data_level1, "/season_date_with_filter.rds"))
normalized_data <- readRDS(paste0(envrmt$path_data_level1, "/normalized_with_filter.rds"))
cor.test(season_dates$SoS_Day_of_Year, season_dates$EoS_Day_of_Year)
colnames(season_dates)[4:5] <- c("Start_of_Season", "End_of_Season")
colnames(season_dates)[13:14] <- c("SoS_Day_of_Year", "EoS_Day_of_Year")
season_dates$SoS_Day_of_Year <- yday(season_dates$Start_of_Season)
season_dates$EoS_Day_of_Year <- yday(season_dates$End_of_Season)
season_dates$Year <- as.numeric(season_dates$Year)


season_dates <- season_dates %>%
  mutate(
    SoS_Month = month(Start_of_Season),
    EoS_Month = month(End_of_Season)
  )

season_dates <- season_dates %>%
  dplyr::filter(!is.na(SoS_Day_of_Year) & !is.na(EoS_Day_of_Year))

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
  slice(c(1:10, 12:n(), 11))

season_dates <- season_dates %>%
  dplyr::filter(SoS_Month %in% 3:5, EoS_Month %in% 9:11)

season_dates <- season_dates %>%
  group_by(x, y) %>%
  dplyr::filter(n_distinct(Year) >= 20) %>%
  ungroup()

model <- glmer(SoS_Day_of_Year ~ Year + (1 | x:y), data = season_dates)
summary(model)

sos_trends <- season_dates %>%
  group_by(x, y) %>%
  do({
    tryCatch({
      mod_sos <- rlm(SoS_Day_of_Year ~ Year, data = .)
      coef_est <- coef(mod_sos)["Year"]
      std_error <- summary(mod_sos)$coefficients["Year", "Std. Error"]
      # Calculate t-value
      t_value <- coef_est / std_error
      # Calculate p-value based on t-distribution
      df <- nrow(.) - 2  # degrees of freedom = number of observations - number of coefficients
      p_value <- 2 * pt(-abs(t_value), df)
      data.frame(Slope = coef_est, PValue = p_value)
    }, error = function(e) {
      data.frame(Slope = NA, PValue = NA)  # Handle potential errors gracefully
    })
  }) %>%
  ungroup()

# Convert to sf object
results_sf_sos <- st_as_sf(sos_trends, coords = c("x", "y"), crs = 4326)
results_significant_sos <- results_sf_sos[results_sf_sos$PValue <= 0.05,]

eos_trends <- season_dates %>%
  group_by(x, y) %>%
  do({
    tryCatch({
      mod_eos <- rlm(EoS_Day_of_Year ~ Year, data = .)
      coef_est <- coef(mod_eos)["Year"]
      std_error <- summary(mod_eos)$coefficients["Year", "Std. Error"]
      # Calculate t-value
      t_value <- coef_est / std_error
      # Calculate p-value based on t-distribution
      df <- nrow(.) - 2  # degrees of freedom = number of observations - number of coefficients
      p_value <- 2 * pt(-abs(t_value), df)
      data.frame(Slope = coef_est, PValue = p_value)
    }, error = function(e) {
      data.frame(Slope = NA, PValue = NA)  # Handle potential errors gracefully
    })
  }) %>%
  ungroup()
results_sf_eos <- st_as_sf(eos_trends, coords = c("x", "y"), crs = 4326)
results_significant_eos <- results_sf_eos[which(results_sf_eos$PValue <= 0.05),]

min_slope <- min(c(min(results_significant_sos$Slope, na.rm = TRUE), min(results_significant_eos$Slope, na.rm = TRUE)), na.rm = TRUE)
max_slope <- max(c(max(results_significant_sos$Slope, na.rm = TRUE), max(results_significant_eos$Slope, na.rm = TRUE)), na.rm = TRUE)

all_data <- rbind(st_coordinates(results_significant_sos), st_coordinates(results_significant_eos))
xlims <- range(all_data[,1], na.rm = TRUE) # Longitude limits
ylims <- range(all_data[,2], na.rm = TRUE) # Latitude limits


g1 <- ggplot(results_significant_sos) +
  geom_sf(aes(color = Slope)) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = c(min_slope, max_slope)) +
  labs(title = "Start of Season (n = 1498)") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, vjust = -2, size = 11)) + # Center title, adjust position, and set size
  coord_sf(xlim = xlims, ylim = ylims)

# Update g2 with centered title and smaller font size
g2 <- ggplot(results_significant_eos) +
  geom_sf(aes(color = Slope)) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = c(min_slope, max_slope)) +
  labs(title = "End of Season (n = 544)", color = "Trend Coefficient") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2, size = 11)) + # Center title, adjust position, and set size
  coord_sf(xlim = xlims, ylim = ylims)

margin_adjustment <- margin(t = 0, r = 5, b = 0, l = 35, unit = "pt")

g1 <- g1 + theme(plot.margin = margin_adjustment)
g2 <- g2 + theme(plot.margin = margin_adjustment)

# Now combine them with a shared legend and a common title
g_combined <- g1 + g2 +
  plot_layout(guides = "collect") 

# Display the adjusted combined plot
print(g_combined)

ggsave(paste0(envrmt$path_docs,"/eos_sos_trends.png"), g_combined, width = 12, height = 6, dpi = 300)





results_significant_sos <- results_significant_sos %>%
  mutate(x = st_coordinates(geometry)[, 1], # Extracting longitude
         y = st_coordinates(geometry)[, 2]) # Extracting latitude
selected_rows <- season_dates %>%
  inner_join(results_significant_sos, by = c("x", "y"))
average_trend1 <- selected_rows %>%
  group_by(Year) %>%
  summarise(Avg_SoS = mean(SoS_Day_of_Year, na.rm = TRUE))

results_significant_eos <- results_significant_eos %>%
  mutate(x = st_coordinates(geometry)[, 1], # Extracting longitude
         y = st_coordinates(geometry)[, 2]) # Extracting latitude
selected_rows <- season_dates %>%
  inner_join(results_significant_eos, by = c("x", "y"))
average_trend2 <- selected_rows %>%
  group_by(Year) %>%
  summarise(Avg_EoS = mean(EoS_Day_of_Year, na.rm = TRUE))

# Create the SoS plot with subtitle
plot_sos <- ggplot() +
  geom_line(data = selected_rows, aes(x = Year, y = SoS_Day_of_Year, group = interaction(x, y)), alpha = 0.2) +
  geom_line(data = average_trend1, aes(x = Year, y = Avg_SoS), color = "red", size = 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2, size = 9)) + # Center title, adjust position, and set size
  labs(y = "SoS (Annual Day Count)", subtitle = "Start of Season (n = 1498)")

# Create the EoS plot with subtitle
plot_eos <- ggplot() +
  geom_line(data = selected_rows, aes(x = Year, y = EoS_Day_of_Year, group = interaction(x, y)), alpha = 0.2) +
  geom_line(data = average_trend2, aes(x = Year, y = Avg_EoS), color = "red", size = 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2, size = 9)) + # Center title, adjust position, and set size
  labs(y = "EoS (Annual Day Count)", subtitle = "End of Season (n = 544)")

# Combine the plots into one grid with a shared title
combined_plot <- plot_sos / plot_eos +
  plot_layout(guides = "collect")

# Display the combined plot
print(combined_plot)
ggsave(paste0(envrmt$path_docs,"/Anual_eos_sos_trends.png"), combined_plot, width = 12, height = 11, dpi = 300)





combined_grid <- (g_combined / combined_plot)

# Then, apply the layout. This step will not create the exact layout you described because
# of the limitations in dynamically configuring rows and columns directly in patchwork without a layout matrix
combined_grid <- combined_grid + 
  plot_layout(heights = c(1,1.5)) +
  plot_annotation(title = "Comprehensive Analysis of SoS and EoS Trends (p <= 0.05)")+
  theme(plot.title = element_text(face = "bold", size = 20)) # Adjust size as needed


print(combined_grid)
ggsave(paste0(envrmt$path_docs,"/full_eos_sos_trends1.png"), combined_grid, width = 12, height = 16, dpi = 300)



cor.test(selected_rows$Year, selected_rows$EoS_Day_of_Year)

results_significant_sos_df <- as.data.frame(results_significant_sos %>% dplyr::select(x, y))
results_significant_eos_df <- as.data.frame(results_significant_eos %>% dplyr::select(x, y))

# Perform the inner join on these data frames
common_pixels <- inner_join(results_significant_sos_df, 
                            results_significant_eos_df, 
                            by = c("x", "y"))

# Count the number of common pixels
nrow(common_pixels)

cor.test(season_dates$SoS_Day_of_Year, season_dates$EoS_Day_of_Year)

replacement_vector <- setNames(c("January", "February", "March", "April", "May", "June", 
                                 "July", "August", "September", "October", "November", "December"), 
                               c("Januar", "Februar", "März", "April", "Mai", "Juni", 
                                 "Juli", "August", "September", "Oktober", "November", "Dezember"))

month_summary$Month <- replacement_vector[month_summary$Month]
colnames(month_summary)[2:3] <- c("SoS", "EoS")
write.csv(month_summary,paste0(envrmt$path_docs, "/month_summery.csv"))


results_significant_eos$geometry <- st_as_text(results_significant_eos$geometry)
results_significant_sos$geometry <- st_as_text(results_significant_sos$geometry)

# Convert sf objects to regular data frames
results_significant_df <- as.data.frame(results_significant_eos)
results_sf_df <- as.data.frame(results_significant_sos)

# Merge based on the geometry column now represented as WKT text
merged_df <- merge(results_significant_df, results_sf_df, by = "geometry")

# Printing the result to check
print(merged_df)



sos_trends <- season_dates %>%
  group_by(x, y) %>%
  # Filter groups with at least 10 distinct years
  do({
    tryCatch({
      mod_sos <- lm(Vegetation_Period_Length ~ Year, data = .)
      data.frame(Slope = coef(mod_sos)["Year"], PValue = summary(mod_sos)$coefficients["Year", "Pr(>|t|)"])
    }, error = function(e) {
      data.frame(Slope = NA, PValue = NA)  # Handle potential errors gracefully
    })
  }) %>%
  ungroup()
results_sf_sos <- st_as_sf(sos_trends, coords = c("x", "y"), crs = 4326)
results_significant_sos <- results_sf_sos[which(results_sf_sos$PValue <= 0.05),]


g1 <- ggplot(results_significant_sos) +
  geom_sf(aes(color = Slope)) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0,) +
  labs(title = "Start of Season (n = 1272)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2, size = 11)) 

results_significant_sos <- results_significant_sos %>%
  mutate(x = st_coordinates(geometry)[, 1], # Extracting longitude
         y = st_coordinates(geometry)[, 2]) # Extracting latitude
selected_rows <- season_dates %>%
  inner_join(results_significant_sos, by = c("x", "y"))
average_trend1 <- selected_rows %>%
  group_by(Year) %>%
  summarise(Avg_SoS = mean(Vegetation_Period_Length, na.rm = TRUE))
plot_sos <- ggplot() +
  geom_line(data = selected_rows, aes(x = Year, y = Vegetation_Period_Length, group = interaction(x, y)), alpha = 0.2) +
  geom_line(data = average_trend1, aes(x = Year, y = Avg_SoS), color = "red", size = 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2, size = 9)) + # Center title, adjust position, and set size
  labs(y = "SoS (Annual Day Count)", subtitle = "Start of Season (n = 1272)")
print(plot_sos)
