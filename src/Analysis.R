source("~/edu/NDVI_germany/src/NDVI_germany_setup.R")
season_dates <- readRDS("D:/Recovery/edu/NDVI_germany/data/raster_data/data_level1/season_date_with_filter.rds")

season_dates <- season_dates %>%
  dplyr::filter(SoS_Month %in% 3:5, EoS_Month %in% 9:11)

yearly_data <- season_dates %>%
  group_by(Year) %>%
  summarise(
    MeanSoS = mean(SoS_Day_of_Year, na.rm = TRUE),
    MeanEoS = mean(EoS_Day_of_Year, na.rm = TRUE)
  )

model <- lm(data = yearly_data, formula = MeanSoS~Year)
summary(model)

plot(x = yearly_data$Year, y = yearly_data$MeanSoS)
abline(model)


library(lme4)

# Assuming 'data' is your dataframe, 'outcome_variable' is the dependent variable you're analyzing,
# 'Year' is your time variable, and each observation is identified by a unique combination of 'x' and 'y' coordinates
# for the pixel it comes from.

# Fit a mixed-effects model
model <- glmer(SoS_Day_of_Year ~ Year + (1 | x:y), data = season_dates)
summary(model)
