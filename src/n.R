results_df <- data.frame(
  Min_Years = integer(),
  SOS_Slope = numeric(),
  PGT_Slope = numeric(),
  EOS_Slope = numeric()
)

# Loop over the range of minimum years
for (min_years in 10:23) {
  # Filter out groups with less than 'min_years' of data
  filtered <- filtered_sos %>%
    group_by(x, y) %>%
    filter(sum(!is.na(SOS)) >= min_years & sum(!is.na(EOS)) >= min_years) %>%
    ungroup() %>%
    mutate(
      Year = as.numeric(Year),  # Ensure Year is numeric
      SOS_julian = yday(SOS),   # Convert SOS to the Julian day
      PGT_julian = yday(PGT),   # Convert PGT to the Julian day
      EOS_julian = yday(EOS)    # Convert EOS to the Julian day
    )
  
  # Only proceed if the filtered dataset is not empty
  if (nrow(filtered) > 0) {
    # Perform linear regression for SOS, PGT, and EOS
    lm_sos <- lm(SOS_julian ~ Year, data = filtered)
    lm_pgt <- lm(PGT_julian ~ Year, data = filtered)
    lm_eos <- lm(EOS_julian ~ Year, data = filtered)
    
    # Extract the slope (coefficient of Year) for each regression
    sos_slope <- coef(lm_sos)["Year"]
    pgt_slope <- coef(lm_pgt)["Year"]
    eos_slope <- coef(lm_eos)["Year"]
    
    # Append results
    results_df <- rbind(results_df, data.frame(Min_Years = min_years, SOS_Slope = sos_slope, PGT_Slope = pgt_slope, EOS_Slope = eos_slope))
  } else {
    # Append NA if no data available for the given minimum number of years
    results_df <- rbind(results_df, data.frame(Min_Years = min_years, SOS_Slope = NA, PGT_Slope = NA, EOS_Slope = NA))
  }
}

# Print the results
print(results_df)



################################################################################
# Function to perform linear regression and return p-value and slope
perform_regression <- function(data) {
  lm_sos <- lm(SOS_julian ~ Year, data = data)
  lm_eos <- lm(EOS_julian ~ Year, data = data)
  
  p_sos <- summary(lm_sos)$coefficients[2,4]  # P-value for SOS
  p_eos <- summary(lm_eos)$coefficients[2,4]  # P-value for EOS
  slope_sos <- coef(lm_sos)["Year"]          # Slope for SOS
  slope_eos <- coef(lm_eos)["Year"]          # Slope for EOS
  
  return(data.frame(p_sos = p_sos, p_eos = p_eos, slope_sos = slope_sos, slope_eos = slope_eos))
}

# Apply the function to each pixel group and add results as new coloumns
reg_results <- filtered_excl %>%
  group_by(x, y) %>%
  do(perform_regression(.)) %>%
  ungroup()

# Filter pixel groups with significant p-values
sig_pix <- reg_results %>%
  filter(p_sos <= 0.05 | p_eos <= 0.05)

# Join with original data to get significant pixel groups only
sig_filtered_excl <- filtered_excl %>%
  inner_join(sig_pix, by = c("x", "y"))

print(paste("Percentage of significant Observations:", nrow(sig_filtered_excl)/nrow(filtered_excl)))


str(sig_pix)
str(sig_filtered_excl)


model <- lm(SOS_julian ~ Year, data = sig_filtered_excl)
summary(model)
plot(sig_filtered_excl$Year, sig_filtered_excl$SOS_julian)
abline(model)
qqnorm(residuals(model))
qqline(residuals(model))
shapiro.test(residuals(model))




glmer