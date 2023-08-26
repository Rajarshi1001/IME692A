library("forecast")
library("ggplot2")

starting_year <- 1956
end_year <- 1956 + 50
data <- read.csv("Data1.csv")

# Create time series object
tsdata <- ts(data$Expendture, start = starting_year, end=end_year-1)
plot(tsdata)
arimamodel <- arima(tsdata, order=c(0,1,0), seasonal=c(0,1,0))
# Forecast for the next 15 years
forecastdata <- forecast(arimamodel, h = 15)

forecast_df <- data.frame(
  Year = seq(max(time(tsdata)) + 1, length.out = 15),
  Forecast = forecastdata$mean
)
# Convert original data to a data frame
data_df <- data.frame(
  Year = time(tsdata),
  Expenditure = tsdata
)

# Combine original and forecast data frames
combined_df <- rbind(data_df, forecast_df)

# plotting the data
ggplot() +
  geom_line(data = data, aes(x = Year, y = Expendture, color="Actual"), lwd = 2) +
  geom_line(data = forecast_df, aes(x = Year, y = Forecast, color = "Forecast"), lwd = 3) +
  labs(title = "Expenditure Over Years",
       x = "Year",
       y = "Expenditure") +
  theme_bw() +
  scale_color_manual(values = c("Actual" = "black", "Forecast" = "red")) +
  theme(plot.title = element_text(hjust = 0.5)) 


