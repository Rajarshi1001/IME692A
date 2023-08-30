library("forecast")
library("dplyr")
library("ggplot2")
library("xts")
library("zoo")
library("rugarch")

q2data <- read.csv("q2_dataset.csv")
print(nrow(q2data))
# fetching data from 1901 to 1990
data <- q2data[1:90,]
# part 1
selected_months <- c("NOV", "DEC", "JAN", "FEB")
filtered_data1 <- data %>%
  select(all_of(selected_months))
# lines(filtered_data1$NOV)
filtered_data <- as.matrix(filtered_data1[,]) # converting to a numeric object

garch_forecasts <- list()
arima_forecasts <- list()
filtered_data <- ts(as.vector(t(filtered_data1)), frequency = 1)
# fitting the ARIMA model
arima_order <- c(1,0,1)
garch_order <- c(2,2)
for (month_idx in seq_along(selected_months)){
  arimamodel <- Arima(filtered_data1[,month_idx], order=arima_order)
  ugarch_spec <- ugarchspec(variance.model = list(model="sGARCH", garchOrder=garch_order))
  arima_residuals <- residuals(arimamodel)
  garch_model <- ugarchfit(data=arima_residuals, spec = ugarch_spec)
  arima_forecasts[[month_idx]] <- forecast(arimamodel, h=nrow(q2data) - nrow(data))
  garch_forecasts[[month_idx]] <- ugarchforecast(garch_model, n.ahead=nrow(q2data) - nrow(data))
}

# Combine ARIMA and GARCH forecasts
combined_forecasts <- lapply(seq_along(selected_months), function(month_idx) {
  arima_mean <- arima_forecasts[[month_idx]]$mean
  garch_volatility <- sigma(garch_forecasts[[month_idx]])
  combined_forecast <- arima_mean + rnorm(n=length(arima_mean), mean = 0, sd = 5.5)
  data.frame(Year = seq(1991, 2011), Forecast = combined_forecast)
})
combined_forecasts[[1]]
# Plot the combined forecasts
# ... (previous code)

# Plot the combined forecasts over original data
par(mfrow = c(2,2))
for (month_idx in seq_along(selected_months)) {
  original_data <- q2data[,selected_months[month_idx]]
  forecast_data <- combined_forecasts[[month_idx]]$Forecast
  
  # Plot the original data
  plot(q2data$YEAR, original_data, type = "l",
       main = paste("Original and Forecasted Rainfall for", selected_months[month_idx]),
       xlab = "Year", ylab = "Rainfall")
  
  # Add the forecast data to the plot
  lines(combined_forecasts[[month_idx]]$Year, forecast_data, col = "red")
}

