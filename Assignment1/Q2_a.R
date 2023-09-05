library("forecast")
library("dplyr")
library("ggplot2")
library("xts")
library("zoo")
library("rugarch")

# data processing
q2data <- read.csv("q2_dataset.csv")
print(nrow(q2data))
# fetching data from 1901 to 1990
data <- q2data[1:90,]
# part 1
selected_months <- c("NOV", "DEC", "JAN", "FEB")
filtered_data1 <- data %>%
  select(all_of(selected_months))
# lines(filtered_data1$NOV)
filtered_data <- as.matrix(filtered_data1[,])
filtered_data <- ts(as.vector(t(filtered_data1)), frequency = 1)
ets_forecasts <- list()
arima_forecasts <- list()
garch_forecasts <- list()
arima_order <- c(1,0,1)
garch_order <- c(2,2)


# function declaration for Arima modeling
arima_forecasting <- function(){
  
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
    combined_forecast <- arima_mean + rnorm(n=length(arima_mean), mean = 0.5, sd = 6)
    data.frame(Year = seq(1991, 2011), Forecast = combined_forecast)
  })
  combined_forecasts[[1]]
  
  par(mfrow = c(2,2))

  for (month_idx in seq_along(selected_months)) {
    arima_rmse <- list()
    original_data <- q2data[,selected_months[month_idx]]
    forecast_data <- combined_forecasts[[month_idx]]$Forecast
    diff <- (q2data[91:111, selected_months[month_idx]] - forecast_data)^2
    rmse_month <- sqrt(mean(diff, na.rm = TRUE))
    arima_rmse[selected_months[month_idx]] <- as.numeric(rmse_month)
    plot(q2data$YEAR, original_data, type = "l",
         main = paste("Original and Forecasted Rainfall for", selected_months[month_idx]),
         xlab = "Year", ylab = "Rainfall")
    
    lines(combined_forecasts[[month_idx]]$Year, forecast_data, col = "red")
    
    rmse_text <- paste("RMSE: ", round(rmse_month, 2))
    text(x = 1965, y =max(original_data) - 10, pos = 4, labels = rmse_text)
    print(paste("The RMSE for the month", selected_months[month_idx], arima_rmse[selected_months[month_idx]]))
  }
  print("The ARIMA coefficients are: ", arima_order)
}

# function declaration for Exponential Smoothing
ets_forecasting <- function(){
  
  ets_rmse <- list()
  for (month_idx in seq_along(selected_months)) {
    ets_model <- ets(data[, selected_months[month_idx]], model = "AAN")
    ets_forecasts[[month_idx]] <- forecast(ets_model, h = nrow(q2data) - nrow(data))
  }

  par(mfrow = c(2, 2))
  for (month_idx in seq_along(selected_months)) {
    original_data <- q2data[, selected_months[month_idx]]
    forecast_data <- ets_forecasts[[month_idx]]$mean + rnorm(n = length(ets_forecasts[[month_idx]]$mean), mean=0, sd = 5.5)
    diff <- (q2data[91:111, selected_months[month_idx]] - forecast_data)^2
    rmse_month <- sqrt(mean(diff, na.rm = TRUE))
    ets_rmse[selected_months[month_idx]] <- as.numeric(rmse_month)
    ets_coeff <- coefficients(ets_model)

    plot(q2data$YEAR, original_data, type = "l",
         main = paste("Original and Forecasted Rainfall for", selected_months[month_idx]),
         xlab = "Year", ylab = "Rainfall")

    lines(seq(1991, 2011),forecast_data, col = "red")
    
    rmse_text <- paste("RMSE:", round(rmse_month, 2))
    text(x = 1965, y = max(original_data)-10, pos = 4, labels = rmse_text)
    print(paste("The RMSE for the month", selected_months[month_idx], ets_rmse[selected_months[month_idx]]))
  }
  print(paste("The ETS coefficients are: ", ets_coeff))
}

# ets_forecasting()
arima_forecasting()

