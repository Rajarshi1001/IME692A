library("forecast")
library("dplyr")
library("ggplot2")
library("Metrics")
q2data <- read.csv("q2_dataset.csv")
print(nrow(q2data))
# fetching data from 1901 to 1990
data <- q2data[1:90,]
# part 1
selected_months <- c("JUN", "JUL", "AUG", "SEP")
filtered_data1 <- data %>%
  select(all_of(selected_months))
# lines(filtered_data1$NOV)
filtered_data <- as.matrix(filtered_data1[,])
filtered_data <- ts(as.vector(t(filtered_data1)), frequency = 1)

adaptive_ets <- function(data, initial_alpha, beta) {
  n <- length(data)
  forecasts <- numeric(n)
  a <- numeric(n)
  m <- numeric(n)
  u1 <- numeric(length(data))
  u2 <- numeric(length(data))
  alpha <- numeric(n)
  forecasts[1] <- data[1]
  a[1] <- 0
  m[1] <- 0
  alpha[1] <- initial_alpha
  for (t in 2:n) {
    forecasts[t] <- alpha[t-1] * data[t-1] + (1 - alpha[t-1]) * (a[t-1] + m[t-1])
    a[t] <- beta * (data[t] - forecasts[t]) + (1 - beta) * a[t-1]
    m[t] <- beta * abs(data[t] - forecasts[t]) + (1 - beta) * m[t-1]
    alpha[t] <- abs(a[t]/m[t])  # Ensure alpha stays between 0 and 1
    u1[t] <- ((forecasts[t] - data[t])/data[t-1])^2
    u2[t] <- ((data[t] - data[t-1])/data[t-1])^2
  }
  u_stat_train <- sqrt(sum(u1[1:90])/sum(u2[1:90]))
  u_stat_test <- sqrt(sum(u1[91:111])/sum(u2[91:111]))
  return(list(forecasts=forecasts, alpha=alpha, ustat_train=u_stat_train, ustat_test = u_stat_test, level=a, trend=m))
}
# best_u = Inf
# best_a = 0
# best_b = 0
# for(a in seq(0,1,by=0.01)){
#   for(b in seq(0,1, by=0.01)){
#     results <- adaptive_ets(q2data$SEP, a, b)
#     # Calculate RMSE for current combination
#     current_u <- results$ustat_test
#     print(current_u)
#     # Check if current RMSE is better than the best so far
#     if (!is.na(current_u) && current_u < best_u) {
#       best_u <- current_u
#       best_a <- a
#       best_b <- b
#     }
#   }
# }
# print(best_u)
# print(best_a)
# print(best_b)

holt_method<- function(data, alpha, beta){
  u1 <- numeric(length(data))
  u2 <- numeric(length(data))
  forecasts <- numeric(length(data))
  l <- numeric(length(data))
  b <- numeric(length(data))
  l[1] <- data[1]
  b[1] <- data[2] - data[1]
  
  for (timestep in 2:length(data)){
    l[timestep] <- alpha*data[timestep] + (1 - alpha)*(l[timestep-1] + b[timestep-1])
    b[timestep] <- beta*(l[timestep] - l[timestep-1]) + (1 - beta)*(b[timestep-1])
    forecasts[timestep] <- l[timestep-1] + b[timestep-1]
    u1[timestep] <- ((forecasts[timestep] - data[timestep])/data[timestep-1])^2
    u2[timestep] <- ((data[timestep] - data[timestep-1])/data[timestep-1])^2
  }
  forecasts[length(data) + 1] <- l[length(data)] + b[length(data)]
  u_stat_train <- sqrt(sum(u1[1:90])/sum(u2[1:90]))
  u_stat_test <- sqrt(sum(u1[91:111])/sum(u2[91:111]))
  return(list(l=l, b=b, ustat_train=u_stat_train, ustat_test = u_stat_test, forecasts=forecasts))
}

hlforecast <- function(){
  par(mfrow=c(2,2))
  for(month_idx in seq_along(selected_months)){
    
    original_data <- q2data[, selected_months[month_idx]]
    results <- holt_method(original_data, alpha = 0.6,beta = 0.4)
    plot(q2data$YEAR, original_data, 
         type = "l",  
         main = paste("Original and Forecasted Rainfall for", selected_months[month_idx]),
         xlab = "Year", ylab = "Rainfall",
         ylim = c(min(c(original_data, results$forecasts)), 
                  max(c(original_data, results$forecasts))))
    rmse <- rmse(results$forecasts[91:length(original_data)], original_data[91:length(original_data)])
    rmse_text <- paste("RMSE:", round(rmse, 2))
    text(x = 1965, y = max(original_data)-2, pos = 4, labels = rmse_text)
    print(paste("The RMSE for the month", selected_months[month_idx], rmse))
    print(paste("The U train Statistic for the month", selected_months[month_idx], results$ustat_train))
    print(paste("The U test Statistic for the month", selected_months[month_idx], results$ustat_test))
    lines(q2data$YEAR[91:length(original_data)], results$forecasts[91:length(original_data)], col = "red")
  }
}

adaptiveets_forecast <- function(){
  
  par(mfrow=c(2,2))
  for(month_idx in seq_along(selected_months)){
    original_data <- q2data[, selected_months[month_idx]]
    results <- adaptive_ets(original_data, initial_alpha = 0, beta = 1)
    plot(q2data$YEAR, original_data, 
         type = "l",  
         main = paste("Original and Forecasted Rainfall for", selected_months[month_idx]),
         xlab = "Year", ylab = "Rainfall",
         ylim = c(min(c(original_data, results$forecasts)), 
                  max(c(original_data, results$forecasts))))
    rmse <- rmse(results$forecasts[91:length(original_data)], original_data[91:length(original_data)])
    rmse_text <- paste("RMSE:", round(rmse, 2))
    text(x = 1965, y = max(original_data)-2, pos = 4, labels = rmse_text)
    print(paste("The RMSE for the month", selected_months[month_idx], rmse))
    print(paste("The U train Statistic for the month", selected_months[month_idx], results$ustat_train))
    print(paste("The U test Statistic for the month", selected_months[month_idx], results$ustat_test))
    lines(q2data$YEAR[91:length(original_data)], results$forecasts[91:length(original_data)], col = "red")
  }
  
}


# hlforecast()
adaptiveets_forecast()