library(tidyverse)
library(dplyr)
library(ggplot2)
library(car)

data <- read.csv("cadata.csv")
init_cols <- colnames(data)

data <- data %>%
  mutate(
    y = log(Median.House.Value, base = exp(1)),
    f1 = Median.Income,
    f2 = Median.Income^2,
    f3 = Median.Income^3,
    f4 = log(Housing.Median.Age, base = exp(1)),
    f5 = log(Total.Rooms/Population, base = exp(1)),
    f6 = log(Total.Bedrooms/Population, base = exp(1)),
    f7 = log(Population/Households, base = exp(1)),
    f8 = log(Households,base = exp(1))
  ) %>%
  select(y, f1:f8)

beta <- as.matrix(c(11.4939, 0.4790, -0.0166, -0.0002, 0.1570, -0.8582, 0.8043, -0.4077, 0.0477))
X <- as.matrix(cbind(1, data[,2:9]))
y <- as.matrix(data[,"y"])
y_hat <-  X %*% beta
errors <- y - y_hat
errors_square <- errors*2
var <- mean(errors_square)
# Plotting the distribution of errors
errors_data <- data.frame(errors = errors)
ggplot(errors_data, aes(x = errors)) +
  geom_histogram(aes(y=..density..), binwidth=0.1, fill="blue", alpha=0.7) + 
  geom_density(color="red") + 
  labs(title = "Distribution of the Errors", y = "Density", x = "Errors") + 
  theme_minimal()

# Checking dependence between the errors
res_data <- data.frame(residuals = errors, fitted = y_hat)
ggplot(res_data, aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  labs(title="Residuals vs Fitted", x="Fitted values", y="Residuals") +
  theme_minimal()

# Checking for multicollinearity
vif_res <- vif(lm(y ~ f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8, data = data))
print(vif_res)


beta_estimator <- solve(t(X) %*% X) %*% t(X) %*% y
estimator_var <- var*(solve(t(X) %*% X))

# part a
sample_size <- 100
num_samples <- 200
predictors <- 8

beta <- as.matrix(c(11.4939, 0.4790, -0.0166, -0.0002, 0.1570, -0.8582, 0.8043, -0.4077, 0.0477))
beta_hats <- matrix(0, nrow = num_samples, ncol = predictors + 1)
beta_probs <- matrix(0, nrow = num_samples, ncol = predictors + 1)
conf_data_1 <- data.frame()
conf_data_2 <- data.frame()

for (i in 1:num_samples) {
  t <- runif(sample_size * predictors)
  X <- cbind(1, matrix(t, sample_size, predictors))
  y <- X %*% beta + rnorm(sample_size)
  error <- y - X %*% beta
  
  var_cov_matrix <- (mean(error^2)) * solve(t(X) %*% X)
  se_beta <- sqrt(diag(var_cov_matrix))
  beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y
  
  
  conf_levels <- c(0.90, 0.95, 0.99)
  z_values <- qnorm(1 - (1 - conf_levels) / 2)
  t_values <- qt(1 - (1 - conf_levels) / 2, df = sample_size-1)
  
  for (k in 1:length(beta)) {
    for (j in 1:length(conf_levels)) {
      lower_bound <- beta_hat[k] - t_values[j] * se_beta[k]
      upper_bound <- beta_hat[k] + t_values[j] * se_beta[k]
      conf_data_2 <- rbind(conf_data_2, data.frame(Beta = k,
                                                   Estimate = beta_hat[k],
                                                   ConfidenceLevel = conf_levels[j] * 100,
                                                   Lower = lower_bound,
                                                   Upper = upper_bound))
    }
  }
  
  # Compute and print confidence intervals
  for (k in 1:length(beta)) {
    for (j in 1:length(conf_levels)) {
      lower_bound <- beta[k] - z_values[j] * se_beta[k]
      upper_bound <- beta[k] + z_values[j] * se_beta[k]
      conf_data_1 <- rbind(conf_data_1, data.frame(beta_index = k,
                                               estimate = beta[k],
                                               confidence_level = conf_levels[j],
                                               lower_int = lower_bound,
                                               upper_int = upper_bound))
    }
  }
  
  beta_hats[i, ] <- beta_hat
  beta_probs[i, ] <- sapply(1:(predictors + 1), function(j) {
    dnorm(x = beta_hat[j], mean = beta[j], sd = sqrt(var_cov_matrix[j, j]))
  })
}

expected_beta_hats <- colMeans(beta_hats)
comp_data <- data.frame(actual = beta, expected_beta_hat = expected_beta_hats, beta_probs = colMeans(beta_probs))

print(comp_data)
# print(head(conf_data))

ggplot(conf_data_1, aes(y = factor(beta_index, levels = predictors:1), x = estimate)) +
  geom_point(aes(color = as.factor(confidence_level)), size = 3) +
  geom_errorbarh(aes(xmin = lower_int, xmax = upper_int, height = 0.2, color = as.factor(confidence_level))) +
  labs(title = "Confidence Intervals for Beta Estimates", x = "Value", y = "Beta Index", color = "Confidence Level") +
  theme_minimal()

# Plotting the confidence intervals
ggplot(conf_data_2, aes(y = factor(Beta, levels = predictors:1), x = Estimate)) +
  geom_point(aes(color = as.factor(ConfidenceLevel)), size = 3) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper, height = 0.2, color = as.factor(ConfidenceLevel))) +
  labs(title = "Confidence Intervals for Beta Estimates", x = "Value", y = "Beta Index", color = "Confidence Level") +
  theme_minimal()

# Return the bar plot
ggplot(comp_data, aes(x = as.factor(round(actual, 4)), y = beta_probs)) + 
  geom_bar(stat = "identity", fill = "black", alpha = 0.5, width = 0.7) +
  geom_text(aes(label=sprintf("%.3f", beta_probs)), vjust=-0.5, size=3.0) + 
  labs(title = "Probability Densities of MLR coefficients", x = "beta's", y = "probability density") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# part b 
a = 0.01
n <- 20640
p <- 9
l <- matrix(c(1,1,1,1,1,1,1,1,1))
# model 1
beta_1 <- solve(t(X) %*% X) %*% t(X) %*% y
b <- t(l) %*% solve(t(X) %*% X) %*% l
ols <- t(X %*% beta - X %*% beta_1) %*% (X %*% beta - X %*% beta_1)
theta <- t(l) %*% beta
theta_hat <- t(l) %*% beta_1
linex <- exp(a*(theta - theta_hat)) - a*(theta - theta_hat) - 1
lambda <- c(0.2, 0.4, 0.6, 0.8)
for(i in 1:4){
  risk_1 <- lambda[i]*(n - p)*var + (1 - lambda[i])*((exp(a^2*var)/2)*(b)-1)
  loss_1 <-lambda[i]*(ols) + (1 - lambda[i])*(linex) 
  print(paste("Risk for Model 1:",risk_1))
  print(paste("Loss for Model 1:",loss_1))
}

# model 2
b <- t(l) %*% solve(t(X) %*% X) %*% l
beta_2 <- solve(t(X) %*% X) %*% t(X) %*% y
theta <- t(l) %*% beta
theta_hat <- t(l) %*% beta_2
ols <- t(beta - beta_2) %*% (beta - beta_2)
linex <- exp(a*(theta - theta_hat)) - a*(theta - theta_hat) - 1
lambda <- c(0.2, 0.4, 0.6, 0.8)
for(i in 1:4){
  risk_1 <- lambda[i]*var + (1 - lambda[i])*((exp(a^2*var)/2)*(b)-1)
  loss_1 <-lambda[i]*(ols) + (1 - lambda[i])*(linex) 
  print(paste("Risk for Model 1:",risk_1))
  print(paste("Loss for Model 1:",loss_1))
}

# model 3
b <- t(l) %*% solve(t(X) %*% X) %*% l
beta_3 <- solve(t(X) %*% X) %*% t(X) %*% y
theta <- t(l) %*% beta
theta_hat <- t(l) %*% beta_3 - (a*var/2)*(t(l) %*% solve(t(X) %*% X) %*% l)
ols <- t(X %*% beta - X %*% beta_3) %*% (X %*% beta - X %*% beta_3)
linex <- exp(a*(theta - theta_hat)) - a*(theta - theta_hat) - 1
lambda <- c(0.2, 0.4, 0.6, 0.8)
for(i in 1:4){
  risk_1 <- lambda[i]*((n-p)*var + a*a*(((n-p)*var)^2)/4) + (1 - lambda[i])*(((a^2*var)/2)*(b)-1)
  loss_1 <-lambda[i]*(ols) + (1 - lambda[i])*(linex) 
  print(paste("Risk for Model 1:",risk_1))
  print(paste("Loss for Model 1:",loss_1))
}

# model 4
b <- t(l) %*% solve(t(X) %*% X) %*% l
beta_4 <- solve(t(X) %*% X) %*% t(X) %*% y
theta <- t(l) %*% beta
theta_hat <- t(l) %*% beta_4 - (a*var/2)*(t(l) %*% solve(t(X) %*% X) %*% l)
ols <- t(beta - beta_4) %*% (beta - beta_4)
linex <- exp(a*(theta - theta_hat)) - a*(theta - theta_hat) - 1
lambda <- c(0.2, 0.4, 0.6, 0.8)
for(i in 1:4){
  risk_1 <- lambda[i]*(var*b + a^2*(var*b)^2/4) + (1 - lambda[i])*(((a^2*var)/2)*(b)-1)
  loss_1 <-lambda[i]*(ols) + (1 - lambda[i])*(linex) 
  print(paste("Risk for Model 1:",risk_1))
  print(paste("Loss for Model 1:",loss_1))
}

# model 5
b <- t(l) %*% solve(t(X) %*% X) %*% l
beta_5 <- solve(t(X) %*% X) %*% t(X) %*% y
theta <- t(l) %*% beta
theta_hat <- t(l) %*% beta_5
ols <- t(theta - theta_hat) %*% (theta - theta_hat)
linex <- exp(a*t(X %*% beta - X %*% beta_5) %*% (X %*% beta - X %*% beta_5)) - a*t(X %*% beta - X %*% beta_5) %*% (X %*% beta - X %*% beta_5) - 1
lambda <- c(0.2, 0.4, 0.6, 0.8)
for(i in 1:4){
  risk_1 <- lambda[i]*(var * b) + (1 - lambda[i])*(exp((n-p)*var*a*a)/2 - 1)
  loss_1 <-lambda[i]*(ols) + (1 - lambda[i])*(linex) 
  print(paste("Risk for Model 1:",risk_1))
  print(paste("Loss for Model 1:",loss_1))
}

# model 6
b <- t(l) %*% solve(t(X) %*% X) %*% l
beta_6 <- solve(t(X) %*% X) %*% t(X) %*% y
theta <- t(l) %*% beta
theta_hat <- t(l) %*% beta_6
ols <- t(theta - theta_hat) %*% (theta - theta_hat)
linex <- exp(a*t(beta - beta_6) %*% (beta - beta_6)) - a*t(beta - beta_6) %*% (beta - beta_6) - 1
lambda <- c(0.2, 0.4, 0.6, 0.8)
for(i in 1:4){
  risk_1 <- lambda[i]*(var * b) + (1 - lambda[i])*(exp(var*a*a*b)/2 - 1)
  loss_1 <-lambda[i]*(ols) + (1 - lambda[i])*(linex) 
  print(paste("Risk for Model 1:",risk_1))
  print(paste("Loss for Model 1:",loss_1))
}

# model 7
b <- t(l) %*% solve(t(X) %*% X) %*% l
beta_7 <- solve(t(X) %*% X) %*% t(X) %*% y - (a*var)/2 * (solve(t(X) %*% X) %*% l)
theta <- t(l) %*% beta
theta_hat <- t(l) %*% beta_7 
ols <- t(theta - theta_hat) %*% (theta - theta_hat)
linex <- exp(a*t(X %*% beta - X %*% beta_7) %*% (X %*% beta - X %*% beta_7)) - a*t(X %*% beta - X %*% beta_7) %*% (X %*% beta - X %*% beta_7) - 1
lambda <- c(0.2, 0.4, 0.6, 0.8)
for(i in 1:4){
  risk_1 <- lambda[i]*(var * a^2*(var*b)^2/4) + (1 - lambda[i])*((a^2*var*(n-p))/2)
  loss_1 <-lambda[i]*(ols) + (1 - lambda[i])*(linex) 
  print(paste("Risk for Model 1:",risk_1))
  print(paste("Loss for Model 1:",loss_1))
}

# model 8
b <- t(l) %*% solve(t(X) %*% X) %*% l
beta_8 <- solve(t(X) %*% X) %*% t(X) %*% y - (a*var)/2 * (solve(t(X) %*% X) %*% l)
theta <- t(l) %*% beta
theta_hat <- t(l) %*% beta_8
ols <- t(theta - theta_hat) %*% (theta - theta_hat)
linex <- exp(a*t(beta - beta_8) %*% (beta - beta_8)) - a*t(beta - beta_8) %*% (beta - beta_8) - 1
lambda <- c(0.2, 0.4, 0.6, 0.8)
for(i in 1:4){
  risk_1 <- lambda[i]*(var*b + a^2*(var*b)^2/4) + (1 - lambda[i])*(a*var*b/2)
  loss_1 <-lambda[i]*(ols) + (1 - lambda[i])*(linex) 
  print(paste("Risk for Model 1:",risk_1))
  print(paste("Loss for Model 1:",loss_1))
}



data <- data[20600:20640,]
X_new <- as.matrix(cbind(1, data[,2:9]))
y <- X_new %*% beta
y_new <- X_new %*% beta_hat             
errors <- (y-y_new)
print(errors)
lines(y, color='red')
lines(y_new, color = 'blue')


plot_data <- data.frame(Index = data[, 1], 
                        y = c(y), 
                        y_new = c(y_new))

plot_data <- data.frame(Index = data[, 1], 
                        y = as.numeric(y), 
                        y_new = as.numeric(y_new))

# Use ggplot to plot y and y_new
gg <- ggplot() +
  geom_line(data = plot_data, aes(x = Index, y = y), color = "blue") +
  geom_line(data = plot_data, aes(x = Index, y = y_new), color = "red") +
  labs(x = "Index", y = "Values", title = "Comparison of y and y_new") +
  theme_minimal()

# Print the plot
print(gg)
