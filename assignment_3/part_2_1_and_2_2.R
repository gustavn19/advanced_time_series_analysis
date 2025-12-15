# Fetching the library and verify that is is loaded properly
# install.packages('ctsmTMB', repos = c('https://phillipbvetter.r-universe.dev'), type="source")
library(ctsmTMB)
############################################################
# Construction of proper dataframe
############################################################
construct_dataframe <- function(data_path){
dat <- read.csv(data_path)

data <- data.frame(
  t = seq(0, length(dat$rainfall) - 1) / 60,
  u = dat$rainfall,
  y = dat$stormwater
)

return(data)
}
############################################################


reservoir_state_model_fit_one_intermediary <- function(n, sigma, data_path){

# Read the data
data <- construct_dataframe(data_path)


############################################################
# Model creation
############################################################

# Create model object
model = ctsmTMB$new()

############################################################
# Test function, looping based
for (i in 1:(n+1)) {
  # build state names
  if (i < n) {
    x_i  <- as.name(paste0("x",  i))
    dx_i <- as.name(paste0("dx", i))
    dw_i <- as.name(paste0("dw", i)) 
    } else if (i == n) {
      x_i  <- as.name("xn")
      dx_i <- as.name("dxn")
      dw_i <- as.name("dwn")
  } else {
  x_i  <- as.name("xn") # Outflow from previous state
  dx_i <- as.name("dxn1") 
  dw_i <- as.name("dwn1")
}
  
  if (i == 1) {
    
    form <- eval(bquote(
      .(dx_i) ~ (A * u - n/k * .(x_i)) * dt + sigma * .(dw_i)
    ))
    
  } else if (i == n + 1) {
    
    form <- eval(bquote(
      .(dx_i) ~ n/k * .(x_i) * dt + sigma * .(dw_i)
    ))
    
} else {
    x_prev <- as.name(paste0("x", i - 1))
    
    form <- eval(bquote(
      .(dx_i) ~ n/k * (.(x_prev) - .(x_i)) * dt + sigma * .(dw_i)
    ))
  }
  
  model$addSystem(form)
}


############################################################
# Add system equations, for only one intermediate
# The following is only going to be for one intermediate
############################################################
# model = ctsmTMB$new()
# model$addSystem(
#   dx1 ~ (A * u - n / k * x1) * dt + sigma * dw1,
#   dx2 ~ n / k * x1 * dt + sigma * dw2
# )
# model$addObs(
#   y ~ x2
# )
# model$setVariance(
#   y ~ sigma_y^2
# )
# model$addInput(u)
# model$setParameter(
#   A   = c(initial = 1, lower=1e-5, upper=200), # Area parameters, estimate
#   k      = c(initial=1.5, lower=0, upper=100000),
#   # Trial for one intermediate
#   sigma = 0.1,
#   n = 1,
#   sigma_y = 0.1 # can make a range
# )
# # Trial for one initial
# model$setInitialState(list(c(0, 0), 1e-1*diag(2)))
# fit <- model$estimate(data, method="ekf")
# fit$cov.fixed
############################################################
# Here the only one intermediate model ends
############################################################

# Add observation equations, should this be something else? The y of the equation
model$addObs(
  y ~ xn1
)

# Set observation equation variances
model$setVariance(
  y ~ sigma_y^2
)

# Add vector input
model$addInput(u) # Controllable disturbance

# Specify parameter initial values and lower/upper bounds in estimation, specify these
model$setParameter(
  A   = c(initial = 1, lower=1e-5, upper=50), # Area parameters, estimate
  k      = c(initial=1.5, lower=0, upper=100),
  sigma = sigma, # Redefine, wrap function
  n = n,
  # Trial for one intermediate
  # sigma = 0.1,
  # n = 1,
  sigma_y = 0.1 # can make a range
)

# Set initial state mean and covariance
initial <- rep(0, times = n+1)
model$setInitialState(list(initial, 1e-1*diag(n+1))) # assumption of initial states



############################################################
# Model estimation
############################################################

# Carry out estimation with default settings (extended kalman filter)
fit <- model$estimate(data, method="ekf")
# Check parameter estimates against truth
print("The model is the following ")
print(model)

############################################################
# Defining the correlation matrix
############################################################
fit$corr_matrix <- matrix(0, nrow = length(fit$par.fixed), ncol = length(fit$par.fixed))

for (row in 1:length(fit$par.fixed)){
  for (col in 1:length(fit$par.fixed)) {
    fit$corr_matrix[row, col] <- fit$cov.fixed[row, col] / (fit$sd.fixed[row] * fit$sd.fixed[col])
  }
}

fit$corr_matrix

return(fit)
}


tested_variances = c(0.1, 0.2, 0.3, 0.4)
extract_stats <- function(n){nll_list <- numeric(length(tested_variances))
A_list <- numeric(length(tested_variances))
k_list <- numeric(length(tested_variances))
for (sigma in 1:length(tested_variances)){
fit <- reservoir_state_model_fit_one_intermediary(n, tested_variances[sigma], "~/Desktop/Desktop/mans_redin_DTU/course_material/advanced_time_series_analysis/exercise3/02427_CE03_Rainfall_runoff_exercise_data/ex1_rainfallrunoff.csv")

# summary(fit,extended=TRUE)

# Plotting the residuals - normally distributed
# plot(fit)
nll_list[sigma] <- fit$nll
A_list[sigma] <- fit$par.fixed[1]
k_list[sigma] <- fit$par.fixed[2]
}

stats_df = data.frame(
  var <- tested_variances,
  nll <- nll_list,
  A <- A_list,
  k <- k_list
)
names(stats_df) <- c("var", "nll", "A", "k")

return(stats_df)
}

# Main function for investigating different amount of intercediates

all_A <- matrix(0, 6, 4)
all_k <- matrix(0, 6, 4)
all_nll <- matrix(0, 6, 4)
for (i in 1:6){
  stats_df_extracted <- extract_stats(i)
  print("Testing number of intermediates:")
  print(i)
  all_A[i, ] <- stats_df_extracted$A
  all_k[i, ] <- stats_df_extracted$k
  all_nll[i, ] <- stats_df_extracted$nll
}

plot_seq <- seq(1, 6)
col_list <- c("red", "green", "blue", "grey")

############################################################
# Plot of all A values
for (p in 1:4){
  if (p==1){
    plot(plot_seq, all_A[, p], col = col_list[p], type="l", xlab="Number of intermediate states (n)", ylab="A", main = "Estimate of A in linear reservoir model")
  } else {
    lines(plot_seq, all_A[, p], col = col_list[p], type="l")
  }
}
legend(x = "topright", legend=c("sigma = 0.1", "sigma = 0.2", "sigma = 0.3", "sigma = 0.4"), fill = c("red", "green", "blue", "grey"))


############################################################
# Plot of all k values
for (p in 1:4){
  if (p==1){
    plot(plot_seq, all_k[, p], col = col_list[p], type="l", xlab="Number of intermediate states (n)", ylab="k", main = "Estimate of k in linear reservoir model")
  } else {
    lines(plot_seq, all_k[, p], col = col_list[p], type="l")
  }
}
legend(x = "topright", legend=c("sigma = 0.1", "sigma = 0.2", "sigma = 0.3", "sigma = 0.4"), fill = c("red", "green", "blue", "grey"))

# Plot of all nll values
for (p in 1:4){
  if (p==1){
    plot(plot_seq, all_nll[, p], col = col_list[p], type="l", xlab="Number of intermediate states (n)", ylab="nll", main = "Estimate of negative log likelihood (nll) in linear reservoir model", ylim=c(0, max(all_nll)))
  } else {
    lines(plot_seq, all_nll[, p], col = col_list[p], type="l")
  }
}
legend(x = "topright", legend=c("sigma = 0.1", "sigma = 0.2", "sigma = 0.3", "sigma = 0.4"), fill = c("red", "green", "blue", "grey"))

# Plot of all nll values, indication of best value

for (p in 1:4){
  if (p==1){
    plot(plot_seq, all_nll[, p], col = col_list[p], type="l", xlab="Number of intermediate states (n)", ylab="nll", main = "Estimate of negative log likelihood (nll) in linear reservoir model", ylim=c(1200, 1800))# ylim=c(0, max(all_nll)))
  } else {
    lines(plot_seq, all_nll[, p], col = col_list[p], type="l")
  }
}
legend(x = "topright", legend=c("sigma = 0.1", "sigma = 0.2", "sigma = 0.3", "sigma = 0.4"), fill = c("red", "green", "blue", "grey"))

############################################################
# The best model, 2.1.5
############################################################


fit <- reservoir_state_model_fit_one_intermediary(3, c(0.1), "~/Desktop/Desktop/mans_redin_DTU/course_material/advanced_time_series_analysis/exercise3/02427_CE03_Rainfall_runoff_exercise_data/ex1_rainfallrunoff.csv")
plot(fit)
names(fit)
fit$corr_matrix

############################################################
# The best model, 2.1.5
############################################################
# Re-define the data 
data <- construct_dataframe("~/Desktop/Desktop/mans_redin_DTU/course_material/advanced_time_series_analysis/exercise3/02427_CE03_Rainfall_runoff_exercise_data/ex1_rainfallrunoff.csv")

model = ctsmTMB$new()
model$addSystem(
  dx1 ~ (A * u - n / k * x1* 1 / (1 + exp(-alpha * (x1 - beta)))) * dt + sigma * dw1,
  dx2 ~ n / k * (x1 * 1 / (1 + exp(-alpha * (x1 - beta))) - x2) * dt + sigma * dw2, 
  dx3 ~ n / k * (x2 - x3) * dt + sigma * dw3, 
  dx4 ~ n / k * x3 * dt + sigma * dw2
)
model$addObs(
  y ~ x4
)
model$setVariance(
  y ~ sigma_y^2
)
model$addInput(u)
model$setParameter(
  A   = c(initial = 1, lower=1e-5, upper=200), # Area parameters, estimate
  k      = c(initial=1.5, lower=0, upper=100000),
  # Trial for one intermediate
  alpha      = c(initial=1.5, lower=0, upper=500),
  beta      = c(initial=1.5, lower=0, upper=40),
  sigma = 0.1,
  n = 1,
  sigma_y = 0.1 # can make a range
)
# Trial for one initial
model$setInitialState(list(c(0, 0, 0, 0), 1e-1*diag(4)))
fit <- model$estimate(data, method="ekf")
fit$nll

############################################################
# Define a function for computing the correlation matrix

compute_correlation_matrix <- function(fit){
fit$corr_matrix <- matrix(0, nrow = length(fit$par.fixed), ncol = length(fit$par.fixed))

for (row in 1:length(fit$par.fixed)){
  for (col in 1:length(fit$par.fixed)) {
    fit$corr_matrix[row, col] <- fit$cov.fixed[row, col] / (fit$sd.fixed[row] * fit$sd.fixed[col])
  }
}
return(fit)
}

fit <- compute_correlation_matrix(fit)
fit$corr_matrix

############################################################
# Simulation of model, 2.1.5
############################################################
# In this part of the exercise, a sample simulation function will be constructed and utilized to analyze model simulation results

simulate_model <- function(data_path, model, nr_simulations = 10){
  
data <- construct_dataframe(data_path)
model_prediction <- model$predict(data, method="ekf")$states
model_simulation <- model$simulate(data, method="ekf", n.sims=nr_simulations)$states$x4$i0
model_simulation.sim <- model_simulation[,(nr_simulations/2 + 1):ncol(model_simulation)]
matplot(model_simulation$t.j, model_simulation.sim, type="l", lty="solid", col="grey70", main="Simulation of model", xlab="Time, t [h]", ylab="Stormwater")
lines(model_prediction$t.j, model_prediction$x4, col="blue", lwd=2)
lines(data$t, data$y, col="red")
legend(x = "topleft", legend=c("Model prediction", "Data", paste0(nr_simulations, " sample simulations")),col=c("blue", "red", "grey70"), lty=c(1, 1, 1))
}
# The following is the data_path
# "~/Desktop/Desktop/mans_redin_DTU/course_material/advanced_time_series_analysis/exercise3/02427_CE03_Rainfall_runoff_exercise_data/ex1_rainfallrunoff.csv"

# Run the function for the best model for the correct data path
simulate_model("~/Desktop/Desktop/mans_redin_DTU/course_material/advanced_time_series_analysis/exercise3/02427_CE03_Rainfall_runoff_exercise_data/ex1_rainfallrunoff.csv", model)

############################################################
# Section 2.2
############################################################
fit <- reservoir_state_model_fit_one_intermediary(3, c(0.1), "~/Desktop/Desktop/mans_redin_DTU/course_material/advanced_time_series_analysis/exercise3/02427_CE03_Rainfall_runoff_exercise_data/ex2_overflow.csv")
plot(fit)
names(fit)
fit$convergence
fit$par.fixed

############################################################
# Section 2.2 - Building of model with activation function
############################################################

# Read the new data
data <- construct_dataframe("~/Desktop/Desktop/mans_redin_DTU/course_material/advanced_time_series_analysis/exercise3/02427_CE03_Rainfall_runoff_exercise_data/ex2_overflow.csv")

activation_model = ctsmTMB$new()
activation_model$addSystem(
  dx1 ~ (A * u - n / k * x1 * 1 / (1 + exp(-alpha * (x1 - beta)))) * dt + sigma * dw1,
  dx2 ~ n / k * (x1 * 1 / (1 + exp(-alpha * (x1 - beta))) - x2) * dt + sigma * dw2, 
  dx3 ~ n / k * (x2 - x3) * dt + sigma * dw3, 
  dx4 ~ n / k * x3 * dt + sigma * dw2
)
activation_model$addObs(
  y ~ x4
)
activation_model$setVariance(
  y ~ sigma_y^2
)
activation_model$addInput(u)
activation_model$setParameter(
  A   = c(initial = 1, lower=1e-5, upper=200), # Area parameters, estimate
  k      = c(initial=1.5, lower=0, upper=100000),
  # Trial for one intermediate
  alpha      = c(initial=1.5, lower=0, upper=500),
  beta      = c(initial=1.5, lower=0, upper=40),
  sigma = 0.1,
  n = 1,
  sigma_y = 0.1 # can make a range
)
# Trial for one initial
activation_model$setInitialState(list(c(0, 0, 0, 0), 1e-1*diag(4)))
fit <- activation_model$estimate(data, method="ekf")
fit
plot(fit)
fit$par.fixed
# fit$cov.fixed
fit$sd.fixed

fit <- compute_correlation_matrix(fit)
fit$corr_matrix
############################################################
# Simulation of data based on model
############################################################

simulate_model("~/Desktop/Desktop/mans_redin_DTU/course_material/advanced_time_series_analysis/exercise3/02427_CE03_Rainfall_runoff_exercise_data/ex2_overflow.csv", activation_model)


