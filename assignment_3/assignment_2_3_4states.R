library(ctsmTMB)
library(ggplot2)
library(readr)


dat <- readr::read_csv("./ex3_largecase.csv")
#dat <- dplyr::filter(dat, Event_ID == 4)
nr_reservoirs <- 4
n <- nr_reservoirs - 1
data <- data.frame(
  t = seq(0, length(dat$Rainfall) - 1) / 12,
  u = dat$Rainfall * 60e-6 * 1e3 , # Rainfall scale to square meter
  y = dat$Volume /1000,
  pumpflow = dat$Pumpflow * 60 /1000 # S
)


# Create model object
model = ctsmTMB$new()


model$addSystem(
  dx1 ~ A * u * dt - n/k_1 * x1 *dt + sigma * dw1, # Ground level
  dx2 ~ n/k_1 * x1*dt-n/k_2 *x2*dt+sigma*dw2, # Cobined sewer system
  dx3 ~ n / k_2 * (1/(1 + exp(-alpha * (x2 - Beta)))) * x2 * dt - n/k_3 * x3 *dt + sigma * dw3, # Stormwater tunnel
  dx4 ~ n / k_3 * x3 * dt - pumpflow * dt + sigma * dw4 # Storage Tower
  
)

# Add observation equations
model$addObs(
  y ~ x4
)

# Set observation equation variances
model$setVariance(
  y ~ sigma_y^2
)

# Add vector input (rain)
model$addInput(u, pumpflow)

# Specify parameter initial values and lower/upper bounds in estimation, specify these
model$setParameter(
  A   = c(initial = 1, lower=1e-5, upper=50), # Area parameters, estimate
  k_1      = c(initial=1.5, lower=0, upper=50),
  k_2      = c(initial=1.5, lower=0, upper=50),
  k_3      = c(initial=1.5, lower=0, upper=50),
  sigma = c(initial=1.5, lower=0, upper=10), 
  n = n,
  sigma_y = c(initial=1.5, lower=0, upper=10),
  alpha = c(initial=5, lower=-10, upper=20),
  Beta = c(initial=2, lower=-4, upper=20)
  #Beta = 2
)


# Set initial state mean and covariance
model$setInitialState(list(c(0, 0, 0, 0), 1e-1*diag(4))) # assumption of initial states - zero water in each reservoir


# Carry out estimation with default settings (extended kalman filter)
fit <- model$estimate(data, method="ekf")
print("The model is the following ")
print(model)
plot(fit)
print("negative log likelihood")
fit$nll

fit$par.fixed

print(summary(fit))

pred_list <-model$predict(data)
pred_df <- pred_list[[1]]

data$predictions <- pred_df$x4


p_event1 <- ggplot(data, aes(x = t)) +
  geom_line(aes(y = y, color = "Observed"), size = 1.2, alpha = 0.8) +
  geom_line(aes(y = predictions, color = "Predicted"), size = 1, linetype = "dashed") +
  labs(title = "ctsmTMB Model Fit for Event 4", 
       subtitle = "4-state Reservoir model with overflow fitted on all data",
       x = "Time", y = "Volume (1000 * m³)") +
  scale_color_manual(values = c("Observed" = "blue", "Predicted" = "red")) +
  theme_bw() + 
  theme(legend.title = element_blank(), legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(-0.5, max(data$y, na.rm = TRUE) * 2.1))

p_event1

