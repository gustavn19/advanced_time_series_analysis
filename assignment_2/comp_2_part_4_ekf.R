##----------------------------------------------------------------
## EKF algorithm for use in Part 4 of computer exercise 2 in
## Advanced Time Series Analysis
##----------------------------------------------------------------
rm(list = ls())
plot.new()
frame()
##----------------------------------------------------------------
## Do the simulation here and keep the values in y
##----------------------------------------------------------------
std_v = 1
std_e = 1


nr_overall_simulation = 20
converged_a = rep(0, nr_overall_simulation)
for (simulation in 1:nr_overall_simulation){
  n_simulations <- 10000
  v <- rnorm(n_simulations, mean = 0, sd = std_v)
  e <- rnorm(n_simulations, mean = 0, sd = std_e)
  x <- rep(0, n_simulations)
  y <- rep(0, n_simulations)
  t_i <- rep(0, n_simulations)
  for (i in 1:n_simulations) {
    t_i[i] = i
    x[i+1] = 0.4 * x[i] + v[i]
    y[i] = x[i] + e[i]
  }

  ##----------------------------------------------------------------
  ## Estimation with the EKF
  ##----------------------------------------------------------------
  ## aInit : The starting guess of the AR coefficient estimate
  aInit <- -0.5
  ## aVarInit : The initial variance for estimation of the AR coefficient
  aVarInit <- 10
  ## sigma.v : Standard deviation of the system noise of x in the filter
  sigma.v <- 10
  
  ## Initialize
  ## Init the state vector estimate
  zt <- c(0,aInit)
  ## Init the variance matrices
  Rv <- matrix(c(sigma.v^2,0,0,0), ncol=2)
  ## sigma.e : Standard deviation of the measurement noise in the filter
  Re <- 1 
  
  ## Init the P matrix, that is the estimate of the state variance
  Pt <- matrix(c(Re,0,0,aVarInit), nrow=2, ncol=2)
  ## The state is [X a] so the differentiated observation function is
  Ht <- t(c(1,0))
  ## Init a vector for keeping the parameter a variance estimates
  aVar <- rep(NA,length(y))
  ## and keeping the states
  Z <- matrix(NA, nrow=length(y), ncol=2)
  Z[,1] <- zt
  
  ## The Kalman filtering
  for(t in 1:(length(y)-1))
  {
    ## Derivatives (Jacobians)
    Ft <- matrix(c(zt[2],0,zt[1],1), ncol=2)  # F_t-1
    # Ht does not change 
    
    ## Prediction step
    zt = c(zt[2]*zt[1],zt[2]) #z_t|t-1 f(z_t-1|t-1)
    Pt = Ft %*% Pt %*% t(Ft) + Rv #P_t|t-1
    
    ## Update step
    res = y[t] - zt[1] # the residual at time t
    St =  Ht %*% Pt %*% t(Ht) + Re # innovation covariance
    Kt = Pt %*% t(Ht) %*% St^-1 # Kalman gain
    zt = zt + Kt * res # z_t|t
    Pt = (diag(2) - Kt%*%Ht)%*%Pt #P_t|t
    
    ## Keep the state estimate
    Z[t+1,] <- zt
    ## Keep the P[2,2], which is the variance of the estimate of a
    aVar[t+1] <- Pt[2,2]
    
  }
  
  
  if (simulation == 1) {
    plot(t_i, Z[, 2], type = "l", col = "blue", ylim = c(-1, 1), xlab="t", ylab="a_t")
    lines(t_i, aVar, col = "orange")
  } else  {
    lines(t_i, Z[, 2], col = "blue")
    lines(t_i, aVar, col = "orange")
  }
  
  converged_a[simulation] = Z[t,2]

}
abline(h = 0.4, col = "red", lty=2)
legend(x = "bottomright", legend = c("a", "a_var", "Actual a"), col = c("blue", "orange", "red"), lty = c(1, 1, 2))
title(main = paste(nr_overall_simulation, "simulations for a_0 = ", aInit, ", a_var_0 = ", aVarInit, "and v_var_0 = ", sigma.v))
