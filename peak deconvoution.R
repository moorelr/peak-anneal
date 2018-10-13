# functions
peak <- function(x, center, imax, beta.sq){
  imax*exp(-pi*(x-center)^2/beta.sq)
}
fit_ys <- function(xs, fit_peaks){
  ys_out <- rep(0, length(xs))
  for(i in 1:nrow(fit_peaks)){
    ys_out <- ys_out + peak(xs, center = fit_peaks$center[i], imax = fit_peaks$imax[i]
                            , beta.sq = fit_peaks$beta.sq[i])
  }
  return(ys_out)
}
cost <- function(data, fit){
  sse <- log10(sum((fit-data)^2))
  return(sse)
}
neighbor <- function(old_soln){
  #old_soln <- old_fit
  new_soln <- old_soln
  # Choose a random parameter and change randomly
  var <- c(sample(1:nrow(old_soln), 1), sample(1:ncol(old_soln), 1))
  change_var <- new_soln[var[1], var[2]]
  new_var <- change_var + (change_var * rnorm(1, mean = 0, sd = 0.01))
  new_var
  if(new_var > limits["max_vals",var[2]]){
    new_var <- limits["max_vals",var[2]]
  }
  else if(new_var < limits["min_vals",var[2]]){
    new_var <- limits["min_vals",var[2]]
  }
  new_soln[var[1], var[2]] <- new_var
  return(new_soln)
}
calc_ap <- function(cost_old, cost_new, temp){
  # cost_new <- new_cost
  # cost_old <- old_cost
  # temp <- 1
  ap <- exp((cost_old - cost_new)/temp)
  if(ap > 1){ap <- 1}
  return(ap)
}
anneal <- function(temp = 1, temp_min = 1E-6, alpha = 0.9, n_wander = 200){
  while(temp > temp_min){
    # random wander
    for(i in 1:n_wander){
      new_fit <- neighbor(old_fit)
      new_cost <- cost(ys, fit_ys(xs, new_fit))
      prob_accept <- calc_ap(cost_old = old_cost, cost_new = new_cost, temp)
      
      random <- runif(n = 1, min = 0, max = 1)
      if(prob_accept > random){
        old_fit <- new_fit
        old_cost <- new_cost
        
        if(FALSE){
          new_row <- data.frame(temp = temp, order = paste(order, collapse = " ")
                                , distance = old_dist, prob_accept = prob_accept)
          save_params <- rbind(save_params, new_row)
        }
        if(FALSE){
          print(data.frame(temp = temp
                           #, order = paste(order, collapse = " ")
                           , distance = old_dist
                           , prob_accept = prob_accept))
        }
      }
    }
    print(paste("T =", temp))
    temp <- temp*alpha
  }
  return(new_fit)
}

# Save .pdf and set up plots
pdf("C:/Users/Lowell Moore/Desktop/peak_anneal.pdf", useDingbats = FALSE
    , width = 8, height = 6
    )

# Generate random spectrum with some convolved peaks
xs <- seq(1, 100, length.out = 1000)
ys <- (rnorm(length(xs))
          + peak(xs, center = 40, imax = 20, beta.sq = 200)
          + peak(xs, center = 50, imax = 40, beta.sq = 200)
          + peak(xs, center = 60, imax = 30, beta.sq = 200)
       )

# Plot and add legend
plot(xs, ys, type = "l", xlab = "Energy", ylab = "Intensity")
legend(x = 0, y = 50, xjust = 0, yjust = 1
       , legend = c("Data", "Guess fit", "Guess residual", "Best Fit", "Peaks", "best residual")
       , pch    = c(NA,     NA,          1,                NA,         NA,      1)
       , col    = c("black", "red",      "blue",           rgb(0, 0.7, 0), "green", "green")
       , lty    = c(1,      1,           NA,               1,          1,       NA)
       , lwd    = c(1,      1,           1,                3,          1,       1)
       )

# Initial guess for fit
n_peaks <- 3
center_pos <- which(ys == max(ys))
center_init <- xs[center_pos]
betasq_init <- 50
imax_init <- ys[center_pos]

# Data frame to hold the fit
old_fit <- data.frame(matrix(0, ncol = 3, nrow = n_peaks))
colnames(old_fit) <- c("center", "imax", "beta.sq")
old_fit$center <- rep(center_init, n_peaks)
old_fit$imax <- rep(imax_init, n_peaks)
old_fit$beta.sq <- rep(betasq_init, n_peaks)

# Limits for fit
min_vals <- c(10, 2, 0.0001)
max_vals <- c(90, 100, 1000)
limits <- rbind(min_vals, max_vals)
colnames(limits) <- colnames(old_fit)

# Compare initial guess with synthetic data
calc_ys <- fit_ys(xs, old_fit)
lines(xs, calc_ys, col = "red")
ys_resid <- calc_ys - ys
points(xs, ys_resid, col = "blue")

# Calculate initial step for annealing "while" loop
old_cost <- cost(ys, fit_ys(xs, old_fit))
new_fit <- neighbor(old_fit)
new_cost <- cost(ys, fit_ys(xs, new_fit))

# Fitting loop
best_fit <- anneal(temp = 1, temp_min = 1E-6, alpha = 0.9, n_wander = 100)

# Draw best fit peaks
for(i in 1:nrow(best_fit)){
  calc_ys <- fit_ys(xs, best_fit[i,])
  lines(xs, calc_ys, col = "green", lwd = 1)
}

# Draw best fit, calculate residual, draw residual
calc_ys <- fit_ys(xs, best_fit)
lines(xs, calc_ys, col = rgb(0, 0.7, 0), lwd = 3)
ys_resid <- calc_ys - ys
points(xs, ys_resid, col = "green", cex = 0.5)

# Are the residuals distributed normally?
shapiro <- shapiro.test(ys_resid)
pstat <- paste("p =", round(shapiro$p.value, 2))
wstat <- paste("W =", round(shapiro$statistic, 2))

# Quantile-quantile plot of residuals
qqplot <- qqnorm(ys_resid)
legend(x = min(qqnorm(ys_resid, plot.it = FALSE)$x), y <- max(qqnorm(ys_resid, plot.it = FALSE)$y)
       , xjust = 0, yjust = 1
       , legend = c(pstat, wstat)
       , bty = "n"
)

# Density plot of residuals
plot(density(ys_resid), main = "Residual density", col = "blue")

# Is the result what it was supposed to be?
print(best_fit)

dev.off()
