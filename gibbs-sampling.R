# Gibbs Sampling for Bivariate Normal Distribution
set.seed(123) # For reproducibility

# Parameters of the bivariate normal
rho <- 0.8 # Correlation coefficient
n_iter <- 10000 # Number of iterations

# Initialize the chain
x1 <- numeric(n_iter)
x2 <- numeric(n_iter)
x1[1] <- 0 # Initial value for x1
x2[1] <- 0 # Initial value for x2

# Gibbs sampling loop
for (t in 2:n_iter) {
  # Sample x1 given x2
  x1[t] <- rnorm(1, mean = rho * x2[t - 1], sd = sqrt(1 - rho^2))
  # Sample x2 given x1
  x2[t] <- rnorm(1, mean = rho * x1[t], sd = sqrt(1 - rho^2))
}

# Combine the samples
samples <- data.frame(x1 = x1, x2 = x2)

# Plot the results
par(mfrow = c(2, 2))
plot(x1, type = "l", main = "Trace Plot of x1", xlab = "Iteration", ylab = "x1")
plot(x2, type = "l", main = "Trace Plot of x2", xlab = "Iteration", ylab = "x2")
plot(samples$x1, samples$x2, pch = 20, main = "Scatter Plot of Samples", xlab = "x1", ylab = "x2")
hist(samples$x1, breaks = 50, probability = TRUE, main = "Histogram of x1", xlab = "x1")
