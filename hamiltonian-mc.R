# Hamiltonian Monte Carlo implementation
set.seed(123) # For reproducibility

# Target log-density and its gradient
log_density <- function(x) {
  -0.5 * x^2 # Standard normal: log-density
}
grad_log_density <- function(x) {
  -x # Gradient of the log-density
}

# HMC function
hmc <- function(iterations, init, epsilon, L, mass = 1) {
  samples <- numeric(iterations)
  samples[1] <- init
  accept_count <- 0
  
  for (t in 2:iterations) {
    # Current state
    x_current <- samples[t - 1]
    p_current <- rnorm(1, mean = 0, sd = sqrt(mass)) # Sample momentum
    
    # Simulate Hamiltonian dynamics
    x <- x_current
    p <- p_current - 0.5 * epsilon * grad_log_density(x) # Half-step for momentum
    
    for (i in 1:L) {
      x <- x + epsilon * p / mass # Full-step for position
      if (i != L) { # Avoid extra gradient calculation at the end
        p <- p - epsilon * grad_log_density(x) # Full-step for momentum
      }
    }
    p <- p - 0.5 * epsilon * grad_log_density(x) # Final half-step for momentum
    
    # Metropolis acceptance step
    current_H <- -log_density(x_current) + 0.5 * (p_current^2 / mass)
    proposed_H <- -log_density(x) + 0.5 * (p^2 / mass)
    acceptance_prob <- exp(current_H - proposed_H)
    
    if (runif(1) < acceptance_prob) {
      samples[t] <- x
      accept_count <- accept_count + 1
    } else {
      samples[t] <- x_current
    }
  }
  
  list(samples = samples, acceptance_rate = accept_count / (iterations - 1))
}

# Run HMC
hmc_result <- hmc(iterations = 10000, init = 0, epsilon = 0.1, L = 10)

# Extract samples and acceptance rate
samples <- hmc_result$samples
acceptance_rate <- hmc_result$acceptance_rate

# Display results
cat("Acceptance rate:", acceptance_rate, "\n")

# Plot the results
par(mfrow = c(2, 1))
hist(samples, breaks = 50, probability = TRUE, main = "Histogram of Samples", xlab = "x")
curve(dnorm(x, mean = 0, sd = 1), add = TRUE, col = "red", lwd = 2)
plot(samples, type = "l", main = "Trace Plot of Samples", xlab = "Iteration", ylab = "x")
