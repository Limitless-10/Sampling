# Metropolis-Hastings algorithm
set.seed(123) # For reproducibility

# Target distribution: Standard normal
target_density <- function(x) {
  dnorm(x, mean = 0, sd = 1)
}

# Metropolis-Hastings function
metropolis_hastings <- function(iterations, init, proposal_sd) {
  # Initialize
  samples <- numeric(iterations)
  samples[1] <- init
  acceptance_count <- 0
  
  for (t in 2:iterations) {
    # Current state
    current <- samples[t - 1]
    
    # Propose a new state
    proposal <- rnorm(1, mean = current, sd = proposal_sd)
    
    # Calculate acceptance ratio
    r <- target_density(proposal) / target_density(current)
    
    # Accept or reject
    if (runif(1) < min(1, r)) {
      samples[t] <- proposal
      acceptance_count <- acceptance_count + 1
    } else {
      samples[t] <- current
    }
  }
  
  # Return samples and acceptance rate
  list(samples = samples, acceptance_rate = acceptance_count / (iterations - 1))
}

# Run the Metropolis-Hastings algorithm
mh_result <- metropolis_hastings(iterations = 10000, init = 0, proposal_sd = 1)

# Extract samples and acceptance rate
samples <- mh_result$samples
acceptance_rate <- mh_result$acceptance_rate

# Display results
cat("Acceptance rate:", acceptance_rate, "\n")

# Plot the results
par(mfrow = c(2, 1))
hist(samples, breaks = 50, probability = TRUE, main = "Histogram of Samples", xlab = "x")
curve(dnorm(x, mean = 0, sd = 1), add = TRUE, col = "red", lwd = 2)
plot(samples, type = "l", main = "Trace Plot of Samples", xlab = "Iteration", ylab = "x")
