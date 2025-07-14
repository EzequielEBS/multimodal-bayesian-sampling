# Simulate data
K <- 4  # number of components
p_0 <- rep(1/K, K)  # weight parameter
mu_0 <- c(-3, 0, 3, 6)  # mean parameter
sigma2_0 <- rep(0.55^2, K)  # variance parameter
N <- 200  # sample size

y <- c()  # observations
for (i in 1:N) {
  z_i <- sample(1:K, size = 1, prob = p_0)  # sample component
  y_i <- rnorm(n = 1, mean = mu_0[z_i], sd = sqrt(sigma2_0[z_i]))  # sample observation
  y <- c(y, y_i)  # append to observations
}

# Create a data frame for the observations
y_df <- data.frame(y = y)

# Save the data frame to a rdata file
save(y_df, file = "data/mix_norm_data.RData")
