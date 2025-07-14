# Run auxiliary functions
source("code/aux_functions_random_beta.R")

#-------------------------------------------------------------------------------
# Sample from the random beta model
#-------------------------------------------------------------------------------

# y_i ~ sum_{j=1}^k w_j * N(mu_j, sigma_j^2)
# w ~ Dirichlet(gamma)
# mu_j ~ N(m, tau^2)
# sigma_j^2 ~ InvGamma(alpha, beta)

# Load the data
load("data/mix_norm_data.RData")
y <- y_df$y  # Extract the observations

# Hyperparameters for p.
gamma_ <- rep(1,K)
# Hyperparameters for mu_k.
m <- 0
tau2 <- 1 
# Hyperparameters for sigma_k^2.
a0 <- 2
b0 <- 1 
# Create a list of hyperparameters
hyper <- list(gamma_ = gamma_, m = m, tau2 = tau2, a0 = a0, b0 = b0)

#-------------------------------------------------------------------------------
# Gibbs sampler
#-------------------------------------------------------------------------------

n_iter <- 150000
burn_in <- 50000
gibbs_samples <- gibbs(y, hyper, n_iter, burn_in)

# Save the samples
save(gibbs_samples, file = "samples/gibbs_samples_random_beta.RData")

#-------------------------------------------------------------------------------
# Parallel tempering
#-------------------------------------------------------------------------------

# Number of temperatures
n_temp <- 10
# Sequence of inverse temperatures
beta_seq <- seq(1,0.1, length.out = n_temp)
n_iter <- 150000
burn_in <- 50000
samples_parallel_temp <- parallel_temp(y, hyper, beta_seq, n_iter, burn_in)

# Save the samples
save(samples_parallel_temp, file = "samples/parallel_temp_samples_random_beta.RData")

#-------------------------------------------------------------------------------
# Tempered transition
#-------------------------------------------------------------------------------

# Number of temperatures
n_temp <- 20
# Sequence of inverse temperatures
beta_seq <- 2^(-(0:(n_temp-1)))
n_iter <- 150000
burn_in <- 50000

# Run the temperated transition
samples_temp_tran <- temp_tran(y, hyper, beta_seq, n_iter, burn_in)

# Save the samples
save(samples_temp_tran, file = "samples/temp_tran_samples_random_beta.RData")

