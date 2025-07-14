# Load necessary libraries
library(ggplot2)
library(invgamma)
library(ggthemes)
library(reshape2)
library(MCMCpack)
library(matrixStats)
library(parallel)
library(dplyr)
library(coda)
library(tidyr)
library(patchwork)
library(scales)

#-------------------------------------------------------------------------------
# Conditional distributions for Gibbs sampling
#-------------------------------------------------------------------------------

# Conditional distribution of z given the rest.
zconditional <- function(y, p, mu, sigma2) {
  N <- length(y)
  K <- length(p)
  next_z <- rep(0, N)
  for (i in 1:N) {
    prob <- p * dnorm(y[i], mean = mu, sd = sqrt(sigma2))
    next_z[i] <- sample(1:K, size = 1, prob = prob)
  }
  return(next_z)
}

# Conditional distribution of p given the rest.
pconditional <- function(y, z, mu, sigma2, gamma_) {
  K <- length(mu)
  n_by_cluster <- sapply(X = 1:K, FUN = function(j) sum(z == j))
  next_p <- rgamma(n = K, shape = gamma_ + n_by_cluster, rate = 1)
  return(next_p/sum(next_p))
}

# Conditional distribution of mu given the rest.
muconditional <- function(y, z, p, sigma2, m, tau2) {
  K <- length(p)
  n_by_cluster <- sapply(X = 1:K, FUN = function(j) sum(z == j))
  sum_y_by_cluster <- sapply(X = 1:K, FUN = function(j) sum(y[z == j]))
  cond_mean <- (sqrt(tau2)^(-2) * m + sum_y_by_cluster * sigma2^(-1))/
    (sqrt(tau2)^(-2) + n_by_cluster * sigma2^(-1))
  cond_var <- 1/(sqrt(tau2)^(-2) + n_by_cluster * sigma2^(-1))
  next_mu <- rnorm(n = K, mean = cond_mean, sd = sqrt(cond_var))
  return(next_mu)
}

# Conditional distribution of sigma^2 given the rest.
sigma2conditional <- function(y, z, p, mu, a0, b0) {
  K <- length(p)
  n_by_cluster <- sapply(X = 1:K, FUN = function(j) sum(z == j))
  sum_square_by_cluster <- sapply(X = 1:K, FUN = function(j) sum((y[z == j] - 
                                                                    mu[j])^2))
  cond_shape <- a0 + n_by_cluster/2
  cond_rate <- b0 + sum_square_by_cluster/2
  next_sigma2 <- 1/rgamma(n = K, shape = cond_shape, rate = cond_rate)
  return(next_sigma2)
}

#-------------------------------------------------------------------------------
# Gibbs sampler
#-------------------------------------------------------------------------------

gibbs <- function(y, hyper, niterations = 10000, burn_in = 2000) {
  # Extract hyperparameters
  gamma_ <- hyper$gamma_
  m <- hyper$m
  tau2 <- hyper$tau2
  a0 <- hyper$a0
  b0 <- hyper$b0
  # Number of components
  K <- length(gamma_)
  # Initialize the Markov chain.
  current_p <- rgamma(n = K, shape = gamma_, scale = 1)
  current_mu <- rnorm(n = K, mean = m, sd = sqrt(tau2))
  current_sigma2 <- 1/rgamma(n = K, shape = a0, rate = b0)
  
  # The following matrices will store the entire trajectory of the Markov
  # chain.
  p_chain <- matrix(rep(current_p, niterations + burn_in), 
                    nrow = niterations + burn_in, 
                    byrow = TRUE)
  mu_chain <- matrix(rep(current_mu, niterations + burn_in), 
                     nrow = niterations + burn_in, 
                     byrow = TRUE)
  sigma2_chain <- matrix(rep(current_sigma2, niterations + burn_in), 
                         nrow = niterations + burn_in, 
                         byrow = TRUE)
  
  # Run the Markov chain.
  for (t in 2:(niterations + burn_in)) {
    current_z <- zconditional(y, current_p, current_mu, current_sigma2)
    current_p <- pconditional(y, current_z, current_mu, current_sigma2, gamma_)
    current_mu <- muconditional(y, current_z, current_p, current_sigma2, m, tau2)
    current_sigma2 <- sigma2conditional(y, current_z, current_p, current_mu, a0, b0)
    p_chain[t, ] <- current_p
    mu_chain[t, ] <- current_mu
    sigma2_chain[t, ] <- current_sigma2
    # Print progress every 100 iterations
    if (t %% 100 == 0) {
      cat("Iteration", t, "/", niterations + burn_in, "\n")
    }
  }
  
  # Remove the burn-in samples and convert to data frames
  p_chain <- data.frame(p_chain[(burn_in + 1):nrow(p_chain), ])
  mu_chain <- data.frame(mu_chain[(burn_in + 1):nrow(mu_chain), ])
  sigma2_chain <- data.frame(sigma2_chain[(burn_in + 1):nrow(sigma2_chain), ])
  
  colnames(p_chain) <- paste0("p", 1:K)
  colnames(mu_chain) <- paste0("mu", 1:K)
  colnames(sigma2_chain) <- paste0("sigma2", 1:K)
  
  chains <- list(p = p_chain, mu = mu_chain, sigma2 = sigma2_chain)
  
  # Acceptance probabilities
  acceptance_probs <- list(
    mu = 1,
    sigma2 = 1,
    p = 1
  )
  # Compute ESS
  ess_p <- effectiveSize(as.mcmc(p_chain))
  ess_mu <- effectiveSize(as.mcmc(mu_chain))
  ess_sigma2 <- effectiveSize(as.mcmc(sigma2_chain))
  ess <- list(p = ess_p, mu = ess_mu, sigma2 = ess_sigma2)
  
  return(list(chains = chains, 
              acceptance_probs = acceptance_probs,
              ess = ess))
}

#-------------------------------------------------------------------------------
# Likelihood and prior functions
#-------------------------------------------------------------------------------

llik <- function(y, p, mu, sigma2) {
  N <- length(y)
  K <- length(p)
  lp <- log(p) - logSumExp(log(p))  
  dens_mat <- sapply(X = 1:K, 
                     FUN = function(k) {
                       dnorm(y, mean = mu[k], sd = sqrt(sigma2[k]), log = TRUE) +
                         lp[k]
                     })
  logL <- sum(apply(dens_mat, 1, logSumExp))  # Log-likelihood
  return(logL)
}

lprior_mu <- function(mu, m, tau2) {
  logL <- sum(dnorm(mu, mean = m, sd = sqrt(tau2), log = TRUE))  # Normal prior for mu
  return(logL)
}

lprior_sigma2 <- function(sigma2, a0, b0) {
  logL <- sum(dgamma(sigma2, shape = a0, scale = b0, log = TRUE)) -
    2 * sum(log(sigma2))  # Inverse Gamma prior for sigma^2
  return(logL)
}

lprior_p <- function(p, gamma_) {
  K <- length(p)
  logL <- dgamma(p, shape = gamma_, scale = 1, log = TRUE)
  return(logL)
}

#-------------------------------------------------------------------------------
# Parallel tempering
#-------------------------------------------------------------------------------

parallel_temp <- function(y, hyper, betas, n_iter = 10000, burn_in = 2000) {
  # Extract hyperparameters
  gamma_ <- hyper$gamma_
  m <- hyper$m
  tau2 <- hyper$tau2
  a0 <- hyper$a0
  b0 <- hyper$b0
  # Number of components
  K <- length(gamma_)
  # Create list to store the chains
  chains <- list()
  n_temp <- length(betas)
  # Acceptance counters
  accept_mu <- rep(0, n_temp)
  accept_sigma2 <- rep(0, n_temp)
  accept_p <- rep(0, n_temp)
  for (j in 1:n_temp) {
    # Initialize the Markov chain.
    beta <- betas[j]  # Current temperature
    chain <- matrix(NA, nrow = n_iter + burn_in, ncol = 3 * K + 1)
    current_p <- rgamma(n = K, shape = gamma_, scale = 1)
    current_mu <- rnorm(n = K, mean = 0, sd = 1)
    current_sigma2 <- 1 / rgamma(n = K, shape = 2, rate = 1)
    # Store the initial values
    chain[1, ] <- c(current_p, current_mu, current_sigma2, beta)
    chains[[j]] <- chain  # Store the chain for this temperature
  }
  
  # Run the Markov chain.
  for (t in 2:(n_iter + burn_in)) {
    for (j in 1:n_temp) {
      beta <- betas[j]  # Current temperature
      # current_chain <- chains[[j]]
      current_p <- as.numeric(chains[[j]][t-1, 1:K])
      current_mu <- as.numeric(chains[[j]][t-1, (K + 1):(2 * K)])
      current_sigma2 <- as.numeric(chains[[j]][t-1, (2 * K + 1):(3 * K)])
      # Propose p from reflective random Walk Metropolis
      proposed_p_raw <- current_p + rnorm(n = K, mean = 0, sd = 0.5)
      proposed_p <- abs(proposed_p_raw)  # Ensure positive p
      lacc_p <- (llik(y, proposed_p, current_mu, current_sigma2) + 
                   lprior_p(proposed_p, gamma_) - 
                   llik(y, current_p, current_mu, current_sigma2) - 
                   lprior_p(current_p, gamma_)) * beta
      lacc_p <- min(0, lacc_p)  # Ensure non-positive acceptance ratio
      if (log(runif(1)) < lacc_p) {
        current_p <- proposed_p  # Accept the proposal
        if (t > burn_in) {
          accept_p[j] <- accept_p[j] + 1  # Count acceptance
        }
      }
      # Propose mu from random Walk Metropolis
      proposed_mu <- current_mu + rnorm(n = K, mean = 0, sd = 0.1)
      lacc_mu <- (llik(y, current_p, proposed_mu, current_sigma2) + 
                    lprior_mu(proposed_mu, m, tau2) - 
                    llik(y, current_p, current_mu, current_sigma2) - 
                    lprior_mu(current_mu, m, tau2)) *
        beta
      lacc_mu <- min(0, lacc_mu)  # Ensure non-positive acceptance ratio
      if (log(runif(1)) < lacc_mu) {
        current_mu <- proposed_mu  # Accept the proposal
        if (t > burn_in) {
          accept_mu[j] <- accept_mu[j] + 1  # Count acceptance
        }
      }
      # Propose sigma2 from reflective random Walk Metropolis
      proposed_sigma2_raw <- current_sigma2 + rnorm(n = K, mean = 0, sd = 0.03)
      proposed_sigma2 <- abs(proposed_sigma2_raw)  # Ensure positive variance
      lacc_sigma2 <- (llik(y, current_p, current_mu, proposed_sigma2) + 
                        lprior_sigma2(proposed_sigma2, a0, b0) - 
                        llik(y, current_p, current_mu, current_sigma2) - 
                        lprior_sigma2(current_sigma2, a0, b0)) * beta 
      lacc_sigma2 <- min(0, lacc_sigma2)  # Ensure non-positive acceptance ratio
      if (log(runif(1)) < lacc_sigma2) {
        current_sigma2 <- proposed_sigma2  # Accept the proposal
        if (t > burn_in) {
          accept_sigma2[j] <- accept_sigma2[j] + 1  # Count acceptance
        }
      }
      # Store the current state of the chain
      chains[[j]][t, ] <- c(current_p, current_mu, current_sigma2, beta)
    }
    
    # # Propose swaps between randomly selected pairs of chains
    if (runif(1) < 0.5) {
      idx <- sample(1:n_temp, size = 2)
      idx1 <- idx[1]
      idx2 <- idx[2]
      
      mu_idx1 <- as.numeric(chains[[idx1]][t, (K + 1):(2 * K)])
      mu_idx2 <- as.numeric(chains[[idx2]][t, (K + 1):(2 * K)])
      sigma2_idx1 <- as.numeric(chains[[idx1]][t, (2 * K + 1):(3 * K)])
      sigma2_idx2 <- as.numeric(chains[[idx2]][t, (2 * K + 1):(3 * K)])
      p_idx1 <- as.numeric(chains[[idx1]][t, 1:K])
      p_idx2 <- as.numeric(chains[[idx2]][t, 1:K])
      # Calculate the log acceptance ratio
      lacc_swap <- betas[idx1] * 
        (llik(y, p_idx2, mu_idx2, sigma2_idx2) +
           lprior_mu(mu_idx2, m, tau2) +
           lprior_sigma2(sigma2_idx2, a0, b0) +
           lprior_p(p_idx2, gamma_) -
           llik(y, p_idx1, mu_idx1, sigma2_idx1) -
           lprior_mu(mu_idx1, m, tau2) -
           lprior_sigma2(sigma2_idx1, a0, b0) -
           lprior_p(p_idx1, gamma_)) + 
        betas[idx2] *
        (llik(y, p_idx1, mu_idx1, sigma2_idx1) +
           lprior_mu(mu_idx1, m, tau2) +
           lprior_sigma2(sigma2_idx1, a0, b0) +
           lprior_p(p_idx1, gamma_) -
           llik(y, p_idx2, mu_idx2, sigma2_idx2) -
           lprior_mu(mu_idx2, m, tau2) -
           lprior_sigma2(sigma2_idx2, a0, b0) -
           lprior_p(p_idx2, gamma_))
      lacc_swap <- min(0, lacc_swap)  # Ensure non-positive acceptance ratio
      if (log(runif(1)) < lacc_swap) {
        # Accept the swap
        chains[[idx1]][t, (K + 1):(2*K)] <- mu_idx2
        chains[[idx1]][t, (2*K + 1):(3*K)] <- sigma2_idx2
        chains[[idx1]][t, 1:K] <- p_idx2
        chains[[idx2]][t, (K+1):(2*K)] <- mu_idx1
        chains[[idx2]][t, (2*K+1):(3*K)] <- sigma2_idx1
        chains[[idx2]][t, 1:K] <- p_idx1
      }
    }
    
    # Print progress every 100 iterations
    if (t %% 100 == 0) {
      cat("Iteration", t, "/", n_iter + burn_in, "\n")
    }
  }
  # Remove the burn-in samples
  chains <- lapply(chains, function(chain) {
    chain[(burn_in + 1):nrow(chain), ]
  })
  # Convert chains to data frames
  chains <- lapply(chains, function(chain) {
    chain <- as.data.frame(chain)
    colnames(chain) <- c(paste0("p", 1:K), 
                         paste0("mu", 1:K), 
                         paste0("sigma2", 1:K), 
                         "beta")
    return(chain)
  })
  # Normalize p to ensure it sums to 1
  chains <- lapply(chains, function(chain) {
    chain_p <- chain[, paste0("p", 1:K)]
    ldenom <- do.call(rbind, 
                      lapply(1:nrow(chain_p), 
                             function(i) logSumExp(log(chain_p[i, ]))))
    chain_p <- exp(log(chain_p) - ldenom)
    chain[, paste0("p", 1:K)] <- chain_p
    return(chain)
  })
  # Acceptance probabilities
  acceptance_probs <- list(
    mu = accept_mu / n_iter,
    sigma2 = accept_sigma2 / n_iter,
    p = accept_p / n_iter
  )
  # Compute ESS
  ess_p <- sapply(chains, function(chain) {
    effectiveSize(as.mcmc(chain[, grep("^p", colnames(chain))]))
  })
  ess_mu <- sapply(chains, function(chain) {
    effectiveSize(as.mcmc(chain[, grep("^mu", colnames(chain))]))
  })
  ess_sigma2 <- sapply(chains, function(chain) {
    effectiveSize(as.mcmc(chain[, grep("^sigma2", colnames(chain))]))
  })
  ess <- list(p = ess_p, mu = ess_mu, sigma2 = ess_sigma2)
  return(list(chains = chains, 
              acceptance_probs = acceptance_probs,
              ess = ess))
}

#-------------------------------------------------------------------------------
# Tempered transition
#-------------------------------------------------------------------------------

temp_tran <- function(y, hyper, betas, n_iter = 10000, burn_in = 2000) {
  # Extract hyperparameters
  gamma_ <- hyper$gamma_
  m <- hyper$m
  tau2 <- hyper$tau2
  a0 <- hyper$a0
  b0 <- hyper$b0
  # Number of components
  K <- length(gamma_)
  chain <- data.frame()
  n_temp <- length(betas)
  current_p <- rgamma(n = K, shape = gamma_, rate = 1)
  current_mu <- rnorm(n = K, mean = 0, sd = 1)
  current_sigma2 <- 1 / rgamma(n = K, shape = 2, rate = 1)
  # Create empty data.frame with 0 rows and the desired columns
  chain <- matrix(NA, nrow = n_iter + burn_in, ncol = 3 * K)
  # Store the initial state of the chain
  chain[1, ] <- c(current_p, current_mu, current_sigma2)
  # Accepted proposals counters
  accept_mu <- 0
  accept_sigma2 <- 0
  accept_p <- 0
  
  for (t in 2:(n_iter + burn_in)) {
    u <- runif(1)
    if (u < 0.5) {
      # Forward pass mu
      mu_up <- list()
      mu_up[[1]] <- current_mu
      for (j in 2:n_temp) {
        prop_mu <- mu_up[[j - 1]] + rnorm(K, 0, 0.08)
        lacc_mu <- (llik(y, current_p, prop_mu, current_sigma2) +
                      lprior_mu(prop_mu, m, tau2) -
                      llik(y, current_p, mu_up[[j - 1]], current_sigma2) -
                      lprior_mu(mu_up[[j - 1]], m, tau2)
                    ) * betas[j]
        lacc_mu <- min(0, lacc_mu)
        if (log(runif(1)) < lacc_mu) {
          mu_up[[j]] <- prop_mu  # Accept the proposal
        } else {
          mu_up[[j]] <- mu_up[[j - 1]]  # Reject the proposal
        }
      }
      
      # Backward pass mu
      mu_down <- list()
      mu_down[[n_temp]] <- mu_up[[n_temp]]
      for (j in (n_temp - 1):1) {
        prop_mu <- mu_down[[j + 1]] + rnorm(K, 0, 0.08)
        lacc_mu <- (llik(y, current_p, prop_mu, current_sigma2) + 
                      lprior_mu(prop_mu, m, tau2) - 
                      llik(y, current_p, mu_down[[j + 1]], current_sigma2) - 
                      lprior_mu(mu_down[[j + 1]], m, tau2)
                      ) * betas[j]
        lacc_mu <- min(0, lacc_mu)
        if (log(runif(1)) < lacc_mu) {
          mu_down[[j]] <- prop_mu  # Accept the proposal
        } else {
          mu_down[[j]] <- mu_down[[j + 1]]  # Reject the proposal
        }
      }
      
      # Accept/reject
      ltarget_current_mu <- llik(y, current_p, current_mu, current_sigma2) + 
        lprior_mu(current_mu, m, tau2)
      ltarget_prop_mu <- llik(y, current_p, mu_down[[1]], current_sigma2) +
        lprior_mu(mu_down[[1]], m, tau2)
      lforward_mu <- sum(sapply(2:n_temp, function(j)
        (llik(y, current_p, mu_up[[j]], current_sigma2) + 
           lprior_mu(mu_up[[j]], m, tau2) - 
           llik(y, current_p, mu_up[[j - 1]], current_sigma2) -
           lprior_mu(mu_up[[j - 1]], m, tau2)
           ) * betas[j]))
      lbackward_mu <- sum(sapply(2:n_temp, function(j)
        (llik(y, current_p, mu_down[[j - 1]], current_sigma2) + 
           lprior_mu(mu_down[[j - 1]], m, tau2) - 
           llik(y, current_p, mu_down[[j]], current_sigma2) -
           lprior_mu(mu_down[[j]], m, tau2)
           ) * betas[j]))
      lacc_mu <- ltarget_prop_mu - ltarget_current_mu + lbackward_mu - lforward_mu
      lacc_mu <- min(0, lacc_mu)
      if (log(runif(1)) < lacc_mu) {
        current_mu <- mu_down[[1]]  # Accept the proposal
        if (t > burn_in) {
          accept_mu <- accept_mu + 1  # Count acceptance
        }
      }
      
      # Forward pass sigma2
      sigma2_up <- list()
      sigma2_up[[1]] <- current_sigma2
      for (j in 2:n_temp) {
        prop_sigma2 <- current_sigma2 + rnorm(K, 0, 0.05)
        prop_sigma2 <- abs(prop_sigma2)  # Ensure positive variance
        lacc_sigma2 <- (llik(y, current_p, current_mu, prop_sigma2) + 
                          lprior_sigma2(prop_sigma2, a0, b0) - 
                          llik(y, current_p, current_mu, sigma2_up[[j - 1]]) - 
                          lprior_sigma2(sigma2_up[[j - 1]], a0, b0)) * betas[j]
        lacc_sigma2 <- min(0, lacc_sigma2)
        if (log(runif(1)) < lacc_sigma2) {
          sigma2_up[[j]] <- prop_sigma2  # Accept the proposal
        } else {
          sigma2_up[[j]] <- sigma2_up[[j - 1]]  # Reject the proposal
        }
      }
      # Backward pass sigma2
      sigma2_down <- list()
      sigma2_down[[n_temp]] <- sigma2_up[[n_temp]]
      for (j in (n_temp - 1):1) {
        prop_sigma2 <- sigma2_down[[j + 1]] + rnorm(K, 0, 0.05)
        prop_sigma2 <- abs(prop_sigma2)  # Ensure positive variance
        lacc_sigma2 <- (llik(y, current_p, current_mu, prop_sigma2) + 
                          lprior_sigma2(prop_sigma2, a0, b0) - 
                          llik(y, current_p, current_mu, sigma2_down[[j + 1]]) - 
                          lprior_sigma2(sigma2_down[[j + 1]], a0, b0)) * betas[j]
        lacc_sigma2 <- min(0, lacc_sigma2)
        if (log(runif(1)) < lacc_sigma2) {
          sigma2_down[[j]] <- prop_sigma2  # Accept the proposal
        } else {
          sigma2_down[[j]] <- sigma2_down[[j + 1]]  # Reject the proposal
        }
      }
      # Accept/reject
      ltarget_current_sigma2 <- llik(y, current_p, current_mu, current_sigma2) + 
        lprior_sigma2(current_sigma2, a0, b0)
      ltarget_prop_sigma2 <- llik(y, current_p, current_mu, sigma2_down[[1]]) +
        lprior_sigma2(sigma2_down[[1]], a0, b0)
      lforward_sigma2 <- sum(sapply(2:n_temp, function(j)
        (llik(y, current_p, current_mu, sigma2_up[[j]]) + 
           lprior_sigma2(sigma2_up[[j]], a0, b0) - 
           llik(y, current_p, current_mu, sigma2_up[[j - 1]]) -
           lprior_sigma2(sigma2_up[[j - 1]], a0, b0)
           ) * betas[j]))
      lbackward_sigma2 <- sum(sapply(2:n_temp, function(j)
        (llik(y, current_p, current_mu, sigma2_down[[j - 1]]) + 
           lprior_sigma2(sigma2_down[[j - 1]], a0, b0) - 
           llik(y, current_p, current_mu, sigma2_down[[j]]) -
           lprior_sigma2(sigma2_down[[j]], a0, b0)
           ) * betas[j]))
      lacc_sigma2 <- ltarget_prop_sigma2 - ltarget_current_sigma2 +
        lbackward_sigma2 - lforward_sigma2
      lacc_sigma2 <- min(0, lacc_sigma2)
      if (log(runif(1)) < lacc_sigma2) {
        current_sigma2 <- sigma2_down[[1]]  # Accept the proposal
        if (t > burn_in) {
          accept_sigma2 <- accept_sigma2 + 1  # Count acceptance
        }
      }
      # Forward pass p
      p_up <- list()
      p_up[[1]] <- current_p
      for (j in 2:n_temp) {
        prop_p <- current_p + rnorm(K, 0, 0.8)
        prop_p <- abs(prop_p)
        lacc_p <- (llik(y, prop_p, current_mu, current_sigma2) + 
                     lprior_p(prop_p, gamma_) - 
                     llik(y, p_up[[j - 1]], current_mu, current_sigma2) - 
                     lprior_p(p_up[[j - 1]], gamma_)) * betas[j]
        lacc_p <- min(0, lacc_p)
        if (log(runif(1)) < lacc_p) {
          p_up[[j]] <- prop_p  # Accept the proposal
        } else {
          p_up[[j]] <- p_up[[j - 1]]  # Reject the proposal
        }
      }
      # Backward pass p
      p_down <- list()
      p_down[[n_temp]] <- p_up[[n_temp]]
      for (j in (n_temp - 1):1) {
        prop_p <- p_down[[j + 1]] + rnorm(K, 0, 0.8)
        prop_p <- abs(prop_p)  # Ensure p is positive 
        lacc_p <- (llik(y, prop_p, current_mu, current_sigma2) + 
                     lprior_p(prop_p, gamma_) - 
                     llik(y, p_down[[j + 1]], current_mu, current_sigma2) - 
                     lprior_p(p_down[[j + 1]], gamma_)) * betas[j]
        lacc_p <- min(0, lacc_p)
        if (log(runif(1)) < lacc_p) {
          p_down[[j]] <- prop_p  # Accept the proposal
        } else {
          p_down[[j]] <- p_down[[j + 1]]  # Reject the proposal
        }
      }
      # Accept/reject
      ltarget_current_p <- llik(y, current_p, current_mu, current_sigma2) + 
        lprior_p(current_p, gamma_)
      ltarget_prop_p <- llik(y, p_down[[1]], current_mu, current_sigma2) +
        lprior_p(p_down[[1]], gamma_)
      lforward_p <- sum(sapply(2:n_temp, function(j)
        (llik(y, p_up[[j]], current_mu, current_sigma2) + 
           lprior_p(p_up[[j]], gamma_) - 
           llik(y, p_up[[j - 1]], current_mu, current_sigma2) -
           lprior_p(p_up[[j - 1]], gamma_)
           ) * betas[j]))
      lbackward_p <- sum(sapply(2:n_temp, function(j)
        (llik(y, p_down[[j - 1]], current_mu, current_sigma2) + 
           lprior_p(p_down[[j - 1]], gamma_) - 
           llik(y, p_down[[j]], current_mu, current_sigma2) -
           lprior_p(p_down[[j]], gamma_)
           ) * betas[j]))
      lacc_p <- ltarget_prop_p - ltarget_current_p + lbackward_p - lforward_p
      lacc_p <- min(0, lacc_p)
      if (log(runif(1)) < lacc_p) {
        current_p <- p_down[[1]]  # Accept the proposal
        if (t > burn_in) {
          accept_p <- accept_p + 1  # Count acceptance
        }
      }
    } else {
      prop_mu <- current_mu + rnorm(K, 0, 0.08)
      lacc_mu <- (llik(y, current_p, prop_mu, current_sigma2) + 
                    lprior_mu(prop_mu, m, tau2) - 
                    llik(y, current_p, current_mu, current_sigma2) - 
                    lprior_mu(current_mu, m, tau2))
      lacc_mu <- min(0, lacc_mu)  # Ensure non-positive acceptance ratio
      if (log(runif(1)) < lacc_mu) {
        current_mu <- prop_mu  # Accept the proposal
        if (t > burn_in) {
          accept_mu <- accept_mu + 1  # Count acceptance
        }
      }
      prop_sigma2 <- current_sigma2 + rnorm(K, 0, 0.05)
      prop_sigma2 <- abs(prop_sigma2)  # Ensure positive variance
      lacc_sigma2 <- (llik(y, current_p, current_mu, prop_sigma2) + 
                        lprior_sigma2(prop_sigma2, a0, b0) - 
                        llik(y, current_p, current_mu, current_sigma2) - 
                        lprior_sigma2(current_sigma2, a0, b0))
      lacc_sigma2 <- min(0, lacc_sigma2)  # Ensure non-positive acceptance ratio
      if (log(runif(1)) < lacc_sigma2) {
        current_sigma2 <- prop_sigma2  # Accept the proposal
        if (t > burn_in) {
          accept_sigma2 <- accept_sigma2 + 1  # Count acceptance
        }
      }
      prop_p_raw <- current_p + rnorm(K, 0, 0.8)
      prop_p <- abs(prop_p_raw)
      lacc_p <- (llik(y, prop_p, current_mu, current_sigma2) + 
                   lprior_p(prop_p, gamma_) - 
                   llik(y, current_p, current_mu, current_sigma2) - 
                   lprior_p(current_p, gamma_))
      lacc_p <- min(0, lacc_p)  # Ensure non-positive acceptance ratio
      if (log(runif(1)) < lacc_p) {
        current_p <- prop_p  # Accept the proposal
        if (t > burn_in) {
          accept_p <- accept_p + 1  # Count acceptance
        }
      }
    }
    # Store the current state of the chain
    chain[t, ] <- c(current_p, current_mu, current_sigma2)
    
    # Print progress every 100 iterations
    if (t %% 100 == 0) {
      cat("Iteration", t, "/", n_iter + burn_in, "\n")
    }
    
  }
  # Remove the burn-in samples
  chain <- chain[(burn_in + 1):nrow(chain), ]
  # Convert chain to data frame
  chain <- as.data.frame(chain)
  colnames(chain) <- c(paste0("p", 1:K), 
                       paste0("mu", 1:K), 
                       paste0("sigma2", 1:K))
  
  # Normalize the p values
  chain_p <- chain %>% 
    select(starts_with("p"))
  lden <- unlist(
    lapply(1:nrow(chain_p), function(i) {
      logSumExp(log(chain_p[i,]))
    })
  )
  chain_p <- exp(log(chain_p) - lden)
  chain[, paste0("p", 1:K)] <- chain_p
  
  # Acceptance rates
  acceptance_probs <- list(
    mu = accept_mu / n_iter,
    sigma2 = accept_sigma2 / n_iter,
    p = accept_p / n_iter
  )
  # Compute ESS
  ess_p <- effectiveSize(as.mcmc(chain[, grep("^p", colnames(chain))]))
  ess_mu <- effectiveSize(as.mcmc(chain[, grep("^mu", colnames(chain))]))
  ess_sigma2 <- effectiveSize(as.mcmc(chain[, grep("^sigma2", colnames(chain))]))
  ess <- list(p = ess_p, mu = ess_mu, sigma2 = ess_sigma2)
  return(list(chains = chain, 
              acceptance_probs = acceptance_probs,
              ess = ess))
}

#-------------------------------------------------------------------------------
# Relabelling function
#------------------------------------------------------------------------------

relabel_by_mu_order <- function(chain_mu, chain_p, chain_sigma2, K) {
  n_iter <- nrow(chain_mu)
  
  relabelled_mu <- matrix(NA, nrow = n_iter, ncol = K)
  relabelled_p <- matrix(NA, nrow = n_iter, ncol = K)
  relabelled_sigma2 <- matrix(NA, nrow = n_iter, ncol = K)
  
  for (i in 1:n_iter) {
    mu_sample <- as.numeric(chain_mu[i, 1:K])
    order_idx <- order(mu_sample)
    
    relabelled_mu[i, ] <- mu_sample[order_idx]
    relabelled_p[i, ] <- as.numeric(chain_p[i, 1:K])[order_idx]
    relabelled_sigma2[i, ] <- as.numeric(chain_sigma2[i, 1:K])[order_idx]
  }
  # Convert matrices to data frames
  relabelled_mu <- data.frame(relabelled_mu)
  colnames(relabelled_mu) <- paste0("mu", 1:K)
  relabelled_p <- data.frame(relabelled_p)
  colnames(relabelled_p) <- paste0("p", 1:K)
  relabelled_sigma2 <- data.frame(relabelled_sigma2)
  colnames(relabelled_sigma2) <- paste0("sigma2", 1:K)
  
  list(
    mu = relabelled_mu,
    p = relabelled_p,
    sigma2 = relabelled_sigma2
  )
}