# Run auxiliary functions
source("code/aux_functions_random_beta.R")

# Load samples
load("samples/gibbs_samples_random_beta.RData")
load("samples/parallel_temp_samples_random_beta.RData")
load("samples/temp_tran_samples_random_beta.RData")

# Simulate data
K <- 4  # number of components
p_0 <- rep(1/K, K)  # weight parameter
mu_0 <- c(-3, 0, 3, 6)  # mean parameter
sigma2_0 <- rep(0.55^2, K)  # variance parameter
N <- 200  # sample size
load("data/mix_norm_data.RData")

# Relabel the samples by the order of mu
relabelled_chains_gibbs_samples <- 
  relabel_by_mu_order(gibbs_samples$chains$mu,
                      gibbs_samples$chains$p,
                      gibbs_samples$chains$sigma2,
                      K)

# Extract the samples from the parallel tempering
chains_par_temp <- samples_parallel_temp$chains
# Extract the cold chain (beta = 1)
chain_par_temp_cold_p <- chains_par_temp[[1]] %>%
  select(starts_with("p"))
chain_par_temp_cold_mu <- chains_par_temp[[1]] %>%
  select(starts_with("mu"))
chain_par_temp_cold_sigma2 <- chains_par_temp[[1]] %>%
  select(starts_with("sigma2"))
# Relabel the chains by the order of mu
relabelled_chains_par_temp <- relabel_by_mu_order(chain_par_temp_cold_mu, 
                                                  chain_par_temp_cold_p, 
                                                  chain_par_temp_cold_sigma2,
                                                  K)

# Extract the chains form the temperated transition
chain_temp_tran <- samples_temp_tran$chain
relabelled_chains_temp_tran <- relabel_by_mu_order(
  chain_temp_tran[, grep("^mu", colnames(chain_temp_tran))],
  chain_temp_tran[, grep("^p", colnames(chain_temp_tran))],
  chain_temp_tran[, grep("^sigma2", colnames(chain_temp_tran))],
  K
)

# Function to plot posterior distributions
plot_post <- function(samples, param_name, true_values) {
  K <- length(true_values)
  # Convert the samples to a long format for ggplot
  samples_long <- samples %>%
    select(starts_with(param_name)) %>%
    mutate(iteration = row_number()) %>%
    pivot_longer(-iteration, names_to = "parameter", values_to = "value")
  # Filter the samples to get the relevant parameter
  samples <- samples %>%
    select(starts_with(param_name))
  name_legend <- ifelse(param_name == "sigma2",
                        expression(sigma^2),
                        ifelse(param_name == "mu", expression(mu),
                               expression(p))
                        )
  dens_plot <- ggplot(samples_long, aes(x = value, color = parameter, fill = parameter)) +
    geom_density(alpha = 0.3) +
    geom_vline(xintercept = true_values, color = "red", linetype = "dashed") +
    labs(title = "",
         x = name_legend,
         y = "Posterior Density") +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_manual(values = hue_pal()(K))
}

# Function to plot trace plots
plot_trace <- function(samples, param_name, true_values) {
  K <- length(true_values)
  # Convert the samples to a long format for ggplot
  samples_long <- samples %>%
    select(starts_with(param_name)) %>%
    mutate(iteration = row_number()) %>%
    pivot_longer(-iteration, names_to = "parameter", values_to = "value")
  names_par <- paste0(param_name, 1:K)
  name_legend <- ifelse(param_name == "sigma2",
                        expression(sigma^2),
                        ifelse(param_name == "mu", expression(mu),
                               expression(p))
  )
  trace_plot <- ggplot(samples_long, aes(x = iteration, y = value, color = parameter)) +
    geom_line() +
    labs(title = "",
         x = "Iteration",
         y = name_legend) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_manual(values = hue_pal()(K),
                       labels = names_par,
                       name = "Parameter") +
    geom_hline(yintercept = true_values, linetype = "dashed", color = "red")
}
  

#-------------------------------------------------------------------------------
# Plot the posterior distributions of mu, sigma^2, and p
#-------------------------------------------------------------------------------

# p plots
post_p_gibbs <- plot_post(relabelled_chains_gibbs_samples$p, "p", p_0)
post_p_par_temp <- plot_post(relabelled_chains_par_temp$p, "p", p_0)
post_p_temp_tran <- plot_post(relabelled_chains_temp_tran$p, "p", p_0)

# mu plots
post_mu_gibbs <- plot_post(relabelled_chains_gibbs_samples$mu, "mu", mu_0)
post_mu_par_temp <- plot_post(relabelled_chains_par_temp$mu, "mu", mu_0)
post_mu_temp_tran <- plot_post(relabelled_chains_temp_tran$mu, "mu", mu_0)

# sigma^2 plots
post_sigma2_gibbs <- plot_post(relabelled_chains_gibbs_samples$sigma2, "sigma2", sigma2_0)
post_sigma2_par_temp <- plot_post(relabelled_chains_par_temp$sigma2, "sigma2", sigma2_0)
post_sigma2_temp_tran <- plot_post(relabelled_chains_temp_tran$sigma2, "sigma2", sigma2_0)

#-------------------------------------------------------------------------------
# Plot the trace plots for p, mu, and sigma^2
#-------------------------------------------------------------------------------

# p trace plots
trace_p_gibbs <- plot_trace(gibbs_samples$chains$p, "p", p_0)
trace_p_par_temp <- plot_trace(chain_par_temp_cold_p, "p", p_0)
trace_p_temp_tran <- plot_trace(chain_temp_tran, "p", p_0)

# mu trace plots
trace_mu_gibbs <- plot_trace(gibbs_samples$chains$mu, "mu", mu_0)
trace_mu_par_temp <- plot_trace(chain_par_temp_cold_mu, "mu", mu_0)
trace_mu_temp_tran <- plot_trace(chain_temp_tran, "mu", mu_0)

# sigma^2 trace plots
trace_sigma2_gibbs <- plot_trace(gibbs_samples$chains$sigma2, "sigma2", sigma2_0)
trace_sigma2_par_temp <- plot_trace(chain_par_temp_cold_sigma2, "sigma2", sigma2_0)
trace_sigma2_temp_tran <- plot_trace(chain_temp_tran, "sigma2", sigma2_0)

#-------------------------------------------------------------------------------
# Plot the observations
#-------------------------------------------------------------------------------

data_plot <- ggplot(y_df, aes(x = y)) +
  geom_density(color = "skyblue", fill = "skyblue", alpha = 0.5) +
  labs(title = "",
       x = "y",
       y = "") +
  theme_bw()

# Save plot
ggsave("figures/data_plot_random_beta.png",
       plot = data_plot,
       width = 8, height = 6, dpi = 300)

#-------------------------------------------------------------------------------
# Combine posterior distributions 
#-------------------------------------------------------------------------------

post_p <- (post_p_gibbs + post_p_par_temp + post_p_temp_tran) +
  plot_layout(ncol = 3, guides = "collect")
post_mu <- (post_mu_gibbs + post_mu_par_temp + post_mu_temp_tran) +
  plot_layout(ncol = 3, guides = "collect")
post_sigma2 <- (post_sigma2_gibbs + post_sigma2_par_temp + post_sigma2_temp_tran) +
  plot_layout(ncol = 3, guides = "collect")

# Save posterior distribution plots
ggsave("figures/post_p_random_beta.png",
       plot = post_p,
       width = 12, height = 4, dpi = 300)
ggsave("figures/post_mu_random_beta.png",
       plot = post_mu,
       width = 12, height = 4, dpi = 300)
ggsave("figures/post_sigma2_random_beta.png",
       plot = post_sigma2,
       width = 12, height = 4, dpi = 300)

#-------------------------------------------------------------------------------
# Combine trace plots
#-------------------------------------------------------------------------------

trace_p <- (trace_p_gibbs + trace_p_par_temp + trace_p_temp_tran) +
  plot_layout(ncol = 1, guides = "collect") &
  theme(text = element_text(size = 16),        # Base text size
        axis.title = element_text(size = 18),  # Axis titles
        axis.text = element_text(size = 16),   # Axis tick labels
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))
trace_mu <- (trace_mu_gibbs + trace_mu_par_temp + trace_mu_temp_tran) +
  plot_layout(ncol = 1, guides = "collect") &
  theme(text = element_text(size = 16),        # Base text size
        axis.title = element_text(size = 18),  # Axis titles
        axis.text = element_text(size = 16),   # Axis tick labels
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))
trace_sigma2 <- (trace_sigma2_gibbs + trace_sigma2_par_temp + trace_sigma2_temp_tran) +
  plot_layout(ncol = 1, guides = "collect") &
  theme(text = element_text(size = 16),        # Base text size
        axis.title = element_text(size = 18),  # Axis titles
        axis.text = element_text(size = 16),   # Axis tick labels
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

# Save trace plots
ggsave("figures/trace_p_random_beta.png",
       plot = trace_p,
       width = 16, height = 12, dpi = 300)
ggsave("figures/trace_mu_random_beta.png",
       plot = trace_mu,
       width = 16, height = 12, dpi = 300)
ggsave("figures/trace_sigma2_random_beta.png",
       plot = trace_sigma2,
       width = 16, height = 12, dpi = 300)
