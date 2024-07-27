# Install and load necessary packages
if(!require(rstan)) install.packages("rstan", dependencies=TRUE)
if(!require(coda)) install.packages("coda", dependencies=TRUE)
if(!require(bayesplot)) install.packages("bayesplot", dependencies=TRUE)
if(!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if(!require(gridExtra)) install.packages("gridExtra", dependencies=TRUE)

library(rstan)
library(coda)
library(bayesplot)
library(ggplot2)
library(gridExtra)

# Set up rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Stan model definition
stan_model_code <- "
data {
  int<lower=0> n_individuals;
  int<lower=0> n_genes;
  int<lower=0, upper=1> X[n_genes, n_individuals];
}
parameters {
  vector<lower=0, upper=1>[n_genes] theta;
}
model {
  for (i in 1:n_genes)
    theta[i] ~ beta(1, 1);
  for (i in 1<n_genes)
    for (j in 1<n_individuals)
      X[i, j] ~ bernoulli(theta[i]);
}
"

# Function to simulate data
simulate_data <- function(n_individuals, n_genes, true_theta) {
  X <- matrix(0, nrow = n_genes, ncol = n_individuals)
  for (i in 1:n_genes) {
    for (j in 1:n_individuals) {
      X[i, j] <- rbinom(1, 1, true_theta[i])
    }
  }
  return(X)
}

# Set seed for reproducibility
set.seed(1234)

# Parameters
n_individuals <- 100
n_genes <- 20
true_theta <- runif(n_genes, 0.1, 0.9)

# Simulate data
X <- simulate_data(n_individuals, n_genes, true_theta)

# Prepare data for Stan
data_stan <- list(X = X, n_genes = n_genes, n_individuals = n_individuals)

# Compile and fit the Stan model
stan_fit <- stan(model_code = stan_model_code, data = data_stan, iter = 6000, warmup = 1000, chains = 3, seed = 1234)

# Extract samples
samples <- extract(stan_fit)
theta_samples <- samples$theta

# Convert stan_fit to an mcmc.list object for coda
stan_mcmc_list <- As.mcmc.list(stan_fit, pars = "theta")

# MCMC diagnostics
gelman_diag <- gelman.diag(stan_mcmc_list)
effective_size <- effectiveSize(stan_mcmc_list)
ess_table <- data.frame(Gene = 1:n_genes, EffectiveSize = effective_size)

# Trace and density plots
trace_plots <- mcmc_trace(stan_fit, pars = paste0("theta[", 1:n_genes, "]"))
density_plots <- mcmc_dens(stan_fit, pars = paste0("theta[", 1:n_genes, "]"))

# Autocorrelation plots
autocorr_plots <- mcmc_acf(stan_fit, pars = paste0("theta[", 1:n_genes, "]"))

# Posterior means
posterior_means <- colMeans(theta_samples)

# Mean Squared Error
mse <- mean((posterior_means - true_theta)^2)

# Plot posterior means vs true values
posterior_vs_true <- data.frame(True = true_theta, Posterior = posterior_means)
plot_posterior_vs_true <- ggplot(posterior_vs_true, aes(x = True, y = Posterior)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Posterior Means vs True Values", x = "True Theta", y = "Posterior Mean")

# Save plots
ggsave("trace_plots.png", plot = trace_plots)
ggsave("density_plots.png", plot = density_plots)
ggsave("posterior_vs_true.png", plot = plot_posterior_vs_true)

# Print MSE
print(paste("Mean Squared Error: ", mse))

# Print results
print(gelman_diag)
print(ess_table)

# Display plots
print(trace_plots)
print(density_plots)
print(plot_posterior_vs_true)
