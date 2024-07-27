# Install and load necessary packages
if(!require(rstan)) install.packages("rstan", dependencies=TRUE)
if(!require(coda)) install.packages("coda", dependencies=TRUE)
if(!require(bayesplot)) install.packages("bayesplot", dependencies=TRUE)
if(!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if(!require(gridExtra)) install.packages("gridExtra", dependencies=TRUE)
if(!require(reshape2)) install.packages("reshape2", dependencies=TRUE)

library(rstan)
library(coda)
library(bayesplot)
library(ggplot2)
library(gridExtra)
library(reshape2)

# Set up rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Corrected Stan model definition
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
  theta ~ beta(1, 1); // Vectorized prior
  for (i in 1:n_genes) {
    for (j in 1:n_individuals) {
      X[i, j] ~ bernoulli(theta[i]);
    }
  }
}
"

# Function to simulate data
simulate_data <- function(n_individuals, n_genes, true_theta) {
  X <- matrix(0, nrow = n_genes, ncol = n_individuals)
  for (i in 1<n_genes) {
    for (j in 1<n_individuals) {
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

#####################

# Compile and fit the Stan model
stan_fit <- stan(model_code = stan_model_code, data = data_stan, iter = 6000, warmup = 1000, chains = 3, seed = 1234)

# Extract samples from each chain without permuting
theta_samples <- extract(stan_fit, pars = "theta", permuted = FALSE)

# Inspect the structure of theta_samples
str(theta_samples)

# Convert samples to mcmc.list
chains <- lapply(1:dim(theta_samples)[2], function(chain) {
  mcmc(as.matrix(theta_samples[, chain, 1:n_genes]))
})
stan_mcmc_list <- mcmc.list(chains)

# Calculate Gelman-Rubin diagnostics
gelman_diag <- gelman.diag(stan_mcmc_list)

# Extract the Point Estimate and Upper CI values
gelman_diag_values <- data.frame(
  Gene = 1:n_genes,
  Point_Estimate = gelman_diag$psrf[, 1],
  Upper_CI = gelman_diag$psrf[, 2]
)

# Save Gelman-Rubin diagnostics to a text file for LaTeX inclusion
write.table(gelman_diag_values, "gelman_diag.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# Calculate effective sample size
effective_size <- effectiveSize(stan_mcmc_list)
effective_size <- as.numeric(effective_size) # Remove row names if necessary

# Create a data frame for the effective sample size
ess_table <- data.frame(Gene = 1:n_genes, EffectiveSize = effective_size)

# Save effective sample size to a text file for LaTeX inclusion
write.table(ess_table, "ess_table.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# Trace and density plots
trace_plots <- mcmc_trace(stan_fit, pars = paste0("theta[", 1:n_genes, "]"))

# Create individual density plots with true values
density_plots_list <- lapply(1:n_genes, function(i) {
  ggplot(data.frame(x = as.vector(theta_samples[, , i])), aes(x = x)) +
    geom_density(fill = "skyblue", alpha = 0.5) +
    geom_vline(xintercept = true_theta[i], linetype = "dotted", color = "red") +
    ggtitle(paste("Density Plot for Theta[", i, "]", sep = "")) +
    labs(x = expression(theta[i]), y = "Density")
})

# Combine all density plots into one
density_plots <- gridExtra::marrangeGrob(density_plots_list, ncol = 4, nrow = 5)

# Autocorrelation plots
autocorr_plots <- mcmc_acf(stan_fit, pars = paste0("theta[", 1:n_genes, "]"))

# Posterior means
posterior_means <- colMeans(as.matrix(theta_samples[, , 1:n_genes]))

# Mean Squared Error
mse <- mean((posterior_means - true_theta)^2)

# Plot posterior means vs true values
posterior_vs_true <- data.frame(True = true_theta, Posterior = posterior_means)
plot_posterior_vs_true <- ggplot(posterior_vs_true, aes(x = True, y = Posterior)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Posterior Means vs True Values", x = "True Theta", y = "Posterior Mean")

# Plot of the raw data
X_melted <- melt(X)
colnames(X_melted) <- c("Gene", "Individual", "Presence")
raw_data_plot <- ggplot(X_melted, aes(x = Individual, y = Gene)) +
  geom_tile(aes(fill = factor(Presence)), color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +
  labs(title = "Raw Data Plot", x = "Individual", y = "Gene") +
  theme_minimal()

# Save plots
ggsave("trace_plots.png", plot = trace_plots)
ggsave("density_plots.png", plot = density_plots, width = 16, height = 12)
ggsave("posterior_vs_true.png", plot = plot_posterior_vs_true)
ggsave("raw_data_plot.png", plot = raw_data_plot)

# Print MSE
print(paste("Mean Squared Error: ", mse))

# Print results
print(gelman_diag)
print(ess_table)

# Display plots
print(trace_plots)
print(density_plots)
print(plot_posterior_vs_true)
print(raw_data_plot)
