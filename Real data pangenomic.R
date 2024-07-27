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
  theta ~ beta(1, 1);  // Prior for theta
  for (i in 1:n_genes) {
    for (j in 1:n_individuals) {
      X[i, j] ~ bernoulli(theta[i]);
    }
  }
}
"

# Install and load necessary packages
if (!require(vegan)) install.packages("vegan", dependencies=TRUE)
if (!require(rstan)) install.packages("rstan", dependencies=TRUE)
if (!require(coda)) install.packages("coda", dependencies=TRUE)
if (!require(bayesplot)) install.packages("bayesplot", dependencies=TRUE)
if (!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if (!require(gridExtra)) install.packages("gridExtra", dependencies=TRUE)
if (!require(reshape2)) install.packages("reshape2", dependencies=TRUE)

library(vegan)
library(rstan)
library(coda)
library(bayesplot)
library(ggplot2)
library(gridExtra)
library(reshape2)

# Load a suitable dataset
data(dune)
X <- as.matrix(dune)

# Binarize the data (convert to presence/absence)
X <- ifelse(X > 0, 1, 0)

# Parameters
n_individuals <- nrow(X)
n_genes <- ncol(X)

# Prepare data for Stan
data_stan <- list(X = t(X), n_genes = n_genes, n_individuals = n_individuals)

# Compile and fit the Stan model
stan_fit <- stan(model_code = stan_model_code, data = data_stan, iter = 6000, warmup = 1000, chains = 3, seed = 1234)

# Extract samples from each chain without permuting
theta_samples <- extract(stan_fit, pars = "theta", permuted = FALSE)

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

# Create individual density plots (no true values for real data)
density_plots_list <- lapply(1:n_genes, function(i) {
  ggplot(data.frame(x = as.vector(theta_samples[, , i])), aes(x = x)) +
    geom_density(fill = "skyblue", alpha = 0.5) +
    ggtitle(paste("Density Plot for Theta[", i, "]", sep = "")) +
    labs(x = expression(theta[i]), y = "Density")
})

# Combine all density plots into one
density_plots <- gridExtra::marrangeGrob(density_plots_list, ncol = 4, nrow = 5)

# Autocorrelation plots
autocorr_plots <- mcmc_acf(stan_fit, pars = paste0("theta[", 1:n_genes, "]"))

# Posterior means
posterior_means <- colMeans(as.matrix(theta_samples[, , 1:n_genes]))

# Plot posterior means
posterior_means_df <- data.frame(Gene = 1:n_genes, PosteriorMean = posterior_means)
plot_posterior_means <- ggplot(posterior_means_df, aes(x = Gene, y = PosteriorMean)) +
  geom_point() +
  labs(title = "Posterior Means for Gene Presence Probabilities", x = "Gene", y = "Posterior Mean")

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
ggsave("posterior_means.png", plot = plot_posterior_means)
ggsave("raw_data_plot.png", plot = raw_data_plot)

# Print results
print(gelman_diag)
print(ess_table)

# Display plots
print(trace_plots)
print(density_plots)
print(plot_posterior_means)
print(raw_data_plot)
