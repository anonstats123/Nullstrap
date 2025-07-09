## Nullstrap
 
**Nullstrap** (Null Bootstrap) is an R package designed for variable selection with rigorous false discovery rate (FDR) control in high-dimensional settings. It generates synthetic null data based on model fits and compares original vs. null statistics through a unified thresholding framework. Nullstrap supports a variety of models, including **linear models**, **generalized linear models (GLM)**, **Cox models**, and **Gaussian graphical models**. Compared to existing FDR methods like knockoffs and data splitting, Nullstrap offers improved power, flexibility, and computational efficiency.


# Table of contents
1. [Installation](#installation-)
2. [Quick Start](#quick-start)
3. [Simulation Examples](#simulation-examples)
4. [Real Data Example](#real-data-examples)


## Installation<a name="installation-"></a>

To install the development version from GitHub, please run:

``` r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("anonstats123/Nullstrap")
```
We are now working on submitting it to Bioconductor and will provide the link once online.



## Quick Start<a name="quick-start"></a>

### Linear Model Example
The following code is a quick example of how to use **Nullstrap** for variable selection in a **linear model** setting. 

``` r
library(Nullstrap)
library(MASS)

n <- 2000
p <- 1000
nonzero_coefs <- 30
Amp <- 0.25
rho <- 0.6
Theta.8 <- toeplitz(rho^(0:(p - 1)))
X <- mvrnorm(n, rep(0, p), Sigma = Theta.8)
X <- scale(X)
beta <- rep(0, p)
beta[1:nonzero_coefs] <- sample(c(-Amp, Amp), nonzero_coefs, replace = TRUE)
true_labels <- beta != 0
Signal_index <- 1:nonzero_coefs
y <- X %*% beta + rnorm(n)
y <- y - mean(y)

result_lm <- nullstrap_filter(
  X, y, fdr_value = 0.1, best_lambda = NULL, B_reps = NULL, dist_type = "normal", 
  model_type = "linear"
)
```

### Generalized Linear Model Example

For a **generalized linear model (GLM)**, you can use the following code:

``` r
library(Nullstrap)
library(MASS)

n <- 3000  
p <- 500  
nonzero_coefs <- 30  
Amp <- 9
rho <- 0.6
Theta.8 <- toeplitz(rho^(0:(p - 1)))
X <- mvrnorm(n, rep(0, p), Sigma = Theta.8)
X <- scale(X)/(sqrt(n))
beta <- rep(0, p)
beta[1:nonzero_coefs] <- sample(c(-Amp, Amp), nonzero_coefs, replace = TRUE)
true_labels <- beta != 0
Signal_index <- 1:nonzero_coefs
linear_predictor <- X %*% beta
prob <- 1 / (1 + exp(-linear_predictor))  # logistic function
y <- rbinom(n, 1, prob)

result_glm <- nullstrap_filter(
  X, y, fdr_value = 0.1, best_lambda = NULL, B_reps = NULL, model_type = "glm"
)
```


### Cox Model Example

For a **Cox model**, you can use the following code:

``` r
library(Nullstrap)
library(MASS)
library(survival)
library(simsurv)

n <- 400
p <- 200
nonzero_coefs <- 30
Amp <- 5
beta <- rep(0, p) 
beta[1:nonzero_coefs] <- Amp
names(beta) <- paste0("X", 1:p)
rho <- 0.6
Theta.8 <- toeplitz(rho^(0:(p - 1)))
X <- mvrnorm(n, rep(0, p), Sigma = Theta.8)
X <- scale(X)/(sqrt(n))
colnames(X) <- paste0("X", 1:p)
Signal_index <- 1:nonzero_coefs
true_labels <- beta != 0
surv_time <- simsurv(lambdas = 1, gammas = 1, betas = beta, x = as.data.frame(X))
surv_time <- surv_time$eventtime
status <- rep(1, n)
time <- surv_time
y <- cbind(time = time, status = status)

result_cox <- nullstrap_filter(
  X, y, fdr_value = 0.1, best_lambda = NULL, B_reps = NULL, model_type = "cox"
)
```


### Gaussian Graphical Model Example

For a **Gaussian graphical model (GGM)**, you can use the following code:

``` r
library(Nullstrap)
library(MASS)

generate_precision_matrix <- function(p, b, graph_type) {
  # Initialize Omega matrix
  Omega_0 <- matrix(0, nrow = p, ncol = p)
  
  if (graph_type == "band") {
    # Band graph
    diag(Omega_0) <- 1
    for (i in 1:(p-1)) {
      for (j in (i+1):min(p, i+10)) {
        Omega_0[i, j] <- sign(b) * abs(b)^(abs(i - j) / 10)
        Omega_0[j, i] <- Omega_0[i, j]  # Symmetric matrix
      }
    }
    
  } else if (graph_type == "block") {
    # Block graph: 10 blocks, each block size 20
    for (k in seq(1, p, by = 20)) {
      block_size <- min(20, p - k + 1)
      Omega_0[k:(k+block_size-1), k:(k+block_size-1)] <- 1
      Omega_0[k:(k+block_size-1), k:(k+block_size-1)] <- b
    }
    
  } else if (graph_type == "erdos") {
    # Erdos-Renyi graph
    diag(Omega_0) <- 1
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        if (rbinom(1, 1, 0.1)) {
          samples_1 <- runif(n, min = -0.6, max = -0.2)  # [-0.6, -0.2]
          samples_2 <- runif(n, min = 0.2, max = 0.6)    # [0.2, 0.6]
          phi_ij = combined_samples <- sample(c(samples_1, samples_2), 1)
          # phi_ij <- sample(c(-0.6, -0.2, 0.2, 0.6), 1)
          Omega_0[i, j] <- Omega_0[j, i] <- phi_ij
        }
      }
    }
    
  } else if (graph_type == "cluster") {
    # Cluster graph: 5 blocks, each block size 40
    for (k in seq(1, p, by = 40)) {
      block_size <- min(40, p - k + 1)
      Omega_0[k:(k+block_size-1), k:(k+block_size-1)] <- 1
      for (i in k:(k+block_size-1)) {
        if (i < (k+block_size-1)) { # in case j exceeds range
          for (j in (i+1):(k+block_size-1)) {
            if (rbinom(1, 1, 0.5)) {
              samples_1 <- runif(n, min = -0.6, max = -0.2)  # [-0.6, -0.2]
              samples_2 <- runif(n, min = 0.2, max = 0.6)    # [0.2, 0.6]
              phi_ij = combined_samples <- sample(c(samples_1, samples_2), 1)
              Omega_0[i, j] <- Omega_0[j, i] <- phi_ij
            }
          }
        }
      }
    }
  }
  
  # Make the precision matrix positive definite
  lambda_min <- eigen(Omega_0)$values[p]
  Omega <- Omega_0 + (abs(lambda_min) + 0.5) * diag(p)
  
  return(Omega)
}

generate_data <- function(n, p, Omega) {
  # Generate the precision matrix Omega
  # Omega <- generate_precision_matrix(p, b, graph_type)
  
  # Generate the covariance matrix Sigma
  Sigma <- solve(Omega)
  
  # Generate n samples from N_p(0, Sigma)
  data <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  return(data)
}

# Function to obtain the true support of the precision matrix
obtain_true_support <- function(Omega, threshold = 1e-10) {
  # Obtain the true support of the precision matrix
  true_support = matrix(0, nrow = ncol(Omega), ncol = ncol(Omega))
  true_support[abs(Omega) > threshold] = 1
  return(true_support)
}

n <- 2000
p <- 200
edge_strength <- -0.8
Omega <- generate_precision_matrix(p, edge_strength, "band")
true_support <- obtain_true_support(Omega)
Signal_index <- which(true_support[upper.tri(true_support)] == 1)
true_labels <- true_support[upper.tri(true_support)]
X <- generate_data(n, p, Omega)
data_scale <- scale(X)

result_ggm <- nullstrap_filter(
  X = data_scale, fdr_value = 0.2, best_lambda = NULL, B_reps = NULL, 
  model_type = "graphical"
)
```

The parameters of `nullstrap_filter()` are:

-   `X`: The design matrix (features).
-   `y`: The response vector.
-   `fdr_value`: The desired false discovery rate level (e.g., 0.1).
-   `best_lambda`: The best lambda value for regularization. If not specified, it will be determined automatically to the default value of each model type.
-   `B_reps`: The number of repetitions for computing the correction factor. If not specified, it will be determined automatically to the default value of each model type.
-   `dist_type`: The type of distribution for the error term. It is only applicable for linear models.
-   `model_type`: The type of model to fit. Options include "linear", "glm", "cox", and "graphical".

The output of `scdesign3()` is a list which includes:

- `threshold`: The threshold value for variable selection.
- `selected`: A vector of selected variable indices.
- `statistics`: The test statistics for each variable.


## Simulation Examples<a name="simulation-examples"></a>

Simulation results for each model type are organized in the [`simulation`](https://github.com/anonstats123/Nullstrap/tree/main/comparison_examples) directory.  
These include the performance of **Nullstrap** as well as the existing methods compared in our manuscript (e.g., knockoff, data splitting).

- Linear model: [`simulation/linear_model/`](https://github.com/anonstats123/Nullstrap/tree/main/comparison_examples/LM)
- Generalized linear model (GLM): [`simulation/glm_model/`](https://github.com/anonstats123/Nullstrap/tree/main/comparison_examples/GLM)
- Cox model: [`simulation/cox_model/`](https://github.com/anonstats123/Nullstrap/tree/main/comparison_examples/Cox)
- Graphical model: [`simulation/graph_model/`](https://github.com/anonstats123/Nullstrap/tree/main/comparison_examples/Graphical)

Each folder contains code to generate synthetic data, apply all methods, and evaluate selection results (e.g., FDR, power).


## Real Data Example<a name="real-data-examples"></a>
A real data analysis is conducted using a longitudinal time-to-labor dataset collected from pregnant women receiving antepartum and postpartum care at Stanford’s Lucile Packard Children’s Hospital ([Stelzer et al., 2021](https://www.science.org/doi/full/10.1126/scitranslmed.abd9898)). The dataset includes 63 participants in their second or third trimester of an uncomplicated pregnancy with a single fetus, each contributing 1 to 3 samples. Each sample comprises 6348 variables, including 3529 metabolites, 1317 plasma proteins, and 1502 single-cell immune variables derived from blood mass cytometry.

The dataset is available at [`real_data_example/Onset of Labor`](https://github.com/anonstats123/Nullstrap/tree/main/real_data_example/Onset%20of%20Labor), and the corresponding analysis results can be found in [`real_data_example`](https://github.com/anonstats123/Nullstrap/tree/main/real_data_example).

