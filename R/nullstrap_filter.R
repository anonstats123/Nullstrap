#' The Nullstrap procedure
#'
#' @description
#' This function runs the full Nullstrap procedure to identify significant variables under a user-specified FDR threshold.
#'
#' @param X n-by-p matrix or data frame of predictors.
#' @param y n-vector of response variables. For cox models, there should be two columns:
#' the first column is the survival time and the second column is the censoring indicator. For graphical models,
#' there is no y.
#' @param fdr_value Numeric value FDR levels, must be between 0 and 1.
#' @param best_lambda Regularization parameter (lambda) used in the model.
#' If not specified by the user (default: \code{NULL}), it will be automatically
#' determined based on the chosen \code{model_type} using cross-validation.
#' For linear models, the default is \code{0.5*fit$lambda.min}.
#' When \eqn{n} is small and \eqn{p > n}, it is recommended to use \code{fit$lambda.min}.
#' For generalized linear models, the default is \code{0.5*fit$lambda.1se}. For cox models, the default is
#' \code{fit$lambda.min}. For graphical models, the default is the one corresponding to the maximum AIC from cross-validation
#' @param B_reps Number of repetitions for computing the correction factor. The correction factor is calculated B_reps times using different random seeds,
#' and the maximum value across all repetitions is used. If not specified by the user (default: \code{NULL}),
#' it will be automatically determined based on the chosen \code{model_type}. For linear models,
#' the default is 5. For generalized linear models, the default is 2. For cox models, the default is 2.
#' For graphical models, the default is 2.
#' @param dist_type Type of distribution for the error term. Applicable only to linear models.
#' The default is \code{"normal"}.
#' @param model_type Type of model, must be one of "linear", "glm", "cox", "graphical".
#'
#' @return Returns a list containing the following elements:
#' \describe{
#'   \item{threshold}{computed selection threshold}
#'   \item{selected}{named vector of selected variables}
#'   \item{statistic}{computed test statistics}
#' }
#'
#' @examples
#' n <- 2000
#' p <- 1000
#' s <- 30
#' Amp <- 0.25
#' rho <- 0.8
#' Theta.8 <- toeplitz(rho^(0:(p - 1)))
#' X <- mvrnorm(n, rep(0, p), Sigma = Theta.8)
#' X <- scale(X)
#' beta <- rep(0, p)
#' beta[1:s] <- sample(c(-Amp, Amp), s, replace = TRUE)
#' true_labels <- beta != 0
#' y <- X %*% beta + rnorm(n)
#' y <- y - mean(y)
#' result_normal <- nullstrap_filter(
#'   X, y, fdr_value = 0.1, best_lambda = NULL,
#'   dist_type = "normal", model_type = "linear"
#' )
#' y <- X %*% beta + rt(n, df = 3)
#' y <- y - mean(y)
#' result_t <- nullstrap_filter(
#'   X, y, fdr_value = 0.1, best_lambda = NULL, B_reps = NULL,
#'   dist_type = "t", model_type = "linear"
#' )
#' result_normal$selected
#' result_t$selected
#'
#' @import glmnet
#' @import MASS
#' @import parallel
#' @import PRROC
#' @import survival
#' @import eha
#' @import simsurv
#' @import huge
#' @import glasso
#'
#' @export nullstrap_filter


nullstrap_filter <- function(X, y, fdr_value, best_lambda = NULL, B_reps = NULL,
                     dist_type = "normal", model_type) {

  library(glmnet)
  library(MASS)
  library(parallel)
  library(PRROC)
  library(survival)
  library(eha)
  library(simsurv)
  library(huge)
  library(glasso)


  if (model_type == "linear") {

    binary_search = function(coef_real, coef_knockoff, q) {
      # coef_knockoff is tilde_beta + gamma = tilde_beta'
      # coef_real is beta_hat from real data
      left = 0
      right = max(coef_real)
      while (abs(left - right) > 1e-8) {
        mid = ((left + right) / 2)
        num_negative = length(which(coef_knockoff >= mid)) + q/2
        num_positive = max(length(which(coef_real >= mid)), 1)
        FDP = num_negative / num_positive
        if (FDP > q) {
          left = mid
        } else {
          right = mid
        }
      }
      return(right)
    }

    # beta_real = coef_corr
    # hat_real = hat_sigma
    correction_factor_estimation = function(X, y, beta_real, hat_real, residuals, best_lambda, q){
      # beta_real: coef of the positive control = beta daggar = beta hat
      n = nrow(X)
      beta_real[which(beta_real>0)] = beta_real[which(beta_real>0)] + best_lambda
      beta_real[which(beta_real<0)] = beta_real[which(beta_real<0)] - best_lambda
      y_correction = X %*% beta_real + sample(residuals, replace = TRUE)
      X = scale(X)
      y_correction = y_correction - mean(y_correction)
      Signal_index_correction = which(beta_real != 0)

      #### correction true ######
      model_correction = glmnet(X, y_correction, alpha = 1, intercept = F, lambda = best_lambda)
      coef_correction = coef(model_correction)
      coef_correction = coef_correction[-1] # remove the intercept
      # hat_sigma_corr = sqrt(sum((y - X %*% as.vector(coef_correction))^2) / max(1,(n - sum(coef_correction != 0))))
      coef_correction = abs(as.vector(coef_correction))
      # coef_correction = coef_correction[-which.max(coef_correction)]
      X_snp = X
      X_snp = scale(X_snp)
      y_snp = sample(residuals, replace = TRUE)
      y_snp = y_snp - mean(y_snp)
      model_snp = glmnet(X_snp, y_snp, alpha = 1, intercept = F, lambda = best_lambda)
      coef_snp = coef(model_snp)
      coef_snp = coef_snp[-1] # remove the intercept
      coef_snp = abs(as.vector(coef_snp))
      # coef_snp = coef_snp[-which.max(coef_snp)]
      left_correction = 0.05
      right_correction = max(coef_correction)/(best_lambda + sqrt(log(p))/sqrt(n))
      while (abs(right_correction - left_correction) > 1e-8) {
        mid = ((right_correction + left_correction) / 2)
        FDR_vector = c()
        # for (b in 1:B){

        # the “inflated” null coefficients after applying the candidate correction.
        corrected = abs(coef_snp) + mid * hat_real * (best_lambda + sqrt(log(p))/sqrt(n))
        # threshold computed for this candidate correction factor
        tau_mid = binary_search(coef_correction, corrected, q)
        # where the “corrected” real coefficients exceed the threshold tau_mid
        rejected_mid = which(coef_correction >= tau_mid)
        # The empirical false discovery proportion at threshold tau_mid
        FDR_mid = length(which(coef_correction[setdiff(1:p, Signal_index_correction)] >= tau_mid)) /  max(length(rejected_mid), 1)
        FDR_vector = c(FDR_vector, FDR_mid)
        # }
        FDR_mean = mean(FDR_vector)
        if (FDR_mean > q) { #
          left_correction = mid
        } else {
          right_correction = mid
        }
      }
      return(right_correction)
    }

    correction_factor_estimation_normal = function(X, y, beta_real, hat_real, best_lambda, q){
      # beta_real: coef of the positive control = beta daggar = beta hat
      n = nrow(X)
      beta_real[which(beta_real>0)] = beta_real[which(beta_real>0)] + best_lambda
      beta_real[which(beta_real<0)] = beta_real[which(beta_real<0)] - best_lambda
      y_correction = X %*% beta_real + rnorm(n, sd = hat_real)
      X = scale(X)
      y_correction = y_correction - mean(y_correction)
      Signal_index_correction = which(beta_real != 0)

      #### correction true ######
      model_correction = glmnet(X, y_correction, alpha = 1, intercept = F, lambda = best_lambda)
      coef_correction = coef(model_correction)
      coef_correction = coef_correction[-1] # remove the intercept
      # hat_sigma_corr = sqrt(sum((y - X %*% as.vector(coef_correction))^2) / max(1,(n - sum(coef_correction != 0))))
      coef_correction = abs(as.vector(coef_correction))
      # coef_correction = coef_correction[-which.max(coef_correction)]
      X_snp = X
      X_snp = scale(X_snp)
      y_snp = rnorm(n, sd = hat_real)
      y_snp = y_snp - mean(y_snp)
      model_snp = glmnet(X_snp, y_snp, alpha = 1, intercept = F, lambda = best_lambda)
      coef_snp = coef(model_snp)
      coef_snp = coef_snp[-1] # remove the intercept
      coef_snp = abs(as.vector(coef_snp))
      # coef_snp = coef_snp[-which.max(coef_snp)]
      left_correction = 0.05
      right_correction = max(coef_correction)/(best_lambda + sqrt(log(p))/sqrt(n))
      while (abs(right_correction - left_correction) > 1e-8) {
        mid = ((right_correction + left_correction) / 2)
        FDR_vector = c()
        # for (b in 1:B){

        # the “inflated” null coefficients after applying the candidate correction.
        corrected = abs(coef_snp) + mid * hat_real * (best_lambda + sqrt(log(p))/sqrt(n))
        # threshold computed for this candidate correction factor
        tau_mid = binary_search(coef_correction, corrected, q)
        # where the “corrected” real coefficients exceed the threshold tau_mid
        rejected_mid = which(coef_correction >= tau_mid)
        # The empirical false discovery proportion at threshold tau_mid
        FDR_mid = length(which(coef_correction[setdiff(1:p, Signal_index_correction)] >= tau_mid)) /  max(length(rejected_mid), 1)
        FDR_vector = c(FDR_vector, FDR_mid)
        # }
        FDR_mean = mean(FDR_vector)
        if (FDR_mean > q) { #
          left_correction = mid
        } else {
          right_correction = mid
        }
      }
      return(right_correction)
    }
    n = nrow(X)
    p = ncol(X)
    X <- scale(X)
    y <- y - mean(y)

    fit <- cv.glmnet(X, y, intercept = FALSE, alpha = 1)
    if (is.null(best_lambda)) {
      best_lambda <- 0.5 * fit$lambda.min
    }
    best_lasso_model <- glmnet(X, y, alpha = 1, intercept = FALSE, lambda = best_lambda)
    coef_real <- coef(best_lasso_model)[-1]
    coef_corr <- coef_real
    hat_sigma = sqrt(sum((y - X %*% as.vector(coef_real))^2) / max(1,(n - sum(coef_real != 0))))
    coef_real = abs(as.vector(coef_real))
    residuals = min(n^(1/4), sqrt(n/(n - sum(coef_real != 0))))*(y - X %*% as.vector(coef(best_lasso_model)[-1]))

    X_knockoff <- X
    X_knockoff <- scale(X_knockoff)
    beta_knockoff <- rep(0, p)

    if (dist_type == "normal" || dist_type == "Normal") {
      y_knockoff = X_knockoff %*% beta_knockoff + rnorm(n, sd = hat_sigma)
      y_knockoff=y_knockoff - mean(y_knockoff)

      ### Lasso: Fit model to synthetic data (X, tilde_y) ###
      knockoff_lasso_model <- glmnet(X_knockoff, y_knockoff, intercept = F,
                                     alpha = 1, lambda = best_lambda)
      coef_knockoff = coef(knockoff_lasso_model)
      coef_knockoff = coef_knockoff[-1] # remove the intercept

      ### Find the data-driven correction factor gamma ###
      nu_vec = c()
      B <- B_reps
      if (is.null(B)) {
        B <- 5
      }
      for (b in 1:B) {
        #print(b)
        set.seed(b)
        nu = correction_factor_estimation_normal(X, y, coef_corr, hat_sigma, best_lambda, fdr_value)
        nu_vec = c(nu_vec, nu)
      }
      fdr_control = max(nu_vec)
    } else {
      y_knockoff = X_knockoff %*% beta_knockoff + sample(residuals, replace = TRUE)
      y_knockoff=y_knockoff - mean(y_knockoff)

      ### Lasso: Fit model to synthetic data (X, tilde_y) ###
      knockoff_lasso_model <- glmnet(X_knockoff, y_knockoff, intercept = F,
                                     alpha = 1, lambda = best_lambda)
      coef_knockoff = coef(knockoff_lasso_model)
      coef_knockoff = coef_knockoff[-1] # remove the intercept

      ### Find the data-driven correction factor gamma ###
      nu_vec = c()

      B <- B_reps
      if (is.null(B)) {
        B <- 5
      }
      for (b in 1:B) {
        #print(b)
        set.seed(b)
        nu = correction_factor_estimation(X, y, coef_corr, hat_sigma, residuals, best_lambda, fdr_value)
        nu_vec = c(nu_vec, nu)
      }
      fdr_control = max(nu_vec)
    }
    coef_knockoff <- abs(coef_knockoff) + fdr_control * hat_sigma * (best_lambda + sqrt(log(p)) / sqrt(n))
    threshold = binary_search(coef_real, coef_knockoff, fdr_value)
    reject = which(coef_real >= threshold)
    model_statistic = coef_real

  } else if (model_type == "glm") {

    binary_search = function(coef_real, coef_knockoff, q) {
      # coef_knockoff is tilde_beta + gamma = tilde_beta'
      # coef_real is beta_hat from real data
      left = 0
      right = max(coef_real)
      while (abs(left - right) > 1e-8) {
        mid = ((left + right) / 2)
        num_negative = length(which(coef_knockoff >= mid)) + q/2
        num_positive = max(length(which(coef_real >= mid)), 1)
        FDP = num_negative / num_positive
        if (FDP > q) {
          left = mid
        } else {
          right = mid
        }
      }
      return(right)
    }
    # Generate design matrix X (i.i.d. rows with AR(1) columns)
    correction_factor_estimation_GLM = function(X, y, beta_real, best_lambda, q){
      # beta_real: coef of the positive control = beta daggar = beta hat
      n = nrow(X)
      beta_real[which(beta_real>0)] = beta_real[which(beta_real>0)] + best_lambda
      beta_real[which(beta_real<0)] = beta_real[which(beta_real<0)] - best_lambda
      # Generate binary response y (from binomial model with logit link function)
      linear_predictor <- X %*% beta_real
      prob <- as.vector(1 / (1 + exp(-linear_predictor)))  # logistic function
      y_correction <- rbinom(n, 1, prob)  # binary response
      Signal_index_correction = which(beta_real != 0)
      # max_coef = max(abs(beta_real[setdiff(1:p, Signal_index_correction)]))
      #### correction true ######
      # model_correction = glmnet(X, y_correction, alpha = 1, intercept = F, lambda = best_lambda)
      model_correction <- glmnet(X, y_correction, alpha = 1, lambda = best_lambda, thresh = 1e-10, family = "binomial")

      coef_correction = coef(model_correction)
      coef_correction = coef_correction[-1] # remove the intercept
      # hat_sigma_corr = sqrt(sum((y - X %*% as.vector(coef_correction))^2) / max(1,(n - sum(coef_correction != 0))))
      coef_correction = abs(as.vector(coef_correction))
      # coef_correction = coef_correction[-which.max(coef_correction)]
      X_snp = X
      # X_snp = scale(X_snp)
      beta_0 <- rep(0, p)
      # X_snp <- matrix(rnorm(n * p, mean = 0, sd = 1 / sqrt(n) ), nrow = n, ncol = p)
      linear_predictor <- X_snp %*% beta_0
      # Generate binary response y (from binomial model with logit link function)
      prob <- as.vector(1 / (1 + exp(-linear_predictor)))  # logistic function
      y_snp <- rbinom(n, 1, prob)  # binary response
      # model_snp = glmnet(X_snp, y_snp, alpha = 1, intercept = F, lambda = best_lambda)
      model_snp <- glmnet(X_snp, y_snp, alpha = 1, lambda = best_lambda, thresh = 1e-10, family = "binomial")
      coef_snp = coef(model_snp)
      coef_snp = coef_snp[-1] # remove the intercept
      coef_snp = abs(as.vector(coef_snp))
      # coef_snp = coef_snp[-which.max(coef_snp)]
      left_correction = 0.00000000000001
      right_correction = max(coef_correction)
      # / (sqrt(log(p)/n))
      # tau = 0
      while (abs(right_correction - left_correction) > 1e-10) {
        # mid = 1
        mid = ((right_correction + left_correction) / 2)
        FDR_vector = c()
        # for (b in 1:B){

        # the “inflated” null coefficients after applying the candidate correction.
        corrected = abs(coef_snp) + mid
        # * (sqrt(log(p)/n))
        # threshold computed for this candidate correction factor
        tau_mid = binary_search(coef_correction, corrected, q)
        tau = tau_mid
        # where the “corrected” real coefficients exceed the threshold tau_mid
        rejected_mid = which(coef_correction >= tau_mid)
        # The empirical false discovery proportion at threshold tau_mid
        FDR_mid = length(which(coef_correction[setdiff(1:p, Signal_index_correction)] >= tau_mid)) /  max(length(rejected_mid), 1)
        FDR_vector = c(FDR_vector, FDR_mid)
        # }
        FDR_mean = mean(FDR_vector)
        if (FDR_mean > q) { #
          left_correction = mid
        } else {
          right_correction = mid
        }
      }
      return(right_correction)
    }

    n = nrow(X)
    p = ncol(X)
    X = scale(X)/(sqrt(n))
    # Generate beta (true coefficients)

    beta_0 <- rep(0, p)
    # X_snp <- matrix(rnorm(n * p, mean = 0, sd = 1 / sqrt(n) ), nrow = n, ncol = p)
    linear_predictor <- X %*% beta_0
    # Generate binary response y (from binomial model with logit link function)
    prob <- 1 / (1 + exp(-linear_predictor))  # logistic function
    y_snp <- rbinom(n, 1, prob)  # binary response
    cv_lasso <- cv.glmnet(X, y, alpha = 1, family = "binomial", measure = "mse")


    if (is.null(best_lambda)) {
      best_lambda <- 0.5 * cv_lasso$lambda.1se
    }
    best_lasso_model <- glmnet(X, y, alpha = 1, lambda = best_lambda,
                               thresh = 1e-10, family = "binomial")
    coef_beta = abs(as.vector(best_lasso_model$beta))
    coef_corr = as.vector(best_lasso_model$beta)

    snp_lasso_model <- glmnet(X, y_snp, alpha = 1, lambda = best_lambda,
                              thresh = 1e-10, family = "binomial")
    nu_vec = c()
    B <- B_reps
    if (is.null(B)) {
      B <- 2
    }
    for (b in 1:B) {
      set.seed(b)
      nu = correction_factor_estimation_GLM(X, y, coef_corr, best_lambda, fdr_value)
      nu_vec = c(nu_vec, nu)
    }
    fdr_control = max(nu_vec)
    coef_snp = abs(as.vector(snp_lasso_model$beta)) + fdr_control

    threshold = binary_search(coef_beta, coef_snp, fdr_value)
    reject = which(coef_beta >= threshold)
    model_statistic = coef_beta

  } else if (model_type == "cox") {

    correction_factor_estimation_Cox = function(X, y, beta_real, baseline_fun, best_lambda, min_time, max_time, q){
      # beta_real: coef of the positive control = beta daggar = beta hat
      n = nrow(X)
      p = ncol(X)
      beta_real[which(beta_real>0)] = beta_real[which(beta_real>0)] +  best_lambda + sqrt(log(p)/n)
      beta_real[which(beta_real<0)] = beta_real[which(beta_real<0)] - best_lambda - sqrt(log(p)/n)

      betas = rep(0, p)
      names(betas) = colnames(X)
      names(beta_real) = colnames(X)
      survival_time_correction <- simsurv(cumhazard = baseline_fun, x = as.data.frame(X), betas = beta_real)
      time_simu = simsurv(cumhazard = baseline_fun, x = as.data.frame(X), betas = betas)
      fit_wei <- fitdistr(time_simu$eventtime, "weibull")
      survival_time_correction <- simsurv(lambdas = fit_wei$estimate[2], gammas = fit_wei$estimate[1], x = as.data.frame(X), betas = beta_real)
      status_correction = rep(1, n)
      time_correction = survival_time_correction$eventtime
      y_correction = cbind(time = time_correction, status = status_correction)
      # y_snp = Surv(time = y_snp[, 1], event = y_snp[, 2])
      fit_correction <- glmnet(X, y_correction, family = "cox",  lambda = best_lambda)
      coef_correction = as.vector(abs(coef(fit_correction)))
      Signal_index_correction = which(beta_real != 0)
      # max_coef = max(abs(beta_real[setdiff(1:p, Signal_index_correction)]))
      #### correction true ######

      snp_survival_time <- simsurv(cumhazard = baseline_fun, x = as.data.frame(X), betas = betas,  interval = c(min_time, max_time))
      # library(MASS)
      # time_weibull = time[status == 0]
      # fit_weibull <- fitdistr(time_weibull, "weibull")
      # View the result
      # sample weibull distribution with parameters from fit_weibull
      # censor_time_weibull_sample = rweibull(n, shape = fit_weibull$estimate[1], scale = fit_weibull$estimate[2])
      # censor_time = rep(max(snp_survival_time$eventtime)+1, n)
      status_snp = rep(1, n)
      time_snp =snp_survival_time$eventtime
      y_snp = cbind(time = time_snp, status = status_snp)
      # y_snp = Surv(time = y_snp[, 1], event = y_snp[, 2])
      fit_snp <- glmnet(X, y_snp, family = "cox",  lambda = best_lambda)
      coef_snp = abs(fit_snp$beta)
      left_correction = 0.000000000001
      right_correction = max(coef_correction) / (best_lambda + sqrt(log(p)/n))
      tau = 0
      while (abs(right_correction - left_correction) > 1e-10) {
        mid = 1
        mid = ((right_correction + left_correction) / 2)
        FDR_vector = c()
        # for (b in 1:B){

        # the “inflated” null coefficients after applying the candidate correction.
        corrected = abs(coef_snp) + mid  * (best_lambda + sqrt(log(p)/n))
        # threshold computed for this candidate correction factor
        tau_mid = binary_search(coef_correction, corrected, q)
        # tau = tau_mid
        # where the “corrected” real coefficients exceed the threshold tau_mid
        rejected_mid = which(coef_correction >= tau_mid)
        # The empirical false discovery proportion at threshold tau_mid
        FDR_mid = length(which(coef_correction[setdiff(1:p, Signal_index_correction)] >= tau_mid)) /  max(length(rejected_mid), 1)
        FDR_vector = c(FDR_vector, FDR_mid)
        # }
        FDR_mean = mean(FDR_vector)
        if (FDR_mean > q) { #
          left_correction = mid
        } else {
          right_correction = mid
        }
      }
      return(right_correction)
    }

    binary_search = function(coef_real, coef_knockoff, q) {
      left = 0
      right = max(coef_real)
      step = 100
      while ((abs(left - right) > 1e-8) & (step > 0)) {
        mid = ((left + right) / 2)
        num_negative = length(which(as.vector(coef_knockoff) >= mid)) + q/2
        num_positive = max(length(which(coef_real >= mid)), 1)
        FDP = num_negative / num_positive
        if (FDP > q) {
          left = mid
        } else {
          right = mid
        }
        step = step - 1
      }
      return(right)
    }

    n = nrow(X)
    p = ncol(X)
    X = scale(X)/(sqrt(n))
    colnames(X) = paste0("X", 1:p)


    cv_lasso <- cv.glmnet(X, y, alpha = 1, family = "cox", measure = "C")
    if (is.null(best_lambda)) {
      best_lambda <- cv_lasso$lambda.min
    }

    fit_cox <- glmnet(X, y, alpha = 1, family = "cox", lambda = best_lambda)
    coef_cox = abs(as.vector(fit_cox$beta))
    coef_corr = as.vector(fit_cox$beta)


    # Step 2: Fit Cox model using glmnet
    y_baseline = Surv(time = y[, 1], event = y[, 2])

    fit_baseline <- coxph(y_baseline ~ X[,which(coef_corr!=0)])
    cum_hazard = basehaz(fit_baseline) # cumulative baseline hazard rate
    cum_hazard_fun = approxfun(cum_hazard$time, cum_hazard$hazard)
    min_time = min(cum_hazard$time)
    max_time = max(cum_hazard$time)

    cum_hazard_fun_extension = function(t) {
      if (t <= min_time) {
        return(0)
      } else if (t >= max_time) {
        return(cum_hazard_fun(max_time-1e-6))
      } else {
        return(cum_hazard_fun(t))
      }
    }
    baseline_fun = function(t, x, betas, ...) {
      return(cum_hazard_fun_extension(t))
    }
    betas = rep(0, p)
    names(betas) = colnames(X)
    snp_survival_time <- simsurv(cumhazard = baseline_fun,
                                 x = as.data.frame(X), betas = betas,
                                 interval = c(min_time, max_time))

    status_snp = rep(1, n)
    time_snp =snp_survival_time$eventtime
    y_snp = cbind(time = time_snp, status = status_snp)
    # y_snp = Surv(time = y_snp[, 1], event = y_snp[, 2])
    fit_snp <- glmnet(X, y_snp, family = "cox",  lambda = best_lambda)
    nu_vec = c()
    B <- B_reps
    if (is.null(B)) {
      B <- 2
    }
    for (b in 1:B) {
      #print(b)
      set.seed(b)
      # correction_factor_estimation_Cox = function(X, y, beta_real, baseline_fun, best_lambda, q)
      nu = correction_factor_estimation_Cox(X, y, coef_corr, baseline_fun,
                                            best_lambda, min_time, max_time,
                                            fdr_value)
      nu_vec = c(nu_vec, nu)
    }
    fdr_control = max(nu_vec)
    coef_snp = abs(as.vector(fit_snp$beta)) + fdr_control * (best_lambda + sqrt(log(p)/n))

    threshold = binary_search(coef_cox, coef_snp, fdr_value)
    reject = which(coef_cox >= threshold)
    model_statistic = coef_cox

  } else if (model_type == "graphical") {

    generate_knockoff_data = function(n, p, diag_vec){
      diag_mat = diag_vec
      data <- mvrnorm(n = n, mu = rep(0, p), Sigma = diag_mat)
      return(data)
    }

    binary_search = function(coef_real, coef_knockoff, q) {
      # coef_knockoff is tilde_beta + gamma = tilde_beta'
      # coef_real is beta_hat from real data
      left = 0
      right = max(coef_real)
      while (abs(left - right) > 1e-8) {
        mid = ((left + right) / 2)
        num_negative = length(which(coef_knockoff >= mid)) + q/2
        num_positive = max(length(which(coef_real >= mid)), 1)
        FDP = num_negative / num_positive
        if (FDP > q) {
          left = mid
        } else {
          right = mid
        }
      }
      return(right)
    }

    correction_factor_estimation_GGM = function(n, p, Omega_dagger, scale, best_lambda, q){
      # beta_real: coef of the positive control = beta daggar = beta hat
      # Omega_dagger[which(Omega_dagger>0)] = Omega_dagger[which(Omega_dagger>0)] + best_lambda
      # Omega_dagger[which(Omega_dagger<0)] = Omega_dagger[which(Omega_dagger<0)] - best_lambda
      ####### Generate true data ######
      data_correction = generate_data(n, p, Omega_dagger)
      true_support_correction = obtain_true_support(Omega_dagger)
      Signal_index_correction = which(true_support_correction[upper.tri(true_support_correction)] == 1)
      data_scale_correction = scale(data_correction)
      glasso_result_correction <- glasso(cov(data_scale_correction), rho = best_lambda)
      omega_hat_correction = glasso_result_correction$wi
      knockoff_data = generate_knockoff_data(n, p, diag_vec = diag(rep(1,p)))
      knockoff_data = scale(knockoff_data)
      glasso_knockoff_result <- glasso(cov(knockoff_data), rho = best_lambda)
      omega_knockoff_hat = glasso_knockoff_result$wi
      upper_tri_correction = abs(omega_hat_correction[upper.tri(omega_hat_correction)])
      upper_tri_knockoff_correction = abs(omega_knockoff_hat[upper.tri(omega_knockoff_hat)])
      left_correction = 0.00000001
      right_correction = max(upper_tri_correction) / scale
      tau = 0
      while (abs(right_correction - left_correction) > 1e-10) {
        # mid = 0.001
        mid = 0.000000001
        # mid = 0.0000001
        mid = ((right_correction + left_correction) / 2)
        FDR_vector = c()
        # for (b in 1:B){

        # the “inflated” null coefficients after applying the candidate correction.
        upper_tri_knockoff_correction_corrected = abs(upper_tri_knockoff_correction) + mid * scale
        # threshold computed for this candidate correction factor
        tau_mid = binary_search(upper_tri_correction, upper_tri_knockoff_correction_corrected, q)
        # tau = tau_mid
        # where the “corrected” real coefficients exceed the threshold tau_mid
        rejected_mid = which(upper_tri_correction >= tau_mid)
        # The empirical false discovery proportion at threshold tau_mid
        FDR_mid = length(which(upper_tri_correction[setdiff(1:length(upper_tri_correction), Signal_index_correction)] >= tau_mid)) /  max(length(rejected_mid), 1)
        FDR_vector = c(FDR_vector, FDR_mid)
        # print(mid)
        # }
        FDR_mean = mean(FDR_vector)
        # print(FDR_mid)
        if (FDR_mid > q) { #
          left_correction = mid
        } else {
          right_correction = mid
        }
      }
      return(right_correction)
    }
    n = nrow(X)
    p = ncol(X)
    data_scale = scale(X)

    cv_result <- huge(data_scale, method = "glasso", nlambda = 100,
                      lambda.min.ratio = 0.05)

    # calculate the AIC
    cv_result$AIC = n/2*2*cv_result$loglik - 2*cv_result$df
    best_lambda_index <- which.max(cv_result$AIC)
    if (is.null(best_lambda)) {
      best_lambda <- cv_result$lambda[best_lambda_index]
    }
    glasso_result <- glasso(cov(data_scale), rho = best_lambda)
    omega_hat = glasso_result$wi

    best_lambda_dagger_index = which.min(abs(cv_result$sparsity - 0.1))
    best_lambda_dagger = cv_result$lambda[best_lambda_dagger_index]
    glasso_result_dagger <- glasso(cov(data_scale), rho = best_lambda_dagger)
    omega_hat_dagger = glasso_result_dagger$wi

    omega_hat_dagger = (omega_hat_dagger + t(omega_hat_dagger))/2
    omega_hat_dagger[which(omega_hat_dagger>0)] = omega_hat_dagger[which(omega_hat_dagger>0)] + best_lambda_dagger
    omega_hat_dagger[which(omega_hat_dagger<0)] = omega_hat_dagger[which(omega_hat_dagger<0)] - best_lambda_dagger
    # omega_hat_dagger = (omega_hat_dagger + t(omega_hat_dagger))/2
    eigen_d = eigen(omega_hat_dagger)$values
    if (min(eigen_d) <= 0) {
      omega_hat_dagger = omega_hat_dagger + diag(rep(1,p))*(0.01 + abs(min(eigen_d)))
    }


    omega_hat_dagger = diag(diag(cov(X))^(-1/2)) %*% omega_hat_dagger %*% diag(diag(cov(X))^(-1/2))
    knockoff_data = generate_knockoff_data(n, p, diag_vec = diag(rep(1,p)))

    knockoff_data = scale(knockoff_data)
    glasso_knockoff_result <- glasso(cov(knockoff_data), rho = best_lambda)
    omega_knockoff_hat = glasso_knockoff_result$wi
    K_gamma1 = cov(data_scale)
    inf_norm = 1/((norm(K_gamma1, type = "1")/p))

    scale = inf_norm * (best_lambda + sqrt(log(p)/n))

    nu_vec = c()
    B <- B_reps
    if (is.null(B)) {
      B <- 2
    }
    for (b in 1:B) {
      set.seed(b)
      nu = correction_factor_estimation_GGM(n, p, omega_hat_dagger, scale, best_lambda, fdr_value)
      nu_vec = c(nu_vec, nu)
    }
    fdr_control = max(nu_vec)


    upper_tri_real_vec <- abs(omega_hat[upper.tri(omega_hat)])
    upper_tri_knockoff_vec <- abs(omega_knockoff_hat[upper.tri(omega_knockoff_hat)])

    upper_tri_knockoff_vec = upper_tri_knockoff_vec +  fdr_control * scale

    threshold = binary_search(upper_tri_real_vec,
                              upper_tri_knockoff_vec, fdr_value)
    reject = which(upper_tri_real_vec >= threshold)
    model_statistic = upper_tri_real_vec
  } else {
    stop("Unknown Model")
  }

  result_list <- list(
    threshold = threshold,
    selected = reject,
    statistic = model_statistic
  )

  return(result_list)
}
