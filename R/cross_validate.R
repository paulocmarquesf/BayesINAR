#' Time series cross-validation for generalized Bayesian INAR models
#' @description Gives predictions for a test set and the mean absolute errors.
#' @param time_series A univariate time series.
#' @param model Model to be used: 'inar', 'adinar', or 'dpinar'
#' @param h Number of steps ahead to be predicted. 
#' @param training_epoch The last observation of the first training set.
#' @param prior List of prior hyperparameters:
#' \describe{
#'   \item{a_alpha}{Hyperparameters of the thinning component.}
#'   \item{a_theta, b_theta}{Hyperparameters of the Geometric innovation parameter.}
#'   \item{a_lambda, b_lambda}{Hyperparameters of the Poisson innovation rate.}
#'   \item{a_w, b_w}{Hyperparameters of the Geometric-Poisson mixture weight.}
#'   \item{a_tau, b_tau}{Hyperparameters of the concentration parameter Gamma prior.}
#'   \item{a0, b0}{Base measure hyperparameters.}
#'   \item{lambda_max}{Hyperparameter of the uniform distribution that minimizes the corresponding D-KL.}
#' }
#' @param burn_in Length of the "burn-in" period.
#' @param chain_length Length of the Markov chain.
#' @param random_seed Value of the random seed.
#' @param verbose  Running status if \code{TRUE}.
#' @return A list with the following elements:
#' \describe{
#'   \item{est}{Predictions for the test set.}
#'   \item{mae}{Mean Absolute Error of the test set predictions.}
#' }
#' @export
cross_validate <- function(time_series,
                           model = 'inar',
                           h = 1,
                           training_epoch = NULL,
                           prior = list(a_alpha = NULL,
                                        a0 = NULL, b0 = NULL,
                                        tau = NULL, k0 = NULL,
                                        a_tau = NULL, b_tau = NULL,
                                        a_w = NULL, b_w = NULL,
                                        a_theta = NULL, b_theta = NULL,
                                        sigma = NULL, 
                                        lambda_max = NULL),
                           burn_in = 10^3,
                           chain_length = 10^4,
                           random_seed = 1761,
                           verbose = TRUE)
{
    if(is.null(training_epoch)) training_epoch <- round(0.7 * length(time_series))
    
    y_hat <- matrix(NA, nrow = length(h), ncol = length(time_series) - (training_epoch + min(h)) + 1)
    diff_y_hat_obs <- matrix(NA, nrow = length(h), ncol = length(time_series) - (training_epoch + min(h)) + 1)
    p <- 1
    if (! model %in% c("inar", "adinar", "dpinar")) stop("Model parameter must be equal to one of the values: ['inar', 'adinar', 'dpinar'].") 
    
    if (model == 'inar') {
        prior <- list(a_alpha = prior[["a_alpha"]],
                      a_lambda = prior[["a_lambda"]], 
                      b_lambda = prior[["b_lambda"]])
        for (j in 1:ncol(y_hat)) {
            if (verbose) cat(sprintf("Training up to epoch %d ...\n", training_epoch + j - 1))
            train_model <- inar(time_series[1:(training_epoch + j - 1)],
                                prior = prior,
                                burn_in = burn_in,
                                chain_length = chain_length,
                                random_seed = random_seed,
                                verbose = verbose)
            for (i in 1:nrow(y_hat)) {
                if (training_epoch + j + i - 1  <= length(time_series)) {
                    if (verbose) cat(sprintf("Predicting epoch %d (h = %d)...\n", training_epoch + j + i - 1, h[i]))
                    y_hat[i, j + i - 1] <- predict(train_model, h = h[i])$est
                    diff_y_hat_obs[i, j + i - 1] = abs(time_series[training_epoch + j + i - 1] - y_hat[i, j + i - 1])
                }
            }
        }
    }
    
    if (model == 'adinar') {
        prior <- list(a_alpha = prior[["a_alpha"]],
                      a_lambda = prior[["a_lambda"]], 
                      b_lambda = prior[["b_lambda"]], 
                      a_w = prior[["a_w"]], 
                      b_w = prior[["b_w"]])
        for (j in 1:ncol(y_hat)) {
            if (verbose) cat(sprintf("Training up to epoch %d ...\n", training_epoch + j - 1))
            train_model <- adinar(time_series[1:(training_epoch + j - 1)],
                                  prior = prior,
                                  burn_in = burn_in,
                                  chain_length = chain_length,
                                  random_seed = random_seed,
                                  verbose = verbose)
            for (i in 1:nrow(y_hat)) {
                if (training_epoch + j + i - 1  <= length(time_series)) {
                    if (verbose) cat(sprintf("Predicting epoch %d (h = %d)...\n", training_epoch + j + i - 1, h[i]))
                    y_hat[i, j + i - 1] <- predict(train_model, h = h[i])$est
                    diff_y_hat_obs[i, j + i - 1] = abs(time_series[training_epoch + j + i - 1] - y_hat[i, j + i - 1])
                }
            }
        }
    }
    
    if (model == 'dpinar') {
        prior <- list(a_alpha = prior[["a_alpha"]],
                      a_tau = prior[["a_tau"]], 
                      b_tau = prior[["b_tau"]], 
                      a0 = prior[["a0"]], 
                      b0 = prior[["b0"]],
                      lambda_max = prior[["lambda_max"]])
        for (j in 1:ncol(y_hat)) {
            if (verbose) cat(sprintf("Training up to epoch %d ...\n", training_epoch + j - 1))
            train_model <- dpinar(time_series[1:(training_epoch + j - 1)],
                                  prior = prior,
                                  burn_in = burn_in,
                                  chain_length = chain_length,
                                  random_seed = random_seed,
                                  verbose = verbose)
            for (i in 1:nrow(y_hat)) {
                if (training_epoch + j + i - 1  <= length(time_series)) {
                    if (verbose) cat(sprintf("Predicting epoch %d (h = %d)...\n", training_epoch + j + i - 1, h[i]))
                    y_hat[i, j + i - 1] <- predict(train_model, h = h[i])$est
                    diff_y_hat_obs[i, j + i - 1] = abs(time_series[training_epoch + j + i - 1] - y_hat[i, j + i - 1])
                }
            }
        }
    }
    
    mae <- numeric(length(h))
    mae <- rowMeans(diff_y_hat_obs, na.rm = TRUE)
    names(mae) <- paste0("MAE (h = ", h, ")")
    
    list(est = y_hat, mae = mae)
}
