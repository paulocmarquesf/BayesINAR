# inar.R

#' INAR model training
#' @description Computes the posterior distribution for the INAR model using a Gibbs sampler.
#' @param time_series A univariate time series.
#' @param prior List of prior hyperparameters:
#' \describe{
#'   \item{a_alpha, b_alpha}{Hyperparameters of the thinning component.}
#'   \item{a_lambda, b_lambda}{Hyperparameters of the Gamma prior for the Poisson innovation rate.}
#' }
#' @param burn_in Length of the "burn-in" period.
#' @param chain_length Length of the Markov chain.
#' @param random_seed Value of the random seed.
#' @param verbose Displays running status if \code{TRUE}.
#' @return \code{inar} returns an object of class "inar".
#' @export
inar <- function(time_series,
                 prior = list(a_alpha = NULL,
                              b_alpha = NULL,
                              a_lambda = NULL,
                              b_lambda = NULL),
                 burn_in = 10^3,
                 chain_length = 10^4,
                 random_seed = 1761,
                 verbose = TRUE)
{
    if (any(time_series %% 1 != 0)) stop("Time series must have only counts.")
    if (any(time_series < 0)) stop("Time series must have only positive counts.")
    if (!length(time_series) > 0) stop("Time series must have positive length.")
    if (!burn_in >= 0) stop("Burn-in must be positive.")
    if (!chain_length >= 0) stop("Chain length must be positive.")
    if (!random_seed >= 0) stop("Random seed must be positive.")
    if (!is.logical(verbose)) stop("Verbose parameter must be TRUE or FALSE.")
    p <- 1
    if (length(time_series) <= p) stop("Time series length must be bigger than p.")
    
    set.seed(random_seed)
    
    begin <- proc.time()
        
    if (is.null(prior[["a_lambda"]]) || is.null(prior[["b_lambda"]])) {
        prior[["a_lambda"]] <- 1
        prior[["b_lambda"]] <- 0.1
    }
    
    if (is.null(prior[["a_alpha"]])) prior[["a_alpha"]] <- 1
    if (is.null(prior[["b_alpha"]])) prior[["b_alpha"]] <- 1
    
    post <- .posterior_inar(time_series,
                            p,
                            prior[c("a_alpha",  "b_alpha", 
                                    "a_lambda", "b_lambda")],
                            burn_in,
                            chain_length,
                            random_seed,
                            verbose)
    
    model <- list(time_series = time_series,
                  prior = prior[c("a_alpha", "b_alpha",
                                  "a_lambda", "b_lambda")],
                  burn_in = list(length = burn_in,
                                 alpha = post$burn_in$alpha,
                                 lambda = post$burn_in$lambda),
                  chain = list(length = chain_length,
                               alpha = post$chain$alpha,
                               lambda = post$chain$lambda),
                  est = list(alpha = apply(post$chain$alpha, 2, mean),
                             lambda = mean(post$chain$lambda)),
                  elapsed = proc.time() - begin)
    
    class(model) <- "inar"
    
    invisible(model)
}

#' Predict method for INAR models
#' @description Predictions and predictive distribution from a trained INAR model object.
#' @param model A trained object of class "inar".
#' @param h Number of steps ahead to be predicted. 
#' @param replications Number of replications for each posterior sample.
#' @return A list with the following elements: 
#' \describe{
#'   \item{est}{The \code{h}-steps-ahead prediction.}
#'   \item{distr}{The \code{h}-steps-ahead predictive distribution.}
#' }
#' @export
predict.inar <- function(model, h = 1, replications = 10^3) {
    pred <- .predictive_distribution_inar(model, h, replications)
    list(est = .generalized_median(pred), distr = pred)
}

#' INAR model summaries
#' Summarizes a trained INAR model.
#' @return A model summary.
#' @export
summary.inar <- function(model) {
    printf <- function(...) cat(sprintf(...))
    printf("\n=====================================================\n")
    printf("                 INAR(1) Model Summary\n")
    printf("=====================================================\n")
    printf("Time series length: %d\n", length(model$time_series))
    printf("Prior parameters:\n")
    printf("  a_alpha = %.2f, b_alpha = %.2f\n", model$prior[["a_alpha"]], model$prior[["b_alpha"]])
    printf("  a_lambda = %.2f, b_lambda = %.2f\n", model$prior[["a_lambda"]], model$prior[["b_lambda"]])
    cat("  ")
    # for (i in 1:model$p) {
    #    printf("a_%d = %.2f", i, model$prior[["a_alpha"]][i])
    #    if (i < model$p) cat(", ")
    #    else cat("\n")
    # }
    printf("Effective Markov chains length: %d (%d burn-in)\n", model$chain$length, model$burn_in$length)
    printf("Posterior means with 95%% credible intervals:\n")
    for (i in 1:1) {
        post_qs <- unname(quantile(model$chain$alpha[, i], probs = c(0.025, 0.975)))
        printf("  alpha = %.2f  [ %.2f ; %.2f ]\n", mean(model$chain$alpha[, i]), post_qs[1], post_qs[2])
    }
    post_qs <- unname(quantile(model$chain$lambda, probs = c(0.025, 0.975)))
    printf("  lambda = %.2f  [ %.2f ; %.2f ]\n", mean(model$chain$lambda), post_qs[1], post_qs[2])
    
    printf("Total simulation time: %.2f seconds\n\n", round(model$elapsed[3]))
}
