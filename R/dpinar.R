# dpinar.R

#' DP-INAR model training
#' @description Computes the posterior distribution for the DP-INAR family using a Gibbs sampler.
#' @param time_series A univariate time series.
#' @param prior List of prior hyperparameters:
#' \describe{
#'   \item{a_alpha, b_alpha}{Hyperparameters of the thinning component.}
#'   \item{a_tau, b_tau}{Hyperparameters of the concentration parameter with Gamma prior.}
#'   \item{a0, b0}{Base measure hyperparameters.}
#'   \item{lambda_max}{Hyperparameter of the uniform distribution that minimizes the corresponding D-KL.}
#' }
#' @param burn_in Length of the "burn-in" period.
#' @param chain_length Length of the Markov chain.
#' @param random_seed Value of the random seed.
#' @param verbose Displays running status if \code{TRUE}.
#' @return \code{dpinar} returns an object of class "dpinar".
#' @export
dpinar <- function(time_series,
                   prior = list(a_alpha = NULL, b_alpha = NULL,
                                a_tau = NULL, b_tau = NULL,
                                a0 = NULL, b0 = NULL,
                                lambda_max = NULL),
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
    p <- 1
    if (length(time_series) <= p) stop("Time series length must be bigger than p.")
    if (!is.logical(verbose)) stop("Verbose parameter must be TRUE or FALSE.")

    set.seed(random_seed)
    
    begin <- proc.time()

    if (is.null(prior[["a0"]]) || is.null(prior[["b0"]])) {
        if (is.null(prior[["lambda_max"]])) prior[["lambda_max"]] <- max(time_series) # heuristic
        opt_G0 <- .KL_base_measure(prior[["lambda_max"]])
        prior[["a0"]] <- opt_G0[1]
        prior[["b0"]] <- opt_G0[2]
    }
    
    if (is.null(prior[["a_tau"]]) || is.null(prior[["b_tau"]])) {
        max_t <- length(time_series)
        if (max_t - p > 170) max_t <- 170 + p
        prior[["a_tau"]] <- prior_tau[max_t - p, 1]
        prior[["b_tau"]] <- prior_tau[max_t - p, 2]
    }

    if (is.null(prior[["a_alpha"]])) prior[["a_alpha"]] <- 1
    if (is.null(prior[["b_alpha"]])) prior[["b_alpha"]] <- 1

    post <- .posterior_dpinar(time_series,
                              p,
                              prior[c("a_alpha", "b_alpha",
                                      "a0", "b0", 
                                      "a_tau", "b_tau")],
                              burn_in,
                              chain_length,
                              random_seed,
                              verbose)

    model <- list(time_series = time_series,
                  prior = prior[c("a_alpha", "b_alpha", 
                                  "a0", "b0", 
                                  "a_tau", "b_tau")],
                  burn_in = list(length = burn_in,
                                 alpha = post$burn_in$alpha,
                                 lambda = post$burn_in$lambda,
                                 tau = post$burn_in$tau,
                                 num_clusters = post$burn_in$num_clusters),
                  chain = list(length = chain_length,
                               alpha = post$chain$alpha,
                               lambda = post$chain$lambda,
                               tau = post$chain$tau,
                               num_clusters = post$chain$num_clusters),
                  est = list(alpha = apply(post$chain$alpha, 2, mean),
                             lambda = apply(post$chain$lambda, 2, mean),
                             tau = mean(post$chain$tau),
                             num_clusters = which.max(tabulate(post$chain$num_clusters))),
                  elapsed = proc.time() - begin)

    class(model) <- "dpinar"

    invisible(model)
}

#' Predict methodfor DP-INAR models
#' @description Predictions and predictive distribution from a trained DP-INAR model.
#' @param model An object of class "dpinar".
#' @param h Number of steps ahead to be predicted. 
#' @param replications Number of replications for each posterior sample.
#' @return A list with the following elements: 
#' \describe{
#'   \item{est}{The \code{h}-steps-ahead prediction.}
#'    \item{distr}{The \code{h}-steps-ahead predictive distribution.}
#' }
#' @export
predict.dpinar <- function(model, h = 1, replications = 10^3) {
    pred <- .predictive_distribution_dpinar(model, h, replications)
    list(est = .generalized_median(pred), distr = pred)
}

#' DP-INAR model summaries
#' Summarizes a trained DP-INAR model
#' @return A model summary.
#' @export
summary.dpinar <- function(model) {
    printf <- function(...) cat(sprintf(...))
    printf("\n=======================================================\n")
    printf("                DP-INAR(1) Model Summary\n")
    printf("=======================================================\n")
    printf("Time series length: %d\n", length(model$time_series))
    printf("Prior parameters:\n")
    printf("  a_alpha = %.2f, b_alpha = %.2f\n", model$prior[["a_alpha"]], model$prior[["b_alpha"]])
    printf("  a0 = %.2f, b0 = %.2f\n", model$prior[["a0"]], model$prior[["b0"]])
    printf("  a_tau = %.2f, b_tau = %.2f\n", model$prior[["a_tau"]], model$prior[["b_tau"]])
    cat("  ")
    # for (i in 1:model$p) {
    #    printf("a_%d = %.2f", i, model$prior[["a_alpha"]][i])
    #    if (i < model$p) cat(", ")
    #    else cat("\n")
    #}
    printf("Effective Markov chains length: %d (%d burn-in)\n", model$chain$length, model$burn_in$length)
    printf("Some posterior means with 95%% credible intervals:\n")
    for (i in 1:1) {
        post_qs <- unname(quantile(model$chain$alpha[, i], probs = c(0.025, 0.975)))
        printf("  alpha = %.2f  [ %.2f ; %.2f ]\n", mean(model$chain$alpha[, i]), post_qs[1], post_qs[2])
    }
    for (i in as.integer(seq(1, length(model$time_series) - 1, length.out = 3))) {
        post_qs <- unname(quantile(model$chain$lambda[, i], probs = c(0.025, 0.975)))
        printf("  lambda_%d = %.2f  [ %.2f ; %.2f ]\n", i, mean(model$chain$lambda[, i]), post_qs[1], post_qs[2])
    }
    post_qs <- unname(quantile(model$chain$tau, probs = c(0.025, 0.975)))
    printf("  tau = %.2f  [ %.2f ; %.2f ]\n", mean(model$chain$tau), post_qs[1], post_qs[2])
    printf("Posterior distribution of number of clusters:")
    print(table(model$chain$num_clusters) / model$chain$length)
    printf("Total simulation time: %.2f seconds\n\n", round(model$elapsed[3]))
}

prior_tau <- structure(c(1.00003047761695, 1.21580601961878, 0.537166013561484, 
0.541601931925979, 0.540860850419251, 0.53832046074914, 0.535102855225786, 
0.53170617454297, 0.528356582933876, 0.525113782940635, 0.522062135576253, 
0.519163339074338, 0.516431921320194, 0.513862906141925, 0.511419415008973, 
0.509132532848726, 0.506991420733701, 0.504909251582994, 0.502676716243083, 
0.501135737647739, 0.499310831073104, 0.497516004467033, 0.495792269510977, 
0.493946227059194, 0.493094309229708, 0.491199719625483, 0.490020907895675, 
0.48825981265592, 0.48751165454133, 0.486512365632474, 0.485280026581696, 
0.48445173470005, 0.483413072303925, 0.482254241446387, 0.481501386912696, 
0.480564397802719, 0.479705891734397, 0.478852766754897, 0.477790535882215, 
0.477406859163749, 0.476693241524884, 0.476083540561368, 0.474691328359971, 
0.474977135379896, 0.472108332458303, 0.472830860431937, 0.473660388278329, 
0.473198090495573, 0.472858827503052, 0.472932456479189, 0.470175014019941, 
0.468800789228387, 0.468617264762362, 0.472143138126425, 0.470147332249766, 
0.470851691222197, 0.471689234419242, 0.46780600908226, 0.471688294171172, 
0.472284462866046, 0.469946654495544, 0.472656875939202, 0.472876677373171, 
0.473423916071545, 0.471969052862078, 0.472284783033178, 0.473408707817921, 
0.474063558528395, 0.472514730563527, 0.47643689817615, 0.474480762483295, 
0.47518406531086, 0.477109954332935, 0.477101691472723, 0.477122398764925, 
0.475584954474563, 0.473694837558509, 0.478393035596908, 0.47650175331117, 
0.479284942618648, 0.479180971896038, 0.475345823336798, 0.476004387081465, 
0.478574427246696, 0.481914874621391, 0.474952079473707, 0.476304562391001, 
0.482823313231777, 0.483798030087496, 0.484402702991, 0.478550152096152, 
0.483108043235274, 0.480816730010043, 0.480044299071912, 0.487797269521275, 
0.487860820994547, 0.480272807270598, 0.490962934218186, 0.481031875402759, 
0.490915815866415, 0.492098211693469, 0.48492704926164, 0.490245371356234, 
0.493450506337682, 0.485368976421061, 0.49233671614621, 0.483774132766943, 
0.497250491732231, 0.496905162404891, 0.490770871631395, 0.49257507321423, 
0.496782329428764, 0.494541176162056, 0.499924525173938, 0.499832973266828, 
0.500054991391687, 0.495364601364249, 0.502612881068253, 0.503304358034546, 
0.504010429784465, 0.505939191544911, 0.506886687126486, 0.507582529683884, 
0.50649869915902, 0.507108000362703, 0.507744536073884, 0.510143727826443, 
0.496555164156655, 0.500516030434019, 0.505459519811782, 0.510830721501961, 
0.504391952178782, 0.500496740079285, 0.498464839176673, 0.49840356903815, 
0.498960351333987, 0.499940114905123, 0.50072907485503, 0.502520307903773, 
0.504505725990444, 0.506525098594307, 0.508586937512437, 0.519467187699007, 
0.517195575600589, 0.510967890393031, 0.509737687949441, 0.508638752687178, 
0.50780683348679, 0.506369864180569, 0.522983355916744, 0.524487519832286, 
0.510352622756595, 0.513971681972424, 0.511625986937666, 0.527046318636252, 
0.517275136846339, 0.521345076826917, 0.520101989028733, 0.521058906174104, 
0.522903467993096, 0.523058252722035, 0.532261728345671, 0.528374202511047, 
0.528538896697241, 0.530231832792565, 0.531824479358988, 0.534522555718601, 
0.533693129237173, 0.525842787873773, 0.536606161256786, 1.00000027798364, 
0.809021581610592, 0.151453922246636, 0.118617019252711, 0.096052506770653, 
0.079962258209283, 0.0680346264119503, 0.0589107722239956, 0.0517500689882478, 
0.0459903285077036, 0.0412987803095725, 0.0373918179354324, 0.034101137134671, 
0.0312971751658861, 0.0288796087957871, 0.0267808449871318, 0.0249422762905001, 
0.023318129578923, 0.0218468050449398, 0.0205906634864367, 0.0194284807647523, 
0.0183838206326149, 0.0174266580731683, 0.0165487210464474, 0.0157925922640236, 
0.0150423744089196, 0.014372090402681, 0.0137380528248571, 0.0132063786924342, 
0.0126997078041488, 0.012207765340201, 0.0117649292010981, 0.0113416778757036, 
0.0109378434976092, 0.0105734758207377, 0.0102288889891131, 0.00990317066984678, 
0.00959650760240495, 0.0092991881890121, 0.0090379274923301, 
0.00878194475454556, 0.00854123212012207, 0.00828975685713939, 
0.00809936863543058, 0.0078014059515205, 0.00763124971906557, 
0.00752164541826327, 0.00734754734512471, 0.00718304970512278, 
0.0070269988431637, 0.00679296095152791, 0.00663803242139151, 
0.00653652681189976, 0.00648426149893308, 0.00631890220159036, 
0.0061570508615102, 0.00610839428059256, 0.00593846922335903, 
0.00585809547402358, 0.00583885388725589, 0.00568385745025897, 
0.00566047580154399, 0.00557701947284373, 0.00549486317848268, 
0.00532199394602384, 0.00524451619537018, 0.0052427508713303, 
0.0051844620339656, 0.00505066093741691, 0.00507104267571499, 
0.005006270942023, 0.00495240110582245, 0.00484814727856456, 
0.00478758110851561, 0.00471856500022095, 0.00467114030186651, 
0.00457659827613521, 0.00464655376901905, 0.00456959769898406, 
0.00455429773706306, 0.00449431467650636, 0.00442727640064569, 
0.00431971207779543, 0.00436365856315858, 0.00435018348590272, 
0.00422789538042264, 0.00412850208357611, 0.00421783735956195, 
0.00417975878393186, 0.00417140695507057, 0.00408687714388045, 
0.00402540309947019, 0.00400471942568552, 0.00390850150248776, 
0.00402212927238359, 0.00399051387090875, 0.00390310026806325, 
0.00393418819633677, 0.00375039414929486, 0.00388828672130386, 
0.00385916840542315, 0.00377685676873691, 0.003787278481519, 
0.00379152021113444, 0.00369542110414178, 0.00367128335759583, 
0.00358946661341874, 0.00365113675829104, 0.0036729499700925, 
0.00351184275889402, 0.00348476005353558, 0.00361427173139389, 
0.00345850752208029, 0.00358381238427815, 0.00355509646443801, 
0.0035413089884019, 0.00348793296835683, 0.00351283318103405, 
0.00349587375699187, 0.00347929570868224, 0.00344204465351712, 
0.00342913260756867, 0.00340402433381922, 0.00341452897102608, 
0.00339859283752854, 0.0033833037907891, 0.00334275017201284, 
0.00318346681653648, 0.00317266319227946, 0.00327442874531207, 
0.00325832205727319, 0.00315568832654417, 0.00311426671871366, 
0.00310089693840792, 0.00310519634195542, 0.00310844288975392, 
0.00310972583369423, 0.00310480202210535, 0.00310352746872411, 
0.00310318663765634, 0.0030980811593062, 0.00309239833881053, 
0.00316863951175012, 0.00313701217245866, 0.00295897748389003, 
0.00295440731390282, 0.00293712684073114, 0.00293152778475587, 
0.00296520506354453, 0.00307991935785115, 0.00307672960296532, 
0.00287435186551589, 0.00291956131355235, 0.0028638898543563, 
0.00303187077094587, 0.00291893018880688, 0.00295867149829054, 
0.00283603698923119, 0.00282462109231953, 0.00283865648080019, 
0.00281164356216553, 0.00298508131680225, 0.002894304738283, 
0.00288533332777725, 0.00288656795800847, 0.00289075317068731, 
0.00293816365203622, 0.00291858915394676, 0.00278528621851238, 
0.00291475161791699), .Dim = c(170L, 2L), .Dimnames = list(NULL, c("a_tau", "b_tau")))

.KL_base_measure <- function(lambda_max) {
    D_KL <- function(a0, b0, lambda_max) lgamma(a0) - a0*log(b0) - (a0 - 1)*(log(lambda_max) - 1) + b0*(lambda_max / 2) - log(lambda_max)
    optim(c(1, 1), function(x) D_KL(x[1], x[2], lambda_max), method = "L-BFGS-B", lower = c(1e-6, 1e-6))$par
}
