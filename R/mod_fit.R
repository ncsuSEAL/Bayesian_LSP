#*******************************************************************************
# Description: Model fit functions.
# Date: 2022-03-29
#*******************************************************************************


#' Get the BLSP fitted curves
#' 
#' The funciton uses the "BlspFit" object returned by the "FitBLSP" function to
#' predict the fitted daily VI curves. 
#' 
#' @param blsp_fit The "BlspFit" object.
#' @param asCI Logical. Default is `TRUE`, which means the fitted curves will be
#' summarized as median and 95% credibel interval.
#'
#' @return The fitted curves, or the summary, for each year.
#'
#' @export
#' @import data.table
BLSPFitted <- function(blsp_fit, asCI = FALSE) {
    if (class(blsp_fit) != "BlspFit") {
        stop("Please provide a BlspFit object!")
    }

    # Unpack data from the object
    date_vec <- blsp_fit$data$date
    vi_vec <- blsp_fit$data$vi
    weights_vec <- blsp_fit$weights
    if (is.null(weights_vec)) {
        weights_vec <- rep(1, length(vi_vec))
    }
    bf_phenos <- blsp_fit$phenos
    yr <- lubridate::year(date_vec) - lubridate::year(date_vec)[1] + 1
    numYears <- length(unique(yr))
    model_str <- blsp_fit$model$model_str

    # ~ Predict fitted value for full dates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bf_pred <- list()
    years <- sort(unique(lubridate::year(date_vec)))
    for (i in 1:numYears) {
        date <- seq(as.Date(paste0(years[i], "-01-01")),
            as.Date(paste0(years[i], "-12-31")),
            by = "day"
        )
        bf_params <- lapply(1:length(blsp_fit$params), function(m) {
            blsp_fit$params[[m]][, i]
        }) %>%
            do.call(cbind, .)

        pred <- matrix(NA, nrow = length(date), ncol = nrow(bf_params))
        for (j in 1:nrow(bf_params)) { 
            # pred based on current parameter samples
            envlist <- lapply(1:ncol(bf_params), function(m) {
                as.numeric(bf_params[j, m])
            })
            names(envlist) <- blsp_fit$model$model_param
            envlist$t <- 1:length(date)
            pred[, j] <- eval(str2expression(model_str), envir = envlist)
        }
        thelist <- list(pred)
        names(thelist) <- years[i]
        bf_pred <- append(bf_pred, thelist)
    }

    # If summarize as CI, the result will just be median, lower, and upper
    if (asCI == TRUE) {
        alpha <- (1 - blsp_fit$cred_int_level) / 2
        
        bf_pred <- lapply(seq_along(years), function(i) {
            date <- seq(as.Date(paste0(years[i], "-01-01")),
                as.Date(paste0(years[i], "-12-31")),
                by = "day"
            )
            theyr <- bf_pred[[as.character(years[i])]]
            predCI <- t(data.table::data.table(
                apply(theyr, 1, function(x) {
                    stats::quantile(x, c(alpha, 1 - alpha))
                })
            ))

            pred <- data.table::data.table(
                apply(theyr, 1, function(x) stats::quantile(x, 0.5))
            )
            return(cbind(date, pred, predCI))
        })
        bf_pred <- do.call(rbind, bf_pred)
        # Make it a data.table
        bf_pred <- data.table::as.data.table(bf_pred)
        colnames(bf_pred) <- c("Date", "Fitted", "Fitted_lower", "Fitted_upper")
        bf_pred$Date <- as.Date(bf_pred$Date, origin = "1970-01-01")
    }

    return(bf_pred)
}



#' Fit a Bayesian mixed hierarchical land urface phenology model.
#' 
#' This function fits a Bayesian mixed hierarchical land surface phenology model 
#' to the supplied data (can be sparse), and returns phenometrics for the 
#' entire time frame. For further explanation, please see the vignette.
#' 
#' @param date_vec The date vector, be sure to convert the vector to "Date" 
#' format or use "yyyy-mm-dd" format string.
#' @param vi_vec The vegetation index vector.
#' @param weights_vec A numeric vector of same length as vi_vec specifying the 
#' weights for the supplied observations. Must be between 0 and 1, inclusive.
#' @param model A string indicating the model name. For now, only support
#' "dblog7" and "dblog6" for the 7- and 6-parameter double-logistic functions.
#' @param init_values Initial values for MCMC sampling. By default, it is 
#' assgined `NULL`. It could also be an object returned from the `FitAvgModel()` 
#' function that fits an averaged model or a numeric vector provided by the user. 
#' @param start_yr The start year of the result. Default is NULL, which means
#' determined by data.
#' @param end_yr The end year of the result. Default is NULL, which means
#' determined by data.
#' @param cred_int_level A scalar value from 0 to 1 (exclusive) that specifies
#' the level for equal-tailed credible intervals of the estimated phenometrics.
#' The default level is 0.9, generating `90%` credible intervals. The end
#' points of these intervals define the upper and lower bounds for the estimated
#' phenometrics.
#' @param opt An option list that contains additional configurations. For now,
#' only support `list(method = "threshold")` to indicate that use the
#' threshold-based method to retrive phenometrics instead of the default SOS and
#' EOS. The threshold-based method will produce 7 phenometrics including
#' Greenup, MidGreenup, Maturity, Peak, Senescence, MidGreendown, and Dormancy
#' using VI amplitude thresholds of 15%, 50%, 90%, and 100%, respectively. 
#' @param verbose logical. If `TRUE`, the progress will be reported.
#' 
#' @return An object of class `BlspFit` will be returned. The object contains the
#' estimated spring and autumn phenometrics for each year, the generated model 
#' parameter samples, and the input data.
#' @examples
#' \dontrun{
#' data(landsatEVI2)
#' blsp_fit <- FitBLSP(date_vec = landsatEVI2$date, vi_vec = landsatEVI2$evi2)
#' }
#' @export 
#' @import data.table
FitBLSP <- function(date_vec, vi_vec, 
    weights_vec = NULL,
    model = "dblog7",
    init_values = NULL,
    start_yr = NULL, 
    end_yr = NULL,
    cred_int_level = 0.9,
    opt = NULL,
    verbose = FALSE
) {
    sc_li <- SanityCheck(date_vec, vi_vec, weights_vec, 
        model, cred_int_level, init_values
    )
    model <- sc_li$model
    cred_int_level <- sc_li$cred_int_level

    # Reorder data to make sure they are sorted by time
    od <- order(date_vec)
    date_vec <- date_vec[od]
    vi_vec <- vi_vec[od]
    weights_vec <- weights_vec[od]

    # Convert data to jags format
    y <- vi_vec
    t <- lubridate::yday(date_vec)
    n <- length(y) # total num of observations
    # year id vector
    if (is.null(start_yr) || is.null(end_yr)) {
        yr <- lubridate::year(date_vec) - lubridate::year(date_vec)[1] + 1
        tmp <- sort(unique(year(date_vec)))
        years <- tmp[1]:tmp[length(tmp)]
        numYears <- length(1:yr[length(yr)])
    } else {
        yr <- lubridate::year(date_vec) - lubridate::year(date_vec)[1] + 1
        years <- start_yr:end_yr
        numYears <- length(years)
    }
    # If user specified weights
    if (is.null(weights_vec)) {
        weights_vec <- rep(1, n)
    }

    # Get model structure
    mod <- GetModel(model)

    # Set initial values
    if (!is.null(init_values) && class(init_values) == "nls") {
        p_pars <- stats::coef(init_values)
    } else if (!is.null(init_values) && class(init_values) == "numeric") {
        p_pars <- init_values
    } else {
        p_pars <- mod$init_val
    }

    data <- list(
        Y = y, t = t, n = n, yr = yr, N = numYears, weights = weights_vec
    )


    inits <- lapply(seq_along(mod$init_param), function(i) {
        rep(p_pars[i], numYears)
    })
    names(inits) <- mod$init_param

    blsp_fit <- tryCatch(
        {
            if (verbose) {
                message("Initialize model...")
            }
            pb_type <- ifelse(verbose, "text", "none")

            model <- rjags::jags.model(textConnection(mod$jags_str),
                data = data, inits = inits,
                n.chains = 3, quiet = TRUE
            )
            stats::update(model, 2000, progress.bar = pb_type)

            if (verbose) {
                message("Sampling (could have multiple chains)...")
            }

            iteration_times <- 0
            repeat {
                samp <- rjags::coda.samples(model,
                    variable.names = mod$out_param,
                    n.iter = 5000,
                    thin = 10,
                    progress.bar = pb_type
                )
                iteration_times <- iteration_times + 5000

                # Try to make it converge
                if (coda::gelman.diag(samp)$mpsrf <= 1.3 |
                    iteration_times > 50000) {
                    break
                }
            }
            
            if (verbose) {
                message("total interation times:", iteration_times)
            }

            # ~ Retrieve parameter estimates
            if (verbose) {
                message("Estimate phenometrics...")
            }
            p_li <- lapply(seq_along(mod$out_param), function(p) {
                p_mat <- lapply(1:numYears, function(i) {
                    c(
                        samp[[1]][, paste0(mod$out_param[p], "[", i, "]")],
                        samp[[2]][, paste0(mod$out_param[p], "[", i, "]")],
                        samp[[3]][, paste0(mod$out_param[p], "[", i, "]")]
                    )
                }) %>%
                    do.call(cbind, .)

                return(p_mat)
            })

            # Construct `blsp_fit` object to return
            if (!is.null(opt) && tolower(opt$method) == "threshold") {
                blsp_fit <- CalPhenoThresh(
                    p_li, mod,
                    years, numYears, 
                    date_vec, vi_vec, weights_vec, 
                    cred_int_level
                )
            } else {
                blsp_fit <- CalPhenoParam(
                    p_li, mod,
                    years, numYears,
                    date_vec, vi_vec, weights_vec,
                    cred_int_level
                )
            }

            if (verbose) {
                message("Done!")
            }

            blsp_fit
        },
        error = function(e) {
            if (verbose) {
                message("Something went wrong!")
            }
            blsp_fit <- EmptyBlspOutput(
                years,
                date_vec,
                vi_vec,
                weights_vec,
                cred_int_level,
                mod,
                method = tolower(opt$method)
            )

            return(blsp_fit)
        }
    )

    return(blsp_fit)
}


#' Fit a Bayesian mixed hierarchical land surface phenology model. Spring only!
#' Note that the result CANNOT be used to plot the fit.
#' 
#' This function fits a Bayesian mixed hierarchical land surface phenology model 
#' to the supplied data (can be sparse), and returns phenometrics for the 
#' entire time frame. For further explanation, please see the vignette.
#' 
#' @param date_vec The date vector, be sure to convert the vector to "Date" 
#' format or use "yyyy-mm-dd" format string.
#' @param vi_vec The vegetation index vector.
#' @param weights_vec A numeric vector of same length as vi_vec specifying the 
#' weights for the supplied observations. Must be between 0 and 1, inclusive.
#' @param initValues Initial values for MCMC sampling. By default, it is 
#' assgined `NULL`. It could also be an object returned from the `FitAvgModel()` 
#' function that fits an averaged model or a numeric vector provided by the user. 
#' @param cred_int_level A scalar value from 0 to 1 (exclusive) that specifies
#' the level for equal-tailed credible intervals of the estimated phenometrics.
#' The default level is 0.9, generating `90%` credible intervals. The end
#' points of these intervals define the upper and lower bounds for the estimated
#' phenometrics.
#' @param verbose logical. If `TRUE`, the progress will be reported.
#' @return An object of class `BlspFit` will be returned. The object contains the
#' estimated spring and autumn phenometrics for each year, the generated model 
#' parameter samples, and the input data.
#' @examples
#' \dontrun{
#' data(landsatEVI2)
#' blsp_fit <- FitBLSP(date_vec = landsatEVI2$date, vi_vec = landsatEVI2$evi2)
#' }
#' @export 
#' @import data.table
FitBLSP_spring <- function(date_vec, vi_vec, 
    weights_vec = NULL, 
    initValues = NULL,
    cred_int_level = 0.9,
    verbose = FALSE
) {
    # Check if date_vec is in Date format
    if (sum(!is.na(lubridate::parse_date_time(date_vec, orders = "ymd"))) != 
        length(date_vec)) {
        stop("There're invalid Date values in the `date_vec`! 
            Be sure to use `yyyy-mm-dd` format.")
    }

    # Check weights to be in the range of [0, 1]
    if (!is.null(weights_vec)) {
        if (min(weights_vec) < 0 | max(weights_vec) > 1) {
            stop("Weights must be within [0, 1].")
        }

        # Check the length of dates, vis, and weights
        if (length(weights_vec) != length(date_vec) |
            length(vi_vec) != length(date_vec)) {
            stop("date_vec, vi_vec, and weights_vec have different lengths.")
        }
    }

    # Check NAs
    if (any(is.na(date_vec)) | any(is.na(vi_vec))) {
        stop("Please remove NAs in the input data.")
    }
    
    # Check if credible interval level is valid
    if (cred_int_level <= 0 | cred_int_level >= 1) {
      stop("Credible interval level must be a value between 0 and 1 (exclusive).")
    }
    # Reorder data to make sure they are sorted by time
    od <- order(date_vec)
    date_vec <- date_vec[od]
    vi_vec <- vi_vec[od]
    weights_vec <- weights_vec[od]

    # Convert data to jags format
    y <- vi_vec
    t <- lubridate::yday(date_vec)
    n <- length(y) # total num of observations
    # year id vector
    yr <- lubridate::year(date_vec) - lubridate::year(date_vec)[1] + 1 
    numYears <- length(unique(yr))

    # If user specified weights
    if (is.null(weights_vec)) {
        weights_vec <- rep(1, n)
    }

    # ~ Format data, inits, and model
    model_string <- "model {
        # Likelihood
        for (i in 1:n) {
            Y[i] ~ dnorm(mu[i], tau_y)
            mu[i] <- weights[i] * (m1[yr[i]] + (m2[yr[i]] - m7[yr[i]] * t[i]) * 
                ((1 / (1 + exp((m3[yr[i]] - t[i]) / m4[yr[i]]))) - 
                (1 / (1 + exp((m5[yr[i]] - t[i]) / m6[yr[i]])))))
        }
    
        # Priors
        for (j in 1:N) {
            M1[j] ~ dnorm(mu_m1, tau[1])
            logit(m1[j]) <- M1[j]
            m2[j] ~ dnorm(mu_m2, tau[2])
            m3[j] ~ dnorm(mu_m3, tau[3])
            m4[j] ~ dnorm(mu_m4, tau[4])
            m5[j] ~ dnorm(mu_m5, tau[5])
            m6[j] ~ dnorm(mu_m6, tau[6])
            M7[j] ~ dbeta(4, 4 * (1 - mu_m7 * 100) / (mu_m7 * 100))
            m7[j] <- M7[j] / 100
        }
    
        mu_m1 ~ dunif(0, 0.3)
        mu_m2 ~ dunif(0.5, 2)
        mu_m3 ~ dunif(0, 185)
        mu_m4 ~ dunif(1, 15)
        mu_m5 ~ dunif(185, 366)
        mu_m6 ~ dunif(1, 15)
        mu_m7 ~ dunif(0, 0.01)
    
        for (k in 1:7) {
            tau[k] ~ dgamma(0.1, 0.1)
        }
        tau_y ~ dgamma(0.1, 0.1)
    }"

    if (!is.null(initValues) && class(initValues) == "nls") {
        p_m1 <- stats::coef(initValues)["m1"]
        p_m2 <- stats::coef(initValues)["m2"]
        p_m3 <- stats::coef(initValues)["m3"]
        p_m4 <- stats::coef(initValues)["m4"]
        p_m5 <- stats::coef(initValues)["m5"]
        p_m6 <- stats::coef(initValues)["m6"]
        p_m7 <- stats::coef(initValues)["m7"]
    } else if (!is.null(initValues) && class(initValues) == "numeric") {
        if (length(initValues) != 7) {
            stop("The length of the initial values does not match",
                "the number of model parameters."
            )
        }
        p_m1 <- initValues[1]
        p_m2 <- initValues[2]
        p_m3 <- initValues[3]
        p_m4 <- initValues[4]
        p_m5 <- initValues[5]
        p_m6 <- initValues[6]
        p_m7 <- initValues[7]
    } else {
        p_m1 <- 0.05
        p_m2 <- 1
        p_m3 <- 120
        p_m4 <- 8
        p_m5 <- 290
        p_m6 <- 8
        p_m7 <- 0.001
    }

    data <- list(
        Y = y, t = t, n = n, yr = yr, N = numYears, weights = weights_vec
    )

    inits <- list(
        M1 = rep(p_m1, numYears),
        m2 = rep(p_m2, numYears), m3 = rep(p_m3, numYears),
        m4 = rep(p_m4, numYears), m5 = rep(p_m5, numYears), 
        m6 = rep(p_m6, numYears)
    )

    tryCatch(
        {
            if (verbose) {
                message("Initialize model...")
            }
            pb_type <- ifelse(verbose, "text", "none")

            model <- rjags::jags.model(textConnection(model_string),
                data = data, inits = inits,
                n.chains = 3, quiet = TRUE
            )
            stats::update(model, 2000, progress.bar = pb_type)

            if (verbose) {
                message("Sampling (could have multiple chains)...")
            }

            iteration_times <- 0
            repeat {
                samp <- rjags::coda.samples(model,
                    variable.names = c("m2", "m3"),
                    n.iter = 5000,
                    thin = 10,
                    progress.bar = pb_type
                )
                iteration_times <- iteration_times + 5000
                
                # Try to make it converge
                if(coda::gelman.diag(samp)$mpsrf <= 1.3 | 
                    iteration_times > 100000) {
                    break
                }
            }
            if (verbose) {
                message("total interation times:", iteration_times)
            }
        },
        error = function(e) {
            years <- sort(unique(year(date_vec)))
            bf_phenos <- NULL
            for (i in 1:numYears) {
                bf_phenos <- rbind(bf_phenos, list(
                    Id = NA, Year = years[i],
                    midgup_lower = NA, midgup = NA, midgup_upper = NA,
                    midgdown_lower = NA, midgdown = NA, midgdown_upper = NA
                ))
            }
            return(list(fitted = NA, phenos = bf_phenos))
        }
    )

    # ~ Retrieve parameter estimates
    if (verbose) {
        message("Estimate phenometrics...")
    }
    m1 <- m2 <- m3 <- m4 <- m5 <- m6 <- m7 <- NULL
    for (i in 1:numYears) {
        m2 <- cbind(m2, c(samp[[1]][, paste0("m2", "[", i, "]")], 
            samp[[2]][, paste0("m2", "[", i, "]")]))
        m3 <- cbind(m3, c(samp[[1]][, paste0("m3", "[", i, "]")], 
            samp[[2]][, paste0("m3", "[", i, "]")]))
    }

    alpha <- (1-cred_int_level)/2
    m2_quan <- data.table::data.table(
        apply(m2, 2, stats::quantile, c(alpha, 0.5, 1-alpha))
    )
    m3_quan <- data.table::data.table(
        apply(m3, 2, stats::quantile, c(alpha, 0.5, 1-alpha))
    )
    
    years <- sort(unique(lubridate::year(date_vec)))
    bf_phenos <- NULL
    for (i in 1:numYears) {
        if (m2_quan[2, ][[i]] > 0.4) { # suppress some amplitude-too-low year
            bf_phenos <- rbind(bf_phenos, data.table::data.table(
                Year = years[i],
                midgup_lower = m3_quan[1, ][[i]], 
                midgup = m3_quan[2, ][[i]], 
                midgup_upper = m3_quan[3, ][[i]]
            ))
        } else {
            bf_phenos <- rbind(bf_phenos, data.table::data.table(
                Year = years[i],
                midgup_lower = NA, 
                midgup = NA, 
                midgup_upper = NA
            ))
        }
    }

    # Construct `blsp_fit` object to return
    blsp_fit <- list(
        phenos = bf_phenos,
        params = list(m2 = m2, m3 = m3),
        data = data.table::data.table(
            date = date_vec, 
            vi = vi_vec, 
            weights = weights_vec
        ),
        cred_int_level = cred_int_level
    )
    class(blsp_fit) <- "BlspFit"

    if (verbose) {
        message("Done!")
    }
    return(blsp_fit)
}


#' Generate pheno from the predicted curve.
#' 
#' Only supports Elmore model (The double-logistic model used in BLSP). 
#' This function is used inside the FitBLSP function.
#' 
#' @param equation The model equation.
#' @param params The Parameter list.
#' @param t Date vector.
#' @return The phenological timing in day of year (DOY)
#' 
#' @noRd
GetPhenosIdx <- function(equation, params, t) {
    y <- eval(equation, envir = list(
        m1 = params[[1]], m2 = params[[2]], m3 = params[[3]], m4 = params[[4]],
        m5 = params[[5]], m6 = params[[6]], m7 = params[[7]], t = t
    ))

    d1 <- stats::D(equation, "t")
    d2 <- stats::D(d1, "t")
    
    y1 <- eval(d1, envir = list(
        m1 = params[[1]], m2 = params[[2]], m3 = params[[3]], m4 = params[[4]],
        m5 = params[[5]], m6 = params[[6]], m7 = params[[7]], t = t
    ))
    y2 <- eval(d2, envir = list(
        m1 = params[[1]], m2 = params[[2]], m3 = params[[3]], m4 = params[[4]],
        m5 = params[[5]], m6 = params[[6]], m7 = params[[7]], t = t
    ))
    k <- abs(y2^2) / ((1 + y1^2)^(3 / 2))

    d <- diff(y)
    d_code <- (d > 0) + (2 * (d < 0)) # 0=no change, 1=inc, 2=dec

    # Find the MidGreenup and MidGreeendown
    midgup <- round(params[[3]])
    midgdown <- round(params[[5]])

    # find peak
    peak <- NULL
    peaks <- unlist(gregexpr("12", paste(d_code, collapse = ""))) # no match is -1
    if (peaks[1] == -1) peaks <- NULL
    # no match is -1
    flat_peaks <- unlist(gregexpr("10+2", paste(d_code, collapse = ""))) 
    if (flat_peaks[1] == -1) flat_peaks <- NULL

    if (is.null(peaks) & is.null(flat_peaks)) {
        print("no peaks found in in GetPhenoIdx function!")
        return(NULL)
    } else {
       peak <- ifelse (!is.null(peaks), peaks[1], flat_peaks[1])
    }

    # the derivative of curvature
    d_k <- c(0, 0, diff(k * 1000, differences = 2))
    # local extremes
    localextr <- which(diff(sign(diff(d_k * 1000))) == -2) + 1
    
    maturity <- localextr[which(localextr < peak & localextr > midgup)]
    maturity <- maturity[length(maturity)]
    
    sene <- localextr[which(localextr > peak & localextr < midgdown)][1]
    
    gup <- localextr[which(localextr < midgup)]
    gup <- gup[length(gup)]
    
    dormancy <- localextr[which(localextr > (midgdown + 7))][1]
    
    return(c(gup = gup, midgup = midgup, maturity = maturity, peak = peak, 
        sene = sene, midgdown = midgdown, dormancy = dormancy))
}


#' Fit the averaged model and get the model parameters.
#' 
#' @param date_vec the date vector, be sure to convert the vector to "Date" 
#' format or use "yyyy-mm-dd" format string.
#' @param vi_vec The vegetation index vector.
#' @param model A string indicating the model name. For now, only support
#' "dblog7" and "dblog6" for the 7- and 6-parameter double-logistic functions.
#' @return Model parameters to be used as MCMC initial parameters in 
#' the FitBLSP function.
#' @export 
FitAvgModel <- function(date_vec, vi_vec, model = "dblog7") {
    # Check the jags model string
    if (!model %in% c("dblog7", "dblog6")) {
        warning("The specified model does not exist, dblog7 will be used.")
        model <- "dblog7"
    }
    
    # Format data
    avg_dt <- FormatAvgData(date_vec, vi_vec)

    # Unpack data
    y <- unlist(avg_dt$evi2)
    t <- unlist(avg_dt$date)
    
    # Fake year information for the averaged year
    cur_start_date <- as.Date("1970-01-01")
    cur_end_date <- as.Date("1970-12-31")
    full_date <- seq(cur_start_date, cur_end_date, by = "day")
    full_t <- as.integer(full_date - cur_start_date + 1)

    # Fit model to get the prior
    avg_fit <- tryCatch(
        {
            mod <- GetModel(model)
            model_str <- mod$model_str
            model_equ <- stats::as.formula(paste("VI", "~", model_str))
            
            start_li <- as.list(mod$init_val)
            names(start_li) <- mod$out_param

            minpack.lm::nlsLM(model_equ, 
                data = list(
                    VI = y,
                    t = as.integer(t)
                ),
                start = start_li,
                lower = mod$init_lwr,
                upper = mod$init_upr,
                control = list(warnOnly = TRUE)
            )
        },
        error = function(e) {
            print(paste("Average fit failed", sep = ":"))
            return(NULL)
        }
    )

    return(avg_fit)
}

