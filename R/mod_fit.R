#************************************************************************************
# Description: Model fit functions.
# Date: 2022-03-29
#************************************************************************************


#' Bayesian mixed hierarchical land surface phenology model.
#' 
#' @param date_vec The date vector, be sure to convert the vector to "Date" format 
#' or use "yyyy-mm-dd" format string.
#' @param vi_vec The vegetation index vector.
#' @param weights_vec For specifying weights to observations. For example, 
#' lower weights can be set to observations with snow. 
#' @param initValues Initial values for MCMC sampling. We get these values from 
#' fitting the averaged model. It could also be `NULL`.
#' @param ifplot: logical. Plot the model fit if TRUE. Note that the fitted curve 
#' with CI will only be returned when `ifplot` is TRUE.
#' @return retrieved phenometrics for each year.
#' @export 
FitBLSP <- function(date_vec, vi_vec, 
    weights_vec = NULL, initValues = NULL, 
    ifplot = FALSE
) {
    # Check if date_vec is in Date format
    if (sum(!is.na(lubridate::parse_date_time(date_vec, orders = "ymd"))) != 
        length(date_vec)) {
        stop("There're invalid Date values in the `date_vec`! 
            Be sure to use `yyyy-mm-dd` format.")
    }
    
    # Convert data to jags format
    y <- vi_vec
    t <- as.numeric(date_vec - as.Date(paste0(year(date_vec), "-01-01"))) + 1
    n <- length(y) # total num of observations
    yr <- year(date_vec) - year(date_vec)[1] + 1 # year id vector
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

    if (!is.null(initValues)) {
        p_m1 <- coef(initValues)["m1"]
        p_m2 <- coef(initValues)["m2"]
        p_m3 <- coef(initValues)["m3"]
        p_m4 <- coef(initValues)["m4"]
        p_m5 <- coef(initValues)["m5"]
        p_m6 <- coef(initValues)["m6"]
        p_m7 <- coef(initValues)["m7"]
    } else {
        p_m1 <- 0.05
        p_m2 <- 1
        p_m3 <- 120
        p_m4 <- 8
        p_m5 <- 290
        p_m6 <- 8
        p_m7 <- 0.001
    }

    data <- list(Y = y, t = t, n = n, yr = yr, N = numYears, weights = weights_vec)

    inits <- list(
        M1 = rep(p_m1, numYears),
        m2 = rep(p_m2, numYears), m3 = rep(p_m3, numYears),
        m4 = rep(p_m4, numYears), m5 = rep(p_m5, numYears), 
        m6 = rep(p_m6, numYears)
    )

    tryCatch(
        {
            model <- rjags::jags.model(textConnection(model_string),
                data = data, inits = inits,
                n.chains = 3, quiet = TRUE
            )
            update(model, 2000, progress.bar = "none")
            iteration_times <- 0
            repeat {
                samp <- rjags::coda.samples(model,
                    variable.names = c("m1", "m2", "m3", "m4", "m5", "m6", "m7"),
                    n.iter = 5000,
                    thin = 10,
                    progress.bar = "none"
                )
                iteration_times <- iteration_times + 5000
                
                # Try to make it converge
                if(coda::gelman.diag(samp)$mpsrf <= 1.3 | 
                    iteration_times > 100000) {
                    break
                }
            }
            print(iteration_times)
        },
        error = function(e) {
            years <- sort(unique(year(landsat$date)))
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
    m1 <- m2 <- m3 <- m4 <- m5 <- m6 <- m7 <- NULL
    for (i in 1:numYears) {
        m1 <- cbind(m1, c(samp[[1]][, paste0("m1", "[", i, "]")], 
            samp[[2]][, paste0("m1", "[", i, "]")]))
        m2 <- cbind(m2, c(samp[[1]][, paste0("m2", "[", i, "]")], 
            samp[[2]][, paste0("m2", "[", i, "]")]))
        m3 <- cbind(m3, c(samp[[1]][, paste0("m3", "[", i, "]")], 
            samp[[2]][, paste0("m3", "[", i, "]")]))
        m4 <- cbind(m4, c(samp[[1]][, paste0("m4", "[", i, "]")], 
            samp[[2]][, paste0("m4", "[", i, "]")]))
        m5 <- cbind(m5, c(samp[[1]][, paste0("m5", "[", i, "]")], 
            samp[[2]][, paste0("m5", "[", i, "]")]))
        m6 <- cbind(m6, c(samp[[1]][, paste0("m6", "[", i, "]")], 
            samp[[2]][, paste0("m6", "[", i, "]")]))
        m7 <- cbind(m7, c(samp[[1]][, paste0("m7", "[", i, "]")], 
            samp[[2]][, paste0("m7", "[", i, "]")]))
    }

    m1_quan <- data.table::data.table(apply(m1, 2, quantile, c(0.05, 0.5, 0.95)))
    m2_quan <- data.table::data.table(apply(m2, 2, quantile, c(0.05, 0.5, 0.95)))
    m3_quan <- data.table::data.table(apply(m3, 2, quantile, c(0.05, 0.5, 0.95)))
    m5_quan <- data.table::data.table(apply(m5, 2, quantile, c(0.05, 0.5, 0.95)))
    
    years <- sort(unique(year(date_vec)))
    bf_phenos <- NULL
    for (i in 1:numYears) {
        if (m2_quan[2, ][[i]] > 0.4) { # suppress some amplitude-too-low year
            bf_phenos <- rbind(bf_phenos, list(
                Year = years[i],
                midgup_lower = m3_quan[1, ][[i]], midgup = m3_quan[2, ][[i]], 
                    midgup_upper = m3_quan[3, ][[i]],
                midgdown_lower = m5_quan[1, ][[i]], midgdown = m5_quan[2, ][[i]], 
                    midgdown_upper = m5_quan[3, ][[i]]
            ))
        } else {
            bf_phenos <- rbind(bf_phenos, list(
                Year = years[i],
                midgup_lower = NA, midgup = NA, midgup_upper = NA,
                midgdown_lower = NA, midgdown = NA, midgdown_upper = NA
            ))
        }
    }


    # ~ If visualize the fit
    # Note: when plotting out the fit, the fitted curve with CI will be returned.
    bf_pred <- NULL
    if (ifplot == TRUE) { # fig: Bayesian Mixed model fit
        #~ Predict fitted value for full dates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        model_str <- "m1 + (m2 - m7 * t) * ((1 / (1 + exp((m3 - t) / m4))) - 
            (1 / (1 + exp((m5 - t) / m6))))"
        years <- sort(unique(year(date_vec)))
        for (i in 1:numYears) { # i = 1
            date <- seq(as.Date(paste0(years[i], "-01-01")), 
                as.Date(paste0(years[i], "-12-31")), by = "day")
            bf_params <- data.table::data.table(m1 = m1[, i], m2 = m2[, i], 
                m3 = m3[, i], m4 = m4[, i], m5 = m5[, i], m6 = m6[, i], 
                m7 = m7[, i])
            phenos_idx <- data.table::data.table(midgup = numeric(nrow(bf_params)), 
                midgdown = numeric(nrow(bf_params)))

            predCI <- NULL
            for (j in 1:nrow(bf_params)) { # j = 1
                # pred based on current parameter samples
                pred <- eval(str2expression(model_str), envir = list(
                    m1 = as.numeric(bf_params[j, 1]), 
                    m2 = as.numeric(bf_params[j, 2]), 
                    m3 = as.numeric(bf_params[j, 3]), 
                    m4 = as.numeric(bf_params[j, 4]),
                    m5 = as.numeric(bf_params[j, 5]), 
                    m6 = as.numeric(bf_params[j, 6]), 
                    m7 = as.numeric(bf_params[j, 7]),
                    t = 1:length(date)
                ))
                predCI <- cbind(predCI, pred)
            }

            predCI <- t(data.table(
                apply(predCI, 1, function(x) {
                    quantile(x, c(0.025, 0.975))
                }
            )))

            pred <- data.table::data.table(
                apply(predCI, 1, function(x) quantile(x, 0.5))
            )
            cur_year_pred <- cbind(date, pred)
            bf_pred <- rbind(bf_pred, cbind(cur_year_pred, predCI))
        }

        bf_pred <- data.table::as.data.table(bf_pred)
        colnames(bf_pred) <- c("Date", "Fitted", "Fitted_lower", "Fitted_upper")
        bf_pred$Date <- as.Date(bf_pred$Date, origin = "1970-01-01")


        # ~ do the plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        bf_phenos <- data.table::data.table(bf_phenos)
        plot(bf_pred$Date, bf_pred$Fitted, cex = 0, ylim = c(-0.1, 1.2), 
            xlab = "Date", ylab = "EVI2")
        polygon(c(bf_pred$Date, rev(bf_pred$Date)), c(bf_pred$Fitted_upper, 
            rev(bf_pred$Fitted_lower)),
            col = Transparent("red", 0.2),
            border = NA
        )
        points(date_vec, vi_vec, pch = 16, col = Transparent(rep("black", 
            length(weights_vec)), weights_vec), cex = 0.5)
        lines(bf_pred$Date, bf_pred$Fitted, type = "l", ylim = c(0, 1), 
            col = "red", lwd = 2)

        pheno_names <- c("midgup", "midgdown")
        pheno_colors <- rev(viridis::viridis(9))
        for (k in 1:length(pheno_names)) { # k = 1
            pheno <- pheno_names[k]
            phn_dates <- bf_phenos[!is.na(get(pheno)), ][[pheno]]
            phn_dates <- as.Date(paste0(years, "-01-01")) + unlist(phn_dates)

            phn_val <- bf_pred[Date %in% as.Date(as.character(phn_dates)), Fitted]

            points(phn_dates, phn_val, pch = 16, col = pheno_colors[k])
            phn_dates_lower <- as.Date(paste0(years, "-01-01")) + 
                unlist(bf_phenos[!is.na(get(pheno)), ][[paste0(pheno, "_lower")]])
            phn_dates_upper <- as.Date(paste0(years, "-01-01")) + 
                unlist(bf_phenos[!is.na(get(pheno)), ][[paste0(pheno, "_upper")]])
            segments(phn_dates_lower, phn_val, phn_dates_upper, phn_val)
        }
        legend("top", ncol = 3, bty = "n", 
            lty = c(NA, 1, rep(NA, 3), 1), 
            pch = c(16, NA, 16, 15, 16, NA),
            col = c("black", "red", pheno_colors[1], 
                Transparent("red", 0.2), pheno_colors[2], "black"),
            legend = c("Observations", "Median Fit", "SOS", "95% C.I. of fit", 
                "EOS", "95% C.I. of phenometrics")
        )
        legend("bottomright", bty = "n", legend = expression(italic(
            "*Observation transparency depends on weight")), cex = 0.8)
    }

    return(list(fitted = bf_pred, phenos = bf_phenos))
}


#' Generate pheno from the predicted curve.
#' Only supports Elmore model (The double-logistic model used in BLSP).
#' 
#' @param equation The model equation.
#' @param params The Parameter list.
#' @param t Date vector.
#' @return The phenological timing.
GetPhenosIdx <- function(equation, params, t) {
    y <- eval(equation, envir = list(
        m1 = params$m1, m2 = params$m2, m3 = params$m3, m4 = params$m4,
        m5 = params$m5, m6 = params$m6, m7 = params$m7, t = t
    ))

    d1 <- D(equation, "t")
    d2 <- D(d1, "t")
    
    y1 <- eval(d1, envir = list(
        m1 = params$m1, m2 = params$m2, m3 = params$m3, m4 = params$m4,
        m5 = params$m5, m6 = params$m6, m7 = params$m7, t = t
    ))
    y2 <- eval(d2, envir = list(
        m1 = params$m1, m2 = params$m2, m3 = params$m3, m4 = params$m4,
        m5 = params$m5, m6 = params$m6, m7 = params$m7, t = t
    ))
    k <- abs(y2^2) / ((1 + y1^2)^(3 / 2))

    d <- diff(y)
    d_code <- (d > 0) + (2 * (d < 0)) # 0=no change, 1=inc, 2=dec

    # Find the MidGreenup and MidGreeendown
    midgup <- round(params$m3)
    midgdown <- round(params$m5)

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
    d_k <- c(0, 0, diff(k*1000, differences = 2))
    # local extremes
    localextr <- which(diff(sign(diff(d_k*1000))) == -2) + 1
    
    maturity <- localextr[which(localextr < peak & localextr > midgup)]
    maturity <- maturity[length(maturity)]
    
    sene <- localextr[which(localextr > peak & localextr < midgdown)][1]
    
    gup <- localextr[which(localextr < midgup)]
    gup <- gup[length(gup)]
    
    dormancy <- localextr[which(localextr > (midgdown + 7))][1]
    
    return(list(gup = gup, midgup = midgup, maturity = maturity, peak = peak, 
        sene = sene, midgdown = midgdown, dormancy = dormancy))
}


#' Fit the averaged model and get the model parameters.
#' 
#' @param date_vec the date vector, be sure to convert the vector to "Date" 
#' format or use "yyyy-mm-dd" format string.
#' @param vi_vec The vegetation index vector.
#' @param ifplot logical. Plot the model fit if TRUE.
#' @export 
FitAvgModel <- function(date_vec, vi_vec, ifplot = FALSE) {
    # check if date_vec is in Date format
    if (sum(!is.na(lubridate::parse_date_time(date_vec, orders = "ymd"))) != 
        length(date_vec)) {
        stop("There're invalid Date values in the `date_vec`! 
            Be sure to use `yyyy-mm-dd` format.")
    }
    vi_dt <- data.table::data.table(date = as.Date(date_vec), evi = vi_vec)
    vi_dt <- na.omit(vi_dt)
    vi_dt <- data.table::setorder(vi_dt, date)

    vi_dt$avg_date <- as.Date(paste0("1970", substr(vi_dt$date, 5, 10)))

    unique_dates <- unique(vi_dt$avg_date)

    merge_dt <- sapply(unique_dates, function(x) {
        # find how many records this day has
        evi <- NA
        find_idx <- which(x == vi_dt$avg_date)
        if (length(find_idx) == 1) {
            evi <- vi_dt[find_idx]$evi
        } else if (length(find_idx) > 1) { # we have multiple values for this date
            # compute the max
            evi <- max(vi_dt[avg_date == x]$evi, na.rm = TRUE)
        }
        return(list(date = x, evi = evi))
    })
    merge_dt <- data.table::as.data.table(t(merge_dt))

    cur_start_date <- as.Date("1970-01-01")
    cur_end_date <- as.Date("1970-12-31")

    y <- unlist(merge_dt$evi)
    t <- unlist(merge_dt$date)
    full_date <- seq(cur_start_date, cur_end_date, by = "day")
    full_t <- as.integer(full_date - cur_start_date + 1)

    # Assign weights
    wgt <- rep(1, length(t))
    # wgt[merge_dt$snow == TRUE] <- 0.5

    # Fit model to get the prior
    avg_fit <- tryCatch(
        {
            model_equ <- as.formula(paste("VI", "~", model_str))
            minpack.lm::nlsLM(model_equ,
                data = list(VI = y, t = as.integer(t)), weights = wgt, start = list(
                    m1 = 0.05, m2 = 1, m3 = 120, m4 = 6, m5 = 290, m6 = 8, 
                    m7 = 0.001), 
                    lower = c(0, 0.1, 1, 0, 1, 0, 0.00001), 
                    upper = c(1, 100, 185, 100, 370, 100, 0.01),
                control = list(warnOnly = TRUE)
            )
        },
        error = function(e) {
            print(paste("Average fit failed", sep = ":"))
            return(NULL)
        }
    )

    

    if (ifplot == TRUE) {
        pred <- NULL
        phenos_idx <- NULL
        phenos <- NULL
        if (!is.null(avg_fit)) {
            pred <- predict(avg_fit, newdata = list(t = full_t))
            phenos_idx <- GetPhenosIdx(str2expression(model_str),
                params = list(
                    m1 = coef(avg_fit)["m1"], 
                    m2 = coef(avg_fit)["m2"],
                    m3 = coef(avg_fit)["m3"], 
                    m4 = coef(avg_fit)["m4"], 
                    m5 = coef(avg_fit)["m5"], 
                    m6 = coef(avg_fit)["m6"], 
                    m7 = coef(avg_fit)["m7"]
                ),
                t = full_t
            )
            phenos <- full_date[unlist(phenos_idx)]
        }
        
        plot(vi_dt[, .(avg_date, evi)], pch = 16, col = "seagreen", 
            xlab = "", ylab = "")
        mtext(text = "Date", side = 1, line = 2)
        mtext(text = "EVI2", side = 2, line = 2)
    
        lines(full_date, pred, col = "orange", lwd = 2)
        points(full_date[phenos_idx$gup], pred[phenos_idx$gup], 
            pch = 21, lwd = 2, cex = 1.5, bg = "red")
        points(full_date[phenos_idx$midgup], pred[phenos_idx$midgup], 
            pch = 21, lwd = 2, cex = 1.5, bg = "red")
        points(full_date[phenos_idx$maturity], pred[phenos_idx$maturity], 
            pch = 21, lwd = 2, cex = 1.5, bg = "red")
        points(full_date[phenos_idx$peak], pred[phenos_idx$peak], 
            pch = 21, lwd = 2, cex = 1.5, bg = "red")
        points(full_date[phenos_idx$sene], pred[phenos_idx$sene], 
            pch = 21, lwd = 2, cex = 1.5, bg = "red")
        points(full_date[phenos_idx$midgdown], pred[phenos_idx$midgdown], 
            pch = 21, lwd = 2, cex = 1.5, bg = "red")
        points(full_date[phenos_idx$dormancy], pred[phenos_idx$dormancy], 
            pch = 21, lwd = 2, cex = 1.5, bg = "red")
    }

    return(avg_fit)
}
