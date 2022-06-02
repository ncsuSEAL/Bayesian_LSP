#************************************************************************************
# Description: Visualize model fit results.
# Date: 2022-05-28
#************************************************************************************


#' Visualize the average model fit result. It will show all points as well as the 
#' averaged model fit curve.
#' 
#' @param date_vec the date vector, be sure to convert the vector to "Date"
#' format or use "yyyy-mm-dd" format string.
#' @param vi_vec The vegetation index vector.
#' @param avg_fit The model fit object returned by `FitAvgModel()`.
#' @return A plot showing the average model fit result.
#' @export
#' 
#' @examples 
#' \dontrun{
#' avg_dt <- FormatAvgData(landsatEVI2$date, landsatEVI2$evi2)
#' avg_fit <- FitAvgModel(landsatEVI2$date, landsatEVI2$evi2)
#' PlotAvg(landsatEVI2$date, landsatEVI2$evi2, avg_fit)
#' }
PlotAvg <- function(date_vec, vi_vec, avg_fit) {
    # Format data
    avg_dt <- FormatAvgData(date_vec, vi_vec)

    # Fake year information for the averaged year
    cur_start_date <- as.Date("1970-01-01")
    cur_end_date <- as.Date("1970-12-31")
    full_date <- seq(cur_start_date, cur_end_date, by = "day")
    full_t <- as.integer(full_date - cur_start_date + 1)

    # Predict from the avg model fit
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

    # Plot the figure
    plot(avg_dt[, .(date, evi2)],
        pch = 16, col = "seagreen",
        xlab = "", ylab = ""
    )
    mtext(text = "DOY", side = 1, line = 2)
    mtext(text = "EVI2", side = 2, line = 2)

    lines(full_date, pred, col = "orange", lwd = 2)
    points(full_date[phenos_idx$midgup], pred[phenos_idx$midgup],
        pch = 21, lwd = 2, cex = 1.5, bg = "red"
    )
    points(full_date[phenos_idx$midgdown], pred[phenos_idx$midgdown],
        pch = 21, lwd = 2, cex = 1.5, bg = "red"
    )
}



#' Plot BLSP model fitting result. 
#' 
#' @param blsp_fit The object of `BlspFit` class returned by `FitBLSP()` function.
#' @param if_return_fit Logic. Determine whether return the fitted values. Default
#' is `FALSE`.
#' @return A plot showing the BLSP model fitting result. If `if_return_fit` is true,
#' the model fitted time series as well as the 95% credible interval will also be
#' returned.
#' @export
#' 
#' @examples 
#' \dontrun{
#' blsp_fit <- FitBLSP(landsatEVI2$date, landsatEVI2$evi2, verbose = TRUE)
#' fitted_dt <- PlotBLSP(blsp_fit, if_return_fit = TRUE)
#' }
#' @import data.table
PlotBLSP <- function(blsp_fit, if_return_fit = FALSE) {
    
    if (class(blsp_fit) != "BlspFit") {
        stop("The input should be the output object of `FitBLSP()` function!")
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

    #~ Predict fitted value for full dates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bf_pred <- NULL
    years <- sort(unique(lubridate::year(date_vec)))
    for (i in 1:numYears) {
        date <- seq(as.Date(paste0(years[i], "-01-01")), 
            as.Date(paste0(years[i], "-12-31")), by = "day")
        bf_params <- data.table::data.table(
            m1 = blsp_fit$params$m1[, i], 
            m2 = blsp_fit$params$m2[, i], 
            m3 = blsp_fit$params$m3[, i], 
            m4 = blsp_fit$params$m4[, i], 
            m5 = blsp_fit$params$m5[, i], 
            m6 = blsp_fit$params$m6[, i], 
            m7 = blsp_fit$params$m7[, i]
        )
        phenos_idx <- data.table::data.table(
            midgup = numeric(nrow(bf_params)), 
            midgdown = numeric(nrow(bf_params))
        )

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

        predCI <- t(data.table::data.table(
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
    # Make it a data.table
    bf_pred <- data.table::as.data.table(bf_pred)
    colnames(bf_pred) <- c("Date", "Fitted", "Fitted_lower", "Fitted_upper")
    bf_pred$Date <- as.Date(bf_pred$Date, origin = "1970-01-01")


    # ~ Do the plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    plot(bf_pred$Date, bf_pred$Fitted, 
        cex = 0, ylim = c(-0.1, 1), 
        xlab = "Date", ylab = "EVI2"
    )
    polygon(c(bf_pred$Date, rev(bf_pred$Date)), 
        c(bf_pred$Fitted_upper, rev(bf_pred$Fitted_lower)),
        col = Transparent("red", 0.2),
        border = NA
    )
    points(date_vec, vi_vec, 
        pch = 16, 
        col = Transparent(rep("black", length(weights_vec)), weights_vec), 
        cex = 0.5
    )
    lines(bf_pred$Date, bf_pred$Fitted, 
        type = "l", ylim = c(0, 1), 
        col = "red", lwd = 2
    )

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
    legend(grconvertX(0.5, "ndc"), grconvertY(0.95, "ndc"), 
        xjust = 0.5,
        ncol = 3, bty = "n", 
        lty = c(NA, 1, rep(NA, 3), 1), 
        pch = c(16, NA, 16, 15, 16, NA),
        col = c("black", "red", pheno_colors[1], 
            Transparent("red", 0.2), pheno_colors[2], "black"),
        legend = c("Observations", "Median Fit", "SOS", "95% C.I. of fit", 
            "EOS", "95% C.I. of phenometrics"),
        xpd = NA
    )
    legend("bottomright", bty = "n", 
        legend = expression(italic("*Observation transparency depends on weight")), 
        cex = 0.8
    )

    if (if_return_fit == TRUE) {
        return(bf_pred)
    }
}

