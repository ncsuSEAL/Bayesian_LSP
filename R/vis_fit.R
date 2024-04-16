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
#' @param model A string indicating the model name. For now, only support
#' "dblog7" and "dblog6" for the 7- and 6-parameter double-logistic functions.
#' 
#' @return A plot showing the average model fit result.
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' avg_dt <- FormatAvgData(landsatEVI2$date, landsatEVI2$evi2)
#' avg_fit <- FitAvgModel(landsatEVI2$date, landsatEVI2$evi2)
#' PlotAvg(landsatEVI2$date, landsatEVI2$evi2, avg_fit)
#' }
#' @import data.table
PlotAvg <- function(date_vec, vi_vec, avg_fit, model = "dblog7") {
    if (!model %in% c("dblog7", "dblog6")) {
        warning("The specified model does not exist, dblog7 will be used.")
        model <- "dblog7"
    }
    mod <- GetModel(model)

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
        pred <- stats::predict(avg_fit, newdata = list(t = full_t))
        phenos_idx <- GetPhenosIdx(str2expression(mod$model_str),
            params = list(
                m1 = stats::coef(avg_fit)["m1"],
                m2 = stats::coef(avg_fit)["m2"],
                m3 = stats::coef(avg_fit)["m3"],
                m4 = stats::coef(avg_fit)["m4"],
                m5 = stats::coef(avg_fit)["m5"],
                m6 = stats::coef(avg_fit)["m6"],
                m7 = stats::coef(avg_fit)["m7"]
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
    graphics::mtext(text = "DOY", side = 1, line = 2)
    graphics::mtext(text = "EVI2", side = 2, line = 2)

    graphics::lines(full_date, pred, col = "orange", lwd = 2)
    
    for (i in 1:length(phenos_idx)) {
        graphics::points(full_date[phenos_idx[i]], pred[phenos_idx[i]],
            pch = 21, lwd = 2, cex = 1.5, bg = "red"
        )
    }
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
    weights_vec <- blsp_fit$data$weights
    if (is.null(weights_vec)) {
        weights_vec <- rep(1, length(vi_vec))
    }
    bf_phenos <- blsp_fit$phenos
    yr <- lubridate::year(date_vec) - lubridate::year(date_vec)[1] + 1
    numYears <- length(unique(yr))
    disp_cred_int_level <- round(blsp_fit$cred_int_level*100)
    
    #~ Predict fitted value for full dates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bf_pred <- BLSPFitted(blsp_fit, asCI = TRUE)
    years <- sort(unique(lubridate::year(date_vec)))


    # ~ Do the plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    plot(bf_pred$Date, bf_pred$Fitted, 
        cex = 0, ylim = c(-0.1, 1), 
        xlab = "Date", ylab = "EVI2",
        bty = "L"
    )
    graphics::polygon(c(bf_pred$Date, rev(bf_pred$Date)), 
        c(bf_pred$Fitted_upper, rev(bf_pred$Fitted_lower)),
        col = adjustcolor("red", 0.2),
        border = NA
    )
    graphics::points(date_vec, vi_vec, 
        pch = 16, 
        col = sapply(weights_vec, function(i) {
            adjustcolor("black", weights_vec[i])
        }),
        cex = 0.5
    )
    graphics::lines(bf_pred$Date, bf_pred$Fitted, 
        type = "l", ylim = c(0, 1), 
        col = "red", lwd = 2
    )

    pheno_names <- colnames(blsp_fit$phenos)[-1]
    pheno_names <- pheno_names[-grep("_", pheno_names)]
    pheno_colors <- rev(viridis::viridis(9))
    for (k in 1:length(pheno_names)) {
        pheno <- pheno_names[k]
        phn_dates <- bf_phenos[!is.na(get(pheno)), ][[pheno]]
        phn_dates <- as.Date(paste0(years, "-01-01")) + unlist(phn_dates)

        phn_val <- bf_pred[Date %in% as.Date(as.character(phn_dates)), Fitted]

        graphics::points(phn_dates, phn_val, pch = 16, col = pheno_colors[k])
        phn_dates_lower <- as.Date(paste0(years, "-01-01")) + 
            unlist(bf_phenos[!is.na(get(pheno)), ][[paste0(pheno, "_lwr")]])
        phn_dates_upper <- as.Date(paste0(years, "-01-01")) + 
            unlist(bf_phenos[!is.na(get(pheno)), ][[paste0(pheno, "_upr")]])
        graphics::segments(phn_dates_lower, phn_val, phn_dates_upper, phn_val)
    }
    graphics::legend(
        graphics::grconvertX(0.5, "ndc"), graphics::grconvertY(0.95, "ndc"), 
        xjust = 0.5, bty = "n", 
        ncol = ifelse(length(pheno_names) == 2, 3, 4), 
        legend = c("Observations", "Median Fit", 
            paste0(disp_cred_int_level, "% C.I. of fit"), 
            paste0(disp_cred_int_level, "% C.I. of phenometrics"),
            pheno_names
        ),
        lty = c(NA, 1, NA, 1, 1, rep(NA, length(pheno_names))), 
        pch = c(16, NA, 15, NA, 16, rep(16, length(pheno_names))),
        col = c("black", "red", 
            adjustcolor("red", 0.2),
            "black",
            pheno_colors[1:length(pheno_names)]
        ),
        xpd = NA
    )
    graphics::legend("bottomright", bty = "n", 
        legend = expression(italic(
            "*Observation transparency depends on weight"
        )), 
        cex = 0.8
    )

    if (if_return_fit == TRUE) {
        return(bf_pred)
    }
}

