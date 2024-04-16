# ******************************************************************************
# Helper functions that are not exposed to users.
# Date: 2024-04-16
# ******************************************************************************


#' Format input date and VI vectors to the structure needed for fitting averaged
#' phenology models such as Fisher et al 2006, Elmore et al 2012.
#'
#' @param date_vec the date vector, be sure to convert the vector to "Date"
#' format or use "yyyy-mm-dd" format string.
#' @param vi_vec The vegetation index vector.
#' @return A list that contains formated data.
#' @import data.table
#'
#' @noRd
FormatAvgData <- function(date_vec, vi_vec) {
    # Check if date_vec is in Date format
    if (sum(!is.na(lubridate::parse_date_time(date_vec, orders = "ymd"))) !=
        length(date_vec)) {
        stop("There're invalid Date values in the `date_vec`!
            Be sure to use `yyyy-mm-dd` format.")
    }

    # Make it a data table
    vi_dt <- data.table::data.table(
        date = as.Date(date_vec),
        evi2 = vi_vec,
        avg_date = ""
    )
    vi_dt[, avg_date := as.Date(paste0("1970", substr(vi_dt$date, 5, 10)))]
    vi_dt <- stats::na.omit(vi_dt)
    data.table::setorder(vi_dt, date)

    # Find unique dates in the averaged year
    unique_dates <- unique(vi_dt$avg_date)

    # Deal with multiple observations on the same date in the averaged year.
    # When that happens, we choose the one whose EVI2 value is the highest.
    merge_dt <- sapply(unique_dates, function(x) {
        # find how many records this day has
        evi2 <- NA
        find_idx <- which(x == vi_dt$avg_date)
        if (length(find_idx) == 1) {
            evi2 <- vi_dt[find_idx]$evi2
        } else if (length(find_idx) > 1) { # we have multiple values for this date
            # compute the max
            evi2 <- max(vi_dt[avg_date == x]$evi2, na.rm = TRUE)
        }
        return(list(date = x, evi2 = evi2))
    })
    merge_dt <- data.table::as.data.table(t(merge_dt))

    return(merge_dt)
}


#' Format empty BLSP output to return
#'
#' @param years The years vector
#' @param date_vec The date vector, be sure to convert the vector to "Date"
#' format or use "yyyy-mm-dd" format string.
#' @param vi_vec The vegetation index vector.
#' @param weights_vec A numeric vector of same length as vi_vec specifying the
#' weights for the supplied observations. Must be between 0 and 1, inclusive.
#' @param cred_int_level A scalar value from 0 to 1 (exclusive) that specifies
#' the level for equal-tailed credible intervals of the estimated phenometrics.
#' The default level is 0.9, generating `90%` credible intervals. The end
#' points of these intervals define the upper and lower bounds for the estimated
#' phenometrics.
#' @param mod The model object returned by the `GetModel()` function.
#' @param method Method used for the phenometric output. For now, only support
#' "default" and "threshold". The "default" method will return `midgup` and
#' `midgdown`, while the "threshold" method will return 7 phenometrics including
#' Greenup, MidGreenup, Maturity, Peak, Senescence, MidGreendown, and Dormancy. 
#'
#' @return An empty BLSP class object.
#'
#' @noRd
EmptyBlspOutput <- function(
    years, date_vec, vi_vec, weights_vec,
    cred_int_level, mod, method = "default"
) {
    if (method == "default") {
        bf_phenos <- data.table::data.table(
            Year = years,
            midgup_lwr = NA, midgup = NA, midgup_upr = NA,
            midgdown_lwr = NA, midgdown = NA, midgdown_upr = NA
        )
    } else if (method == "threshold") {
        bf_phenos <- data.table::data.table(
            Year = years,
            Greenup = NA, Greenup_lwr = NA, Greenup_upr = NA,
            MidGreenup = NA, MidGreenup_lwr = NA, MidGreenup_upr = NA,
            Maturity = NA, Maturity_lwr = NA, Maturity_upr = NA,
            Peak = NA, Peak_lwr = NA, Peak_upr = NA,
            Senescence = NA, Senescence_lwr = NA, Senescence_upr = NA,
            MidGreendown = NA, MidGreendown_lwr = NA, MidGreendown_upr = NA,
            Dormancy = NA, Dormancy_lwr = NA, Dormancy_upr = NA
        )
    }

    bf_phenos[,
        colnames(bf_phenos) := lapply(.SD, as.numeric),
        .SDcols = colnames(bf_phenos)
    ]

    blsp_fit <- list(
        phenos = bf_phenos,
        model = list(
            model_str = mod$model_str,
            model_param = mod$out_param
        ),
        params = NULL,
        data = data.table::data.table(
            date = date_vec,
            vi = vi_vec,
            weights = weights_vec
        ),
        cred_int_level = cred_int_level
    )
    class(blsp_fit) <- "BlspFit"

    return(blsp_fit)
}


#' Check if user inputs are reasonable
#'
#' @param date_vec The date vector, be sure to convert the vector to "Date"
#' format or use "yyyy-mm-dd" format string.
#' @param vi_vec The vegetation index vector.
#' @param weights_vec A numeric vector of same length as vi_vec specifying the
#' weights for the supplied observations. Must be between 0 and 1, inclusive.
#' @param model The model string.
#' @param cred_int_level A scalar value from 0 to 1 (exclusive) that specifies
#' the level for equal-tailed credible intervals of the estimated phenometrics.
#' The default level is 0.9, generating `90%` credible intervals. The end
#' points of these intervals define the upper and lower bounds for the estimated
#' phenometrics.
#' @param init_values Initial values for MCMC sampling. By default, it is
#' assgined `NULL`. It could also be an object returned from the `FitAvgModel()`
#' function that fits an averaged model or a numeric vector provided by the
#' user. 
#'
#' @import data.table
SanityCheck <- function(
    date_vec, vi_vec, weights_vec,
    model, cred_int_level, init_values
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

    # Check the jags model string
    if (!model %in% c("dblog7", "dblog6")) {
        warning("The specified model does not exist, dblog7 will be used.")
        model <- "dblog7"
    }

    # Check init value
    if (!is.null(init_values)) {
        # Check if the length matches w/ model
        model_init_vals <- GetModel(model)$init_val
        if (class(init_values) == "nls") {
            if (length(stats::coef(init_values)) != length(model_init_vals)) {
                stop("Init values does not match w/ the model!")
            }
        } else if (class(init_values) == "numeric") {
            if (length(init_values) != length(model_init_vals)) {
                stop("Init values does not match w/ the model!")
            }
        }
    }

    if (cred_int_level >= 1 || cred_int_level <= 0) {
        warning("`cred_int_level` is not in (0, 1), 0.9 will be used.")
        cred_int_level <- 0.9
    }

    return(list(
        model = model, cred_int_level = cred_int_level
    ))
}


#' Calculate phenometrics based on the model parameters
#' 
#' @param p_li A list containing the model parameters.
#' @param mod The model object returned by the `GetModel()` function.
#' @param years The years vector
#' @param numYears Number of years.
#' @param date_vec The date vector, be sure to convert the vector to "Date"
#' format or use "yyyy-mm-dd" format string.
#' @param vi_vec The vegetation index vector.
#' @param weights_vec A numeric vector of same length as vi_vec specifying the
#' weights for the supplied observations. Must be between 0 and 1, inclusive.
#' @param cred_int_level A scalar value from 0 to 1 (exclusive) that specifies
#' the level for equal-tailed credible intervals of the estimated phenometrics.
#' The default level is 0.9, generating `90%` credible intervals. The end
#' points of these intervals define the upper and lower bounds for the estimated
#' phenometrics.
#' 
#' @return An "BlspFit" object filled with retrieved phenology and parameters.
#' 
#' @import data.table
#' @noRd
CalPhenoParam <- function(
    p_li, mod,
    years, numYears,
    date_vec, vi_vec, weights_vec,
    cred_int_level
) {
    blsp_fit <- EmptyBlspOutput(
        years,
        date_vec,
        vi_vec,
        weights_vec,
        cred_int_level,
        mod,
        method = "default"
    )
    p_qt <- lapply(p_li, function(li) {
        apply(li, 2, stats::quantile, c(0.05, 0.5, 0.95))
    })

    phenos <- lapply(1:numYears, function(i) {
        # suppress some amplitude-too-low year
        amp <- p_qt[[2]][2, i]
        if (mod$model_name == "dblog7" & amp > 0.4) {
            c(
                years[i],
                p_qt[[3]][1, i], p_qt[[3]][2, i], p_qt[[3]][3, i],
                p_qt[[5]][1, i], p_qt[[5]][2, i], p_qt[[5]][3, i]
            )
        } else if (mod$model_name == "dblog6" & amp > 0.1) {
            c(
                years[i],
                p_qt[[3]][1, i], p_qt[[3]][2, i], p_qt[[3]][3, i],
                p_qt[[5]][1, i], p_qt[[5]][2, i], p_qt[[5]][3, i]
            )
        } else {
            c(years[i], rep(NA, 6))
        }
    }) %>%
        do.call(rbind, .) %>%
        data.table::as.data.table() %>%
        setnames(c(
            "Year",
            "midgup_lwr", "midgup", "midgup_upr",
            "midgdown_lwr", "midgdown", "midgdown_upr"
        ))

    blsp_fit$phenos <- phenos
    blsp_fit$params <- p_li

    return(blsp_fit)
}


#' Calculate phenometrics using the threshold-based method
#'
#' @param p_li A list containing the model parameters.
#' @param mod The model object returned by the `GetModel()` function.
#' @param years The years vector
#' @param numYears Number of years.
#' @param date_vec The date vector, be sure to convert the vector to "Date"
#' format or use "yyyy-mm-dd" format string.
#' @param vi_vec The vegetation index vector.
#' @param weights_vec A numeric vector of same length as vi_vec specifying the
#' weights for the supplied observations. Must be between 0 and 1, inclusive.
#' @param cred_int_level A scalar value from 0 to 1 (exclusive) that specifies
#' the level for equal-tailed credible intervals of the estimated phenometrics.
#' The default level is 0.9, generating `90%` credible intervals. The end
#' points of these intervals define the upper and lower bounds for the estimated
#' phenometrics.
#' 
#' #' @return An "BlspFit" object filled with retrieved phenology and parameters.
#' 
#' @import data.table
#' @noRd
CalPhenoThresh <- function(
    p_li, mod,
    years, numYears,
    date_vec, vi_vec, weights_vec,
    cred_int_level
) {
    # Format MCD12Q2-like phenometrics
    blsp_fit <- EmptyBlspOutput(
        years,
        date_vec,
        vi_vec,
        weights_vec,
        cred_int_level,
        mod,
        method = "threshold"
    )
    blsp_fit$params <- p_li

    bf_pred <- BLSPFitted(blsp_fit)

    yrs <- names(bf_pred)
    phenos_dt <- NULL
    for (yr in yrs) {
        pred <- bf_pred[eval(yr)][[1]]
        phenos <- lapply(1:ncol(pred), function(i) {
            p <- pred[, i]
            peakdate <- which.max(p)
            peak <- max(p)

            # Spring amp
            spring_min <- min(p[1:peakdate])
            spring_amp <- peak - spring_min
            # Autumn amp
            autumn_min <- min(p[peakdate:length(p)])
            autumn_amp <- peak - autumn_min

            if (spring_amp < 0.2 | autumn_amp < 0.2) {
                return(rep(NA, 7))
            }

            gup <- which(
                p[1:peakdate] > (spring_amp * 0.15 + min(p[1:peakdate]))
            )[1]
            midgup <- which(
                p[1:peakdate] > (spring_amp * 0.5 + min(p[1:peakdate]))
            )[1]
            mat <- which(
                p[1:peakdate] > (spring_amp * 0.90 + min(p[1:peakdate]))
            )[1]

            sens <- which(
                p[peakdate:length(p)] < autumn_amp * 0.90 + autumn_min
            )[1] + peakdate
            midgdown <- which(
                p[peakdate:length(p)] < autumn_amp * 0.5 + autumn_min
            )[1] + peakdate
            dorm <- which(
                p[peakdate:length(p)] < autumn_amp * 0.15 + autumn_min
            )[1] + peakdate

            return(c(gup, midgup, mat, peakdate, sens, midgdown, dorm))
        })
        phenos <- do.call(cbind, phenos)

        # Get CIs
        phenos <- apply(phenos, 1, quantile,
            c(0.025, 0.5, 0.975),
            na.rm = TRUE
        )
        colnames(phenos) <- c(
            "gup", "midgup", "mat", "peakdate", "sens", "midgdown", "dorm"
        )

        therow <- data.table(
            Year = yr,
            Greenup = round(phenos[2, "gup"]),
            Greenup_lwr = round(phenos[1, "gup"]),
            Greenup_upr = round(phenos[3, "gup"]),
            MidGreenup = round(phenos[2, "midgup"]),
            MidGreenup_lwr = round(phenos[1, "midgup"]),
            MidGreenup_upr = round(phenos[3, "midgup"]),
            Maturity = round(phenos[2, "mat"]),
            Maturity_lwr = round(phenos[1, "mat"]),
            Maturity_upr = round(phenos[3, "mat"]),
            Peak = round(phenos[2, "peakdate"]),
            Peak_lwr = round(phenos[1, "peakdate"]),
            Peak_upr = round(phenos[3, "peakdate"]),
            Senescence = round(phenos[2, "sens"]),
            Senescence_lwr = round(phenos[1, "sens"]),
            Senescence_upr = round(phenos[3, "sens"]),
            MidGreendown = round(phenos[2, "midgdown"]),
            MidGreendown_lwr = round(phenos[1, "midgdown"]),
            MidGreendown_upr = round(phenos[3, "midgdown"]),
            Dormancy = round(phenos[2, "dorm"]),
            Dormancy_lwr = round(phenos[1, "dorm"]),
            Dormancy_upr = round(phenos[3, "dorm"])
        )

        phenos_dt <- rbind(phenos_dt, therow)
    }

    blsp_fit$phenos <- phenos_dt
    blsp_fit$params <- p_li

    return(blsp_fit)
}
