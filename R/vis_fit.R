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
