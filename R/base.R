#*******************************************************************************
# Description: Needed library and pre-defined functions
# Date: 2020-11-22
#*******************************************************************************

.datatable.aware <- TRUE

usethis::use_pipe(export = TRUE)

# The double-logistic function equation
model_str <- "m1 + (m2 - m7 * t) * ((1 / (1 + exp((m3 - t) / m4))) -
            (1 / (1 + exp((m5 - t) / m6))))"


#' Make a standard color transparent.
#' This function is borrowed from 'yarrr' package, but I changed the trans.val
#' to use alpha value directly.
#' @param orig.col: the original color, can be a color name, a hexadecimal code,
#' or a rgb vector.
#' @param alpha: define the transparent level.
#' @param maxColorValue: used to convert the color to rgb format before making
#' it transparent.
#' @return color code.
#' 
#' @noRd
Transparent <- function(orig.col, alpha = 1, maxColorValue = 255) {
    n.cols <- length(orig.col)
    orig.col <- grDevices::col2rgb(orig.col)
    final.col <- rep(NA, n.cols)
    for (i in 1:n.cols) {
        final.col[i] <- grDevices::rgb(
            orig.col[1, i], orig.col[2, i], orig.col[3, i],
            alpha = alpha[i] * 255,
            maxColorValue = maxColorValue
        )
    }
    return(final.col)
}


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