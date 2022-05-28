#************************************************************************************
# Description: Needed library and pre-defined functions
# Date: 2020-11-22
#************************************************************************************

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
Transparent <- function(orig.col, alpha = 1, maxColorValue = 255) {
    n.cols <- length(orig.col)
    orig.col <- col2rgb(orig.col)
    final.col <- rep(NA, n.cols)
    for (i in 1:n.cols) {
        final.col[i] <- rgb(orig.col[1, i], orig.col[2, i], orig.col[3, i],
            alpha = alpha[i] * 255,
            maxColorValue = maxColorValue
        )
    }
    return(final.col)
}



