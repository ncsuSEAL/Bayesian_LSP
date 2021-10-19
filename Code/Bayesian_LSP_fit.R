
#! Set working directory
setwd("../Bayesian_LSP/Code")

source("base.R")


# First, need to load the Landsat time series
landsat <- readRDS("test_ts.Rds")

# We can plot the time series
{
    plot(landsat[, .(date, evi_L5)], pch = 16, col = "seagreen", xlab = "", ylab = "")
    mtext(text = "Date", side = 1, line = 2)
    mtext(text = "EVI2", side = 2, line = 2)
    # Landsat 5
    lines(landsat[!is.na(evi_L5), .(date, evi_L5)], lwd = 2, col = Transparent("seagreen", 0.3))
    # Landsat 7
    points(landsat[!is.na(evi_L7), .(date, evi_L7)], pch = 16, col = "purple")
    lines(landsat[!is.na(evi_L7), .(date, evi_L7)], lwd = 2, col = Transparent("purple", 0.3))
    # Landsat 8
    points(landsat[!is.na(evi_L8), .(date, evi_L8)], pch = 16, col = "red")
    lines(landsat[!is.na(evi_L8), .(date, evi_L8)], lwd = 2, col = Transparent("red", 0.3))
    # snow
    points(landsat[snow == TRUE, .(date, all_evi)], lwd = 2)

    legend("bottom",
        pch = c(16, 16, 16, 1), col = c("seagreen", "purple", "red", "black"), bg = "white",
        legend = c("Landsat 5", "Landsat 7", "Landsat 8", "snow"), ncol = 4
    )
}

# From the plot, we can see the sparsity of the Landsat time series. 
# There are few observations for some years, especially before the launch of Landsat 7 and 8.

# Before fitting the Bayesian model, we'd like to get some estimates of the model parameters as the initial 
#   values for the Monte Calro Markov Chain (MCMC) sampling. We pool all years of observations to fit an 
#   averaged model to get initial values for model parameters. 
# This is not a requirement and the initial values would have little impact on the final parameter values, 
#   but good initial values can reduce the computational time needed for MCMC. 

model_init <- FitAvgModel(landsat, model_str, ifplot = TRUE)
if (is.null(model_init)) model_init <- NULL


# Now, we can fit the Bayesian model and plot it out. 
# The code line below would take ~5 min to run (1 min for the model, 4 min for the plot).
bf_fit <- FitBayesianModel(model_str, landsat, initValues = model_init, ifplot = FALSE)

# The retrieved phenometrics and their 95% credible intervals are:
colnames(bf_fit$phenos) <- c("Id", "Year", "SOS_lower", "SOS", "SOS_upper", "EOS_lower", "EOS", "EOS_upper")
apply(bf_fit$phenos, 2, as.integer)

