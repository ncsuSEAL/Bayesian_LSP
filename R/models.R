# ******************************************************************************
# This file defines functions to fit the VI curves.
# Date: 2024-04-15
# ******************************************************************************



#' Get the user selected model structure and parameters
#' 
#' @param model_name The name of the model. For now only support `dblog7` and
#' `dblog6` corresponding to the 7- and 6-parameter double-logistic functions.
#'
#' @return A list containing model structure and parameters.
GetModel <- function(model_name) {
    switch(tolower(model_name), 
        "dblog7" = dblog7,
        "dblog6" = dblog6,
        dblog7
    )
}


# Double logistic function w/ the "greendown" parameter
dblog7 <- list(
    model_name = "dblog7",

    model_str = "m1 + (m2 - m7 * t) * ((1 / (1 + exp((m3 - t) / m4))) -
        (1 / (1 + exp((m5 - t) / m6))))",
    
    jags_str = "model {
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
    }", 

    out_param = paste0("m", 1:7),
    init_param = c("M1", paste0("m", 2:6)),
    init_val = c(0.05, 1, 120, 8, 290, 8, 0.001),
    init_lwr = c(0, 0.1, 1, 0, 1, 0, 0.00001),
    init_upr = c(1, 100, 185, 100, 370, 100, 0.01)
)


# Double logistic function w/o the "greendown" parameter
dblog6 <- list(
    model_name = "dblog6",

    model_str = "m1 + m2 * ((1 / (1 + exp((m3 - t) / m4))) -
        (1 / (1 + exp((m5 - t) / m6))))",
    
    jags_str = "model {
        # Likelihood
        for (i in 1:n) {
            Y[i] ~ dnorm(mu[i], tau_y)
            mu[i] <- weights[i] * (m1[yr[i]] + m2[yr[i]] *
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
        }

        mu_m1 ~ dunif(0, 0.3)
        mu_m2 ~ dunif(0.5, 2)
        mu_m3 ~ dunif(0, 185)
        mu_m4 ~ dunif(1, 15)
        mu_m5 ~ dunif(185, 366)
        mu_m6 ~ dunif(1, 15)

        for (k in 1:6) {
            tau[k] ~ dgamma(0.1, 0.1)
        }
        tau_y ~ dgamma(0.1, 0.1)
    }",

    out_param = paste0("m", 1:6),
    init_param = c("M1", paste0("m", 2:6)),
    init_val = c(0.05, 1, 120, 8, 290, 8),
    lwr = c(0, 0.1, 1, 0, 1, 0),
    upr = c(1, 100, 185, 100, 370, 100)
)
