#************************************************************************************
# Description: Test model fitting related functions
# Date: 2022-04-08
#************************************************************************************

# Test `FitAvgModel`
test_that("Fitting the average model works", {
    data(landsatEVI2)
    avg_mod <- FitAvgModel(landsatEVI2$date, landsatEVI2$evi2)
    expect_type(avg_mod, "list")

    est_coef <- coef(avg_mod)
    
    # to allow some tolerance
    est_coef <- round(est_coef, 3)
    ept_coef <- c(
       m1 = 1.74e-01,
       m2 = 7.26e-01,
       m3 = 1.35672e+02,
       m4 = 8.403e+00,
       m5 = 2.93534e+02,
       m6 = 1.1898e+01,
       m7 = 1.0e-03
    )
    expect_equal(est_coef, ept_coef)
})

# Test `GetPhenosIdx`
test_that("Phenometrics retrieval works", {
    skip("Skip GetPhenosIdx")
})

# Test `FitBLSP`
test_that("BLSP model works", {
    skip_on_cran()

    data(landsatEVI2)
    model_init <- FitAvgModel(landsatEVI2$date, landsatEVI2$evi2)
    bf_fit <- FitBLSP(
        date_vec = landsatEVI2$date,
        vi_vec = landsatEVI2$evi2,
        weights_vec = ifelse(landsatEVI2$snow == TRUE, 0.1, 1),
        initValues = model_init,
        ifplot = FALSE, 
        verbose = FALSE
    )

    expect_type(bf_fit, "list")
    expect_equal(bf_fit$fitted, NULL)

    # The output of FitBLSP function is a bit weird, should revise it and chagne 
    #   here later.
    est_lsp <- data.frame(bf_fit$phenos)
    est_lsp <- lapply(est_lsp, function(x) { unlist(x) })
    est_lsp <- data.frame(est_lsp)
    est_lsp <- round(est_lsp, 0)

    # Check the values
    expect_equal(est_lsp$Year, 1984:2019)

    # Allow 5 days deviance
    CheckDeviance <- function(val, ept) {
        if (all(abs(val - ept) < 5)) {
            return(TRUE)
        } else {
           return(FALSE)
        }
    }
    
    expect_true(CheckDeviance(est_lsp$midgup_lower, c(
        132, 133, 133, 137, 136, 134, 133, 131, 137, 133, 138, 136, 133, 137,
        132, 135, 135, 131, 135, 134, 132, 137, 136, 135, 136, 132, 126, 138,
        133, 131, 137, 134, 138, 134, 135, 137
    )))
    expect_true(CheckDeviance(est_lsp$midgup, c(
        139, 139, 139, 141, 140, 140, 140, 135, 141, 138, 142, 141, 140, 144, 
        137, 140, 139, 136, 140, 140, 136, 142, 141, 139, 140, 136, 132, 143, 
        139, 136, 141, 139, 142, 139, 139, 141
    )))
    expect_true(CheckDeviance(est_lsp$midgdown_lower, c(
        279, 276, 277, 280, 276, 278, 279, 277, 283, 280, 279, 279, 281, 279,
        278, 279, 279, 279, 281, 280, 281, 282, 279, 282, 282, 279, 281, 277,
        275, 275, 282, 284, 283, 283, 282, 281
    )))
    expect_true(CheckDeviance(est_lsp$midgdown, c(
        285, 284, 285, 285, 284, 285, 285, 284, 287, 285, 285, 285, 286, 285,
        285, 285, 285, 285, 286, 285, 286, 286, 285, 286, 286, 285, 286, 285,
        284, 283, 286, 288, 287, 287, 286, 286
    )))
    expect_true(CheckDeviance(est_lsp$midgdown_upper, c(
        290, 289, 290, 291, 289, 289, 290, 288, 295, 290, 290, 291, 292, 291,
        291, 290, 291, 290, 295, 291, 292, 292, 289, 294, 292, 289, 291, 290,
        288, 288, 293, 298, 298, 295, 295, 293
    )))


})

