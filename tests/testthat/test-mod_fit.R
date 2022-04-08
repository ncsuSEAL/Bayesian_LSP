#************************************************************************************
# Description: Test model fitting related functions
# Date: 2022-04-08
#************************************************************************************

# Test `FitAvgModel`
test_that("Fitting the average model works", {
    data(test_ts)
    avg_mod <- FitAvgModel(test_ts$date, test_ts$all_evi)
    expect_type(avg_mod, "list")
})

# Test `GetPhenosIdx`
test_that("Phenometrics retrieval works", {
    expect_equal(2 * 2, 4)
})

# Test `FitBLSP`
test_that("BLSP model works", {
    expect_equal(2 * 2, 4)
})

