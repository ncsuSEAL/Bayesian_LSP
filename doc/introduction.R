## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(blsp)

## ----set-up 1-----------------------------------------------------------------
data(landsatEVI2)

## ----code 1-------------------------------------------------------------------
sub_dt <- landsatEVI2[lubridate::year(date) %in% 1984:1994, ]
results <- FitBLSP(
  date_vec = sub_dt$date, 
  vi_vec = sub_dt$evi2, 
  verbose = TRUE
)

## ----code 2-------------------------------------------------------------------
str(results)

## ----code 3-------------------------------------------------------------------
fitted_dt <- PlotBLSP(results, if_return_fit = TRUE)

## ----code 4-------------------------------------------------------------------
head(fitted_dt)

## ----code 5-------------------------------------------------------------------
sub_dt <- landsatEVI2[lubridate::year(date) %in% 1984:1994, 
  .(date, evi2, weights = ifelse(snow, 0.1, 1), snow)
]
head(sub_dt)

## ----code 6-------------------------------------------------------------------
results <- FitBLSP(
  date_vec = sub_dt$date, 
  vi_vec = sub_dt$evi2, 
  weights_vec = sub_dt$weights,
  verbose = TRUE
)
head(results$phenos)

