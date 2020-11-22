# Bayesian_LSP
A Bayesian hierarchical model that quantifies long-term annual land surface phenology from sparse 30 m Landsat time series (well, it doesn't have to be Landsat).

## Installation
The scripts are written in R programming language and use JAGS software to conduct the MCMC sampling for the Bayesian model. To run the scripts, users should install R along with some packages and JAGS. 

### R 
We use R v3.6.2. Although it should not be limited to this R version, but all of the scripts were tested under v3.6.2. 
Needed R packages are:
* data.table (most of the data in the scripts are processed by functions of data.table. Well, I like data.table!)
* rjags (for communicating with JAGS software)
* RColorBrewer (for plotting the result)
* minpack.lm (it provides functions for non-linear least square fit)

### JAGS
We use JAGS v4.3.0. Again, it should not be limited to this JAGS version, but all of the scripts were tested under v4.3.0.
Users can download and install JAGS from this webpage. To use our scripts, there's no need to know how to use JAGS, actually after installation, users don't even need to open JAGS, we'll use R code to communitate with it.

## Test data
We provide test data in the repository for users to quickly test their development environment as well as the Bayesian model. The test data is a cached R dataset file named "test_ts.Rds". It contains a Landsat EVI2 time series with columns including date and EVI2 value. 

## To run the model



