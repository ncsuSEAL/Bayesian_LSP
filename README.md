# blsp: Bayesian land surface phenology model
Welcome to the blsp R package for creating a Bayesian land surface phenology model. This is a hierarchical model that quantifies long-term annual land surface phenology from temporally sparse optical remote sensing time series (originally developed for 30 m Landsat time series).

For a more in-depth description, please read the paper: [Long-term, medium spatial resolution annual land surface phenology with a Bayesian hierarchical model](https://doi.org/10.1016/j.rse.2021.112484), with the citation:
> Gao, X., Gray, J. M., & Reich, B. J. (2021). Long-term, medium spatial resolution annual land surface phenology with a Bayesian hierarchical model. Remote Sensing of Environment, 261, 112484. https://doi.org/10.1016/j.rse.2021.112484

```
For the exact version in the paper, please go to the `release/reproduce_paper` branch.
```
## Updating the package (move this before fully finishing)
Once you have updated the `.R` files with documentation or new code, run `devtools::document()`, which will update the readonly files. Push the changes to github, then re-install the package.

## How to install blsp package in R
To install this package, please run `devtools::install_github("ncsuSEAL/Bayesian_LSP")` in an R session. Afterwards, you can call the package using `library(blsp)`. 

## Note:
We are currently (as of Apr 2022) improving the computing speed of the BLSP algorithm, thanks to Matt Shisler and Dr. Brian Reich's help. Be sure to watch or star this repo to keep up with our updates.

## Installation
The scripts are written in R programming language and use JAGS software to conduct the MCMC sampling for the Bayesian model. To run the scripts, users need to install certain dependencies, the JAGS software, and the blsp package itself.

### R dependencies
The original code was developed using R v3.6.2, and has been tested on 4.1.1.
Required R packages are:
* `data.table` - most of the data in the scripts are processed by functions of data.table. Well, I like data.table!
* `rjags` - for communicating with JAGS software.
* `minpack.lm` - it provides functions for non-linear least square fit.
* `RColorBrewer` - for plotting the result.
* `viridis` - for plotting the result.
* `lubridate` - for easy parsing date strings.

### JAGS
Please visit the [JAGS website](http://mcmc-jags.sourceforge.net/) for installation. Specifically, we tested the code using JAGS v4.3.0. Once it is installed, users do not need to interact with the software as all communication will be conducted via R.

## Run the model
There are 3 files in the respository's `Code` folder:
* base.R
* Bayesian_LSP_fit.R
* test_ts.Rds

`base.R` contains the needed libraries and pre-defined functions, it'll be sourced in `Bayesian_LSP_fit.R`, which runs the model.

`test_ts.Rds` is a cached R dataset file. It contains a Landsat EVI2 time series with columns including date and EVI2 value. Users can use the test data to quickly run the Bayesian model.

`Bayesian_LSP_fit.R` is easy to understand and run. The only thing need to do before running the script is changing the working directory specified in `setwd()` function, make sure the script can find the `base.R` file and `test_ts.Rds` file. After running the script, there will be a plot showing the result of LSP fit and a table that contains all retreived phenometrics.

Detail information about the model can be found in the paper.

Program result - Model fit:

![](img/model_fit_plot.png)

Program result - Retrieved phenos:

<img src="img/model_fit_phenos.png" alt="" width="500"/>


# Docker
The BLSP docker container installs the required R packages and JAGS.
Assuming that both docker and docker-compose are installed, the container can
be built and started from the CLI by running
```bash
cd scripts && ./build.sh
```
You can enter the container directly by running
```bash
docker exec -ti blsp bash
```

# Acknowledgments
We thank the following people for their assistance with the creation of this package:
- Matt Shisler: optimization of the MCMC code
- [Isabella Hinks](https://github.com/iHinks): translation of MCMC code to C++
- [Owen Smith](https://github.com/ocsmit): development of Docker container
- [Ian McGregor](https://github.com/mcgregorian1): package development, documentation, and vignette
