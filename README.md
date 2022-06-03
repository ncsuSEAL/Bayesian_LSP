# blsp: Bayesian land surface phenology model
Welcome to the `blsp` R package for creating a Bayesian land surface phenology model. This is a hierarchical model that quantifies long-term annual land surface phenology from temporally sparse optical remote sensing time series (originally developed for 30 m Landsat time series).

For a more in-depth description, please read the paper: [Long-term, medium spatial resolution annual land surface phenology with a Bayesian hierarchical model](https://doi.org/10.1016/j.rse.2021.112484), with the citation:
> Gao, X., Gray, J. M., & Reich, B. J. (2021). Long-term, medium spatial resolution annual land surface phenology with a Bayesian hierarchical model. Remote Sensing of Environment, 261, 112484. https://doi.org/10.1016/j.rse.2021.112484

```
For the exact version in the paper, please go to the `release/reproduce_paper` branch.
```

## How to install blsp package in R

We use JAGS (Just Another Gibbs Sampler) software to conduct Markov Chain Monte Carlo (MCMC) sampling for the Bayesian model. Please install JAGS software before installing the `blsp` package. Please visit the [JAGS website](http://mcmc-jags.sourceforge.net/) for installation. Don't worry if you know nothing about JAGS, you don't even need to open it after installing. We use R to communicate with it.

Next, in the R terminal, run:
```
devtools::install_github("ncsuSEAL/Bayesian_LSP", build_vignettes = TRUE)
```
Afterwards, you can call the package using `library(blsp)`. Run `??blsp` to see the vignettes and functions available in the package. 

## Note:
We are currently (as of Apr 2022) improving the computing speed of the BLSP algorithm, thanks to Matt Shisler and Dr. Brian Reich's help. Be sure to watch or star this repo to keep up with our updates.

## The package functionality 
The package takes sparse vegetation index observations from the entirety of the Landsat time series (for example), and create a continuous estimate of annual land surface phenology. In addition to calculating start of season (SOS) and end of season (EOS) dates, the model also calculated pixel-wise uncertainty estimates for each of these phenometrics. 

The model fit is shown in the below figure:

![](img/model_fit_plot.png)

And, the estimated phenometrics and their 95% credible intervals are stored in a table returned by the `FitBLSP()` function of the `blsp` package:

<img src="img/model_fit_phenos.png" alt="" width="500"/>

For detailed introduction of the package usage, please use `??blsp` to see the vignettes.


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
- [Matt Shisler](https://github.com/mattshisler): optimization of the MCMC code
- [Isabella Hinks](https://github.com/iHinks): translation of MCMC code to C++ (the work hasn't been merged to this package, but will come soon.)
- [Owen Smith](https://github.com/ocsmit): development of Docker container
- [Ian McGregor](https://github.com/mcgregorian1): package development, documentation, and vignette
