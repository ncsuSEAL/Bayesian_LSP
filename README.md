<img src="man/figures/logo.png" align="right" height="138" />

# blsp: Bayesian land surface phenology model
Welcome to the `blsp` R package for creating a Bayesian land surface phenology model. This is a hierarchical model that quantifies long-term annual land surface phenology from temporally sparse optical remote sensing time series (originally developed for 30 m Landsat time series).

For a more in-depth description, please read the paper: [Long-term, medium spatial resolution annual land surface phenology with a Bayesian hierarchical model](https://doi.org/10.1016/j.rse.2021.112484), with the citation:
> Gao, X., Gray, J. M., & Reich, B. J. (2021). Long-term, medium spatial resolution annual land surface phenology with a Bayesian hierarchical model. Remote Sensing of Environment, 261, 112484. https://doi.org/10.1016/j.rse.2021.112484

> **Note**
> 
> For the exact version in the paper, please go to the `release/reproduce_paper` branch.

To cite the package, please use:

> Xiaojie Gao, Ian R. McGregor, Owen Smith, Isabella Hinks, & Matt Shisler. (2022). The blsp R package with a Bayesian land surface phenology model (1.0). Zenodo. https://doi.org/10.5281/zenodo.6824017

## How to install blsp package in R

We use JAGS (Just Another Gibbs Sampler) software to conduct Markov Chain Monte Carlo (MCMC) sampling for the Bayesian model. Please install JAGS software before installing the `blsp` package. Please visit the [JAGS website](http://mcmc-jags.sourceforge.net/) for installation. Don't worry if you know nothing about JAGS, you don't even need to open it after installing. We use R to communicate with it.

Next, in the R terminal, run:
```r
devtools::install_github("ncsuSEAL/Bayesian_LSP", build_vignettes = TRUE)
```
Afterwards, you can call the package using `library(blsp)`. Run `help(package = "blsp")` to see the vignette and functions available in the package. 

## Note
We are currently (as of June 2022) improving the computing speed of the BLSP algorithm, thanks to Matt Shisler and Dr. Brian Reich's help. Be sure to watch or star this repo to keep up with our updates.

## The package functionality 
The package takes sparse vegetation index observations from the entirety of the Landsat time series (for example), and create a continuous estimate of annual land surface phenology. In addition to calculating start of season (SOS) and end of season (EOS) dates, the model also calculated pixel-wise uncertainty estimates for each of these phenometrics. 

The model fit is shown in the below figure:

![](img/model_fit_plot.png)

And, the estimated phenometrics and their 95% credible intervals are stored in a table returned by the `FitBLSP()` function of the `blsp` package:

| Year | midgup_lwr | midgup | midgup_upr | midgdown_lwr | midgdown | midgdown_upr |
| :--: | :--------: | :----: | :--------: | :----------: | :------: | :----------: |
| 1984 |    130     |  139   |    146     |     277      |   284    |     291      |
| 1985 |    132     |  138   |    140     |     272      |   281    |     288      |
| ...  |    ...     |  ...   |    ...     |     ...      |   ...    |     ...      |
| 2023 |    ...     |  ...   |    ...     |     ...      |   ...    |      ..      |

Starting from v1.5, in addition to `midgup` (SOS) and `midgdown` (EOS), we also support getting more detailed phenometrics using a threshold-based method. The method can be configured when using `FitBLSP(.., opt = list(method = "threshold"))`. The detailed phenometrics and their amplitude threshold are shown in the following table and figure:

| Phenometric  | Threshold               |
| :----------: | ----------------------- |
|   Greenup    | 15% amplitude in spring |
|  MidGreenup  | 50% amplitude in spring |
|   Maturity   | 90% amplitude in spring |
|     Peak     | 100% amplitude          |
|  Senescence  | 90% amplitude in autumn |
| MidGreendown | 50% amplitude in autumn |
|   Dormancy   | 15% amplitude in autumn |

![](img/model_fit_more_phenos.png)

Also from v1.5, we support both 6- and 7-parameter double-logistic functions. To specify which function to use, pass a `model` string to the `FitBLSP()` function, e.g.,  `FitBLSP(..., model = "dblog6")`. To use the 6-parameter model, do `model = "dblog6"`; while `model = "dblog7"` will use the 7-parameter model, which is the default value.

Starting from v1.7, we added a `greendown_aware` parameter to account for the summer EVI2 greendown phenomenon when using the threshold-based phenometrics. Specifically, when the greendown phenomenon is substantial (e.g., in PhenoCam data), the `Senescence` metric, which is defined as 90% amplitude in autumn, can be biased early. However, this EVI2 decrease is not necessarily `Senescence` (the exact mechanisms that induced this summer greendown are still unclear. Some previous studies have attributed it to shadow and/or leaf angle). So, by using `FitBLSP(.., opt = list(method = "threshold", greendown_aware = TRUE))`, the `Senescence` metric will be retrieved as the end date of summer greendown (or, the edge of the curve), and the `MidGreendown` will be the date with the mininum first derivative of the autumn EVI2 curve after `Senescence`. 

For detailed introduction of the package usage, please use `help(package = "blsp")` to see the vignettes. We also provide Google Earth Engine javascript script and Microsoft Planetary Computer R functions to help users get Landsat time series for any latitude and longitude points so that users can try the `blsp` package with minimal effort in preparing data (see the vignettes).

> **Note** 
> 
> Unlike other land surface phenology products, we don't have QA/QC flags. The reason is, from our current experience, that the quality of the retrieved phenometrics can be indicated from the uncertainty. For example, if the uncertainty for a phenometric is very large, it indicates that the phenometric might be of low quality; otherwise, the pheometirc is trustable. This strategy may be changed based on future experience with the BLSP model, though.
> 
> Some data pre-processing such as filling in the extremly low values in the winter period using 2th percentile and removing abnormal low values in the summer period can help fitting the model better. Those abnormal observations should be captured by cloud detection but sometimes it fails. 

# Known limitations
Here are some limitations users frequently asked, we appreciate the feedback and want to notify future users to be aware of them. 

- **Computing speed**. As the method uses Markov Chain Monte Carlo (MCMC) sampling to estimate model parameters, it requires some computing power and can be time consuming. This is the reason we do not provide functions for image processing in the package. For image processing, users need to run `blsp` on a super computer using parallel processing.
- **Interannual variability**. As the algorithm computes a particular year's LSP by considering data available within the current year and using information from other years as prior, when data in the current year are very limited, especially for the seasonal transition periods, the prior information would get more weight in the final LSP calculation and thus the overall LSP time series may lack interannual variability. We encourage users to check the uncertainty of phenometrics in the `blsp` result as well as the fitted time series.


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

---

_Graphs used in the icon are created by <a href="https://www.flaticon.com/free-icons" title="seal icons">Freepik - Flaticon</a>_
