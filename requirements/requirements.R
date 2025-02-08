#!/usr/bin/env Rscript

pkgs <- c(
          "data.table",
          "rjags",
          "minpack.lm",
          "RColorBrewer",
          "viridis",
          "lubridate",
          "pkgdown",
          "profvis",
          "bench",
          "miniUI",
          "pkgdown",
          "devtools",
          "minpack.lm",
          "rstac",
          "terra"
          )

install.packages(pkgs, repos='http://cran.us.r-project.org', dependencies=TRUE)
