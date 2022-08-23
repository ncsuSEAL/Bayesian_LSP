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
          "devtools",
          "geojsonio",
          "geojsonR",
          "getPass",
          "minpack.lm"
          )

install.packages(pkgs, repos='http://cran.us.r-project.org', dependencies=TRUE)
