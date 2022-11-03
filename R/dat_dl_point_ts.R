#******************************************************************************
# Description: Download Landsat point time series using Microsoft Planetary
# Computer with STAC API.
# 
# Date: 2022-11-02
#******************************************************************************


#' Calculate EVI2 values.
#' 
#' @param nir_band NIR band.
#' @param red_band RED band.
#' @return EVI2 values.
#' @noRd
CalEVI2 <- function(nir_band, red_band) {
    red <- red_band * 0.0000275 - 0.2
    nir <- nir_band * 0.0000275 - 0.2

    evi2 <- 2.5 * ((nir - red) / (1 + nir + 2.4 * red))

    return(as.numeric(evi2))
}


#' Parse Landsat cloud and snow QA values.
#' 
#' @param x The QA value.
#' @return A list with logic values indicating `fill`, `cloud`, `cloudShadow`,
#' and `snow`.
#' @noRd
LandsatCloudSnowQA <- function(x) {
    ## Bit 0 - if pixel is fill, then true
    fill <- ifelse(bitwAnd(x, 1), TRUE, FALSE)
    ## Bit 3 - if cloud, then true
    cloud <- ifelse(bitwAnd(bitwShiftR(x, 3), 1), TRUE, FALSE)
    ## Bit 4 - if cloud shadow, then true
    cloudShadow <- ifelse(bitwAnd(bitwShiftR(x, 4), 1), TRUE, FALSE)
    ## Bit 5 - if snow, then true
    snow <- ifelse(bitwAnd(bitwShiftR(x, 5), 1), TRUE, FALSE)

    return(list(fill = fill, 
        cloud = cloud, cloudShadow = cloudShadow, 
        snow = snow
    ))
}

#' Use Microsoft Planetary Computer with STAC API to get Landsat EVI2 time 
#' series for any point location specified by longitude and latitude.
#' 
#' @param pt_coords Point location. Longitude and latitude.
#' @param focalDates Temporal period. 
#' @param ncore Number of cores used to parallel the process.
#' @return A data.table containing EVI2 time series along with QA values. 
#' @export
#' @import data.table
#' 
#' @examples
#' \dontrun{
#'   pt_coords <- data.table::data.table(x = -71.700975, y = 43.945733)
#'   focalDates <- "1984-01-01/1989-06-04"
#'   ncores <- 5
#'   val_dt <- GetEvi2PointTs(pt_coords, ncore = ncores)
#'   val_dt <- data.table::setorder(val_dt, date)
#'   plot(val_dt[qa == 0, .(date, evi2)], pch = 16, type = "b")
#' }
GetEvi2PointTs <- function(pt_coords, focalDates = "1984-01-01/2022-12-31", 
    ncores = 1
) {
    # Make the bounding box for each point
    pt <- terra::vect(cbind(pt_coords$x, pt_coords$y), 
        crs = "EPSG:4326"
    )
    pt_buf <- terra::buffer(pt, width = 1)
    pt_ext <- terra::ext(pt_buf)

    # Split the entire time period
    start_date <- strsplit(focalDates, "/")[[1]][1] %>% as.Date()
    end_date <- strsplit(focalDates, "/")[[1]][2] %>% as.Date()

    if ((data.table::year(end_date) - data.table::year(start_date)) <= 5) {
        start_end_dates <- focalDates
    } else {
        yrs <- seq(data.table::year(start_date), 
            data.table::year(end_date), 
            by = 5
        )

        start_end_dates <- NULL
        for (i in seq_along(yrs)) {
            if (i == 1) {
                sd <- start_date
                ed <- paste0(yrs[i + 1] - 1, "-12-31") %>% as.Date()
            } else if (i == length(yrs)) {
                sd <- paste0(yrs[i], "-01-01") %>% as.Date()
                ed <- end_date
            } else {
                sd <- paste0(yrs[i], "-01-01") %>% as.Date()
                ed <- paste0(yrs[i + 1] - 1, "-12-31") %>% as.Date()
            }

            start_end_dates <- c(start_end_dates, paste(sd, ed, sep = "/"))
        }
    }
    

    # Request images from the server
    
    # Planetary Computer API
    s_obj <- rstac::stac("https://planetarycomputer.microsoft.com/api/stac/v1/")

    it_obj <- lapply(start_end_dates, function(focal_dates) {
        obj <- rstac::stac_search(s_obj,
            collections = "landsat-c2-l2",
            ids = NULL, # could specify this if wanted
            bbox = pt_ext[c(1, 3, 2, 4)],
            datetime = focal_dates,
            limit = 1000
        ) %>%
            rstac::get_request() %>%
            rstac::items_sign(sign_fn = rstac::sign_planetary_computer()) %>%
            suppressWarnings()
        return(obj)
    })


    # Iterate temporal periods
    res_dt <- lapply(it_obj, function(it) {
        # Parallel processing
        ncores <- ifelse(ncores > parallel::detectCores() - 1, 
            parallel::detectCores() - 1,
            ncores
        )
        cl <- parallel::makeCluster(ncores)
        calls <- parallel::clusterCall(cl, function() {
            suppressWarnings({
                require(terra)
                require(magrittr)
            })
        })
        parallel::clusterExport(cl,
            c("CalEVI2", "pt_coords"),
            envir = environment()
        )

        # Iterate features
        val_dt <- parallel::parLapply(cl, X = it$features, function(img) {
        # val_dt <- lapply(it$features, function(img) { browser() # For debug
            img_id <- img$id
            date <- as.Date(gsub("T.*", "", img$properties$datetime))
            epsg <- img$properties$`proj:epsg`

            # Project the point buffer to the epsg
            # Has to create the points again as terra doesn't serielize spatial
            # objects.
            pt <- terra::vect(cbind(pt_coords$x, pt_coords$y),
                crs = "EPSG:4326"
            ) %>%
                terra::project(paste0("EPSG:", epsg))

            row <- tryCatch({
                # Red
                red_band <- paste0("/vsicurl/", img$assets$red$href) %>%
                    terra::rast() %>%
                    terra::extract(pt) %>%
                    as.numeric()

                # Nir
                nir_band <- paste0("/vsicurl/", img$assets$nir08$href) %>%
                    terra::rast() %>%
                    terra::extract(pt) %>%
                    as.numeric()

                # QA
                qa_band <- paste0("/vsicurl/", img$assets$qa_pixel$href) %>%
                    terra::rast() %>%
                    terra::extract(pt) %>%
                    as.numeric()

                # Calculate EVI2
                evi2 <- CalEVI2(nir_band, red_band)

                return(data.table::data.table(
                    img_id = img_id,
                    lon = pt_coords$x,
                    lat = pt_coords$y,
                    evi2 = evi2,
                    date = date,
                    qa = qa_band
                ))

            }, error = function(e) {
                return(NULL)
            })

            return(row)

        }) %>%
            do.call(rbind, .)

        parallel::stopCluster(cl)

        return(val_dt)
    }) %>%
        do.call(rbind, .)
    
    # Parse the QA band
    qa_parse <- LandsatCloudSnowQA(res_dt$qa)
    res_dt$snow <- qa_parse$snow
    res_dt <- res_dt[
        qa_parse$fill == FALSE &
        qa_parse$cloud == FALSE & 
        qa_parse$cloudShadow == FALSE,
    ]

    res_dt <- data.table::setorder(res_dt, date)

    return(res_dt)
}

