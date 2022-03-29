#************************************************************************************
# Description: Working with HLS data API, including image query and download.
#     Note the `.netrc` file login of earthdata doesn't work for me, so I figured a 
#     way around. Currently, the `DownloadQueriedHLS` function would download the 
#     whole HLS tile and if an roi is provided, it would crop and overwrite the 
#     downloaded image file.
# Author: Xiaojie(J) Gao
# Date: 2021-12-19
#************************************************************************************


library(rgdal)
library(raster)
library(httr)
library(data.table)
library(jsonlite)
library(tools)
library(geojsonio)
library(geojsonR)

#' @examples
#' \dontrun{
#'
#' # read a study region
#' roi_geojson <- readOGR("/Volumes/GoogleDrive/My Drive/
#'   Research/urban_pheno/Data/boston_bbox.geojson")
#' # view the study region
#' mapview::mapview(roi_geojson, col.regions = NA)
#'
#' # search for images within a time period
#' search_req <- QueryHLS(roi = extent(roi_geojson), 
#'    start_date = "2021-08-01", end_date = "2021-08-15")
#' # Format the query result as a data frame
#' img_df <- FormatReq2DT(search_req)
#' DT::datatable(img_df)
#'
#' # Can take a look at the first requested image
#' BriefView(search_req$features[1, ]$assets$browse$href)
#'
#'
#' # Download images
#' DownloadQueriedHLS(
#'     img_urls = img_df$Asset_Link[1:3],
#'     out_dir = "/Users/xgao26/Desktop/Local_doc/Temp/HLS",
#'     username = "[Earthdata username]",
#'     password = "[Earthdata password]",
#'     roi = roi_geojson,
#'     crop = TRUE
#' )
#' }
#'
search_URL <- "https://cmr.earthdata.nasa.gov/stac/LPCLOUD/search"


# start_date and end_date should be YYYY-MM-DD strings
QueryHLS <- function(roi, pt = NULL, start_date, end_date,
                     prod = list("HLSS30.v2.0", "HLSL30.v2.0")) {
    bbox_str <- paste(roi[1], roi[3], roi[2], roi[4], sep = ",")
    datetime <- paste0(start_date, "T00:00:00Z/", end_date, "T23:59:59Z")

    # Submit a query
    search_body <- list(
        limit = 100,
        datetime = datetime,
        bbox = bbox_str,
        collections = prod
    )

    search_req <- POST(search_URL, body = search_body, encode = "json") %>%
        content(as = "text") %>%
        fromJSON()

    return(search_req)
}


# Format a search request to a data frame
FormatReq2DT <- function(search_req,
                         s_bands = c("B8A", "B04", "Fmask"),
                         l_bands = c("B05", "B04", "Fmask")) {
    search_features <- search_req$features
    # A searchable data table
    granule_list <- list()
    n <- 1
    # Get the NIR, Red, and Quality band layer names
    for (item in row.names(search_features)) { 
        if (search_features[item, ]$collection == "HLSS30.v2.0") {
            bands <- s_bands
        } else {
            bands <- l_bands
        }
        for (b in bands) {
            f <- search_features[item, ]
            b_assets <- f$assets[[b]]$href

            # Make a data frame including links and other info
            df <- data.frame(
                Collection = f$collection, 
                Granule_ID = f$id,
                Cloud_Cover = f$properties$`eo:cloud_cover`,
                band = b,
                Asset_Link = b_assets, stringsAsFactors = FALSE
            )
            granule_list[[n]] <- df
            n <- n + 1
        }
    }
    search_df <- do.call(rbind, granule_list)

    return(search_df)
}


# Quickly see an requested online image
BriefView <- function(browse_img_url) {
    browse_req <- GET(browse_img_url) %>%
        httr::content()

    plot(0:1, 0:1, type = "n", ann = FALSE, axes = FALSE)
    rasterImage(browse_req, 0, 0, 1, 1)
}


# This function doesn't work, but I keep it here for a future fix
CropCloudHLS <- function(roi, img_urls) {
    rgdal::setCPLConfigOption(ConfigOption = "GDAL_HTTP_UNSAFESSL", 
        value = "YES")
    rgdal::setCPLConfigOption(ConfigOption = "GDAL_HTTP_COOKIEFILE", 
        value = ".rcookies")
    rgdal::setCPLConfigOption(ConfigOption = "GDAL_HTTP_COOKIEJAR", 
        value = ".rcookies")
    rgdal::setCPLConfigOption(ConfigOption = "GDAL_DISABLE_READDIR_ON_OPEN", 
        value = "YES")
    rgdal::setCPLConfigOption(ConfigOption = "CPL_VSIL_CURL_ALLOWED_EXTENSIONS", 
        value = "TIF")

    # make sure the spatial reference is consistent
    r1 <- raster(img_urls[1])
    if (crs(roi) != crs(r1)) {
        roi <- spTransform(roi, crs(r1))
    }
    remove(r1)

    # now corp images
    img_stack <- list()
    pb <- txtProgressBar(min = 0, max = length(img_urls), initial = 0, style = 3)
    for (i in seq(length(img_urls))) {
        setTxtProgressBar(pb, i)
        r <- raster(img_urls[i])
        r_crop <- mask(crop(r, extent(roi)), roi)
    }
}


# Download queried HLS tiles
# If an roi is provided and `crop` is set to TRUE, this function would crop 
# and overwrite the downloaded images
DownloadQueriedHLS <- function(img_urls, out_dir, username, password, roi = NULL, 
    crop = TRUE) {
    
    if (dir.exists(out_dir) == FALSE) dir.create(file.path(out_dir))

    for (i in 1:length(img_urls)) {
        print(paste("Downloading...", i, "/", length(img_urls)))
        filename <- file.path(out_dir, basename(img_urls[i]))
        response <- GET(
            img_urls[i],
            write_disk(filename, overwrite = TRUE),
            progress(),
            authenticate(username, password)
        )
    }

    # Crop the images using the roi
    if (crop == TRUE) {
        if (is.null(roi)) stop("Must provide an roi if cropping!")

        imgfiles <- list.files(out_dir, pattern = "*.tif$", full.names = TRUE)

        # make sure the spatial reference is consistent

        r1 <- raster(imgfiles[1])
        if (!identical(crs(roi), crs(r1))) {
            roi <- spTransform(roi, crs(r1))
        }
        remove(r1)

        print("Cropping...")
        pb <- txtProgressBar(0, length(imgfiles), initial = 0, style = 3)
        for (i in 1:length(imgfiles)) {
            setTxtProgressBar(pb, i)
            r <- raster(imgfiles[i])
            r_crop <- mask(crop(r, extent(roi)), roi)
            writeRaster(r_crop, filename = imgfiles[i], overwrite = TRUE)
        }
    }
}
