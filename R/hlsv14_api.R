#************************************************************************
# Description: Download Harmonized Landsat and Sentinel-2 imagery. See https://hls.gsfc.nasa.gov/
# Author: Xiaojie Gao
# Date: 2020-10-14
#************************************************************************
library(utils)
library(httr)
library(XML)
library(raster)
library(rgdal)
library(tools)
library(gdalUtils)


#' Download the HLS v1.4 image product. 
#' @param out_dir Data save directory.
#' @param tiles Tile names needed.
#' @param yrs Years needed.
#' @param crop_file A study region vector file path. If provided, the image will be croped to the specified region and saved to GTiff format.
#' @return NULL
#' @export
#' @example 
#' /dontrun {
#' Download_HLS(
#'     out_dir = "D:/Temp",
#'     tiles = "15TYH", 
#'     yrs = 2019, 
#'     crop_file = "Q:/My Drive/Research/Urban_pheno/Data/ARD/madison_bbox.geojson")
#' }
Download_HLS <- function(out_dir, tiles, yrs, crop_file = NULL) {
    if(!dir.exists(out_dir)) dir.create(out_dir)

    root_url <- "https://hls.gsfc.nasa.gov/data/v1.4"

    for(prod in c("S30", "L30")) {
        if (!dir.exists(file.path(out_dir, prod))) dir.create(file.path(out_dir, prod))

        for(yr in yrs) {
            for(tile in tiles) {
                # construct data file urls
                tile_folder_url <- file.path(root_url, prod, yr, substr(tile, 1, 2), substr(tile, 3, 3), substr(tile, 4, 4), substr(tile, 5, 5))
                # request for contents
                r <- GET(tile_folder_url)
                doc <- htmlParse(r)
                links <- xpathSApply(doc, "//a/@href")
                data_links <- links[grepl(".*.hdf", links)]
                # construct data links
                data_urls <- paste0(tile_folder_url, "/", data_links)
                # download
                sapply(seq_along(data_urls), function(i) {
                    file_name <- file.path(out_dir, prod, data_links[i])
                    if (!file.exists(file_name)) {
                        download.file(data_urls[i], file_name, method = "curl")
                        # crop file to study region
                        if (!is.null(crop_file)) {
                            # find out if the .hdf and .hdr files are all downloaded
                            if (grepl(".hdr", file_name)) {
                                hdr_file <- file_name
                                hdf_file <- tools::file_path_sans_ext(file_name)
                            } else {
                                hdf_file <- file_name
                                hdr_file <- paste0(file_name, ".hdr")
                            }
                            if (file.exists(hdf_file) & file.exists(hdr_file)) {
                                sds <- get_subdatasets(hdf_file)
                                r <- stack(sds)
                                crop_extent <- readOGR(crop_file, verbos = FALSE)
                                crop_extent_reproj <- spTransform(crop_extent, crs(r))
                                r_crop <- raster::crop(r, crop_extent_reproj, snap = "in")
                                writeRaster(r_crop, hdf_file, format = "GTiff", overwrite = TRUE)
                                message(paste("croped file and saved it to GTiff: ", hdf_file))
                                # remove the original hdf file
                                file.remove(hdf_file)
                                file.remove(hdr_file)
                            }
                        }
                    }
                })
            }
        }
    }
}























