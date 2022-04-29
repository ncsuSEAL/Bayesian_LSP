#************************************************************************************
# Description: Download Landsat time series using AppEEARS.
# Date: 2022-03-27
#************************************************************************************

#' Log in to AppEEARS system.
#' 
#' @param usr,pwd Username and password of AppEEARS system.
#' @return Return the token string
#' @export
Login <- function(usr, pwd) {
    # if username and password are not provided, input them in the ternimal window
    if (is.null(usr) | is.null(pwd)) {
        message("username and password must be provided.
            You can create one on the EarthData website.")
        # Enter NASA Earthdata Login Username
        usr <- getPass::getPass(msg = "Enter NASA Earthdata Login Username: ")
        # Enter NASA Earthdata Login Password
        pwd <- getPass::getPass(msg = "Enter NASA Earthdata Login Password: ")
    }

    if (usr == "" | pwd == "") {
        message("username and password must be provided.
            You can create one on the EarthData website.")
        return(NULL)
    }

    secret <- jsonlite::base64_enc(paste(usr, pwd, sep = ":"))

    response <- httr::POST(
        paste0(API_URL, "login"),
        httr::add_headers(
            "Authorization" = paste("Basic", gsub("\n", "", secret)),
            "Content-Type" = "application/x-www-form-urlencoded;charset=UTF-8"
        ),
        body = "grant_type=client_credentials"
    )
    response_content <- httr::content(response)
    token <- paste("Bearer", response_content$token)

    return(token)
}


#' log out from AppEEARS.
#' 
#' @param token AppEEARS login token. Use `Login()` to get a token.
#' @return 
#' @export
Logout <- function(token) {
    req <- httr::POST(
        paste0(API_URL, "logout"),
        httr::add_headers(Authorization = token)
    )
}


#' Submit a point task. `task_type` can either be 'point' or 'area', but relative 
#' params need to be provided.
#' 
#' @param token AppEEARS login token. Use `Login()` to get a token.
#' @param task_name The name of the task.
#' @param task_type The AppEEARS task type, either "point" or "area".
#' @param start_date,end_date The start and end dates of the requested time series.
#' @param recursive If the specified start and end dates are recursive. See AppEEARS
#' documentation.
#' @param layers The interested layers.
#' @param out_format The desired image format, default is `geotiff`.
#' @param out_proj The desired image projection.
#' 
#' @return 
#' 
#' @examples 
#' \dontrun{
#' 
#' }
SubmitTask <- function(token, task_name, task_type = "point", 
    start_date = NULL, end_date = NULL, recursive = FALSE, layers = NULL, 
    point_df = NULL, polygon_file = NULL, 
    out_format = "geotiff", out_proj = NULL
) {
    # check arguments
    if (sum(is.null(start_date), is.null(end_date), is.null(layers)) > 0) {
        stop("Please specify parameters!")
    }
    if (tolower(task_type) == "point" && is.null(point_df)) {
        stop("Point task but no points provided!")
    }
    if (tolower(task_type) == "Area" && 
        (is.null(polygon_file) | is.null(out_format) | is.null(out_proj))) {
        stop("Area task but paramter(s) not specified!")
    }
    

    # Format Dates
    dates <- data.frame(startDate = start_date, endDate = end_date)

    if (tolower(task_type) == "point") { # ~ Point tasks
        task_info <- list(dates, layers, point_df) # Create a list of data frames
        names(task_info) <- c("dates", "layers", "coordinates") # Assign names

        task <- list(task_info, task_name, task_type) # Create a nested list
        names(task) <- c("params", "task_name", "task_type") # Assign names

        # Convert to JSON object
        task_json <- jsonlite::toJSON(task, auto_unbox = TRUE) 

        response <- httr::POST(paste0(API_URL, "task"),
            body = task_json,
            encode = "json",
            httr::add_headers(Authorization = token, "Content-Type" = "application/json")
        )

        task_content <- httr::content(response) # Retrieve content of the request
        # Convert the content to JSON object
        task_response <- jsonlite::prettify(
            jsonlite::toJSON(task_content, auto_unbox = TRUE)
        ) 
        
        return(task_response)

    } else if (tolower(task_type) == "area") { # ~ Area tasks
        # read the polygon file
        if (file_ext(polygon_file) == "geojson") {
            polygon_f <- rgdal::readOGR(polygon_file)
        } else if (tools::file_ext(polygon_file) == "shp") {
            polygon_f <- rgdal::readOGR(dsn = dirname(polygon_file), 
                layer = tools::file_path_sans_ext(basename(polygon_file)))
        } else {
            stop("Please provide a valid shp or geojson file!")
        }

        # Convert the data frame to GeoJSON
        gc_json <- geojsonio::geojson_json(polygon_f, geometry = "polygon")
        gc_js <- geojsonR::FROM_GeoJson(gc_json) # Read the GeoJSON
        gc_js$features[[1]]$geometry$coordinates <- list(
            gc_js$features[[1]]$geometry$coordinates)

        out <- list(out_proj)
        names(out) <- c("projection")
        out$format$type <- out_format

        task_info <- list(dates, layers, out, gc_js) # Create a list of data frames
        names(task_info) <- c("dates", "layers", "output", "geo") # Assign names
        task <- list(task_info, task_name, task_type) # Create a nested list
        names(task) <- c("params", "task_name", "task_type") # Assign names
        task_json <- jsonlite::toJSON(task, auto_unbox = TRUE, digits = 10)

        response <- httr::POST(paste0(API_URL, "task"),
            body = task_json, encode = "json",
            add_headers(Authorization = token, "Content-Type" = "application/json")
        )

        task_content <- httr::content(response) # Retrieve content of the request
        # Convert the content to JSON and prettify it
        task_response <- jsonlite::toJSON(task_content, auto_unbox = TRUE) 
        jsonlite::prettify(task_response)
    }
}


#' Check the status of current tasks
#' 
#' @param token AppEEARS login token. Use `Login()` to get a token.
#' @param limit The query limit amount.
#' @param task_name The name of the task.
#' @param brief Logical. `TRUE` to show brief information, `FALSE` to show the full
#'  information (default).
#' @return Job status information
#' @export
#' 
#' @examples 
#' \dontrun{
#' 
#' }
CheckTaskStatus <- function(token, limit, task_name = NULL, brief = FALSE) {
    params <- list(limit = limit, pretty = TRUE)
    response_req <- httr::GET(paste0(API_URL, "task"), 
        query = params, 
        httr::add_headers(Authorization = token)
    )

    # Retrieve content of the request
    response_content <- httr::content(response_req) 
    # Convert the content to JSON object
    status_response <- jsonlite::toJSON(response_content, auto_unbox = TRUE) 
    response_df <- jsonlite::fromJSON(status_response)
    names(response_content) <- response_df$task_name

    if(is.null(task_name) == FALSE) {
        my_task <- response_content[grepl(task_name, response_content)]
        my_task_json <- jsonlite::toJSON(my_task, auto_unbox = TRUE)

        if (brief == TRUE) {
            my_task_df <- jsonlite::fromJSON(my_task_json)[[task_name]]
            my_task_df_b <- data.frame(
                task_name = my_task_df$task_name, 
                status = my_task_df$status, 
                id = my_task_df$task_id
            )
            return(my_task_df_b)
        } else {
           return(jsonlite::prettify(my_task_json))
        }
    }

    if (brief == TRUE) {
        return(response_df[, c("task_name", "status", "task_id")])
    } else {
       return(jsonlite::prettify(jsonlite::toJSON(response_content)))
    }
}


#' Download the job from AppEEARS.
#' 
#' @param token AppEEARS login token. Use `Login()` to get a token.
#' @param task_id The ID of the task. Use `CheckTaskStatus()` to get the IDs.
#' @param out_dir The target directory to store downloaded files.
#' @return 
#' @export
#' 
#' @examples 
#' \dontrun{
#' 
#' }
DownloadTask <- function(token, task_id, out_dir) {
    # Request the task bundle info from API bundle URL
    response <- httr::GET(paste0(API_URL, "bundle/", task_id), 
        httr::add_headers(Authorization = token))
    response_content <- httr::content(response) # Retrieve content of the request
    # Convert the content to JSON object
    bundle_response <- jsonlite::toJSON(response_content, auto_unbox = TRUE) 
    jsonlite::prettify(bundle_response)

    bundle <- jsonlite::fromJSON(bundle_response)$files
    for (id in bundle$file_id) {
        # retrieve the filename from the file_id
        filename <- bundle[bundle$file_id == id, ]$file_name
        # create a destination directory to store the file in
        filepath <- paste(out_dir, filename, sep = "/")
        suppressWarnings(dir.create(dirname(filepath)))
        # write the file to disk using the destination directory and file name
        response <- httr::GET(
            paste0(API_URL, "bundle/", task_id, "/", id),
            httr::write_disk(filepath, overwrite = TRUE),
            httr::progress(),
            httr::add_headers(Authorization = token)
        )
    }
}


#' Request Landsat ARD point time series data using AppEEARS.
#' 
#' @param pts_df A data frame containing point information, including 4 columns:
#'  `id`, `longitude`, `latitude`, `category`.
#' @param task_name Task name.
#' @param token AppEEARS login token. Use `Login()` to get a token.
#' @return TRUE if the job is successfully submited.
#' @export
#' 
#' @examples 
#' \dontrun{
#' lat <- c(36.206228, 37.289327) # Latitude of the point sites
#' lon <- c(-112.127134, -112.973760) # Longitude of the point sites
#' id <- c("0", "1") # ID for the point sites
#' pts_df <- data.frame(id = id, longitude = lon, latitude = lat, 
#'     category = category)
#' }
ReqArdPoint <- function(pts_df, task_name, token,
    start_date = "01-01-1984", end_date = "12-31-2022"
) {
    # Target layers
    landsatARD_layers <- data.frame(
        product = c(
            rep("CU_LT05.001", 3), 
            rep("CU_LE07.001", 3), 
            rep("CU_LC08.001", 3)),
        layer = c(
            "SRB3", "SRB4", "PIXELQA", 
            "SRB3", "SRB4", "PIXELQA", 
            "SRB4", "SRB5", "PIXELQA"
        )
    )

    # Submit point task
    SubmitTask(
        token = token, task_name = task_name, task_type = "point",
        start_date = start_date, end_date = end_date, layers = landsatARD_layers,
        point_df = pts_df
    )

    return(TRUE)
}


