# Not-exported nhSDM functions

# tospf
#'
#' Return a list with sp-flag and standardized sf object, given either sf or sp object
#' 
#' @param spf input sf or sp object
#' @param rastproj raster dataset with desired output projection
#' 
#' @return sf object
#' 
#' @importFrom sf st_as_sf st_geometry st_zm st_crs st_transform st_crs<-
#' @importFrom methods as
#' @importFrom raster projection
#' 
#' @keywords internal

tospf <- function(spf, rastproj) {
  if (grepl("^Spatial*", class(spf)[1])) {
    sp <- TRUE
    spf <- st_as_sf(spf)
  } else if (grepl("sf", class(spf)[1])) {
    sp <- FALSE
  } else {
    stop("Must provide either 'sp' or 'sf'-class spatial object.")
  }
  # standardize geometry column name
  spf$geometry <- st_geometry(spf) 
  spf <- st_zm(st_set_geometry(spf, "geometry"))
  spf$geom <- NULL
  
  # spf <- st_zm(spf)

  # transform if necessary
  if (!missing(rastproj)) {
    if (!is.na(st_crs(spf)$proj4string)) {
      if (st_crs(spf)$proj4string != st_crs(rastproj)$proj4string) spf <- st_transform(spf, crs = st_crs(rastproj)$proj4string)
    } else {
      message("No projection on input features. Assuming features are using raster's projection...")
      st_crs(spf) <- projection(rastproj)
    }
  }
  
  return(list(sp,spf))
}

# gRasterize
#'
#' Rasterize with temp file parameters and all polygon values = 1
#' 
#' A temporary raster and shapefile are created in the raster temp directory, and 
#' deleted. 
#' 
#' This is used internally in nh_burn instead of raster::rasterize, which was found to
#' be very slow for large and/or many feature rasterizing. In some cases a 
#' 'striping' was noticed in the resulting raster, in areas where there were no features.
#' 
#' Previous versions used GDAL, now uses R package fasterize.
#' 
#' @param spf input sf or sp object
#' @param rast raster dataset with desired output projection, extent, cell size
#' @param value integer value to apply to areas covered by spf
#' @param background value to apply to areas not covered by spf
#' 
#' @return raster object
#' 
#' @import raster
#' @importFrom sf st_buffer st_cast
#' @importFrom methods as
#' @importFrom fasterize fasterize
#' 
#' @keywords internal

gRasterize <- function(spf, rast, value = 1, background = NA) {
  
  # handle sp/sf class
  spf <- tospf(spf, rast)[[2]]
  
  # fasterize needs polygons
  if (grepl("POINT|LINESTRING", st_geometry_type(spf)[1]))
    spf <- st_cast(st_buffer(spf, res(rast)[1]*sqrt(2) / 2), "MULTIPOLYGON")
  if (is.numeric(value)) {
      spf$burnval <- value
    } else {
      spf$burnval <- as.data.frame(spf)[,value]
    }
  rast <- fasterize(spf, rast, field = "burnval", fun = "max")
  return(rast)
}

# nodes
#'
#' Calculate nodes in a line network, using line direction
#' 
#' Returns same object with added columns startNode and endNode
#' 
#' @param spf input sf object
#' 
#' @importFrom lwgeom st_snap_to_grid st_startpoint st_endpoint
#' @import sf
#' 
#' @return sf object
#' 
#' @keywords internal

nodes <- function(spf) {
  message("Calculating line connections...")
  st <- st_snap_to_grid(st_transform(st_startpoint(spf), 3857), 0.01) # need to transform from latlong to use snap to grid
  spN <- unlist(lapply(st_equals(st), min))
  ep <- st_snap_to_grid(st_transform(st_endpoint(spf), 3857), 0.01)
  epN <- unlist(lapply(st_equals(ep), min))
  epN <- epN + max(spN)
  
  spf$startNode <- spN
  spf$endNode <- epN
  
  # re-assign endnode nums to match startnodes
  suppressWarnings({
    match <- unlist(lapply(st_equals(ep, st), min))
    match[match==Inf] <- NA
    spf$endNode[!is.na(match)] <- spN[match[!is.na(match)]]
  })
  return(spf)
}
