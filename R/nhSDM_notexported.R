# Not-exported nhSDM functions

# tospf
#'
#' Return a list with sp-flag and sf object, given either sf or sp object
#' 
#' @param spf input sf or sp object
#' @param rastproj raster dataset with desired output projection
#' 
#' @return sf object
#' 
#' @import sf
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
    names(spf)[length(names(spf))] <- "geometry"
  } else {
    stop("Must provide either 'sp' or 'sf'-class spatial object.")
  }
  spf <- st_zm(spf)
  
  # transform if necessary
  if (!missing(rastproj)) {
    if (!is.na(st_crs(spf)$proj4string)) {
      if (st_crs(spf)$proj4string != projection(rastproj)) spf <- st_transform(spf, crs = projection(rastproj))
    } else {
      message("No projection on input features. Assuming features are using raster's projection...")
      st_crs(spf) <- projection(rastproj)
    }
  }
  
  return(list(sp,spf))
}

# gRasterize
#'
#' GDAL Rasterize with temp file parameters and all polygon values = 1
#' 
#' A temporary raster and shapefile are created in the raster temp directory, and 
#' deleted. 
#' 
#' This is used internally in nh_burn instead of raster::rasterize, which was found to
#' be very slow for large and/or many feature rasterizing. In some cases a 
#' 'striping' was noticed in the resulting raster, in areas where there were no features.
#' 
#' @param spf input sf or sp object
#' @param rast raster dataset with desired output projection, extent, cell size
#' @param value integer value to apply to areas covered by spf
#' @param background value to apply to areas not covered by spf
#' 
#' @return raster object
#' 
#' @import sf
#' @import raster
#' @importFrom gdalUtils gdal_rasterize
#' @importFrom methods as
#' 
#' @keywords internal

gRasterize <- function(spf, rast, value = 1, background = NA) {
  
  # handle sp/sf class
  spf <- tospf(spf, rast)[[2]]
  try(spf <- st_union(spf), silent = T)
  # temp names
  tmp <- gsub(".grd", "", rasterTmpFile())
  tmpr <- paste0(tmp, ".tif")
  tmpshp <- paste0(tmp, ".shp")
  values(rast) <- background
  
  # try to rasterize, remove temp files in any case
  tryCatch({
    # gdal rasterize the mask area (fixed paths, just overwrites layer each time)
    writeRaster(rast, tmpr, datatype = "INT2U")
    st_write(obj = spf, dsn = tmpshp, driver = "ESRI Shapefile", quiet = T)
    gdal_rasterize(tmpshp, dst_filename = tmpr, burn = value)
    r1 <- raster(tmpr)
    values(rast) <- values(r1)
    names(rast) <- "gRasterize"
  }, finally = {
    unlink(paste0(tmp,"*"))
  })
  return(rast)
}
