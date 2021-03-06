# nh_burn

#' Add areas represented by features ('burn-in') to a binary raster SDM output
#' 
#' Takes spatial features, a continuous raster (values between 0 and 1), and 
#' returns a classified binary (0/1) raster. Areas greater than \code{orig.thresh} are
#' set to 1, all others 0. Areas intersecting \code{spf} features
#' are assigned a value of 1 in the returned classified raster. See details
#' for default calculation of \code{orig.thresh} and usage of \code{buffer}.
#' 
#' When \code{buffer} is used or \code{orig.thresh} is not provided,
#' a minimum cell value (min.cell) across all \code{spf} is calculated.
#' The min.cell value will be used as a threshold
#' for the full raster, or just areas within \code{buffer} distance of \code{spf},
#' if \code{buffer} is greater than 0.
#' 
#' If \code{buffer} is greater than 0 and \code{orig.thresh} is not NULL, \code{orig.thresh}
#' is used as the threshold for the full raster. Additionally,
#' cells within \code{buffer} distance of \code{spf} which are greater than the
#' calculated min.cell value are also set to 1. If a given \code{orig.thresh} is lower 
#' than the calculated min.cell value, \code{buffer} will have no impact on the output.
#' 
#' The default is to return the raster only. If \code{return.thresh = TRUE}, the function will
#' return a list with 3 named objects: \code{rast}, the output raster; \code{orig.thresh}, 
#' the global threshold (if used); \code{min.cell}, the minimum cell value threshold (if used).
#'
#' @param spf input spatial features (sp or sf spatial object)
#' @param rast input raster model output (values between 0 and 1)
#' @param orig.thresh numeric between 0 and 1; 'global' threshold value to apply to raster
#' @param buffer numeric; spatial buffer around spf to include in burn-in
#' @param return.thresh logical; whether to return thresholds along with raster in a list
#' @param ... Other arguments as to \code{raster::writeRaster}
#' 
#' @return RasterLayer
#' 
#' @author David Bucklin
#' 
#' @import raster
#' @importFrom methods as
#' @importFrom sf st_buffer
#'
#' @export
#'
#' @examples
#' \dontrun{
#' spf <- st_read("ambymabe/polygon_data/ambymabe.shp")
#' rast <- raster("ambymabe/grids/ambymabe_20171018_130837.tif")
#' class_burn <- nh_burn(spf, rast, 0.75, 250, filename = "model_classified.tif", datatype = "INT2U")
#' }

nh_burn <- function(spf, rast, orig.thresh = NULL, buffer = 0, return.thresh = FALSE, ...) {
  
  if (!is.null(orig.thresh) && (orig.thresh > 1 | orig.thresh < 0)) stop("orig.thresh value must be between 0 and 1.")
  if (is.null(orig.thresh)) omiss <- TRUE else omiss <- FALSE
  # not testing rast since it takes a while
  
  # handle sp/sf class
  spf1 <- tospf(spf, rast)
  sp <- spf1[[1]]
  spf <- spf1[[2]]
  rm(spf1)
  
  # cell size buffer (use/don't use?)
  csz <- max(res(rast))
  # spf <- st_buffer(spf, sqrt(2)*(csz/2)) # use this to add cellsize buffer always ensure capture of all intersecting cells (not just cell centers)
  
  # crop, extract minimum value in spf
  r1 <- crop(rast, extent(extent(spf)[1]-(csz+buffer), extent(spf)[2]+(csz+buffer), extent(spf)[3]-(csz+buffer), extent(spf)[4]+(csz+buffer)))
  # update area within features
  message("Updating feature intersection areas...")
  r2 <- gRasterize(spf, r1, value = 1) # sets = 1 areas within original features
  
  if (buffer == 0 & !is.null(orig.thresh)) {
    # case when buffer is 0 and threshold is given - simple burn in
    message("Burning in features...")
    rburn <- max(rast, extend(r2,rast), na.rm = T)
    rcl <- data.frame(from = c(NA, -1, orig.thresh-1e-10), to = c(NA, orig.thresh-1e-10, 1.1), becomes = c(NA,0,1)) # subtract very small number from threshold
    message("Reclassifying raster...")
    rast <- reclassify(rburn, rcl, ...)
    if (!return.thresh) return(rast) else return(list(rast = rast, orig.thresh = orig.thresh, min.cell = NA))
  } else {
    if (is.null(orig.thresh)) {
      orig.thresh <- NA
      # other cases where orig.thresh not given, calculate the minimum value within features
      message("Calculating minimum value in features...")
      mval <- min(unlist(extract(r1, spf)), na.rm = T)
      message("Calculated minimum value threshold = ", mval)
    } else {
      mval <- NA
    }
  }
  
  # set threshold to use
  if(!omiss) tuse <- orig.thresh else tuse <- mval
    
  if (buffer > 0) {
    message("Reclassifying buffer areas using threshold...")
    spfb <- st_buffer(spf, buffer)
    r3 <- gRasterize(spfb, r2, value = 1, background = NA) # buffer area
    r2 <- sum(r2, r3, na.rm = T) # features (2) + buffer area (1)
    # find areas in buffer area above mval, set to 1
    values(r2) <- ifelse(test = values(r1+r2) > (tuse+1), yes = 1, no = NA)
    # r2 is now areas of original features + buffer areas > tuse
    r2[is.na(r2)] <- 0
    rast <- mask(extend(r2, rast, value = 0), rast, ...)
    if (!return.thresh) return(rast) else return(list(rast = rast, orig.thresh = orig.thresh, min.cell = mval))
  } else {
    message("Reclassifying raster...")
    # reclassify
    rcl <- data.frame(from = c(NA, -1, tuse-1e-10), to = c(NA, tuse-1e-10, 1.1), becomes = c(NA,0,1)) # subtract very small number from threshold
    rburn <- max(rast, extend(r2,rast), na.rm = T)
    rast <- reclassify(rburn, rcl, ...)
    if (!return.thresh) return(rast) else return(list(rast = rast, orig.thresh = orig.thresh, min.cell = mval))
  }

  # old code  
  # if (orig.thresh < mval) {
  #   message("Original threshold is lower than minimum feature value. Reclassifying with original threshold...")
  #   # return thresholded rast, no changes needed
  #   # all areas covered by mval are already covered by orig.thresh
  #   rcl <- data.frame(from = c(NA, -1, orig.thresh-1e-10), to = c(NA, orig.thresh-1e-10, 1.1), becomes = c(NA,0,1)) # subtract very small number from threshold
  #   rast <- reclassify(rast, rcl)
  #   if (!return.thresh) return(rast) else return(list(rast = rast, orig.thresh = orig.thresh, min.cell = NA))
  # }
  # 
  # # update within buffer
  # if (buffer > 0 & orig.thresh >= mval) {
  #   # message("Thresholding buffer areas using minimum feature value...")
  #   spfb <- st_buffer(spf, buffer)
  #   r3 <- rasterize(spfb, r2, field = 1, background = NA) # buffer area
  #   r2 <- sum(r2, r3, na.rm = T) # features (2) + buffer area (1)
  #   # find areas in buffer area above mval, set to 1
  #   values(r2) <- ifelse(test = values(r1+r2) > (mval+1), yes = 1, no = NA)
  #   # r2 is now areas of original features + buffer areas > mval
  #   if (omiss) {
  #     r2[is.na(r2)] <- 0
  #     rast <- mask(extend(r2, rast, value = 0), rast)
  #     if (!return.thresh) return(rast) else return(list(rast = rast, orig.thresh = NA, min.cell = mval))
  #   }
  #   message("Burning in features...")
  #   rburn <- max(rast, extend(r2,rast), na.rm = T)
  # } else {
  #   rburn <- rast
  # }
  # 
  # message("Reclassifying raster...")
  # # reclassify
  # rcl <- data.frame(from = c(NA, -1, orig.thresh-1e-10), to = c(NA, orig.thresh-1e-10, 1.1), becomes = c(NA,0,1)) # subtract very small number from threshold
  # rast <- reclassify(rburn, rcl)
  # if (!return.thresh) return(rast) else return(list(rast = rast, orig.thresh = orig.thresh, min.cell = mval))
  
}
