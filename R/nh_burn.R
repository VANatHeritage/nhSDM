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
#' If \code{buffer} is greater than 0 \code{orig.thresh} is not NULL, \code{orig.thresh}
#' is used as the threshold for the full raster. Additionally,
#' cells within \code{buffer} distance of \code{spf} which are greater than min.cell
#' are also set to 1. If a given \code{orig.thresh} is lower 
#' than the calculated min.cell value, \code{buffer}
#' distance will have no impact on the output.
#'
#' @param spf input spatial features (sp or sf spatial object)
#' @param rast input raster model output (values between 0 and 1)
#' @param orig.thresh numeric between 0 and 1; threshold value to apply to full raster
#' @param buffer numeric; spatial buffer around spf to include in burn-in
#' @param ... additional parameters to \code{raster::writeRaster}
#' 
#' @return RasterLayer
#' 
#' @author David Bucklin
#'
#' @importFrom methods as
#' @import sf
#' @import raster
#'
#' @export
#'
#' @examples
#' \dontrun{
#' spf <- st_read("ambymabe/polygon_data/ambymabe.shp")
#' rast <- raster("ambymabe/grids/ambymabe_20171018_130837.tif")
#' class_burn <- nh_burn(spf, rast, 0.75, 250, filename = "model_classified.tif", datatype = "INT2U")
#' }

nh_burn <- function(spf, rast, orig.thresh = NULL, buffer = 0, ...) {
  
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
  r2 <- rasterize(spf, r1, field = 1) # sets = 1 areas within original features
  
  if (buffer == 0 & !is.null(orig.thresh)) {
    # case when buffer is 0 and threshold is given - simple burn in
    message("Burning in features...")
    rburn <- max(rast, extend(r2,rast), na.rm = T)
    rcl <- data.frame(from = c(NA, -1, orig.thresh-1e-10), to = c(NA, orig.thresh-1e-10, 1.1), becomes = c(NA,0,1)) # subtract very small number from threshold
    message("Reclassifying raster...")
    rast <- reclassify(rburn, rcl, ...)
    return(rast)
  } else {
    # other cases, calculated the minimum value within features
    message("Calculating minimum value in features...")
    mval <- min(unlist(extract(r1, spf)), na.rm = T)
    message("Calculated minimum value threshold = ", mval)
  }
  
  if (is.null(orig.thresh)) {
    orig.thresh <- mval
    # this continues with orig.thresh == mval. Buffer not used.
  } else if (orig.thresh < mval) {
    message("Original threshold is lower than minimum feature value. Reclassifying with original threshold...")
    # return thresholded rast, no changes needed
    # all areas covered by mval are already covered by orig.thresh
    rcl <- data.frame(from = c(NA, -1, orig.thresh-1e-10), to = c(NA, orig.thresh-1e-10, 1.1), becomes = c(NA,0,1)) # subtract very small number from threshold
    rast <- reclassify(rast, rcl, ...)
    return(rast)
  }
  
  # update within buffer
  if (buffer > 0 & orig.thresh >= mval) {
    message("Thresholding buffer areas using minimum feature value...")
    spfb <- st_buffer(spf, buffer)
    r3 <- rasterize(spfb, r2, field = 1, background = NA)
    r2 <- sum(r2, r3, na.rm = T)
    # find areas in buffer area above mval, set to 1
    values(r2) <- ifelse(test = values(r1+r3) > (mval+1), yes = 1, no = NA)
    # r2 is now areas of original features + buffer areas > mval
    if (omiss) {
      r2[is.na(r2)] <- 0
      rast <- mask(extend(r2, rast, value = 0), rast)
      return(rast)
    }
    message("Burning in features...")
    rburn <- max(rast, extend(r2,rast), na.rm = T)
  } else {
    rburn <- rast
  }
  
  message("Reclassifying raster...")
  # reclassify
  rcl <- data.frame(from = c(NA, -1, orig.thresh-1e-10), to = c(NA, orig.thresh-1e-10, 1.1), becomes = c(NA,0,1)) # subtract very small number from threshold
  rast <- reclassify(rburn, rcl, ...)
  return(rast)
  
}