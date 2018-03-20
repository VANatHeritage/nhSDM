# nh_burn

#' Add areas represented by features ('burn-in') to a binary raster SDM output
#' 
#' Takes spatial features, a continuous raster (values between 0 and 1), and 
#' returns a classified binary (0/1) raster. Areas greater than \code{orig.thresh} are
#' set to 1, all others 0. Areas intersecting \code{spf} features
#' are assigned a value of 1 in the returned classified raster. See details
#' for default calculation of \code{orig.thresh} and usage of \code{buffer}.
#' 
#' A minimum cell value (min.cell) across all \code{spf} is calculated.
#' If \code{orig.thresh} is NULL, the min.cell value will be used as a threshold
#' for the full raster.
#' 
#' If \code{buffer} is greater than 0 \code{orig.thresh} is not NULL, \code{orig.thresh}
#' is used as the threshold for the full raster. Additionally,
#' cells within \code{buffer} distance of \code{spf} which are greater than min.cell
#' are also set to 1. As a result, if a given \code{orig.thresh} is lower 
#' than the calculated min.cell value, \code{buffer}
#' distance will have no impact on the output.

#' @param spf input spatial features (sp or sf spatial object)
#' @param rast input raster model output (values between 0 and 1)
#' @param orig.thresh numeric between 0 and 1; threshold value to apply to raster
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
  r2 <- rasterize(as(spf, "Spatial"), r1, field = 1, update = TRUE) # sets = 1 areas within original features + cellsz buffer
  
  if (buffer == 0 & !is.null(orig.thresh)) {
    message("Burning in features...")
    # case when buffer is 0 and threshold is given - simple burn in
    r2[r2 != 1] <- NA
    rout <- max(rast, extend(r2,rast), na.rm = T)
    rout[rout >= orig.thresh] <- 1
    rout[rout != 1] <- 0
    return(rout)
  } else {
    # other cases, calculated the minimum value within features
    message("Calculating minimum value in features...")
    mval <- min(unlist(extract(r1, spf)), na.rm = T)
  }
  
  if (is.null(orig.thresh)) {
    orig.thresh <- mval
  } else if (orig.thresh < mval) {
    message("Original threshold is lower than minimum feature value. Returning classified raster...")
    # return thresholded rast, no changes needed
    # all areas covered by mval are already covered by orig.thresh
    rast[rast >= orig.thresh] <- 1
    rast[rast != 1] <- 0
    return(rast)
  }
  
  # update within buffer
  if (buffer > 0) {
    message("Thresholding buffer areas using value = ", mval)
    spfb <- st_buffer(spf, buffer)
    r3 <- rasterize(as(spfb, "Spatial"), r2, field = 1, background = NA)
    # find areas in buffer area above mval, set to 1
    r4 <- r2 + r3
    r4[r4 >= (mval+1)] <- 1
    r4[r4 != 1] <- NA
    r2 <- max(r2,r4, na.rm = T)
  }
  
  message("Combining thresholded and burn-in areas...")
  rout <- max(rast, extend(r2,rast), na.rm = T)
  
  message("Reclassifying raster...")
  # reclassify
  rout[rout >= orig.thresh] <- 1
  rout[rout != 1] <- 0
  return(rout)
  
}