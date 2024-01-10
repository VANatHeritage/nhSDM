# nh_burn

#' Add areas represented by features ('burn-in') to a binary raster SDM output, 
#' plus areas within a buffer distance above a threshold
#' 
#' Takes spatial features, a continuous raster (values between 0 and 1), and 
#' returns a classified binary (0/1) raster. Areas intersecting \code{spf} features
#' are assigned a value of 1 in the returned classified raster. Additionally, areas
#' within \code{buffer} distance of \code{spf} and above the threshold are set to 1. 
#' See details for default calculation of \code{orig.thresh} and usage of \code{buffer}.
#' 
#' When \code{buffer} is used or \code{orig.thresh} is not provided,
#' a minimum cell value (min.cell) of values within \code{spf} is calculated and 
#' used as the threshold.
#' 
#' The default is to return the raster only. If \code{return.thresh = TRUE}, the function will
#' return a list with 3 named objects: \code{rast}, the output raster; \code{orig.thresh}, 
#' the global threshold (if used); \code{min.cell}, the minimum cell value threshold (if used).
#'
#' @param spf input spatial features (sp or sf spatial object)
#' @param rast input raster model output (values between 0 and 1)
#' @param orig.thresh numeric between 0 and 1; threshold value to apply to raster
#' @param buffer numeric; spatial buffer distance around spf where threshold will be applied.
#' @param return.thresh logical; whether to return thresholds along with raster in a list
#' 
#' @return SpatRaster
#' 
#' @author David Bucklin
#' 
#' @import terra
#' @importFrom methods as
#' @importFrom sf st_buffer
#'
#' @export
#'
#' @examples
#' \dontrun{
#' spf <- sf::st_read("_data/occurrence/ambymabe.shp")
#' rast <- terra::rast("_data/species/ambymabe/outputs/model_predictions/ambymabe_20171018_130837.tif")
#' full_burn <- nh_burn(spf, rast, buffer = NA)
#' buff_burn <- nh_burn(spf, rast, buffer = 10000)
#' }

nh_burn <- function(spf, rast, orig.thresh = NULL, buffer = NA, return.thresh = FALSE) {
  
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
  if (!is.na(buffer)) {
    r1 <- crop(rast, ext(ext(spf)[1]-(csz+buffer), ext(spf)[2]+(csz+buffer), ext(spf)[3]-(csz+buffer), ext(spf)[4]+(csz+buffer)))
  } else { 
    r1 <- rast
  }
  
  # update area within features
  message("Updating feature intersection areas...")
  r2 <- gRasterize(spf, r1, value = 1) # sets = 1 areas within original features


  if (is.na(buffer) & !is.null(orig.thresh)) {
    # case when buffer is NA and threshold is given simple burn in
    message("Burning in features...")
    rburn <- max(rast, extend(r2,rast), na.rm = T)
    rcl <- data.frame(from = c(NA, -1, orig.thresh-1e-10), to = c(NA, orig.thresh-1e-10, 1.1), becomes = c(NA,0,1)) # subtract very small number from threshold
    message("Reclassifying raster...")
    rast <- classify(rburn, rcl)
    if (!return.thresh) return(rast) else return(list(rast = rast, orig.thresh = orig.thresh, min.cell = NA))
  } else {
    if (is.null(orig.thresh)) {
      orig.thresh <- NA
      # other cases where orig.thresh not given, calculate the minimum value within features
      message("Calculating minimum value in features...")
      mval <- min(extract(r1, vect(spf))[,2], na.rm = T)
      message("Calculated minimum value threshold = ", mval)
    } else {
      mval <- NA
    }
  }
  
  # set threshold to use
  if(!omiss) tuse <- orig.thresh else tuse <- mval
    
  if (!is.na(buffer)) {
    message("Reclassifying buffer areas using threshold...")
    spfb <- st_buffer(spf, buffer)
    r3 <- gRasterize(spfb, r2, value = 1, background = NA) # buffer area
    r2 <- sum(r2, r3, na.rm = T) # features (2) + buffer area (1)
    # find areas in buffer area above mval, set to 1
    values(r2) <- ifelse(test = values(r1+r2) > (tuse+1), yes = 1, no = NA)
    # r2 is now areas of original features + buffer areas > tuse
    r3 <- extend(r2, rast)
    r3[is.na(r3)] <- 0
    rast <- mask(r3, rast)
    if (!return.thresh) return(rast) else return(list(rast = rast, orig.thresh = orig.thresh, min.cell = mval))
  } else {
    message("Reclassifying raster...")
    # reclassify
    rcl <- data.frame(from = c(NA, -1, tuse-1e-10), to = c(NA, tuse-1e-10, 1.1), becomes = c(NA,0,1)) # subtract very small number from threshold
    rburn <- max(rast, extend(r2,rast), na.rm = T)
    rast <- classify(rburn, rcl)
    if (!return.thresh) return(rast) else return(list(rast = rast, orig.thresh = orig.thresh, min.cell = mval))
  }
  
}
