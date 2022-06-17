# nh_map

#' Make a suitable/unsuitable binary map, given a continuous SDM output or raster template
#' 
#' This function performs a reclassification of a continuous map to a binary one,
#' and then performs additional changes to the binary map, depending on inputs. If
#' no continuous map is to be used, then a template raster should be given to \code{rast}.
#' 
#' Reclassification (when \code{thresh} is given) is always performed first. 
#' After this, optional masks are applied where areas within features (\code{feature.mask}) 
#' and not NoData (\code{raster.mask}) are kept (value = 1). All other cells are given (value = 0).
#' 
#' Patches (contiguous groups of cells = 1) are then be dropped by specifying a 
#' minimum patch size to \code{patch.drop} in area units of the raster's CRS 
#' (see nhSDM::nh_patchdrop). 
#' 
#' The final step is to 'burn in' (rasterize with value = 1) features 
#' from \code{feature.occ} and \code{feature.burn}. Because it is the last step, it is not
#' affected by masks or patch dropping.
#' 
#' @param rast input raster model output (values between 0 and 1); if no model, a template raster must be provided
#' @param thresh numeric between 0 and 1; 'global' threshold value to apply to rast
#' @param feature.occ input feature occurrences, to burn-in to map (sp or sf spatial object)
#' @param feature.burn input feature(s), to burn-in to map (sp or sf spatial object)
#' @param feature.mask input feature mask to apply to map (sp or sf spatial object)
#' @param raster.mask input raster mask to apply to map (should match \code{rast} resolution)
#' @param patch.drop patch size (in \code{rast} area units) to remove from map
#' 
#' @return SpatRaster
#' 
#' @author David Bucklin
#' 
#' @import terra
#' @import sf
#'
#' @export
#'
#' @examples
#' \dontrun{
#' spf <- sf::st_read("_data/occurrence/ambymabe.shp")
#' rast <- terra::rast("_data/species/ambymabe/outputs/model_predictions/ambymabe_20171018_130837.tif")
#' feature.burn <- st_buffer(spf, 250)
#' feature.mask <- st_buffer(spf, 50000)
#' raster.mask <- rast("D:/PSH/inputs/species_masks/no_med_hi_dev/no_med_hi_dev.tif")
#' map <- nh_map(rast, thresh=0.75, feature.occ=spf, feature.burn=feature.burn, 
#'    feature.mask=feature.mask, raster.mask=raster.mask, NULL)
#' }

nh_map <- function(rast, thresh = NULL, feature.occ = NULL, feature.burn = NULL, feature.mask = NULL, 
                   raster.mask = NULL, patch.drop = NULL) {
  run.pd <- TRUE
  # r is working raster (values always 1/0; NAs reclassed)
  if (!is.null(thresh)) {
    message("Reclassifying raster...")
    if (thresh > 1 | thresh < 0) stop("thresh value must be between 0 and 1.")
    r <- rast
    rcl <- data.frame(from = c(NA, -1, thresh-1e-10), to = c(NA, thresh-1e-10, 1.1), becomes = c(0,0,1)) # subtract very small number from threshold
    r <- classify(r, rcl)
  } else {
    r <- rast(rast)
    values(r) <- 0
    # no need to run patch drop if this is starting point
    run.pd <- FALSE
  }

  # feature mask
  if (!is.null(feature.mask)) {
    message("Masking by feature mask...")
    feature.mask <- tospf(feature.mask, r)[[2]]
    r <- mask(r, vect(feature.mask), updatevalue=0)
  }

  # raster mask
  if (!is.null(raster.mask)) {
    message("Masking by raster mask...")
    r <- crop(r, raster.mask)
    msk <- crop(raster.mask, r)
    r <- mask(r, msk, updatevalue = 0)
    # this would re-trigger patch drop
    run.pd <- TRUE
  }

  # patch size filter
  if (!is.null(patch.drop) & run.pd) {
    px <- res(r)[1] * res(r)[2]
    if (!patch.drop < px) {
      r <- nh_patchdrop(rast = r, min.patch = patch.drop, updatevalue = 0)
    } else {
      message("patch.drop size smaller than raster cell resolution, skipping patch removal...")
    }
  }
  
  # burn-in features
  # this collects all features into one sf-object (does not union). 
  ls.burn <- list()
  if (!is.null(feature.occ)) ls.burn[[(length(ls.burn)+1)]] <- st_sf(a = 1, geometry = tospf(feature.occ, rast)[[2]]$geometry)
  if (!is.null(feature.burn)) ls.burn[[(length(ls.burn)+1)]] <- st_sf(a = 1, geometry = tospf(feature.burn, rast)[[2]]$geometry)
  if (length(ls.burn) > 0) {
    message("Burning in features...")
    burn <- st_cast(do.call(rbind, ls.burn), "MULTIPOLYGON")
    r <- extend(r, burn)
    burnr <- rasterize(vect(burn), r, background = 0)
    r <- max(r, burnr)
  }
  
  return(r)
}
