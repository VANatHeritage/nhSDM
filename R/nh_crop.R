# nh_crop

#' Crop extra rows/columns from all sides of a raster
#' 
#' This function will reduce the extent of a raster, by removing rows/columns
#' from all sides if they do not have any non-NA cells with a value greater than zero.
#' One buffer row/column is left on each side.
#' 
#' @param rast input raster
#' 
#' @return SpatRaster
#' 
#' @author David Bucklin
#' 
#' @import terra
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rast <- terra::rast("_data/species/ambymabe/outputs/model_predictions/ambymabe_20171018_130837.tif")
#' values(rast) <- ifelse(values(rast) > 0.3, 1, NA)
#' rast.crop <- nh_crop(rast)
#' 
#' # should be TRUE
#' sum(values(rast), na.rm=T) == sum(values(rast.crop), na.rm=T)
#' }

nh_crop <- function(rast) {
  message("Cropping extra rows/columns...")
  system.time(v <- values(rast))
  vl <- length(v)
  row.all <- seq(1, vl, by=ncol(rast))
  # top
  for (i in 1:length(row.all)) {
    x <- row.all[i]
    xe <- x + (ncol(rast)-1)
    if (sum(v[x:xe], na.rm = T) > 0) break
  }
  if (i == 1) y0 <- 1 else y0 <- i-1
  # bottom
  for (i in length(row.all):1) {
    x <- row.all[i]
    xe <- x + (ncol(rast)-1)
    if (sum(v[x:xe], na.rm = T) > 0) break
  }
  if (i == length(row.all)) y1 <- i else y1 <- i+1
  # left
  col.all <- 1:ncol(rast)
  for (i in 1:length(col.all)) {
    x <- col.all[i]
    x <- seq(x, vl, by = ncol(rast))
    if (sum(v[x], na.rm = T) > 0) break
  }
  if (i == 1) x0 <- 1 else x0 <- i-1
  # right
  for (i in length(col.all):1) {
    x <- col.all[i]
    x <- seq(x, vl, by = ncol(rast))
    if (sum(v[x], na.rm = T) > 0) break
  }
  if (i == length(col.all)) x1 <- i else x1 <- i+1
  
  # cr <- extent(rast, y0, y1, x0, x1)
  rout <- rast[y0:y1, x0:x1, drop=FALSE]
  # rout <- crop(rast, cr)
  return(rout)
}
