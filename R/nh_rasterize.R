# nh_rasterize

#' Convert vector format SDM predictions to raster format
#' 
#' Requires spatial features with a prediction values attribute, and a 
#' template raster. The buffer is optional; when provided, the
#' buffer distance can be provided as a single numeric value or a 
#' vector matching length of \code{spf}, for a variable buffer by feature.
#' 
#' Buffer units should be given in the units of \code{rast}. Make 
#' sure to keep in mind the resolution of the raster when choosing 
#' a buffer. If rasterizing lines, and a one-cell width is desired,
#' do not use a buffer.
#' 
#' A vector can be given to `priority` for sorting prior to rasterization, 
#' where higher values have priority. When \code{priority = NULL} (default),
#' priority will be defined as the \code{pred.vals}.
#' 
#' If \code{rast.out} is not specified, the raster will remain in temp folder.
#'
#' @param spf input spatial features (model predictions; sp or sf spatial object)
#' @param rast input raster template
#' @param pred.vals column name in spf holding feature prediction values
#' @param buffer numeric (single value or vector matching length of spf); spatial buffer around spf to include in burn-in
#' @param priority column name in spf holding priority values for feature rasterization (higher values have priority)
#' @param touches from \code{terra::rasterize}: If TRUE, all cells touched by lines or polygons are affected, not just those on the line render path, 
#' or whose center point is within the polygon.
#' @param rast.out Output raster file name (with file extension)
#' @param ... Additional arguments to writeRaster (e.g. overwrite)
#' 
#' @return SpatRaster
#' 
#' @author David Bucklin
#'
#' @import terra
#' @importFrom methods as
#' @importFrom sf st_buffer st_cast
#'
#' @export
#'
#' @examples
#' \dontrun{
#' spf <- st_read("_data/species/acipoxyr/outputs/model_predictions/acipoxyr_20180105_133929_results.shp")
#' rast <- terra::rast("_data/species/ambymabe/outputs/model_predictions/ambymabe_20171018_130837.tif")
#' 
#' # rasterize
#' rast_poly <- nh_rasterize(spf, rast, pred.vals = "prbblty", priority = "strord", buffer = spf$strord*15, touches=F, 
#'    rast.out = "C:/David/scratch/nh_rasterize_poly.tif", overwrite = T)
#' rast_line <- nh_rasterize(spf, rast, pred.vals = "prbblty", priority = "strord", touches=T, 
#'    rast.out = "C:/David/scratch/nh_rasterize_line.tif", overwrite = T)
#' }

nh_rasterize <- function(spf, rast, pred.vals, buffer = NULL, priority = NULL, touches = TRUE, rast.out = NULL, ...) {
  
  if (!pred.vals %in% names(spf)) stop(paste0("Column `", pred.vals, "` not found."))
  if (!is.null(priority) && !priority %in% names(spf)) stop(paste0("Column `", priority, "` not found."))
  
  # handle sp/sf class
  spf1 <- tospf(spf, rastproj = rast)
  sp <- spf1[[1]]
  spf <- spf1[[2]]
  rm(spf1)
  
  if (!is.null(buffer)) {
    message("Buffering input features...")
    spf <- st_cast(st_buffer(spf, buffer), "MULTIPOLYGON")
  }
  
  message("Rasterizing...")
  rast <- extend(crop(rast, ext(spf)), ext(spf))
  if (!is.null(priority)) {
    spf <- spf[order(as.data.frame(spf[priority])[,1], decreasing = F),]
    rout <- rasterize(vect(spf), rast, field=pred.vals, fun="last", background = NA, touches=touches)
  } else {
    rout <- rasterize(vect(spf), rast, field=pred.vals, fun="max", background = NA, touches=touches)
  }
  
  if (!is.null(rast.out)) {
    message("Writing raster...")
    writeRaster(rout, rast.out, ...)
    rout <- rast(rast.out)
  }
  return(rout)
}
