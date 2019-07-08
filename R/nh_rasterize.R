# nh_rasterize

#' Convert vector format SDM prections to raster format
#' 
#' Requires spatial features, a template raster, and prediction
#' values vector. The buffer is optional; when provided, the
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
#' @param rast.out Output raster file name (with file extension)
#' @param ... Additional arguments to writeRaster (e.g. overwrite)
#' 
#' @return RasterLayer
#' 
#' @author David Bucklin
#'
#' @import raster
#' @importFrom methods as
#' @importFrom sf st_buffer st_cast
#' @importFrom fasterize fasterize
#'
#' @export
#'
#' @examples
#' \dontrun{
#' spf <- st_read("acipoxyr/shapefiles/acipoxyr_20180105_133929_results.shp")
#' rast <- raster("template.tif")
#' 
#' # rasterize
#' bla <- nh_rasterize(spf, rast, pred.vals = "prbblty", buffer = spf$strord*15)
#' }

nh_rasterize <- function(spf, rast, pred.vals, buffer = NULL, priority = NULL, rast.out = NULL, ...) {
  
  if (!pred.vals %in% names(spf)) stop(paste0("Column `", pred.vals, "` not found."))
  if (!is.null(priority) && !priority %in% names(spf)) stop(paste0("Column `", priority, "` not found."))
  
  # handle sp/sf class
  spf1 <- tospf(spf, rastproj = rast)
  sp <- spf1[[1]]
  spf <- spf1[[2]]
  rm(spf1)
  
  if (!is.null(buffer)) {
    spf <- st_cast(st_buffer(spf, buffer), "MULTIPOLYGON")
  } else {
    # buffer linstring to at least one cell
    if (grepl("POINT|LINESTRING", st_geometry_type(spf)[1]))
      spf <- st_cast(st_buffer(spf, res(rast)[1]*sqrt(2) / 2), "MULTIPOLYGON")
  }
  
  message("Rasterizing...")
  rast <- extend(crop(rast, extent(spf)), extent(spf))
  if (!is.null(priority)) {
    spf <- spf[order(as.data.frame(spf[priority])[,1], decreasing = F),]
    rout <- fasterize(spf, rast, field = pred.vals, fun = "last", background = NA)
  } else {
    rout <- fasterize(spf, rast, field = pred.vals, fun = "max", background = NA)
  }
  
  if (!is.null(rast.out)) {
    writeRaster(rout, rast.out, ...)
    rout <- raster(rast.out)
  }
  return(rout)
}
