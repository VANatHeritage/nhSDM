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
#' a buffer. If rasterizing lines and a one-cell width is desired,
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
#' @param pred.vals prediction values for spf features
#' @param buffer numeric; spatial buffer around spf to include in burn-in
#' @param priority numeric; vector of priority values for feature rasterization (higher values have priority)
#' @param rast.out Output raster file name (with file extension)
#' @param ... Additional arguments to writeRaster (e.g. overwrite)
#' 
#' @return RasterLayer
#' 
#' @author David Bucklin
#'
#' @import raster
#' @importFrom methods as
#' @importFrom sf st_buffer st_write
#' @importFrom gdalUtils gdal_rasterize
#'
#' @export
#'
#' @examples
#' \dontrun{
#' spf <- st_read("acipoxyr/shapefiles/acipoxyr_20180105_133929_results.shp")
#' rast <- raster("template.tif")
#' 
#' # rasterize
#' bla <- nh_rasterize(spf, rast, pred.vals = spf$prbblty, buffer = spf$strord*15)
#' }

nh_rasterize <- function(spf, rast, pred.vals, buffer = 0, priority = NULL, rast.out = NULL, ...) {
  
  # handle sp/sf class
  spf1 <- tospf(spf, rastproj = rast)
  sp <- spf1[[1]]
  spf <- spf1[[2]]
  rm(spf1)
  
  # assign pred vals, buffer
  spf$pred <- pred.vals
  
  if (length(buffer) == 1 && buffer != 0) {
    spf$buffer <- rep(buffer, nrow(spf))
    buffer <- TRUE
  } else if (length(buffer) > 1) {
    spf$buffer <- buffer
    buffer <- TRUE
  } else {
    buffer <- FALSE
  }
  
  # order 
  if (!is.null(priority)) {
    # order using priority vector
    spf <- spf[order(priority, decreasing = F, na.last = F),]
  } else {
    # order using pred
    spf <- spf[order(spf$pred, decreasing = F, na.last = F),]
  }
  row.names(spf) <- 1:nrow(spf)
  spf$orig_rn <- row.names(spf)
  
  message("Prepping raster...")
  # temp names
  tmp <- gsub(".grd", "", rasterTmpFile())
  if (!is.null(rast.out)) tmpr <- rast.out else tmpr <- paste0(tmp, ".tif")
  tmpshp <- paste0(tmp, ".shp")
  values(rast) <- NA
  writeRaster(rast, tmpr, ...)
  
  # buffer
  if (buffer) {
    message("Buffering features...")
    s0 <- spf[c("orig_rn", "pred")][0,]
    for (so in sort(unique(spf$buffer))) {
      s2 <- st_buffer(spf[c("orig_rn", "pred")][spf$buffer == so,], so)
      s0 <- rbind(s0, s2)
    }
  } else {
    s0 <- spf[c("orig_rn", "pred")]
  }
  # order by original row names
  s1 <- s0[order(s0$orig_rn, decreasing = F),]
  row.names(s1) <- s1$orig_rn           
  
  message("Rasterizing...")
  tryCatch({
    st_write(s1, tmpshp, delete_layer = T, quiet = T)
    gdalUtils::gdal_rasterize(tmpshp, tmpr, a = "pred")
  }, finally = {
    del <- paste(dirname(tmp), list.files(dirname(tmp), pattern = paste0(basename(tmp),"*")), sep = "/")
    unlink(del[!grepl(".tif$", del)])
  })
  rout <- raster(tmpr)
  return(rout)
}
