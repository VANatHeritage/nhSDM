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
#' If \code{rast.out} is not specified, the raster will remain in temp folder.
#'
#' @param spf input spatial features (model predictions; sp or sf spatial object)
#' @param rast input raster template
#' @param pred.vals prediction values for spf features
#' @param buffer numeric; spatial buffer around spf to include in burn-in
#' @param rast.out Output raster file name (with file extension)
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
#' # rasterize order is taken from row.names - higher row.names values have preference
#' spf <- spf[order(spf$strord, decreasing = F),]
#' row.names(spf) <- 1:length(spf$geometry)
#' rast <- raster("template.tif")
#' 
#' # rasterize
#' bla <- nh_rasterize(spf, rast, pred.vals = spf$prbblty, buffer = spf$strord*15)
#' }

nh_rasterize <- function(spf, rast, pred.vals, buffer = 0, rast.out = NULL) {
  
  spf$orig_rn <- row.names(spf)
  
  # handle sp/sf class
  spf1 <- tospf(spf, rastproj = rast)
  sp <- spf1[[1]]
  spf <- spf1[[2]]
  rm(spf1)
  
  # assign pred vals, buffer
  spf$pred <- pred.vals
  
  if (length(buffer) == 1 && buffer != 0) {
    spf$buffer <- rep(buffer, length(spf$geometry))
    buffer <- TRUE
  } else if (length(buffer) > 1) {
    spf$buffer <- buffer
    buffer <- TRUE
  } else {
    buffer <- FALSE
  }
  
  message("Prepping raster...")
  # temp names
  tmp <- gsub(".grd", "", rasterTmpFile())
  if (!is.null(rast.out)) tmpr <- rast.out else tmpr <- paste0(tmp, ".tif")
  tmpshp <- paste0(tmp, ".shp")
  values(rast) <- NA
  writeRaster(rast, tmpr, overwrite = TRUE)
  
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
    gdal_rasterize(tmpshp, tmpr, a = "pred")
  }, finally = {
    unlink(tmpshp)
  })
  rout <- raster(tmpr)
  
  if (!is.null(rast.out)) unlink(paste0(tmp,"*"))
  
  return(rout)
}