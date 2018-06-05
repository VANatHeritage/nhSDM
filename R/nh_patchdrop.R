# nh_patchdrop

#' Remove patches (contiguous cells representing habitat area) smaller than a given patch size from thresholded SDM output
#' 
#' Takes a binary/thresholded raster (values either NA, 0, or 1), and 
#' returns a binary raster. Clumps of contiguous cells with the value (1) are
#' given a value of 0 (if 0 values are present in the input 'rast'); otherwise, they
#' are given an NA value.
#' 
#' If 'spf' is given, the smallest feature's area will be used to derive
#' a 'min.patch' value, and any given 'min.patch' is ignored. If 'spf' is 
#' not given, a 'min.patch' value must be given, in units of the input raster.
#'
#' @param spf input spatial features (sp or sf spatial object)
#' @param rast input binary raster output (values either 0 or 1)
#' @param min.patch area of minimum patch size, in area units used in \code{rast}
#' @param directions Integer. Which cells are considered adjacent? Should be 8 (default; Queen's case) or 4 (Rook's case)
#' @param ... Other arguments as to \code{raster::writeRaster}
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
#' spf <- st_read("D:/SDM/Tobacco/inputs/species/ambymabe/polygon_data/ambymabe.shp")
#' rast <- raster("E:/projects/SDM_ancilliary/Hypergrid/VA_methods/raster/screen_lpsh/thumb/ambymabe.tif")
#' 
#' # use minimum patch size from presence features
#' rast_contig_minpres <- nh_patchdrop(spf, rast, filename = "bla.tif", datatype = "INT2U")
#' 
#' # use a minimum patch size of 10000 (in units from 'rast')
#' rast_contig_10km <- nh_patchdrop(rast = rast, min.patch = 10000, filename = "bla.tif", datatype = "INT2U")
#' }

nh_patchdrop <- function(spf = NULL, rast, min.patch = NULL, directions = 8, ...) {
  
  if (!"igraph" %in% installed.packages()) stop("The package 'igraph' is required for this function.")
  if (is.null(spf) & is.null(min.patch)) stop("Must provide either spatial features (spf) or min.patch size.")

  # handle sp/sf class
  if (!is.null(spf)) {
    spf1 <- tospf(spf, rast)
    sp <- spf1[[1]]
    spf <- spf1[[2]]
    rm(spf1)
    
    min.patch <- min(as.numeric(st_area(spf)))
  }
  
  # minimum number of cells
  cell1 <- min.patch /(res(rast)[1] * res(rast)[2])
  if (as.integer(cell1) == cell1) cells <- cell1 else cells <- floor(cell1) + 1

  if (cells == 1) {
    message("Minimum patch size = 1 cell. Returning original raster...")
    return(rast)
  } else {
    message("Minimum patch size = ", cells, " cells...")
  }
  
  # new value for excluded patches
  if (0 %in% unique(rast)) upd <- 0 else upd <- NA
  
  # run clump (extend to break grouping from one x edge to another)
  r2 <- crop(clump(extend(rast, y=c(1,1)), directions = directions), rast)

  # figure out what to drop
  drop <- as.data.frame(freq(r2))
  if (min(drop$count) >= cells) {
    message("No patches smaller than minimum patch size. Returning original raster...")
    return(rast)
  }
  drop <- drop[drop$count < cells,]
  
  r2[r2 %in% drop$value] <- NA

  rout <- mask(rast, r2, updatevalue = upd, ...)
  return(rout)
}