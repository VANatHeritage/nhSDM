# nh_patchdrop

#' Remove contiguous patches smaller than a given patch size from binary output
#' 
#' Takes a binary/thresholded raster (values either NA, 0, or 1), and 
#' returns a binary raster. Clumps of contiguous cells with the value (1) 
#' that are smaller than the \code{min.patch} size are given \code{updatevalue} (default = 0).
#' 
#' If 'spf' is given, the smallest feature's area will be used to derive
#' a \code{min.patch} value, and any given 'min.patch' is ignored. If 'spf' is 
#' not given, a 'min.patch' value must be given, in area units of the input raster.
#'
#' @param spf input spatial features (sp or sf spatial object)
#' @param rast input binary raster output (values either NA/0 or 1)
#' @param min.patch area of minimum patch size, in area units used in \code{rast}
#' @param directions Integer. Which cells are considered adjacent? Should be 8 (default; Queen's case) or 4 (Rook's case). From \code{terra::patches}
#' @param updatevalue Integer or NA. Value to apply to cells which do not meet the min.patch size. Default = 0.
#' 
#' @return SpatRaster
#' 
#' @author David Bucklin
#'
#' @importFrom methods as
#' @importFrom sf st_area
#' @import terra
#'
#' @export
#'
#' @examples
#' \dontrun{
#' spf <- sf::st_read("_data/occurrence/ambymabe.shp")
#' rast <- terra::rast("_data/species/ambymabe/outputs/model_predictions/ambymabe_20171018_130837.tif")
#' values(rast) <- ifelse(values(rast) > 0.9, 1, NA)
#' 
#' # use minimum patch size from presence features
#' rast_contig_minpres <- nh_patchdrop(spf, rast)
#' 
#' # use a minimum patch size of 10000 (in units from 'rast')
#' rast_contig_10km <- nh_patchdrop(rast = rast, min.patch = 10000)
#' }

nh_patchdrop <- function(spf = NULL, rast, min.patch = NULL, directions = 8, updatevalue=0) {
  # ISSUE: clump/patches functions are incredibly slow. Likely not using.
  
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
  
  # run clump
  message("Finding contiguous patches...")
  if (sum(rast[1:nrow(rast), 1], na.rm=T) > 0) {
    # extend to break grouping from one x-edge to another, when there are some values 
    # r2 <- crop(clump(extend(rast, y=c(1,1)), directions = directions), rast) 
    r2 <- crop(patches(extend(rast, y=c(1,1)), directions = directions), rast) 
  } else {
    r2 <- patches(rast, directions = directions)
  }
  
  # figure out what to drop
  drop <- as.data.frame(freq(r2))
  if (min(drop$count) >= cells) {
    message("No patches smaller than minimum patch size. Returning original raster...")
    return(rast)
  }
  drop <- drop[drop$count < cells,]
  
  r2[r2 %in% drop$value] <- NA

  rout <- mask(rast, r2, updatevalue = updatevalue)
  return(rout)
}
