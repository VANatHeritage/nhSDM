# nh_best

#' Extract areas with best (highest value) as polygons from an SDM prediction raster, with optional feature mask
#' 
#' If \code{spf} is given, areas intersecting these features (plus a buffer, if \code{min.dist} is specified) 
#' will not be included in the returned polygons.
#' 
#' \code{top.percent} and \code{min.size} will both be derived from \code{spf} if they are null
#' and \code{spf} is given. In this case, \code{top.percent} will be set to be equal to \code{spf}'s cell coverage 
#' (prevalence relative to non-NA areas in raster). \code{min.size} will be set to the area 
#' of the smallest feature in \code{spf}. If \code{spf} are line or point features, this will be zero. 
#' 
#' \code{rank.by} only subsets the outputs if \code{num.patches} is set. \code{rank.by}
#' must be set to \code{"value"} to return the mean prediction value of the patch.
#' 
#' For better performance, crop the raster prior to running, and/or use a very low 
#' \code{top.percent} value, to return a smaller proportion of the cells as polygons.
#' 
#' @param rast input raster model output with continuous values
#' @param spf input spatial features (sp or sf spatial object); if given, it can be used to modify areas selected from \code{rast}
#' @param top.percent numeric; percent (e.g.; 0.01 = 0.01\%) of highest cell values in raster to extract
#' @param min.size numeric; optional minimal area of an extracted polygon
#' @param min.dist numeric; optional minimal distance from spf for patches. Default is 0 (can be adjacent to 'spf')
#' @param num.patches numeric; optional number of patches to return, using rank.by criteria
#' @param rank.by character; used for ranking patches when num.patches is not null; either 'area' (default) or 'value'
#' 
#' @return sp or sf object (polygons)
#' 
#' @author David Bucklin
#'
#' @importFrom methods as
#' @importFrom sf st_buffer st_area st_length st_cast st_union st_as_sf st_sfc st_sf
#' @import terra
#'
#' @export
#'
#' @examples
#' \dontrun{
#' spf <- sf::st_read("_data/occurrence/ambymabe.shp")
#' rast <- terra::rast("_data/species/ambymabe/outputs/model_predictions/ambymabe_20171018_130837.tif")
#' rast <- terra::crop(rast, spf)
#' best_polys <- nh_best(rast, spf, top.percent = 0.01, min.size = 10000, min.dist <- 10000, 
#'    num.patches = 100, rank.by = "value")
#' }

nh_best <- function(rast, spf = NULL, top.percent = NULL, min.size = NULL, min.dist = NULL, num.patches = NULL, rank.by = 'area') {
  
  if (!rank.by %in% c("area","value")) {
    stop("rank.by must be either 'area' or 'value'.")
  }
  
  if (!is.null(spf)) {
    # handle sp/sf class
    spf1 <- tospf(spf, rast)
    sp <- spf1[[1]]
    spf <- spf1[[2]]
    rm(spf1)
    
    message("Removing areas covered by 'spf'...")
    if (!is.null(min.dist)) {
      spf.buff <- st_buffer(spf, min.dist)
      r <- mask(rast, vect(spf.buff), inverse = TRUE)
    } else {
      r <- mask(rast, vect(spf), inverse = TRUE)
    }
  } else {
    if (!is.null(min.dist)) message("Ignoring 'min.dist', as 'spf' is not specified.")
    r <- rast
  }
  
  vals <- as.matrix(r[])
  if (!is.null(top.percent)) {
    top.percent <- top.percent/100
  } else {
    if (!is.null(spf)) {
      ar <- sum(as.numeric(st_area(spf)))
      notna <- sum(!is.na(vals))
      if (ar > 0) {
          # ncells
          ncell.spf <- ar / (res(rast)[1] * res(rast)[2])
        } else {
          al <- sum(as.numeric(st_length(spf)))
          if (al > 0) ncell.spf <- al/res(rast)[1] else ncell.spf <- length(spf[[1]]) # lines/points
        }
      top.percent <- ncell.spf / notna
    } else {
      top.percent <- 0.01
    }
  }
  message("Using top.percent value of [", round(top.percent * 100, 3), " %].")
  
  # use smallest patch size from spf as min.size
  if (is.null(min.size) & !is.null(spf)) {
    min.size <- as.numeric(min(st_area(spf)))
  }
  message("Using min.size of [", min.size, "].")
  
  # calc cutoff value
  best <- quantile(vals, 1-top.percent, na.rm = T)[1]
  
  # set values below cutoff to 0, convert to polys
  r[r < best] <- NA
  # spf.poly <- rasterToPolygons(r)
  spf.poly <- st_as_sf(as.polygons(r, dissolve=F, values=T))
  names(spf.poly)[1] <- "mean_value"
  
  # dissolve adjacent, summarize value
  if (rank.by == "value") {
    spf.gr <- st_cast(st_union(spf.poly, by_feature = FALSE),"POLYGON") # dump to intersecting polygons
    spf.out <- aggregate(spf.poly, by = spf.gr, mean)
  } else {
    spf.out <- st_sfc(st_cast(st_union(spf.poly, by_feature = FALSE),"POLYGON")) # dump to intersecting polygons
    spf.out <- st_sf(geometry=spf.out)
  }
  # add columns
  spf.out$id <- row.names(spf.out)
  spf.out$area <- as.numeric(st_area(spf.out))
  
  # remove smaller than minimum size
  if (!is.null(min.size)) spf.out <- spf.out[spf.out$area >= min.size,]

  if (!is.null(num.patches)) {
    if (rank.by == "area") {
      patch.area <- sort(spf.out$area, decreasing = T)
      if (length(patch.area) > num.patches) {
        min.size <- patch.area[num.patches]
        message("Returning patches with area >= than ", num.patches, "-ranked largest patch...")
        spf.out <- spf.out[spf.out$area >= min.size,]
      }
    } else {
      mean.val <- sort(spf.out$mean_value, decreasing = T)
      if (length(mean.val) > num.patches) {
        min.val <- mean.val[num.patches]
        message("Returning patches with mean value >= than ", num.patches, "-ranked highest patch...")
        spf.out <- spf.out[spf.out$mean_value >= min.val,]
      }
    }
  }
  
  if (!is.null(spf) && sp) return(as(spf.out,Class = "Spatial")) else return(spf.out)
}
