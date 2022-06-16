# nh_sample

#' Create points in features in reference raster cells
#'
#' For each spatial feature, a given number (\code{num_samps})
#' of points are created in cells that the feature intersects.
#'
#' \code{num.samps} can be a a proportion (a decimal value < 1), single integer,
#' or vector of integers equal to length of spf indicated the number of samples
#' to take from each feature. If left NULL, \code{num.samps} will
#' be set to the number of cells [n] intersecting the feature. If a proportion is given
#' (e.g., 0.5), than [n * num.samps] will be returned. If a single integer is given,
#' \code{num.samps} points will be sampled in each feature.
#'
#' When \code{replace = FALSE} and \code{force.min = FALSE} (defaults),
#' each cell can only contain one point (across
#' each feature, and the entire returned set of points). In this case, when the number
#' of samples points to create exceeds the number of unique cells intersected by a given
#' feature, the number of samples for that feature equals the number of cells.
#' If \code{replace = TRUE}, sampling is done with replacement and duplicates
#' may be taken. The special case \code{replace = FALSE} and \code{force.min = TRUE}
#' will always return \code{num.samps} per feature. It only produces duplicates
#' if \code{num.samps} exceeds the number of cells intersecting the feature,
#' in which case it will replicate the samples until \code{num.samps} is reached.
#'
#' If CRS do not match, the features will be transformed to the CRS
#' of the raster.
#' 
#' A column `feat.id` is added to the output point features to indicate the row number
#' of the feature that the point was generated within.

#' @param spf input spatial features (sp or sf spatial object)
#' @param rast raster dataset with extent overlapping spf
#' @param num.samps number of samples to create in each feature (see details)
#' @param replace whether to sample with or without replacement
#' @param force.min whether to force \code{num.samps} points in features, even if they are duplicates
#' 
#' @return sp or sf object (points)
#' 
#' @author David Bucklin
#'
#' @importFrom sf st_crs st_transform st_cast st_buffer st_as_sf st_intersection st_geometry_type st_length st_area st_centroid st_crs<-
#' @importFrom methods as
#' @import terra
#' @export
#'
#' @examples
#' \dontrun{
#' spf <- sf::st_read("_data/occurrence/ambymabe.shp")
#' rast <- terra::rast("_data/species/ambymabe/outputs/model_predictions/ambymabe_20171018_130837.tif")
#' spf$num.samps <- sample(1:5, nrow(spf), replace = T)
#' spf.samps <- nh_sample(spf, rast, num.samps =spf$num.samps, replace = F, force.min = T)
#' nrow(spf.samps) == sum(spf$num.samps)  # Should be TRUE when force.min = T.
#' }

nh_sample <- function(spf, rast, num.samps = NULL, replace = FALSE, force.min = FALSE) {
  
  message("Generating samples...")
  # handle sp/sf class
  spf1 <- nhSDM:::tospf(spf, rast)
  sp <- spf1[[1]]
  spf <- spf1[[2]]
  rm(spf1)

  if (!is.null(num.samps)) {
    if (length(num.samps) == 1) {
      num.samps <- rep(num.samps, length(spf$geometry))
    } else if (length(num.samps) != length(spf$geometry)) {
      stop("num.samps must a vector of length one, or length equal to spf.")
    }
  }

  # make seq raster
  r1 <- crop(rast, st_sf(st_buffer(st_as_sfc(st_bbox(spf)), dist = res(rast)[1]))) # crop to 1-cell buffered extent of spf
  r2 <- r1
  values(r2) <- 1:ncell(r1)
  r2 <- mask(r2, r1)
  names(r2) <- "seqid"
  rm(r1)

  # loop over polygons
  for (i in 1:nrow(spf)) {
    suppressWarnings(rm(a2s))
    p <- spf[i,]
    # Dump to single-part features
    pall <- st_cast(p)
  
    # Loop over each feature
    for (xp in 1:length(pall$geometry)) {
      pxp <- pall[xp,]
      gb <- st_buffer(pxp, dist = res(rast)[1])
      try(a2 <- as.polygons(crop(r2, gb)), silent = T)
      if (!exists("a2")) message("Polygon at index value [", i , ".", xp , "] does not intersect raster. Skipping...")
      if (!exists("a2")) next
      a2 <- st_as_sf(a2)
      row.names(a2) <- paste0(row.names(a2), ".", xp)
      a2s1 <- suppressWarnings(st_intersection(a2, pxp))
      rm(a2)
      # rbind by original feature
      if (!exists("a2s")) a2s <- a2s1 else a2s <- rbind(a2s, a2s1)
    }
    if (!exists("a2s")) next

    # score features relative to length/area (used for sampling probability)
    if (all(st_geometry_type(a2s) == "LINESTRING")) a2s$score <- as.numeric(st_length(a2s)) else
      a2s$score <- as.numeric(st_area(a2s))
    suppressWarnings(pts <- st_centroid(a2s)) #not sure about weird polygon centroids...
    # select samples
    if (is.null(num.samps)) {
      pts.s <- pts
    } else {
      ns <- num.samps[i]
      if (ns < 1) ns <- ceiling(length(pts$geometry) * ns) # case of fraction
      if (!replace) {
        o.ns <- ns
        if (ns > length(pts$geometry)) ns <- length(pts$geometry) # case of more samples than available
        pts.s <- pts[row.names(pts) %in% sample(row.names(pts), size = ns, prob = a2s$score, replace = FALSE),]
        if (length(pts.s$geometry) < o.ns & force.min) {
          pts.o <- pts.s
          samp.rm <- sample(1:length(pts.s$geometry))
          while(length(pts.s$geometry) < o.ns) pts.s <- rbind(pts.s, pts.o)
          row.names(pts.s) <- 1:length(pts.s$geometry)
          if (length(pts.s$geometry) != o.ns) pts.s <- pts.s[!row.names(pts.s) %in% samp.rm[1:(length(pts.s$geometry)-o.ns)],]
        }
      } else {
        pts.s <- pts[c(sample(row.names(pts), size = ns, prob = pts$score , replace = TRUE)),]
      }
      row.names(pts.s) <- 1:length(pts.s$geometry)
    }
    pts.out <- pts.s
    row.names(pts.out) <- paste0(row.names(p),".",row.names(pts.s))
    pts.out$feat.id <- rep(as.character(row.names(p)), length(pts.s$geometry))

    # rbind pts
    if (!exists("sp.out")) sp.out <- pts.out else sp.out <- rbind(sp.out, pts.out)
  }
  # remove duplicates
  if (!replace & !force.min) sp.out <- sp.out[!duplicated(extract(r2, sp.out)),]
  if (sp) return(as(sp.out,Class = "Spatial")) else return(sp.out)
}
