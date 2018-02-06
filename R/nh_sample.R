# nh_sample

#' Create points in polygons in reference raster cells
#'
#' For each spatial polygon, a given number (\code{num_samps})
#' of points are created in cells that the polygon intersects.
#'
#' \code{num.samps} can be a a proportion (a decimal value < 1), single integer,
#' or vector of integers equal to length of spP indicated the number of samples
#' to take from each polygon. If left NULL, \code{num.samps} will
#' be set to the number of cells [n] intersecting the polygon. If a proportion is given
#' (e.g., 0.5), than [n * num.samps] will be returned. If a single integer is given,
#' \code{num.samps} points will be sampled in each polygon.
#'
#' When \code{replace = FALSE} and {force.min = FALSE} (defaults),
#' each cell can only contain one point (across
#' each polygon, and the entire returned set of points). In this case, when the number
#' of samples points to create exceeds the number of unique cells intersected by a given
#' polygon, the number of samples for that polygon is reduced to the number of cells.
#' If \code{replace = TRUE}, sampling is done with replacement and duplicates
#' may be taken. The special case \code{replace = FALSE} and \code{force.min = TRUE}
#' will always return \code{num.samps} per polygon. It only produces duplicates
#' if \code{num.samps} exceeds the number of cells intersecting the polygon,
#' in which case it will replicate the samples until \code{num.samps} is reached.
#'
#' If CRS do not match, the SpatialPolygons will be transformed to the CRS
#' of the raster.

#' @param spP SpatialPolygons[DataFrame]
#' @param rast raster dataset with extent overlapping spP
#' @param num.samps number of samples to create in each polygon (see details)
#' @param replace whether to sample with or without replacement
#' @param force.min whether to always create \code{num.samps} points, even if they are duplicates
#'
#' @author David Bucklin
#'
#' @importFrom raster crs crop extent ncell rasterToPolygons res values mask extract area values<-
#' @importFrom sp SpatialPointsDataFrame spTransform proj4string disaggregate
#' @importFrom rgeos gBuffer gIntersection gPointOnSurface
#' @importFrom utils installed.packages
#'
#' @export
#'
#' @examples
#' \dontrun{
#' r<-raster::raster("D:/SDM/Tobacco/env_vars/Tobacco/AnnMnTemp.tif")
#' spP <- rgdal::readOGR("D:/SDM/Tobacco/inputs/species/ambymabe/polygon_data", "ambymabe")
#' spg <- nh_sample(spP, 1000, TRUE)
#' rast <- r
#' num.samps <- 10
#' replace <- FALSE
#' force.min <- FALSE
#' sp.sing <- nh_sample(spP, r, num.samps)
#' sp.mult <- nh_sample(spg, r, num.samps)
#' }

nh_sample <- function(spP, rast, num.samps = NULL, replace = FALSE, force.min = FALSE) {

  # fixed number for num.samps
  if (!is.null(num.samps)) {
    if (length(num.samps) == 1) {
      num.samps <- rep(num.samps, length(spP))
    } else if (length(num.samps) != length(spP)) {
      stop("num.samps must a vector of length one, or length equal to spP.")
    }
  }
  # transform if necessary
  if (proj4string(spP) != proj4string(rast)) spP <- spTransform(spP, CRSobj = crs(rast))
  # seq raster
  r1 <- crop(rast, extent(spP))
  r2 <- r1
  values(r2) <- 1:ncell(r1)
  r2 <- mask(r2, r1)
  names(r2) <- "seqid"
  rm(r1)

  # loop over polygons
  for (i in 1:length(spP)) {
    suppressWarnings(rm(a2s))
    p <- spP[i,]
    pall <- disaggregate(p)

    for (xp in 1:length(pall)) {
      pxp <- pall[xp,]
      gb <- gBuffer(pxp, byid = FALSE, width = res(r2)[1])
      try(a2 <- rasterToPolygons(crop(r2, gb)), silent = TRUE)
      if (!exists("a2")) message("Polygon at index value [", i , ".", xp , "] does not intersect raster. Skipping...")
      if (!exists("a2")) next
      row.names(a2) <- paste0(row.names(a2), ".", xp)
      a2s1 <- gIntersection(a2, pxp, byid = TRUE, id = row.names(a2))
      rm(a2)
      # rbind polygons
      if (!exists("a2s")) a2s <- a2s1 else a2s <- rbind(a2s, a2s1)
    }
    if (!exists("a2s")) next
    pts <- gPointOnSurface(a2s, byid = TRUE, id = row.names(a2s))

    if (is.null(num.samps)) {
      pts.s <- pts
    } else {
      ns <- num.samps[i]
      if (ns < 1) ns <- ceiling(length(pts) * ns) # case of fraction
      if (!replace) {
        o.ns <- ns
        if (ns > length(pts)) ns <- length(pts) # case of more samples than available
        pts.s <- sample(pts, size = ns, prob = area(a2s)^2, replace = FALSE)
        if (length(pts.s) < o.ns & force.min) {
          pts.o <- pts.s
          samp.rm <- sample(1:length(pts.s))
          while(length(pts.s) < o.ns) pts.s <- rbind(pts.s, pts.o)
          row.names(pts.s) <- 1:length(pts.s)
          if (length(pts.s) != o.ns) pts.s <- pts.s[!row.names(pts.s) %in% samp.rm[1:(length(pts.s)-o.ns)]]
        }
      } else {
        pts.s <- sample(pts, size = ns, prob = area(a2s)^2, replace = TRUE)
      }
      row.names(pts.s) <- 1:length(pts.s)
    }
    pts.out <- SpatialPointsDataFrame(pts.s,
                                      data.frame(poly.id = rep(as.character(row.names(p)), length(pts.s)),
                                                 row.names = paste0(row.names(p),".",row.names(pts.s))))
    # rbind pts
    if (!exists("sp.out")) sp.out <- pts.out else sp.out <- rbind(sp.out, pts.out)
  }
  # remove duplicates
  if (!replace & !force.min) sp.out <- sp.out[!duplicated(extract(r2, sp.out)),]
  return(sp.out)
}
