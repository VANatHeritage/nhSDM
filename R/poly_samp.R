# poly_samp

#' Create points in polygons in unique raster cells
#'
#' For each spatial polygon, a given number (\code{num_samps})
#' of points are created in cells that the polygon intersects.
#'
#' \code{num.samps} can be a a proportion (decimal value < 1), single integer,
#' or vector of integers equal to length of spP. If left NULL, \code{num.samps} will
#' be set to the number of cells [n] intersecting the polygon. If a proportion is given
#' (e.g., 0.5), than [n * num.samps] will be returned. If a single integer is given,
#' \code{num.samps} points will be sampled in each polygon.
#'
#' If \code{replace = FALSE} (default), each cell can only contain one point (across
#' each polygon, and the entire returned set of points). If the number
#' of samples points to create exceeds the number of unique cells in a given
#' polygon, the number of samples for that polygon is set to the number of cells.
#' If \code{replace = TRUE}, sampling is done with replacement and duplicates
#' will be allowed. This is the only way to return a number of points
#' greater than the number of cells within a polygon.
#'
#' If CRS do not match, the SpatialPolygons will be transformed to the CRS
#' of the raster.

#' @param spP SpatialPolygons[DataFrame]
#' @param rast raster dataset with extent overlapping spP
#' @param num.samps number of samples to create in each polygon (see details)
#'
#' @author David Bucklin
#'
#' @importFrom raster crs crop extent ncell rasterToPolygons res values mask extract area values<-
#' @importFrom sp SpatialPointsDataFrame spTransform
#' @importFrom rgeos gBuffer gIntersection gPointOnSurface
#'
#' @export
#'
#' @examples
#' \dontrun{
#' r<-raster("D:/SDM/Tobacco/env_vars/Tobacco/AnnMnTemp.tif")
#' spP <- readOGR("D:/SDM/Tobacco/inputs/species/ambymabe/polygon_data", "ambymabe_expl")
#' s <- c(4,3,5)[spP$RA]
#' sp <- poly_samp(spP, r, s)
#' }

poly_samp <- function(spP, rast, num.samps = NULL, replace = FALSE) {

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
    p <- spP[i,]

    gb <- gBuffer(p, byid = FALSE, width = res(r2)[1])
    try(a2 <- rasterToPolygons(crop(r2, gb)), silent = TRUE)
    if (!exists("a2")) message("Polygon at index value [", i , "] does not intersect raster. Skipping...")
    if (!exists("a2")) next
    a2s <- gIntersection(a2, p, byid = TRUE, id = row.names(a2))
    rm(a2)
    pts <- gPointOnSurface(a2s, byid = TRUE, id = row.names(a2s))

    if (is.null(num.samps)) {
      pts.s <- pts
    } else {
      ns <- num.samps[i]
      if (ns < 1) ns <- ceiling(length(pts) * ns) # case of fraction
      if (!replace) {
        if (ns > length(pts)) ns <- length(pts) # case of more samples than available
        pts.s <- sample(pts, size = ns, prob = area(a2s)^2, replace = FALSE)
      } else {
        pts.s <- sample(pts, size = ns, prob = area(a2s)^2, replace = TRUE)
      }
      row.names(pts.s) <- 1:length(pts.s)
    }
    pts.out <- SpatialPointsDataFrame(pts.s,
                                      data.frame(poly.id = rep(as.character(row.names(p)), length(pts.s)),
                                                 row.names = paste0(row.names(p),".",row.names(pts.s))))

    if (!exists("sp.out")) sp.out <- pts.out else sp.out <- rbind(sp.out, pts.out)
  }
  # remove duplicates
  if (!replace) sp.out <- sp.out[!duplicated(extract(r2, sp.out)),]
  return(sp.out)

}
