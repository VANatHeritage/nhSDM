# poly_samp

#' Create points in polygons in unique raster cells
#'
#' For each spatial polygon, a given number (\code{num_samps})
#' of points are created in cells that the polygon intersects,
#' with only one point allowed per cell (in the polygon and the
#' entire dataset).
#'
#' \code{num.samps} can be a a proportion (decimal value < 1), single integer,
#' or vector of integers equal to length of spP. If left NULL, all
#' cells [n] intersecting the polygon will contain a point. If a proportion is given
#' (e.g., 0.5), than [n * num.samps] will be returned. If a single integer is given,
#' \code{num.samps} points will be sampled in each polygon. If the number of samples to
#' create exceeds the number of unique cells in a given
#' polygon, the number of samples for that polygon is set to the number of cells (n).

#' @param spP SpatialPolygons[DataFrame]
#' @param rast raster dataset with extent greater or equal to spP
#' @param num.samps number of samples to create in each polygon (see details)
#'
#' @author David Bucklin
#'
#' @importFrom raster crop extent ncell rasterToPolygons res values mask extract area values<-
#' @importFrom sp SpatialPointsDataFrame
#' @importFrom rgeos gBuffer gIntersection gPointOnSurface
#'
#' @export
#'
#' @examples
#' \dontrun{
#' a<-raster("D:/SDM/Tobacco/env_vars/Tobacco/AnnMnTemp.tif")
#' spP <- readOGR("D:/SDM/Tobacco/inputs/species/ambymabe/polygon_data", "ambymabe_expl")
#' s <- c(4,3,5)[spP$RA]
#' sp <- poly_samp(spP, a, s)
#' }

poly_samp <- function(spP, rast, num.samps = NULL) {

  # fixed number for num.samps
  if (!is.null(num.samps)) {
    if (length(num.samps) == 1) {
      num.samps <- rep(num.samps, length(spP))
    } else if (length(num.samps) != length(spP)) {
      stop("num.samps must a vector of length one, or length equal to spP.")
    }
  }

  # loop over polygons
  for (i in 1:length(spP)) {

      p <- spP[i,]

      # seq raster
      r1 <- crop(rast, extent(p))
      r2 <- r1
      values(r2) <- 1:ncell(r1)
      r2 <- mask(r2, r1)
      names(r2) <- "seqid"
      rm(r1)

      gb <- gBuffer(p, byid = FALSE, width = res(r2)[1])
      a2 <- rasterToPolygons(crop(r2, gb))
      a2s <- gIntersection(a2, p, byid = TRUE, id = row.names(a2))
      pts <- gPointOnSurface(a2s, byid = TRUE, id = row.names(a2s))

      if (is.null(num.samps)) {
        pts.s <- pts
      } else {
        ns <- num.samps[i]
        if (ns < 1) ns <- ceiling(length(pts) * ns) # case of fraction
        if (ns > length(pts)) ns <- length(pts) # case of more samples than available
        pts.s <- sample(pts, size = ns, prob = area(a2s)^2, replace = FALSE)
      }
      pts.out <- SpatialPointsDataFrame(pts.s,
                                        data.frame(poly.id = rep(row.names(p), length(pts.s)),
                                                   row.names = paste0(row.names(p),".",row.names(pts.s))))

      if (i == 1) sp.out <- pts.out else sp.out <- rbind(sp.out, pts.out)
  }

  # remove duplicates
  sp.out <- sp.out[!duplicated(extract(r2, sp.out)),]
  return(sp.out)

}
