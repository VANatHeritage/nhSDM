# nh_sample

#' Create points in polygons in reference raster cells
#'
#' For each spatial polygon, a given number (\code{num_samps})
#' of points are created in cells that the polygon intersects.
#'
#' \code{num.samps} can be a a proportion (a decimal value < 1), single integer,
#' or vector of integers equal to length of spf indicated the number of samples
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

#' @param spf input spatial features (sp or sf spatial object)
#' @param rast raster dataset with extent overlapping spf
#' @param num.samps number of samples to create in each polygon (see details)
#' @param replace whether to sample with or without replacement
#' @param force.min whether to always create \code{num.samps} points, even if they are duplicates
#'
#' @author David Bucklin
#'
#' @import sf
#' @importFrom methods as
#' @import raster
#' @export
#'
#' @examples
#' \dontrun{
#' r<-raster::raster("AnnMnTemp.tif")
#' spf <- rgdal::readOGR("ambymabe/polygon_data", "ambymabe")
#' # can also use sf: spf <- sf::st_read("ambymabe/polygon_data/ambymabe.shp")
#' spf.samps <- nh_sample(spf, r, num.samps)
#' }

nh_sample <- function(spf, rast, num.samps = NULL, replace = FALSE, force.min = FALSE) {

  if (grepl("^Spatial*", class(spf)[1])) {
    sp <- TRUE
    spf <- st_as_sf(spf)
  } else if (grepl("sf", class(spf)[1])) {
    sp <- FALSE
  } else {
    stop("Must provide either 'sp' or 'sf'-class spatial object.")
  }
  spf <- st_zm(spf)

  if (!is.null(num.samps)) {
    if (length(num.samps) == 1) {
      num.samps <- rep(num.samps, length(spf$geometry))
    } else if (length(num.samps) != length(spf$geometry)) {
      stop("num.samps must a vector of length one, or length equal to spf.")
    }
  }
  # transform if necessary
  if (!is.na(st_crs(spf)$proj4string)) {
    if (st_crs(spf)$proj4string != projection(rast)) spf <- st_transform(spf, crs = projection(rast))
  } else {
    message("No projection on input features. Assuming features are using raster's projection...")
    st_crs(spf) <- projection(rast)
  }
  # seq raster
  r1 <- crop(rast, extent(extent(spf)[1], extent(spf)[2], extent(spf)[3], extent(spf)[4])) # bug when applying extent object from spf directly
  r2 <- r1
  values(r2) <- 1:ncell(r1)
  r2 <- mask(r2, r1)
  names(r2) <- "seqid"
  rm(r1)

  # loop over polygons
  for (i in 1:length(spf$geometry)) {
    suppressWarnings(rm(a2s))
    p <- spf[i,]
    pall <- st_cast(p)

    for (xp in 1:length(pall$geometry)) {
      pxp <- pall[xp,]
      gb <- st_buffer(pxp, res(r2)[1])
      try(a2 <- rasterToPolygons(crop(r2, extent(extent(gb)[1], extent(gb)[2], extent(gb)[3], extent(gb)[4]))), silent = TRUE)
      if (!exists("a2")) message("Polygon at index value [", i , ".", xp , "] does not intersect raster. Skipping...")
      if (!exists("a2")) next
      a2 <- st_as_sf(a2)
      row.names(a2) <- paste0(row.names(a2), ".", xp)
      a2s1 <- suppressWarnings(st_intersection(a2, pxp))
      rm(a2)
      # rbind polygons
      if (!exists("a2s")) a2s <- a2s1 else a2s <- rbind(a2s, a2s1)
    }
    if (!exists("a2s")) next
    # pts <- gPointOnSurface(a2s, byid = TRUE, id = row.names(a2s))
    pts <- st_centroid(a2s) #not sure about weird polygon cetroids...

    if (is.null(num.samps)) {
      pts.s <- pts
    } else {
      ns <- num.samps[i]
      if (ns < 1) ns <- ceiling(length(pts$geometry) * ns) # case of fraction
      if (!replace) {
        o.ns <- ns
        if (ns > length(pts$geometry)) ns <- length(pts$geometry) # case of more samples than available
        pts.s <- pts[row.names(pts) %in% sample(row.names(pts), size = ns, prob = as.numeric(st_area(a2s)), replace = FALSE),]
        if (length(pts.s$geometry) < o.ns & force.min) {
          pts.o <- pts.s
          samp.rm <- sample(1:length(pts.s$geometry))
          while(length(pts.s$geometry) < o.ns) pts.s <- rbind(pts.s, pts.o)
          row.names(pts.s) <- 1:length(pts.s$geometry)
          if (length(pts.s$geometry) != o.ns) pts.s <- pts.s[!row.names(pts.s) %in% samp.rm[1:(length(pts.s$geometry)-o.ns)],]
        }
      } else {
        pts.s <- pts[c(sample(row.names(pts), size = ns, prob = st_area(a2s), replace = TRUE)),]
      }
      row.names(pts.s) <- 1:length(pts.s$geometry)
    }
    pts.out <- pts.s
    row.names(pts.out) <- paste0(row.names(p),".",row.names(pts.s))
    pts.out$poly.id <- rep(as.character(row.names(p)), length(pts.s$geometry))

    # rbind pts
    if (!exists("sp.out")) sp.out <- pts.out else sp.out <- rbind(sp.out, pts.out)
  }
  # remove duplicates
  if (!replace & !force.min) sp.out <- sp.out[!duplicated(extract(r2, sp.out)),]
  if (sp) return(as(sp.out,Class = "Spatial")) else return(sp.out)
}
