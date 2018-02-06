# nh_group

#' Group spatial features using a defined separation distance
#'
#' Input features are allocated into new groups depending
#' on whether they lie within the seperation distance of
#' another input feature.
#'
#' The grouping is done by input ID, so multi-features
#' are considered one input feature (they will not be
#' be ungrouped).
#'
#' The seperation distance \code{sep.dist} is numeric and in the
#' units of \code{spP}'s coordinate system, unless the
#' coordinate system uses latitude/longitude as the unit (e.g. WGS 84). In
#' these cases, geodesic distances will be used and \code{sep.dist}
#' should be specified in meters.
#'
#' A column 'group' will be added to the output features.
#' Specifying \code{union = TRUE} will output one (multi)feature per group,
#' meaning original attributes will be discarded - only `group` and `count`
#' (the number of original features in the group) will be returned. This feature
#' requires the package \code{dplyr} to be installed.
#'
#' @param spP input sp or sf spatial object
#' @param sep.dist separation distance with which to define groups (see description)
#' @param union whether to union output groups into multi-features
#'
#' @author David Bucklin
#'
#' @importFrom methods as
#' @importFrom sf st_as_sf st_distance st_is_longlat
#'
#' @export
#'
#' @examples
#' \dontrun{
#' spP <- rgdal::readOGR("D:/SDM/Tobacco/inputs/species/ambymabe/polygon_data", "ambymabe")
#' spg <- nh_group(spP, 1000)
#' }

nh_group <- function(spP, sep.dist = 0, union = FALSE) {

  # make sure spdf

  # handle name conflicts...
  # dt <- gsub("\\s|:|-", "", as.character(Sys.time()))
  #  if ("group" %in% names(spP))

  if (grepl("^Spatial*", class(spP)[1])) {
    sp <- TRUE
    spP <- st_as_sf(spP)
  } else if (grepl("sf", class(spP)[1])) {
    sp <- FALSE
  } else {
    stop("Must provide either 'sp' or 'sf'-class spatial object.")
  }

  row.names(spP) <- 1:length(spP$geometry)
  spP$TEMPRN <- row.names(spP)

  if (isTRUE(st_is_longlat(spP))) {
      if (!"lwgeom" %in% installed.packages()) stop("Need to install package 'lwgeom' for lat/lon distance calculations.")
      message("Using geodesic distance calculation (in meters)...")
  } else {
      message("Distances calculated in coordinate system units...")
  }

  gd <- as.data.frame(st_distance(spP))
  for (i in names(gd)) gd[i] <- as.numeric(gd[,i]) # remove units
  names(gd) <- row.names(gd)

  # set up group list
  grps <- as.list(rep(NA, length(colnames(gd))))
  names(grps) <- colnames(gd)

  # group
  for (n in colnames(gd)) {
    if (is.na(grps[n])) grps[n] <- n else next
    unass <- names(grps)[is.na(grps)]
    if (length(unass) == 0) break
    d <- gd[n,]
    subg <- names(d[unass])[d[unass] <= sep.dist]
    lb <- 0
    la <- 1
    done <- character(0)
    while (lb != la) {
      lb <- length(subg)
      for (i in subg) {
        if (!i %in% done) {
        d <- gd[i,]
        add <- names(d[unass])[d[unass] <= sep.dist]
        subg <- c(subg, add[!add %in% subg])
        done <- c(done, i)
        }
      }
      la <- length(subg)
    }
    grps[subg] <- n
  }

  polyg <- data.frame(TEMPRN = names(grps), group = as.numeric(unlist(grps)))
  spP <-merge(spP, polyg, by = "TEMPRN")
  spP$TEMPRN <- NULL

  if (union) {
    if (!"dplyr" %in% installed.packages()) {
      message("Install 'dplyr' to enable union. Returning non-unioned features...")
      break
    }
    spP <- spP %>% dplyr::group_by(group) %>% dplyr::summarize(count = n())
    # spP <- unionSpatialPolygons(spP, spP$group)
    # spP <- SpatialPolygonsDataFrame(spP, data = data.frame(group=row.names(spP), row.names = row.names(spP)), match.ID = TRUE)
  }

  if (sp) return(as(spP,Class = "Spatial")) else return(spP)
}
