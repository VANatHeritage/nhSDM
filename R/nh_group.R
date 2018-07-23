# nh_group

#' Group spatial features using a defined separation distance
#'
#' Input features are allocated into new groups depending
#' on whether they lie within the separation distance of
#' another input feature.
#'
#' The grouping is done by input ID, so multi-features
#' are considered one input feature (they will not be
#' be ungrouped).
#'
#' The separation distance \code{sep.dist} is numeric and in the
#' units of \code{spf}'s coordinate system, unless the
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
#' @param spf input spatial features (sp or sf spatial object)
#' @param sep.dist separation distance with which to define groups (see description)
#' @param union whether to union output groups into multi-features
#' 
#' @return sp or sf object (same as input)
#'
#' @author David Bucklin
#'
#' @importFrom methods as
#' @importFrom sf st_as_sf st_distance st_is_longlat
#' @importFrom utils installed.packages
#'
#' @export
#'
#' @examples
#' \dontrun{
#' spf <- rgdal::readOGR("ambymabe/polygon_data", "ambymabe")
#' spg <- nh_group(spf, 1000)
#' }

nh_group <- function(spf, sep.dist = 0, union = FALSE) {

  # make sure spdf

  # handle name conflicts...
  # dt <- gsub("\\s|:|-", "", as.character(Sys.time()))
  #  if ("group" %in% names(spf))

  # handle sp/sf class
  spf1 <- tospf(spf)
  sp <- spf1[[1]]
  spf <- spf1[[2]]
  rm(spf1)

  row.names(spf) <- 1:length(spf$geometry)
  spf$TEMPRN <- row.names(spf)

  if (isTRUE(st_is_longlat(spf))) {
      if (!"lwgeom" %in% installed.packages()) stop("Need to install package 'lwgeom' for lat/lon distance calculations.")
      message("Using geodesic distance calculation (in meters)...")
  } else {
      message("Distances calculated in coordinate system units...")
  }

  gd <- as.data.frame(st_distance(spf))
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
  spf <-merge(spf, polyg, by = "TEMPRN")
  spf$TEMPRN <- NULL

  if (union) {
    if (!"dplyr" %in% installed.packages()) {
      message("Install package 'dplyr' to enable union. Returning non-unioned features...")
    } else {
      group <- NULL # just here to avoid check() notes (undefined global variable)
      spf <- spf %>% dplyr::group_by(group) %>% dplyr::summarize(count = n())
    }
  }

  if (sp) return(as(spf,Class = "Spatial")) else return(spf)
}
