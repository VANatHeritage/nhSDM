# Not-exported nhSDM functions

# tospf
#'
#' Return a list with sp-flag and sf object, given either sf or sp object
#' 
#' @param spf input sf or sp object
#' @param rastproj raster dataset with desired output projection
#' 
#' @return sf object
#' 
#' @import sf
#' @importFrom methods as
#' @importFrom raster projection
#' 
#' @keywords internal

tospf <- function(spf, rastproj) {
  if (grepl("^Spatial*", class(spf)[1])) {
    sp <- TRUE
    spf <- st_as_sf(spf)
  } else if (grepl("sf", class(spf)[1])) {
    sp <- FALSE
    names(spf)[length(names(spf))] <- "geometry"
  } else {
    stop("Must provide either 'sp' or 'sf'-class spatial object.")
  }
  spf <- st_zm(spf)
  
  # transform if necessary
  if (!missing(rastproj)) {
    if (!is.na(st_crs(spf)$proj4string)) {
      if (st_crs(spf)$proj4string != projection(rastproj)) spf <- st_transform(spf, crs = projection(rastproj))
    } else {
      message("No projection on input features. Assuming features are using raster's projection...")
      st_crs(spf) <- projection(rastproj)
    }
  }
  
  return(list(sp,spf))
}

