# nh_stack

#' Stack multiple binary SDM rasters into one raster layer.
#' 
#' Takes a \code{rastfiles} of binary raster's filenames (values 0/1), 
#' a template raster covering the extent desired for the stack, 
#' and optionally the \code{codes} (names) to use
#' for the rasters. Returns a single raster layer (factor type with an
#' attribute table). 
#' 
#' The raster attribute table (accessed using \code{levels()}) 
#' has four columns:  `ID`: the unique integer value in the raster; `VALUE`:
#' the internal nh_stack unique value for that `ID`, `ALLCODES`: the
#' identity of species codes, pasted in a character vector seperated by `;`,
#' `ALLCODES_CT`: the number of unique codes for that value.
#' 
#' If \code{codes} is not given, the raster layer name will be used as the 
#' layer's code. If these are not unique, the internal nh_stack unique value will
#' be pasted to the end of the original code, and a message will be printed.
#' 
#' All rasters must have the same projection and resolution, though
#' they can have different extents - the processing extent is defined by
#' \code{rast}. To summarize all rasters across their
#' entire extents, \code{template} should essentially be the union (mosaic) 
#' of all raster extents.
#' 
#' When \code{return.table} is TRUE, a list with two objects (1), the 
#' stack raster, and (2) a summary table of included rasters is returned (with
#' internal nh_stack unique values, species codes, and file names).
#' 
#' @param rastfiles raster file names; character vector
#' @param rast raster template dataset
#' @param codes species codes; character vector. If given, must match length of \code{rastfiles}
#' @param return.table Whether to return a table with nh_stack unique values, species codes, and filenames
#' 
#' @return RasterLayer
#' 
#' @author David Bucklin
#'
#' @import raster
#' @importFrom stringi stri_rand_strings stri_length stri_sub
#' @export
#'
#' @examples
#' \dontrun{
#' rast<-raster::raster("project_mask.tif")
#' list <- sort(list.files(paste0(rastout, t, "/thumb"), full.names = T,
#'  pattern = "^[[:lower:]]{5,12}\\.tif$"))
#'
#' stack <- nh_stack(list, rast)
#' # view raster attribute table
#' levels(stack)[[1]]
#' }

nh_stack <- function(rastfiles, rast, codes = NULL, return.table = FALSE) {
  
  list <- rastfiles
  if (!is.null(codes) & length(codes) != length(list)) stop("`codes` length must match `rastfiles` length.")
  
  message("Prepping template raster...")
  rcont <- setValues(raster(rast), 1:ncell(rast))
  dataType(rcont) <- "INT4S"
  
  # make unique codes
  sz <- length(list)
  
  len <- 0
  cd <- character(0)
  while (length(cd) < sz) {
    len <- len + 1
    cd <- unique(stri_rand_strings(sz * 10, len, pattern = "[a-z0-9]"))  
  }
  # get codes for length you need
  stk_order <- data.frame(uval = sample(cd, size = sz, replace = F), code_nm = NA, file = list)
  
  bigd <- rep("",ncell(rcont))
  # loop over layers
  message("Stacking layers...")
  for (i in 1:length(stk_order$uval)) { 
    # get nh_stack uval
    spcd <- stk_order$uval[i]
    # get/crop raster
    r <- raster(stk_order$file[i])
    r <- crop(r, rcont)
    # assign long code name
    if (!is.null(codes)) names(r) <- codes[i]
    if (names(r) %in% stk_order$code_nm) {
      message("Non-unique code '", names(r), "' changed to '",paste0(names(r), "_", spcd),"'.")
      names(r) <- paste0(names(r), "_", spcd)
    }
    message(paste0("Working on raster ", names(r), "..."))
    stk_order$code_nm[stk_order$uval == spcd] <- names(r)
    
    # crop/mask raster, get values
    cr <- crop(rcont, r, datatype = "INT4S")
    rval <- mask(cr, r, maskvalue = 1, inverse = T)
    v <- values(rval)
    v <- v[!is.na(v)]
    
    v2 <- rval[]
    v2 <- v2[!is.na(v2)]
    # paste to bigd
    if (length(v2) > 0) bigd[v] <- paste0(bigd[v], spcd) 
  }
  
  message("Calculating stack attributes...")
  r2 <- setValues(raster(rcont), as.factor(bigd))

  # get unique values
  uvals <- levels(r2)[[1]]
  
  # parse
  pars <- lapply(uvals$VALUE, 
                 function(x) if (stri_length(x) > len) stri_sub(x, from = seq(1, stri_length(x), by = len), length = len) else x)
  parssp <- unlist(lapply(pars, function(x) paste(stk_order$code_nm[stk_order$uval %in% x], collapse = ";")))
  parsct <- unlist(lapply(pars, length))
  
  uvals$ALLCODES <- parssp
  uvals$ALLCODES_CT <- parsct
  uvals$ALLCODES_CT[uvals$VALUE == ""] <- 0
  
  levels(r2) <- uvals
  names(r2) <- "nh_stack"
  
  if (return.table) {
    names(stk_order) <- c("nh_stack_uval","CODE","FILENAME")
    return(list(r2, stk_order))
  } else {
    return(r2)
  }
}
