# nh_proj

#' Initiate new SDM project
#' 
#' Creates directory structure for a new SDM project and 
#' new SDM project SQLite database. New project folder is 
#' created in working directory, unless \code{folder} is specified.
#' 
#' @param proj.name The name of the project. New folder created and name appended to a new project database file.
#' @param folder The SDM folder in which to initiate the project
#' 
#' @return nothing
#'
#' @author David Bucklin
#'
#' @import RSQLite
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' setwd("D:/testing_SDM")
#' nh_proj("new_proj")
#' }

nh_proj <- function(proj.name, folder = ".") {
  
  s <- paste0(folder , "/", proj.name)
  if (!dir.exists(s)) {
    dir.create(s)
  } else {
    stop("Project folder at that location ('",s,"') already exists, choose a different name.")
  }
  
  # create folders
  dbfile <- paste0(s, "/databases/", "sdm_tracking_", proj.name, ".sqlite")
  if (file.exists(dbfile)) stop("Database file '", dbfile, "' already exists. Choose a different name or use nh_updateproj() to update database.")
  
  folders <- c("databases", "env_vars/background", "env_vars/raster", "env_vars/tabular",
               "other_spatial/raster", "other_spatial/feature",
               "species/_template/inputs/presence", "species/_template/inputs/scripts","species/_template/inputs/model_input",
               "species/_template/outputs/rdata", "species/_template/outputs/metadata", "species/_template/outputs/model_predictions")
  
  for (i in folders) {
    if (!dir.exists(paste0(s, "/" , i))) dir.create(paste0(s, "/" , i), recursive = T)
  }
  
  insert_file <- system.file("sql/sqlite_template_db_nodata.sql", package = "nhSDM")
  sql <- paste(readLines(insert_file, encoding = "UTF-8"), collapse = "")
  sql <- paste(strsplit(sql, ";")[[1]], ";",sep = "")
  tryCatch(expr = {
      conn <- dbConnect(SQLite(), dbfile)
      dbBegin(conn)
      for (i in sql) dbExecute(conn, i)
    },
    finally = {
      dbCommit(conn)
      dbDisconnect(conn)
    })
  message("New project and database '", dbfile , "' created.")
}
