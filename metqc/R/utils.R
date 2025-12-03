# ---------------------------------------------------------------------------- #
#' @title Pre-process dataset
#'
#' @description Map the name of metabolites to the existing ID files.
#'
#' @param dt_file Data table file
#'
#' @return Numeric vector
#' @importFrom data.table data.table
#'
#' @export
map_hmdb_id <- function(dt_file) {
  dt <- fread(dt_file)

  map <- MetInfo$HMDBID
  names(map) <- MetInfo$metabolites

  dt[is.na(HMDBID)]$HMDBID <- map[dt[is.na(HMDBID)]$metabolites]
  if (sum(is.na(dt$HMDBID))) {
    dt[is.na(HMDBID)]$HMDBID <- sprintf("Unknown%04s", 171:(170 + nrow(dt[is.na(HMDBID)])))
  }
  return(dt)
}
