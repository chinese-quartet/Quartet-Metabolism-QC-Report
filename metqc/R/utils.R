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

  # --- [新增修改] 表达矩阵前两列标准化 ---
  # 目的：兼容 Metabolites/metabolites, HMDBID/hmdbid 等写法
  # 必须标准化为 'metabolites' 和 'HMDBID'，因为后续代码严格依赖这两个名字

  # 1. 第一列: 只要拼写像 metabolites，统一为 metabolites (小写)
  if (ncol(dt) >= 1 && tolower(colnames(dt)[1]) == "metabolites") {
    colnames(dt)[1] <- "metabolites"
  }

  # 2. 第二列: 只要拼写像 hmdbid，统一为 HMDBID (大写)
  if (ncol(dt) >= 2 && tolower(colnames(dt)[2]) == "hmdbid") {
    colnames(dt)[2] <- "HMDBID"
  }
  # ------------------------------------

  map <- MetInfo$HMDBID
  names(map) <- MetInfo$metabolites

  dt[is.na(HMDBID)]$HMDBID <- map[dt[is.na(HMDBID)]$metabolites]
  if (sum(is.na(dt$HMDBID))) {
    dt[is.na(HMDBID)]$HMDBID <- sprintf("Unknown%04s", 171:(170 + nrow(dt[is.na(HMDBID)])))
  }
  return(dt)
}
