# ---------------------------------------------------------------------------- #
#' @title Get output table
#'
#' @description Get information of overall performance with historical data.
#'
#' @param dt_file Data table file
#' @param metadata_file Metadata table file
#' @param output_dir Output dir
#'
#' @return Numeric vector
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#'
#' @examples
#' get_performance(dt_file = sample_data, metadata_file = sample_metadata)
#'
#' @export
get_performance <- function(dt_file, metadata_file, output_dir = NULL) {
  dt <- map_hmdb_id(dt_file)
  metadata <- fread(metadata_file)

  metadata <- metadata[metadata$sample %in% c("D5", "D6", "F7", "M8"), ]
  cols <- c("metabolites", "HMDBID", metadata$col_names)

  ### check input
  # 检查列名是否存在于dt中
  if (!all(cols %in% colnames(dt))) {
    missing_cols <- cols[!cols %in% colnames(dt)]
    stop("The following columns are missing in 'dt': ", paste(missing_cols, collapse = ", "), ". Please check your data format and keep it identical with our example data.")
  }

  # 尝试选择列，并在错误时停止执行
  dt <- tryCatch(
    {
      dt[, ..cols]
    },
    error = function(e) {
      stop("Error in selecting columns: ", e$message)
    }
  )

  SNR <- count_snr(dt_file = dt_file, metadata_file = metadata_file, output_dir = output_dir)
  RC <- count_rc(dt_file = dt_file, metadata_file = metadata_file, output_dir = output_dir)
  Recall <- count_recall(dt_file = dt_file, metadata_file = metadata_file)

  # dt.result <- data.table("QUERIED DATA","SNR"=SNR,"RC"=RC,"Recall"=Recall)
  dt.result <- data.table("QUERIED DATA", "SNR" = SNR$SNR, "RC" = RC$RC, "Recall" = Recall)
  names(dt.result)[1] <- names(HistoricalData)[1]
  percentile <- ecdf(HistoricalData$Total_normalized)
  HistoricalData$perc <- round(percentile(HistoricalData$Total_normalized) * 100, 0)
  HistoricalData[Total_normalized == 1]$perc <- 0
  HistoricalData[Total_normalized == 10]$perc <- 100

  dt.overall <- rbindlist(list(HistoricalData, dt.result), fill = T)

  fn.query <- function(x) {
    xMin <- min(x[1:(length(x) - 1)], na.rm = TRUE)
    xMax <- max(x[1:(length(x) - 1)], na.rm = TRUE)

    if (x[length(x)] <= xMin) {
      return(1)
    }
    if (x[length(x)] >= xMax) {
      return(10)
    }

    x.norm <- round((10 - 1) * (x[length(x)] - xMin) / (xMax - xMin) + 1, 3)
    return(x.norm)
  }
  dt.overall[batchcode == "QUERIED DATA"]$SNR_normalized <- fn.query(dt.overall$SNR)
  dt.overall[batchcode == "QUERIED DATA"]$RC_normalized <- fn.query(dt.overall$RC)
  dt.overall[batchcode == "QUERIED DATA"]$Recall_normalized <- fn.query(dt.overall$Recall)
  cols <- c("SNR_normalized", "RC_normalized", "Recall_normalized")
  dt.overall$Total <- apply(dt.overall[, ..cols], 1, function(x) exp(mean(log(x), na.rm = T)))
  dt.overall[batchcode == "QUERIED DATA"]$Total_normalized <- fn.query(dt.overall$Total)

  perc.tt <- HistoricalData[Total_normalized == dt.overall[batchcode == "QUERIED DATA"]$Total_normalized]
  percentile2 <- ecdf(dt.overall$Total_normalized)
  dt.overall[batchcode == "QUERIED DATA"]$perc <- ifelse(nrow(perc.tt) > 0, perc.tt$perc,
    round(percentile2(dt.overall[batchcode == "QUERIED DATA"]$Total_normalized) * 100, 0)
  )

  get_class <- function(aList) {
    quantiles <- quantile(aList, probs = seq(0, 1, .2), na.rm = T)
    quantiles2 <- quantile(aList, probs = seq(0, 1, .25), na.rm = T)
    cutoffs <- c(quantiles[1], quantiles[2], quantiles2[3], quantiles[5], quantiles[6])
    result <- as.character(cut(aList,
      breaks = cutoffs,
      include.lowest = T, labels = c("Bad", "Fair", "Good", "Great")
    ))
    return(result)
  }

  dt.overall$Class <- get_class(dt.overall$Total_normalized)
  names(dt.overall)[1] <- "batch"
  setkey(dt.overall, batch)
  if (!is.null(output_dir)) {
    write.csv(x = dt.overall, file = file.path(output_dir, "rank_table.csv"), row.names = F, quote = T)
  }

  cols <- c("SNR", "RC", "Recall", "Total")
  dt.overall.hist <- dt.overall[batch != "QUERIED DATA"]
  dt.overall.hist.stat <- apply(
    dt.overall.hist[batch != "QUERIED DATA", ..cols], 2,
    function(x) paste(round(mean(x, na.rm = T), 2), "±", SD = round(sd(x, na.rm = T), 2))
  )

  dt.overall$SNR_Class <- get_class(dt.overall$SNR_normalized)
  dt.overall$RC_Class <- get_class(dt.overall$RC_normalized)
  dt.overall$Recall_Class <- get_class(dt.overall$Recall_normalized)
  dt.overall$SNR_rank <- paste(rank(-dt.overall$SNR_normalized, ties.method = "min"), "/", nrow(dt.overall))
  dt.overall$RC_rank <- paste(rank(-dt.overall$RC_normalized, ties.method = "min"), "/", nrow(dt.overall))
  dt.overall$Recall_rank <- paste(rank(-dt.overall$Recall_normalized, ties.method = "min"), "/", nrow(dt.overall))
  dt.overall$Total_rank <- paste(rank(-dt.overall$Total_normalized, ties.method = "min"), "/", nrow(dt.overall))

  outputdt <- data.table(
    "Quality metrics" = c("Signal-to-Noise Ratio (SNR)", "Relative Correlation with Reference Datasets (RC)", "Recall of DAMs in Reference Datasets (Recall)", "Total Score"),
    "Value" = as.numeric(round(dt.overall["QUERIED DATA", ..cols], 3)),
    "Historical value (mean ± SD)" = as.character(dt.overall.hist.stat),
    "Rank" = as.character(dt.overall[batch == "QUERIED DATA", c("SNR_rank", "RC_rank", "Recall_rank", "Total_rank")]),
    "Performance" = as.character(dt.overall[batch == "QUERIED DATA", c("SNR_Class", "RC_Class", "Recall_Class", "Class")])
  )

  quantiles <- quantile(dt.overall$Total_normalized, probs = seq(0, 1, .2), na.rm = T)
  quantiles2 <- quantile(dt.overall$Total_normalized, probs = seq(0, 1, .25), na.rm = T)
  cutoffs <- data.frame(c(quantiles[1], quantiles[2], quantiles2[3], quantiles[5], quantiles[6]))
  colnames(cutoffs)[1] <- "Percentile"
  cutoffs$"Cut-off" <- rownames(cutoffs)
  dt.overall.sub <- dt.overall[batch == "QUERIED DATA", c("Total_normalized", "perc")]
  dt.overall.sub$perc <- paste0(dt.overall.sub$perc, "%")
  colnames(dt.overall.sub) <- c("Percentile", "Cut-off")
  cutoffs.all <- rbind(cutoffs, dt.overall.sub)

  if (!is.null(output_dir)) {
    write.csv(x = outputdt, file = file.path(output_dir, "conclusion_table.csv"), row.names = F, quote = T)
    write.csv(x = cutoffs.all, file = file.path(output_dir, "cutoff_table.csv"), row.names = F, quote = T)
  }

  return(list(
    "rank_table" = dt.overall,
    "conclusion_table" = outputdt,
    "cutoff_table" = cutoffs.all,
    "scplot" = RC$scplot,
    "pcaplot" = SNR$pcaplot
  ))
}

# ---------------------------------------------------------------------------------- #
#' @title Calculate SNR and plot PCA
#'
#' @description Plot PCA and calculate SNR based on the first two principal componants of PCA.
#'
#' @param dt_file Data table file
#' @param metadata_file Metadata table file
#' @param output_dir Output dir
#'
#' @return Numeric vector
#' @importFrom data.table setkey
#' @importFrom data.table fread
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @importFrom reshape2 melt
#' @importFrom data.table setDT
#' @importFrom data.table :=
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 element_text
#' @examples
#' count_snr(dt = sample_data, metadata = sample_metadata)
#'
#' @export
count_snr <- function(dt_file, metadata_file, output_dir = NULL) {
  dt <- map_hmdb_id(dt_file)
  metadata <- fread(metadata_file)

  cols <- metadata$col_names
  setkey(setDT(metadata), col_names)
  dt.num.0 <- dt[, ..cols]
  dt.num.0[dt.num.0 == 0] <- NA
  dt.num.log2 <- apply(dt.num.0, 2, function(x) log2(x))
  dt.num.t.cp <- t(data.frame(dt.num.log2[complete.cases(dt.num.log2), ]))

  pca_prcomp <- prcomp(x = dt.num.t.cp, scale. = T)

  pcs <- as.data.frame(predict(pca_prcomp))
  dt.perc.pcs <- data.table(
    PCX = 1:ncol(pcs),
    Percent = summary(pca_prcomp)$importance[2, ],
    AccumPercent = summary(pca_prcomp)$importance[3, ]
  )

  dt.dist <- data.table(
    ID.A = rep(rownames(pcs), each = nrow(pcs)),
    ID.B = rep(rownames(pcs), time = nrow(pcs))
  )

  map <- metadata$sample
  names(map) <- metadata$col_names
  dt.dist$Group.A <- map[dt.dist$ID.A]
  dt.dist$Group.B <- map[dt.dist$ID.B]

  dt.dist[, Type := ifelse(ID.A == ID.B, "Same",
    ifelse(Group.A == Group.B, "Intra", "Inter")
  )]
  dt.dist[, Dist := (dt.perc.pcs[1]$Percent * (pcs[ID.A, 1] - pcs[ID.B, 1])^2 + dt.perc.pcs[2]$Percent * (pcs[ID.A, 2] - pcs[ID.B, 2])^2)]

  dt.dist.stats <- dt.dist[, .(Avg.Dist = mean(Dist)), by = Type]
  setkey(dt.dist.stats, Type)
  signoise <- dt.dist.stats["Inter"]$Avg.Dist / dt.dist.stats["Intra"]$Avg.Dist

  signoise_db <- 10 * log10(signoise)

  pcs$col_names <- rownames(pcs)
  dt.forPlot <- merge(pcs, metadata, by = "col_names")
  colors.Quartet <- c(D5 = "#4CC3D9", D6 = "#7BC8A4", F7 = "#FFC65D", M8 = "#F16745")

  pcaplot <- ggplot(dt.forPlot, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample),
      size = 3
    ) +
    theme_bw() +
    labs(
      x = sprintf("PC1: %.1f%%", summary(pca_prcomp)$importance[2, 1] * 100),
      y = sprintf("PC2: %.1f%%", summary(pca_prcomp)$importance[2, 2] * 100),
      color = "sample",
      title = sprintf("SNR = %.2f", signoise_db)
    ) +
    theme(plot.title = element_text(hjust = 0.5, size = 12)) +
    scale_color_manual(values = colors.Quartet) +
    theme(legend.position = "right")


  if (!is.null(output_dir)) {
    write.csv(x = dt.forPlot, file = file.path(output_dir, "PCAtable.csv"), row.names = F)
    ggsave(
      filename = file.path(output_dir, "PCA_withSNR.png"), pcaplot,
      device = "png", width = 10, height = 8, units = c("cm"), dpi = 300
    )
  }

  ### 返回snr，pcaplot
  # return(signoise_db)
  return(list(
    "SNR" = signoise_db,
    "pcaplot" = pcaplot
  ))
}

# ---------------------------------------------------------------------------------- #
#' @title Calculate and plot RC
#'
#' @description Calculate correlation to reference datasets and plot related scatter plot
#'
#' @param dt_file Data table file
#' @param metadata_file Metadata table file
#' @param output_dir Output dir
#'
#' @return Numeric vector
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @importFrom reshape2 melt
#' @importFrom data.table copy
#' @importFrom data.table setDT
#' @importFrom data.table :=
#' @importFrom data.table dcast.data.table
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_alpha_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 geom_abline
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 element_text
#' @importFrom stringr str_split_fixed
#'
#' @examples
#' count_rc(dt = sample_data, metadata = sample_metadata)
#'
#' @export
count_rc <- function(dt_file, metadata_file, output_dir = NULL) {
  dt <- map_hmdb_id(dt_file)
  metadata <- fread(metadata_file)

  RefDataset <- MetReference
  metsinRef <- unique(RefDataset$HMDBID)

  metadata$platform <- gsub("Untargeted", "U", metadata$strategy)
  metadata$platform <- gsub("Targeted", "T", metadata$platform)
  metadata$batchcode <- paste(metadata$platform, metadata$lab, metadata$batch, sep = "_")
  metadata$samplecode <- paste(metadata$batchcode, metadata$sample, sep = "_")

  cols <- metadata$col_names
  setkey(setDT(metadata), col_names)
  # remove metabolites without full record
  # 0 recode as NA
  dt.num.0 <- dt[, ..cols]
  dt.num.0[dt.num.0 == 0] <- NA
  dt.rmNA <- dt[complete.cases(dt.num.0)]

  # log2
  cols <- c("HMDBID", metadata$col_names)
  dt.long <- melt(dt.rmNA[, ..cols], id.vars = "HMDBID", variable.name = "col_names")
  cols <- c("col_names", "sample")
  dt.long.info <- data.table(merge(dt.long, metadata[, ..cols], by = "col_names"))
  dt.long.info.log2 <- copy(dt.long.info)
  dt.long.info.log2$value <- log2(dt.long.info.log2$value)

  # get relative abundance to D6
  togetRe.log2 <- dt.long.info.log2[HMDBID %in% metsinRef, ]
  cols <- c("HMDBID", "sample", "value")
  togetRe.ave <- togetRe.log2[, .(log2mean = mean(value)), by = c("HMDBID", "sample")]
  togetRe.ave.wide <- dcast.data.table(togetRe.ave, HMDBID ~ sample, value.var = "log2mean")
  togetRe.ave.wide$D5toD6 <- togetRe.ave.wide$D5 - togetRe.ave.wide$D6
  togetRe.ave.wide$F7toD6 <- togetRe.ave.wide$F7 - togetRe.ave.wide$D6
  togetRe.ave.wide$M8toD6 <- togetRe.ave.wide$M8 - togetRe.ave.wide$D6

  cols <- c("HMDBID", "D5toD6", "F7toD6", "M8toD6")
  togetRe.ave.long <- melt(togetRe.ave.wide[, ..cols],
    id.vars = "HMDBID",
    variable.name = "dataset", value.name = "log2FC"
  )

  # subset DAMs in RDs
  togetRe.withRef <- merge(togetRe.ave.long, RefDataset, by = c("dataset", "HMDBID"))

  # Pearson correlation to RDs
  RC <- cor(togetRe.withRef$Mean, togetRe.withRef$log2FC,
    use = "pairwise.complete.obs", method = "pearson"
  )

  # Scatter plot
  scplot <- ggplot(togetRe.withRef, aes(x = Mean, y = log2FC, color = dataset)) +
    geom_point() +
    scale_color_manual(values = c(D5toD6 = "#4CC3D9", F7toD6 = "#FFC65D", M8toD6 = "#F16745")) +
    theme(legend.position = "right") +
    scale_y_continuous(limits = c(-3, 3)) +
    scale_x_continuous(limits = c(-3, 3)) +
    labs(
      x = "Ratios in reference datasets",
      y = "Measured ratios",
      color = "Sample pair",
      title = sprintf("RC = %.3f", RC)
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 12))

  if (!is.null(output_dir)) {
    write.csv(x = togetRe.withRef, file = file.path(output_dir, "RCtable.csv"), row.names = F)
    ggsave(
      filename = file.path(output_dir, "ScatterPlot_withRC.png"), scplot, device = "png",
      width = 10, height = 8, units = c("cm"), dpi = 300
    )
  }

  ### 返回RC和scplot
  # return(RC)
  return(list(
    "RC" = RC,
    "scplot" = scplot
  ))
}

# ---------------------------------------------------------------------------------- #
#' @title count_recall
#'
#' @description Calculate Recall of DAMs
#'
#' @param dt_file Data table file
#' @param metadata_file Metadata table file
#'
#' @return Numeric vector
#' @importFrom data.table data.table
#' @importFrom reshape2 melt
#' @importFrom data.table copy
#' @importFrom data.table :=
#' @importFrom data.table setDT
#'
#' @examples
#' count_recall(dt = sample_data, metadata = sample_metadata)
#'
#' @export
count_recall <- function(dt_file, metadata_file) {
  dt <- map_hmdb_id(dt_file)
  metadata <- fread(metadata_file)

  cols <- metadata$col_names
  setkey(setDT(metadata), col_names)

  # remove metabolites without full record
  # 0 recode as NA
  dt.num.0 <- dt[, ..cols]
  dt.num.0[dt.num.0 == 0] <- NA
  dt.rmNA <- dt[complete.cases(dt.num.0)]

  # log2
  cols <- c("HMDBID", metadata$col_names)
  dt.long <- melt(dt.rmNA[, ..cols], id.vars = "HMDBID", variable.name = "col_names")
  cols <- c("col_names", "sample")
  dt.long.info <- data.table(merge(dt.long, metadata[, ..cols], by = "col_names"))
  dt.long.info.log2 <- copy(dt.long.info)
  dt.long.info.log2$value <- log2(dt.long.info.log2$value)

  # pairwise t test (default)
  mets.hmdb <- unique(dt.long.info.log2$HMDBID)
  dt.long.info.log2$sample <- factor(dt.long.info.log2$sample, levels = c("D6", "D5", "F7", "M8"))
  dt.log2.stat <- rbindlist(lapply(mets.hmdb, function(xMet) {
    dt.sub <- dt.long.info.log2[HMDBID == xMet, ]

    if (sd(dt.sub$value) == 0) {
      return(NULL)
    }

    dt.ttest <- pairwise.t.test(x = dt.sub$value, g = dt.sub$sample, p.adjust.method = "none")
    dt.ttest.df <- reshape2::melt(dt.ttest$p.value)
    dt.ttest.df.rmNA <- dt.ttest.df[complete.cases(dt.ttest.df), ]
    colnames(dt.ttest.df.rmNA) <- c("group1", "group2", "p")
    dt.ttest.df.rmNA$p.adj <- p.adjust(dt.ttest.df.rmNA$p, method = "holm")
    dt.ttest.df.rmNA$HMDBID <- xMet
    return(dt.ttest.df.rmNA)
  }))
  dt.log2.stat$p.adj2 <- p.adjust(dt.log2.stat$p, method = "fdr")

  # get adjusted p of specific sample pairs (relative to D6)
  dt.log2.stat.D6 <- dt.log2.stat[group2 == "D6", ]
  dt.log2.stat.D6$dataset <- paste0(dt.log2.stat.D6$group1, "to", dt.log2.stat.D6$group2)
  dt.log2.stat.D6.sig <- dt.log2.stat.D6[p.adj2 < 0.05, ]

  # detected DAMs in RDs
  dt.log2.stat.D6.inRD <- merge(dt.log2.stat.D6, MetReference, by = c("HMDBID", "dataset"))
  # identified as DAMs
  dt.log2.stat.D6.inRD.sig <- dt.log2.stat.D6.inRD[p.adj2 < 0.05, ]

  result <- data.table(
    N.detected.pairs = nrow(dt.log2.stat.D6.inRD),
    N.detected.sigpairs = nrow(dt.log2.stat.D6.inRD.sig)
  )

  Recall <- result$N.detected.sigpairs / result$N.detected.pairs

  return(Recall)
}
