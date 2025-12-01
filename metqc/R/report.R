# ---------------------------------------------------------------------------------- #
#' @title Generate Quartet Metabolomics report
#'
#' @description Use calculated Met result to generate report
#'
#' @param qc_result list
#' @param report_template character
#' @param report_dir character
#' @param report_name character
#'
#' @importFrom dplyr %>%
#' @importFrom flextable flextable
#' @importFrom flextable theme_vanilla
#' @importFrom flextable color
#' @importFrom flextable set_caption
#' @importFrom flextable align
#' @importFrom flextable width
#' @importFrom flextable bold
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 theme
#' @importFrom officer body_add_par
#' @importFrom flextable body_add_flextable
#' @importFrom officer body_add_gg
#' @importFrom officer body_add_break
#' @importFrom officer read_docx
#'
#' @examples
#' # 加载示例 qc_result 对象
#' qc_result_path <- system.file("extdata", "qc_result_example.RData", package = "metqc")
#' load(qc_result_path)
#'
#' # 指定包内文档的路径
#' report_template <- system.file("extdata", "quartet_template.docx", package = "metqc")
#'
#' # 运行函数
#' generate_met_report(qc_result = met_result, report_template = report_template)
#'
#' @export
generate_met_report <- function(qc_result,
                                report_template,
                                report_dir = NULL,
                                report_name = NULL) {
  if (is.null(qc_result) || is.null(report_template)) {
    stop("All arguments (qc_result, report_template) are required.")
  }

  if (is.null(report_dir)) {
    path <- getwd()
    report_dir <- file.path(path, "output")
    dir.create(report_dir, showWarnings = FALSE)
  }

  ### 读取quarter报告模板并生成报告
  if (is.null(report_name)) {
    report_name <- "Quartet_met_report.docx"
  }
  output_file <- file.path(report_dir, report_name)

  ### 创建Evaluate Metrics 表格

  ft1 <- flextable(qc_result$conclusion_table)
  ft1 <- ft1 %>%
    color(~ Performance == "Bad", color = "#B80D0D", ~Performance) %>%
    color(~ Performance == "Fair", color = "#D97C11", ~Performance) %>%
    color(~ Performance == "Good", color = "#70C404", ~Performance) %>%
    color(~ Performance == "Great", color = "#0F9115", ~Performance) %>%
    width(width = 1.25) %>%
    align(align = "center", part = "all") %>%
    bold(i = 4, part = "body")


  ### 绘制Total score 历史分数排名散点图
  p_rank_scatter_plot <- ggplot(data = qc_result$rank_table) +
    # 添加四个区域
    geom_rect(aes(xmin = 1, xmax = 5.64, ymin = -Inf, ymax = Inf), fill = "#B80D0D", alpha = 0.08) +
    geom_rect(aes(xmin = 5.64, xmax = 7.13, ymin = -Inf, ymax = Inf), fill = "#D97C11", alpha = 0.08) +
    geom_rect(aes(xmin = 7.13, xmax = 7.95, ymin = -Inf, ymax = Inf), fill = "#70C404", alpha = 0.08) +
    geom_rect(aes(xmin = 7.95, xmax = 10, ymin = -Inf, ymax = Inf), fill = "#0F9115", alpha = 0.08) +
    # 添加基础点图层
    geom_point(aes(x = Total, y = reorder(batch, Total))) +
    # 突出显示 "QUERIED DATA" 对应的点
    geom_point(
      data = subset(qc_result$rank_table, batch == "QUERIED DATA"),
      aes(x = Total, y = reorder(batch, Total)),
      color = "orange", size = 3
    ) +
    # 自定义x轴刻度
    scale_x_continuous(breaks = c(1, 5.64, 7.13, 7.95, 10)) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    ) +
    labs(
      x = " ",
      y = " ",
      title = "Total Score"
    )

  #### 设置输出文本

  ###### 第一部分
  ### Assessment Summary

  text_1 <- "The performance of the submitted data will be graded as Bad, Fair, Good, or Great based on the ranking by comparing the total score with the historical datasets.The total score is the geometric mean of the scaled values of the number of Signal-to-Noise Ratio (SNR), relative correlation with reference datasets (RC), and recall of DAMs in Reference Datasets (Recall)."
  ### Four levels of performance
  text_1_sup_1 <- "Based on the scaled total score, the submitted data will be ranked together with all Quartet historical datasets. The higher the score, the higher the ranking. After this, the performance levels will be assigned based on their ranking ranges."
  text_1_sup_2 <- "· Bad - the bottom 20%."
  text_1_sup_3 <- "· Fair - between bottom 20% and median 50%."
  text_1_sup_4 <- "· Good - between median 50% and top 20%."
  text_1_sup_5 <- "· Great - the top 20%."

  #### 第二部分 Quality control metric

  ### Performance Score
  text_2 <- "Scores of evaluation metrics for the current batch and all historical batches assessed.Please note that the results shown here are scaled values for all batches in each metric. The name of your data is Queried_Data."
  ### Signal-to-Noise Ratio
  text_3 <- "SNR is established to characterize the power in discriminating multiple groups. The PCA plot is used to visualise the metric."
  ### Correlation with Reference Datasets
  text_4 <- "Relative correlation with reference datasets metric which was representing the numerical consistency of the relative expression profiles."

  ### Method
  supplementary_info_1_1 <- "We apply SNR in the reliability assessment of metabolome data based on the built-in biological differences between Quartet samples. SNR is the fraction of distances between different Quartet samples ('signal') and distances between technical replicates ('noise') on 2D-PCA scatter plot, where a high SNR indicates the tight clustering of technical and wide dispersion of different Quartet samples replicates, as well as good reproducibility and discriminability overall the batch level."
  supplementary_info_1_2 <- "RC is used for assessment of quantitative consistency with the reference datasets (RDs) at relative levels. To evaluate the performance of both targeted and untargeted metabolomics, the RDs were established with historical datasets of high quality by benchmarking the relative abundance values for each sample pair (D5/D6, F7/D6, M8/D6) at metabolite abundance level. We calculate relative abundance values (ratios to D6) of the queried data for metabolites overlapped with the RDs. Then we calculate the Pearson correlation as RC of measured relative abundance values and those in the RDs."
  supplementary_info_1_3 <- "Recall is used for qualitative assessment of the accuracy of biological difference detecting, as the fraction of the differential abundancial metabolites (DAMs) in RDs that are successfully retrieved. Here recall is the number of measured DAMs (p < 0.05, t test) divided by the number of DAMs should be identified as RDs."

  ### Reference
  supplementary_info_ref1 <- "1. Zheng, Y. et al. Multi-omics data integration using ratio-based quantitative profiling with Quartet reference materials. Nature Biotechnology 1–17 (2023)."
  supplementary_info_ref2 <- "2. Zhang, N. et al. Quartet metabolite reference materials for assessing inter-laboratory reliability and data integration of metabolomic profiling. bioRxiv 2022.11. 01.514762 (2022)."

  ### Contact us
  supplementary_info_2_1 <- "Fudan University Pharmacogenomics Research Center"
  supplementary_info_2_2 <- "Project manager: Quartet Team"
  supplementary_info_2_3 <- "Email: quartet@fudan.edu.cn"

  ### Disclaimer
  supplementary_info_3 <- "This quality control report is only for this specific test data set and doesn’t represent an evaluation of the business level of the sequencing company. This report is only used for scientific research, not for clinical or commercial use. We don’t bear any economic and legal liabilities for any benefits or losses (direct or indirect) from using the results of this report."

  read_docx(report_template) %>%
    ## 添加报告标题
    body_add_par(value = "Quartet Report for Metabolomics", style = "heading 1") %>%
    ## 第一部分，Assessment Summary
    body_add_par(value = "Assessment Summary", style = "heading 2") %>%
    body_add_flextable(ft1) %>%
    body_add_break() %>%
    ### 第二部分 Quality control metric
    body_add_par(value = "Quality Control Metric", style = "heading 2") %>%
    body_add_par(value = "Signal-to-Noise Ratio (SNR):", style = "heading 3") %>%
    body_add_par(value = supplementary_info_1_1, style = "Normal") %>%
    body_add_par(value = "Relative Correlation with Reference Datasets (RC):", style = "heading 3") %>%
    body_add_par(value = supplementary_info_1_2, style = "Normal") %>%
    body_add_par(value = "Recall of DAMs in Reference Datasets (Recall):", style = "heading 3") %>%
    body_add_par(value = supplementary_info_1_3, style = "Normal") %>%
    body_add_par(value = "Total Score:", style = "heading 3") %>%
    body_add_par(value = text_1, style = "Normal") %>%
    body_add_par(value = "Performance Category:", style = "heading 3") %>%
    body_add_par(value = text_1_sup_1, style = "Normal") %>%
    body_add_par(value = text_1_sup_2, style = "Normal") %>%
    body_add_par(value = text_1_sup_3, style = "Normal") %>%
    body_add_par(value = text_1_sup_4, style = "Normal") %>%
    body_add_par(value = text_1_sup_5, style = "Normal") %>%
    ### 排名散点图
    body_add_par(value = "Performance Score", style = "heading 2") %>%
    body_add_gg(value = p_rank_scatter_plot, style = "centered") %>%
    body_add_par(value = text_2, style = "Normal") %>%
    ## 分页
    body_add_break() %>%
    ## 信噪比
    body_add_par(value = "Signal-to-Noise Ratio", style = "heading 2") %>%
    body_add_gg(qc_result$pcaplot, style = "centered") %>%
    body_add_par(value = text_3, style = "Normal") %>%
    ## 分页
    body_add_break() %>%
    ## RC
    body_add_par(value = "Correlation with Reference Datasets", style = "heading 2") %>%
    body_add_gg(qc_result$scplot, style = "centered") %>%
    body_add_par(value = text_4, style = "Normal") %>%
    ## 分页
    body_add_break() %>%
    ### 附加信息
    body_add_par(value = "Supplementary Information", style = "heading 2") %>%
    # body_add_par(value = "Method", style = "heading 3") %>%
    # body_add_par(value = supplementary_info_1_1, style = "Normal") %>%
    # body_add_par(value = supplementary_info_1_2, style = "Normal") %>%
    # body_add_par(value = supplementary_info_1_3, style = "Normal") %>%

    body_add_par(value = "Reference", style = "heading 3") %>%
    body_add_par(value = supplementary_info_ref1, style = "Normal") %>%
    body_add_par(value = supplementary_info_ref2, style = "Normal") %>%
    body_add_par(value = "Contact us", style = "heading 3") %>%
    body_add_par(value = supplementary_info_2_1, style = "Normal") %>%
    body_add_par(value = supplementary_info_2_2, style = "Normal") %>%
    body_add_par(value = supplementary_info_2_3, style = "Normal") %>%
    body_add_par(value = "Disclaimer", style = "heading 3") %>%
    body_add_par(value = supplementary_info_3, style = "Normal") %>%
    ## 输出文件
    print(target = output_file)
}
