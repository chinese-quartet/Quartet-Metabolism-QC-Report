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
#' @importFrom flextable theme_box
#' @importFrom flextable color
#' @importFrom flextable set_caption
#' @importFrom flextable align
#' @importFrom flextable width
#' @importFrom flextable bold
#' @importFrom flextable bg
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

  # --- 1. 定义中文文本内容 ---
  
  # 摘要
  text_sum_intro <- "本报告基于多项组学关键质量控制指标，总结了 Quartet Metabolite 参考物质所生成数据的质量情况。质量控制流程从用户输入代谢物表达矩阵开始，分别计算外部质控品的信噪比（Signal-to-Noise Ratio, SNR）、与参考数据集的Pearson相关系数 (Pearson correlation coefficient, PCC) 及整体质量判断。"
  
  # 质量控制指标定义
  # SNR
  text_snr_title <- "信噪比（Signal-to-Noise Ratio, SNR）"
  text_snr_desc <- "在代谢组数据的可靠性评估中，SNR 基于 Quartet 样本之间内在的生物学差异进行计算。具体而言，SNR 定义为二维 PCA 散点图中，不同 Quartet 样本之间距离（“信号”）与技术重复之间距离（“噪声”）的比值。较高的 SNR 表明技术重复聚类更紧密、不同 Quartet 样本之间分离度更高，反映出该批次在整体层面具有良好的重复性和区分能力。"
  
  # RC
  text_rc_title <- "与参考数据集的Pearson相关系数 (Pearson correlation coefficient, PCC) "
  text_rc_desc <- "用于评估测试数据在相对定量层面与参考数据集（Reference Datasets, RDs）之间的一致性。为同时适用于靶向与非靶向代谢组学分析，参考数据集基于历史高质量数据构建，通过比较各样本对（D5/D6、F7/D6、M8/D6）在代谢物丰度水平上的相对丰度值。评估时，首先计算测试数据中与参考数据集重叠代谢物相对于 D6 的丰度比值，然后计算这些相对丰度值与参考数据集中对应数值之间的 Pearson 相关系数。"
  
  # 参考文献
  text_ref_title <- "参考文献"
  text_ref_1 <- "1. Zheng Y, Liu Y, Yang J, et al. Multi-omics data integration using ratio-based quantitative profiling with Quartet reference materials. Nature Biotechnology, 2024."
  text_ref_2 <- "2. Zhang N, Chen Q, Zhang P, et al. Quartet metabolite reference materials for inter-laboratory proficiency test and data integration of metabolomics profiling. Genome Biology, 2024."
  text_ref_3 <- "3. 上海临床队列组学检测工作指引（征求意见稿）, 2025/11/26."
  
  # 免责声明
  text_disclaimer_title <- "免责声明"
  text_disclaimer_content <- "本数据质量报告仅针对所评估的特定数据集提供分析结果，仅供信息参考之用。尽管已尽最大努力确保分析结果的准确性和可靠性，但本报告按“现状（AS IS）”提供，不附带任何形式的明示或暗示担保。报告作者及发布方不对基于本报告内容所采取的任何行动承担责任。本报告中的结论不应被视为对任何产品或流程质量的最终判定，也不应用于关键应用场景、商业决策或法规合规用途，除非经过专业核查和独立验证。对于分析结果的正确性、准确性、可靠性或适用性，不作任何明示或暗示的保证。"
  
  # --- 2. 创建符合新格式的表格 ---
  
  # 从 conclusion_table 中提取数据
  # 注意：performance.R 生成的 conclusion_table 包含 "Quality metrics" 和 "Value" 列
  raw_table <- qc_result$conclusion_table
  
  # 提取 SNR 和 RC 的数值
  # 假设 conclusion_table 的第一列是指标名称，第二列是数值
  # 使用模糊匹配确保鲁棒性
  snr_row <- raw_table[grep("Signal-to-Noise Ratio", raw_table$`Quality metrics`), ]
  rc_row  <- raw_table[grep("Relative Correlation", raw_table$`Quality metrics`), ]
  
  snr_val <- as.numeric(snr_row$Value)
  rc_val  <- as.numeric(rc_row$Value)
  
  # 定义判断逻辑 (参考标准: SNR >= 10, RC >= 0.80)
  batch_name_str <- "Queried_Data"
  
  # 格式化显示字符串 (如果未达标，添加向下箭头 ↓)
  snr_str <- sprintf("%.2f", snr_val)
  if (!is.na(snr_val) && snr_val < 10) {
    snr_str <- paste0(snr_str, " ↓")
  }
  
  rc_str <- sprintf("%.3f", rc_val)
  if (!is.na(rc_val) && rc_val < 0.80) {
    rc_str <- paste0(rc_str, " ↓")
  }
  
  # 整体质量判断
  # is_pass <- (!is.na(snr_val) && snr_val >= 10) && (!is.na(rc_val) && rc_val >= 0.80)
  is_pass <- (!is.na(snr_val) && snr_val >= 10)
  quality_str <- ifelse(is_pass, "全部通过", "No")
  
  # 构建新数据框
  new_df <- data.frame(
    "样本组" = c("推荐质量标准", batch_name_str),
    "信噪比" = c("≥10", snr_str),
    "Pearson相关系数" = c("≥0.80", rc_str),
    "整体质量" = c("全部通过", quality_str),
    check.names = FALSE
  )
  
  # 生成 Flextable 样式
  ft1 <- flextable(new_df) %>%
    theme_box() %>%                                   # 基础边框主题
    align(align = "center", part = "all") %>%         # 全局居中
    width(width = 1.5) %>%                            # 列宽
    bold(part = "header") %>%                         # 表头加粗
    bold(i = 1, part = "body") %>%                    # 第一行(推荐标准)加粗
    bg(part = "header", bg = "#EFEFEF") %>%           # 表头背景色
    # 动态上色：如果整体质量是 No，标红
    color(i = 2, j = "整体质量", color = ifelse(quality_str == "No", "#B80D0D", "black")) %>%
    # 动态上色：如果指标未达标，标红
    color(i = 2, j = "Pearson相关系数", color = ifelse(rc_val < 0.80, "#B80D0D", "black")) %>%
    color(i = 2, j = "信噪比", color = ifelse(snr_val < 10, "#B80D0D", "black"))
  
  
  # --- 3. 生成报告文档 ---
  
  read_docx(report_template) %>%
    # 标题
    body_add_par(value = "Quartet代谢组质量报告", style = "heading 1") %>%
    
    # 摘要
    body_add_par(value = "摘要", style = "heading 2") %>%
    body_add_par(value = text_sum_intro, style = "Normal") %>%
    body_add_par(value = " ", style = "Normal") %>% # 空行
    
    # 插入表格
    body_add_flextable(ft1) %>%
    # body_add_break() %>%
    
    # 质量控制指标说明
    body_add_par(value = "质量控制指标", style = "heading 2") %>%
    
    # SNR 定义
    body_add_par(value = text_snr_title, style = "heading 3") %>%
    body_add_par(value = text_snr_desc, style = "Normal") %>%
    
    # RC 定义
    body_add_par(value = text_rc_title, style = "heading 3") %>%
    body_add_par(value = text_rc_desc, style = "Normal") %>%
    
    # 参考文献
    body_add_par(value = text_ref_title, style = "heading 2") %>%
    body_add_par(value = text_ref_1, style = "Normal") %>%
    body_add_par(value = text_ref_2, style = "Normal") %>%
    body_add_par(value = text_ref_3, style = "Normal") %>%
    body_add_par(value = " ", style = "Normal") %>%
    
    # 免责声明
    body_add_par(value = text_disclaimer_title, style = "heading 3") %>%
    body_add_par(value = text_disclaimer_content, style = "Normal") %>%
    body_add_break() %>%
    
    # 插入图片
    
    # SNR Plot
    body_add_par(value = "Signal-to-Noise Ratio", style = "heading 2") %>%
    body_add_gg(value = qc_result$pcaplot, style = "centered") %>%
    
    # RC Plot
    body_add_break() %>%
    body_add_par(value = "Correlation with Reference Datasets", style = "heading 2") %>%
    body_add_gg(value = qc_result$scplot, style = "centered") %>%
    
    # 输出文件
    print(target = output_file)
}
