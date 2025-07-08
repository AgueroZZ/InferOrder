subsample_cell_types <- function (x, n = 1000) {
  cells <- NULL
  groups <- levels(x)
  for (g in groups) {
    i  <-  which(x == g)
    n0 <- min(n,length(i))
    i  <- sample(i,n0)
    cells <- c(cells,i)
  }
  return(sort(cells))
}

plot_structure <- function(loadings,
                           order   = NULL,
                           palette = NULL) {
  # 检查
  if (is.null(rownames(loadings))) {
    stop("loadings 矩阵必须带有行名，作为 Sample ID。")
  }

  # 转成 data.frame 并加上 Sample 列
  df <- as.data.frame(loadings)
  df <- tibble::rownames_to_column(df, var = "Sample")

  # 给列命名为 Factor1…FactorK
  K <- ncol(loadings)
  colnames(df)[-1] <- paste0("Factor", seq_len(K))

  # reshape 为 long 格式
  df_long <- tidyr::pivot_longer(
    df,
    cols      = paste0("Factor", seq_len(K)),
    names_to  = "Factor",
    values_to = "Loading"
  )

  # 如果用户指定了 order，就把 Sample 转为因子并设定水平
  if (!is.null(order)) {
    if (!all(order %in% df_long$Sample)) {
      stop("提供的 order 必须包含所有 Sample 名称。")
    }
    df_long$Sample <- factor(df_long$Sample, levels = order)
  }

  # 绘图
  p <- ggplot(df_long, aes(x = Sample, y = Loading, fill = Factor)) +
    geom_bar(stat = "identity", width = 1) +
    scale_y_continuous(expand = expansion(mult = c(0, .05))) +
    labs(
      x     = NULL,
      y     = "Loading (membership)",
      title = "Structure–style Plot of Loadings"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  if (!is.null(palette)) {
    p <- p + scale_fill_manual(values = palette)
  }

  return(p)
}

plot_highlight_types <- function(type_vec,
                                 subset_types,
                                 ordering_metric = NULL,
                                 order_vec       = NULL,
                                 other_color     = "grey80") {

  # 样本数
  n <- length(type_vec)

  # 检查 type_vec
  if (is.null(names(type_vec))) {
    stop("type_vec 必须是命名向量，names(type_vec)=样本 ID。")
  }

  if (!is.null(ordering_metric)) {
    # 用 ordering_metric 定义样本和位置
    if (is.null(names(ordering_metric))) {
      stop("ordering_metric 必须是命名向量，names(ordering_metric)=样本 ID。")
    }
    metric <- ordering_metric
    # 只保留在 type_vec 中的样本
    common <- intersect(names(metric), names(type_vec))
    metric <- metric[common]
    types  <- as.character(type_vec[common])
    df <- data.frame(
      Sample = names(metric),
      Metric = as.numeric(metric),
      Type   = types,
      stringsAsFactors = FALSE
    )
    # 样本按 Metric 排序
    df <- df[order(df$Metric), ]
  } else {
    # 没有 ordering_metric，使用 order_vec 定义顺序和等间距度量
    if (is.null(order_vec)) {
      stop("必须提供 order_vec（或 ordering_metric）。")
    }
    if (!all(order_vec %in% names(type_vec))) {
      stop("order_vec 中的所有样本 ID 必须出现在 type_vec 中。")
    }
    df <- data.frame(
      Sample = order_vec,
      Metric = seq_len(n) / n,
      Type   = as.character(type_vec[order_vec]),
      stringsAsFactors = FALSE
    )
  }

  # 标记高亮
  df$Highlight <- ifelse(df$Type %in% subset_types, df$Type, "Other")
  df$Highlight <- factor(df$Highlight, levels = c(subset_types, "Other"))

  # 配色
  cols <- c(scales::hue_pal()(length(subset_types)), other_color)
  names(cols) <- levels(df$Highlight)

  # 绘图
  ggplot(df, aes(x = Metric, y = 1, fill = Highlight)) +
    geom_col(width = if (is.null(ordering_metric)) 1/n else diff(range(df$Metric)) / n) +
    scale_fill_manual(values = cols) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Ordering Metric", y = NULL,
         title = "Samples Highlighted by Specified Types") +
    theme_minimal() +
    theme(
      axis.text.x      = element_text(size = 10),
      axis.ticks.x     = element_line(),
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      panel.grid       = element_blank(),
      legend.title     = element_blank()
    )
}

library(ggplot2)

#’ 绘制各类型在 Ordering Metric 或 order_vec 上的分布（核密度或小提琴图）
#’
#’ @param type_vec        命名字符向量，names(type_vec)=样本 ID，values=对应的类别。
#’ @param subset_types    要可视化分布的类别向量。
#’ @param ordering_metric 可选命名数值向量，names=样本 ID，对应的度量值。
#’ @param order_vec       可选字符向量，指定样本 ID 顺序；当 ordering_metric=NULL 且提供了 order_vec 时，会用等间距 seq(1/n,2/n,…,1) 代替 metric。
#’ @param density         逻辑值，TRUE 时画叠加的核密度曲线，FALSE 时画小提琴图。
#’ @param other_color     “Other” 类别在图中的颜色（默认灰色），此版本中不直接绘制 Other。
#’ @return ggplot2 对象
distribution_highlight_types <- function(type_vec,
                                         subset_types,
                                         ordering_metric = NULL,
                                         order_vec       = NULL,
                                         density         = TRUE,
                                         other_color     = "grey80") {
  # 检查 type_vec
  if (is.null(names(type_vec))) {
    stop("type_vec 必须是命名向量，names(type_vec)=样本 ID。")
  }

  # 如果有 ordering_metric，优先使用；否则尝试用 order_vec
  if (!is.null(ordering_metric)) {
    if (is.null(names(ordering_metric))) {
      stop("ordering_metric 必须是命名数值向量，names=样本 ID。")
    }
    common <- intersect(names(ordering_metric), names(type_vec))
    metric <- ordering_metric[common]
    types  <- as.character(type_vec[common])
  } else {
    if (is.null(order_vec)) {
      stop("必须提供 ordering_metric 或 order_vec。")
    }
    if (!all(order_vec %in% names(type_vec))) {
      stop("order_vec 中的所有样本 ID 必须出现在 type_vec 中。")
    }
    n      <- length(order_vec)
    common <- order_vec
    metric <- seq_len(n) / n
    names(metric) <- order_vec
    types  <- as.character(type_vec[order_vec])
  }

  # 构造数据框，且只保留 subset_types
  df <- data.frame(
    Metric = metric[common],
    Type   = factor(types, levels = subset_types),
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$Type), , drop=FALSE]

  # 配色
  cols <- scales::hue_pal()(length(subset_types))
  names(cols) <- subset_types

  if (density) {
    # 叠加核密度曲线
    p <- ggplot(df, aes(x = Metric, color = Type, fill = Type)) +
      geom_density(alpha = 0.3, size = 1) +
      scale_color_manual(values = cols) +
      scale_fill_manual(values = cols) +
      labs(x = "Ordering Metric", y = "Density",
           title = "Density of Ordering Metric by Type") +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid  = element_blank()
      )
  } else {
    # 小提琴图
    p <- ggplot(df, aes(x = Type, y = Metric, fill = Type)) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_jitter(width = 0.1, size = 0.5, alpha = 0.4) +
      scale_fill_manual(values = cols) +
      labs(x = "Type", y = "Ordering Metric",
           title = "Violin Plot of Ordering Metric by Type") +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid  = element_blank(),
        legend.position = "none"
      )
  }

  return(p)
}

