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
                           plot_name = NULL,
                           order   = NULL,
                           palette = NULL) {
  if (is.null(rownames(loadings))) {
    stop("loadings must have row names representing Sample IDs.")
  }

  df <- as.data.frame(loadings)
  df <- tibble::rownames_to_column(df, var = "Sample")

  K <- ncol(loadings)
  colnames(df)[-1] <- paste0("Factor", seq_len(K))

  df_long <- tidyr::pivot_longer(
    df,
    cols      = paste0("Factor", seq_len(K)),
    names_to  = "Factor",
    values_to = "Loading"
  )

  if (!is.null(order)) {
    if (!all(order %in% df_long$Sample)) {
      stop("Order must contain all Sample names in the data.")
    }
    df_long$Sample <- factor(df_long$Sample, levels = order)
  }

  p <- ggplot(df_long, aes(x = Sample, y = Loading, fill = Factor)) +
    geom_bar(stat = "identity", width = 1) +
    scale_y_continuous(expand = expansion(mult = c(0, .05))) +
    labs(
      x     = NULL,
      y     = "Loading (membership)",
      title = plot_name
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

  n <- length(type_vec)

  if (is.null(names(type_vec))) {
    stop("type_vec must be a named vector with names as Sample IDs.")
  }

  if (!is.null(ordering_metric)) {
    if (is.null(names(ordering_metric))) {
      stop("ordering_metric must be a named vector with names as Sample IDs.")
    }
    metric <- ordering_metric
    common <- intersect(names(metric), names(type_vec))
    metric <- metric[common]
    types  <- as.character(type_vec[common])
    df <- data.frame(
      Sample = names(metric),
      Metric = as.numeric(metric),
      Type   = types,
      stringsAsFactors = FALSE
    )
    df <- df[order(df$Metric), ]
  } else {
    if (is.null(order_vec)) {
      stop("Must provide order_vec (or ordering_metric).")
    }
    if (!all(order_vec %in% names(type_vec))) {
      stop("order_vec must contain all Sample IDs in type_vec.")
    }
    df <- data.frame(
      Sample = order_vec,
      Metric = seq_len(n) / n,
      Type   = as.character(type_vec[order_vec]),
      stringsAsFactors = FALSE
    )
  }

  df$Highlight <- ifelse(df$Type %in% subset_types, df$Type, "Other")
  df$Highlight <- factor(df$Highlight, levels = c(subset_types, "Other"))

  cols <- c(scales::hue_pal()(length(subset_types)), other_color)
  names(cols) <- levels(df$Highlight)

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
distribution_highlight_types <- function(type_vec,
                                         subset_types,
                                         ordering_metric = NULL,
                                         order_vec       = NULL,
                                         density         = TRUE,
                                         other_color     = "grey80") {
  if (is.null(names(type_vec))) {
    stop("type_vec must be a named vector with names as Sample IDs.")
  }

  if (!is.null(ordering_metric)) {
    if (is.null(names(ordering_metric))) {
      stop("ordering_metric must be a named vector with names as Sample IDs.")
    }
    common <- intersect(names(ordering_metric), names(type_vec))
    metric <- ordering_metric[common]
    types  <- as.character(type_vec[common])
  } else {
    if (is.null(order_vec)) {
      stop("Must provide ordering_metric or order_vec.")
    }
    if (!all(order_vec %in% names(type_vec))) {
      stop("order_vec must contain all Sample IDs in type_vec.")
    }
    n      <- length(order_vec)
    common <- order_vec
    metric <- seq_len(n) / n
    names(metric) <- order_vec
    types  <- as.character(type_vec[order_vec])
  }

  df <- data.frame(
    Metric = metric[common],
    Type   = factor(types, levels = subset_types),
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$Type), , drop=FALSE]

  cols <- scales::hue_pal()(length(subset_types))
  names(cols) <- subset_types

  if (density) {
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

