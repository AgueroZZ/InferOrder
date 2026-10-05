# Show paired noise conditions side by side without changing saved summaries.
grouped_noise_table <- function(data, keys, values, key_labels, value_labels,
                                caption, digits = 3) {
  stopifnot(length(keys) == length(key_labels),
            length(values) == length(value_labels),
            all(data$noise %in% c('low', 'high')))
  low <- data[data$noise == 'low', c(keys, values), drop = FALSE]
  high <- data[data$noise == 'high', c(keys, values), drop = FALSE]
  row_key <- function(x) do.call(paste, c(x[keys], sep = '\r'))
  stopifnot(!anyDuplicated(row_key(low)), !anyDuplicated(row_key(high)),
            setequal(row_key(low), row_key(high)))
  high <- high[match(row_key(low), row_key(high)), , drop = FALSE]
  wide <- cbind(low, high[values])
  html <- as.character(knitr::kable(
    wide, format = 'html', row.names = FALSE, digits = digits,
    col.names = c(key_labels, value_labels, value_labels), caption = caption,
    table.attr = 'class="table grouped-noise-table"'))
  group_row <- sprintf(
    '<tr><th colspan="%d"></th><th colspan="%d" scope="colgroup" style="text-align:center">Low noise (SNR 16)</th><th colspan="%d" scope="colgroup" style="text-align:center">High noise (SNR 1)</th></tr>',
    length(keys), length(values), length(values))
  knitr::asis_output(paste0('<div style="overflow-x:auto">',
                          sub('<thead>', paste0('<thead>', group_row), html,
                              fixed = TRUE), '</div>'))
}
