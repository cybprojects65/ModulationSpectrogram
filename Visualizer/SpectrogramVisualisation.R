library(ggplot2)
library(data.table)

# Run in R/RStudio after setting this path, or use:
# Rscript SpectrogramVisualisation.R input.csv [output.png]
args <- commandArgs(trailingOnly = TRUE)
file <- if (length(args) >= 1L) args[1] else "../samples/PS14Audio_noisy_speech_modspec.csv"

# These settings must match the Java extraction parameters.
nfeatures <- 8L
minfreq <- 165.4
maxfreq <- 3000
feature_names <- paste0("F", 0:(nfeatures - 1L))

# New layout: one frame per row, two timestamps, then base/delta/double-delta.
data <- fread(file)
required <- c("Time_start_s", "Time_end_s", feature_names)
if (!all(required %in% names(data))) {
  stop("Missing CSV columns: ", paste(setdiff(required, names(data)), collapse = ", "))
}
if (nrow(data) == 0L) stop("The CSV contains no frames.")
if (!all(vapply(data[, ..required], is.numeric, logical(1)))) {
  stop("Timestamps and MS features must be numeric.")
}
mat <- as.matrix(data[, ..feature_names])  # frames x bands; exclude derivatives
start_s <- data$Time_start_s
end_s <- data$Time_end_s
if (any(!is.finite(mat)) || any(!is.finite(start_s)) || any(!is.finite(end_s))) {
  stop("Timestamps and MS features must contain only finite values.")
}
if (any(end_s <= start_s) || any(diff(start_s) <= 0)) {
  stop("Frame starts must increase and each end must follow its start.")
}

# Preserve the previous palette: lower three quartiles white; divide the
# upper-quartile VALUE RANGE into four equal intervals, light to dark red.
qs <- quantile(mat, probs = c(0, .25, .5, .75, 1), names = FALSE)
upper_seq <- seq(qs[4], qs[5], length.out = 5)
colors <- c("white", "white", "white", "#ffcccc", "#ff6666", "#ff1a1a", "#990000")
cat("Quantiles:", qs, "\nUpper-range boundaries:", upper_seq, "\n")
colour_values <- function(values) {
  result <- rep("white", length(values))
  # A constant upper range has no contrast: leave it white.
  if (qs[5] > qs[4]) {
    for (j in 1:4) result[values >= upper_seq[j]] <- colors[j + 3L]
  }
  result
}

# Java exports F0 as the highest band. Reverse bands for a low-to-high y axis.
# Labels reflect the Mel-spaced centers for corrected Java indices 1..N.
freq_to_mel <- function(f) 2595 * log10(1 + f / 700)
mel_to_freq <- function(m) 700 * (10^(m / 2595) - 1)
centers_hz <- mel_to_freq(seq(freq_to_mel(minfreq), freq_to_mel(maxfreq),
                              length.out = nfeatures + 2L)[2:(nfeatures + 1L)])
mat_low_to_high <- mat[, nfeatures:1, drop = FALSE]

# Display each analysis value over one hop, not over its overlapping 250 ms
# analysis window. X coordinates refer to window START times, as in the old plot.
windowshift <- if (length(start_s) > 1L) median(diff(start_s)) else 0.0125
right_s <- c(start_s[-1L], tail(start_s, 1) + windowshift)
n_y <- nfeatures + 2L  # preserve the old white boundary rows
plot_values <- as.vector(t(cbind(NA_real_, mat_low_to_high, NA_real_)))
fill <- rep("white", length(plot_values))
valid <- is.finite(plot_values)
fill[valid] <- colour_values(plot_values[valid])

dt <- data.table(
  xmin = rep(start_s, each = n_y),
  xmax = rep(right_s, each = n_y),
  y = rep(seq_len(n_y), times = nrow(data)),
  fill = fill
)
freq_hz_labels <- round(c(minfreq, centers_hz, maxfreq), 1)
tick_step <- 0.25
x_breaks <- seq(ceiling(min(start_s) / tick_step) * tick_step,
                max(right_s), by = tick_step)

p <- ggplot(dt) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = y - .5, ymax = y + .5,
                fill = fill), colour = NA) +
  scale_fill_identity() +
  scale_x_continuous(name = "Window start time (s)", breaks = x_breaks,
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Acoustic band frequency (Hz)",
                     breaks = seq_len(n_y), labels = freq_hz_labels,
                     expand = c(0, 0)) +
  labs(title = "Modulation Spectrogram (4 Hz)",
       subtitle = sprintf("%.1f ms analysis windows; %.1f ms hop",
                          1000 * median(end_s - start_s), 1000 * windowshift)) +
  theme_minimal(base_size = 14) +
  theme(panel.ontop = TRUE,
        panel.grid.major.x = element_line(color = "grey50", linewidth = 0.4),
        panel.grid.major.y = element_line(color = "grey70", linewidth = 0.3),
        panel.grid.minor = element_blank())

print(p)
if (length(args) >= 2L) {
  ggsave(args[2], plot = p, width = 11, height = 5, dpi = 300, bg = "white")
}
