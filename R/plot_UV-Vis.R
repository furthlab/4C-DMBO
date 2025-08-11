#' Plot and Fit UV-Vis Absorbance Decay Curves from CSV Files
#'
#' This function reads all `.csv` files in a given folder containing UV-Vis
#' absorption data (with time and absorbance columns), plots the absorbance
#' over time, and fits an exponential decay model to estimate the rate constant \code{k}.
#'
#' @param folder_path Character string. Path to the folder containing CSV files.
#'        Each CSV is expected to have a header and contain absorbance vs time data,
#'        with absorbance at 307 nm in the column \code{Absorbance} and time in \code{Time.s.}.
#' @param drift Logical. If \code{TRUE}, a linear drift term is included in the model
#'        (\code{A0 * exp(-k * t) + Ainf + m * t}); otherwise, the drift term is excluded.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Automatically lays out subplots based on the number of CSV files.
#'   \item Fits an exponential decay model to each absorbance trace.
#'   \item Adds a fitted curve (in red) and the estimated decay constant \code{k} to each plot.
#'   \item Displays a message if the nonlinear fit fails for a given file.
#' }
#'
#' The function assumes data starts at row 16 (i.e., skips 15 header lines).
#'
#' @return This function is called for its side effect: it creates a series of base R plots.
#' No value is returned.
#'
#' @examples
#' \dontrun{
#' plot_all_uvvis("data/uvvis_measurements/", drift = FALSE)
#' }
#'
#' @export
plot_all_uvvis <- function(folder_path, drift = TRUE) {
  # Get all CSV files in the folder
  files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)

  # Decide grid layout (square-ish)
  n <- length(files)
  ncol <- ceiling(sqrt(n))
  nrow <- ceiling(n / ncol)

  # Set up plotting window
  par(mfrow = c(nrow, ncol), mar = c(4, 4, 2, 1)) # adjust margins if needed


  for (f in files) {
    data <- read.table(f, header = TRUE, skip = 15, sep = ",")

    plot(data$Time.s., data$Absorbance,
         type = 'l', ylab = 'Absorbance (307 nm)', xlab = 'Time (sec.)',
         las = 1, lwd = 2, main = basename(f))

    # Estimate initial parameters
    A0_start <- max(data$Absorbance) - min(data$Absorbance)
    k_start <- 0.001
    Ainf_start <- min(data$Absorbance)
    m_start <- 0

    # Choose formula based on drift argument
    if (drift) {
      formula_fit <- Absorbance ~ A0 * exp(-k * Time.s.) + Ainf + m * Time.s.
      start_list <- list(A0 = A0_start, k = k_start, Ainf = Ainf_start, m = m_start)
    } else {
      formula_fit <- Absorbance ~ A0 * exp(-k * Time.s.) + Ainf
      start_list <- list(A0 = A0_start, k = k_start, Ainf = Ainf_start)
    }

    # Try fitting, handle possible errors
    fit_success <- try({
      fit <- nls(formula_fit,
                 data = data,
                 start = start_list)

      # Add fitted line
      lines(data$Time.s., predict(fit), col = 'red', lwd = 2)

      # Extract fitted k value
      k_val <- coef(fit)["k"]

      # Add legend with k value
      legend("topright",
             legend = c("Data", "Fit", sprintf("k = %.4f", k_val)),
             col = c("black", "red", NA),
             lwd = c(2, 2, NA),
             bty = "n")

    }, silent = TRUE)

    if (inherits(fit_success, "try-error")) {
      mtext("Fit failed", side = 3, line = -1, col = "red")
    }
  }

  # Reset plotting layout
  par(mfrow = c(1, 1))
}

plot_all_uvvis('data/kinetics/DBCO-amine')
