# Assuming you have a folder with all the data.


#' Extract Pseudo-First-Order Rate Constants from UV-Vis Data
#'
#' This function loops through all `.csv` files in a given folder,
#' fits a pseudo-first-order exponential decay model to UV-Vis absorbance data,
#' and extracts the fitted rate constant \eqn{k} from each file.
#'
#' The model fitted is: \deqn{A(t) = A_0 \cdot e^{-k t} + A_{\infty}}
#'
#' The data files are assumed to have a common format:
#' - CSV file
#' - First 15 lines skipped
#' - Columns named `Time.s.` and `Absorbance`
#'
#' @param folder_path Character string. Path to the folder containing the `.csv` data files.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{File}{Name of the file}
#'   \item{k}{Estimated pseudo-first-order rate constant}
#' }
#'
#' @examples
#' \dontrun{
#' # Extract rate constants from all UV-Vis traces in folder
#' k_results <- extract_k_values("raw-data/UV-Vis/")
#' print(k_results)
#' }
#'
#' @importFrom stats nls coef
#' @export
extract_k_values <- function(folder_path, drift = FALSE, plot_results = TRUE) {
  # List all CSV files in the folder
  files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)

  # Initialize results list
  results <- data.frame(File = character(), k = numeric(), stringsAsFactors = FALSE)

  for (file in files) {
    # Read the data (skip header lines and assume same format)
    data <- tryCatch({
      read.table(file, header = TRUE, skip = 15, sep = ",")
    }, error = function(e) return(NULL))

    if (is.null(data) || !"Time.s." %in% names(data) || !"Absorbance" %in% names(data)) next

    # Estimate initial parameters
    A0_start <- max(data$Absorbance) - min(data$Absorbance)
    k_start <- 0.005
    Ainf_start <- min(data$Absorbance)
    m_start <- 0  # baseline slope for drift model

    # Choose formula based on drift argument
    if (drift) {
      formula_fit <- Absorbance ~ A0 * exp(-k * Time.s.) + Ainf + m * Time.s.
      start_list <- list(A0 = A0_start, k = k_start, Ainf = Ainf_start, m = m_start)
    } else {
      formula_fit <- Absorbance ~ A0 * exp(-k * Time.s.) + Ainf
      start_list <- list(A0 = A0_start, k = k_start, Ainf = Ainf_start)
    }

    # Fit the model
    fit <- tryCatch({
      nls(formula_fit,
          data = data,
          start = start_list,
          control = nls.control(warnOnly = TRUE))
    }, error = function(e) return(NULL))

    # Extract k value if fitting succeeded
    if (!is.null(fit)) {
      k_value <- summary(fit)$coefficients["k", "Estimate"]
      results <- rbind(results, data.frame(File = basename(file), k = k_value))

      # Plot if requested
      if (plot_results) {
        plot(data$Time.s., data$Absorbance, type = 'l',
             ylab = 'Absorbance (307 nm)', xlab = 'Time (sec.)',
             las = 1, lwd = 2,
             main = paste("Fit for", basename(file)))
        lines(data$Time.s., predict(fit), col = 'red', lwd = 2)
        legend("topright", legend = c("Data", "Fit"), col = c("black", "red"), lwd = 2)


      }
    }
  }

  return(results)
}
