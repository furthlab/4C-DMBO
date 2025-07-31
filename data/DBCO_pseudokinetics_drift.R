data<-read.table('kinetics/DBCO-amine/1_BD1-5mM.csv', header=T, skip=15, sep=',')

plot(data$Time.s., data$Absorbance, type='l', ylab='Absorbance (307 nm)', xlab='Time (sec.)', las=1,lwd=2)


# Estimate initial parameters
A0_start <- max(data$Absorbance) - min(data$Absorbance)
k_start <- 0.001  # guess based on timescale
Ainf_start <- min(data$Absorbance)

# Fit pseudo-first-order kinetics
fit <- nls(Absorbance ~ A0 * exp(-k * Time.s.) + Ainf,
           data = data,
           start = list(A0 = A0_start, k = k_start, Ainf = Ainf_start))

# Summary of fit
summary(fit)

# Add fitted curve to the plot
lines(data$Time.s., predict(fit), col='red', lwd=2)
legend("topright", legend=c("Data", "Fit"), col=c("black", "red"), lwd=2)

summary(fit)$coefficients

k <- summary(fit)$coefficients["k", "Estimate"]
print(k)



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


par(mfrow=c(3,2))
k_results <- extract_k_values("kinetics/DBCO-amine", drift=TRUE, plot_results = TRUE)

print(k_results)

# Molarity of azide across conditions
molarity <- c(0.005, 0.015, 0.030, 0.050, 0.1)  # in M

# Ensure k_results is a numeric vector of the same length
k_values <- k_results$k

# Simulate pseudo-first-order rate constants with some noise
#set.seed(42)  # for reproducibility
#true_k2 <- 25  # M^-1 s^-1
#k_values <- true_k2 * molarity + rnorm(length(molarity), mean = 0, sd = 0.5)


# Plot pseudo-first-order rate constants vs azide concentration
plot(molarity, k_values,
     xlab = "[Azide] (M)", ylab = expression(k[obs]~"(" * s^{-1} * ")"),
     pch = 16, las = 1, col = "blue", main = "Second-order kinetics fit")

# Fit linear model: k_obs = k2 * [azide]
fit_lm <- lm(k_values ~ molarity)

# Add regression line
abline(fit_lm, col = "red", lwd = 2)

# Extract second-order rate constant (slope)
k2 <- coef(fit_lm)[["molarity"]]
print(paste("Second-order rate constant k2 =", signif(k2, 4), "M^-1 s^-1"))

# Optional: Add R^2 to plot
r2 <- summary(fit_lm)$r.squared
# Legend with math-style units
legend("topleft", legend = c(
  bquote(k[2] == .(signif(k2, 4)) ~ M^{-1} %.% s^{-1}),
  bquote(R^2 == .(signif(r2, 3)))
), bty = "n", text.col = "darkred")







#_________________

plot_all_uvvis <- function(folder_path) {
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
    
    # Try fitting, handle possible errors
    fit_success <- try({
      fit <- nls(Absorbance ~ A0 * exp(-k * Time.s.) + Ainf,
                 data = data,
                 start = list(A0 = A0_start, k = k_start, Ainf = Ainf_start))
      
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

plot_all_uvvis('kinetics/DBCO-amine')


