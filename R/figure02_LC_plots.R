# Load LC-MS data from text files, skipping instrument metadata lines (first 43 lines)

read_lcms_data <- function(file_path) {
  data <- read.table(file_path, skip = 43)
  names(data) <- c("Time.min", "Step.s", "Value.mAU")
  data$Value.mAU <- gsub(",", "", data$Value.mAU)
  return(data)
}

# Read individual datasets
dmbo_photo <- read_lcms_data("data/reactivity_assay/photo-4C-DMBO.txt")
dmbo_alkyne <- read_lcms_data("data/reactivity_assay/4C-DMBO_alkyne.txt")
spaac_paz <- read_lcms_data("data/reactivity_assay/4C-DMBO_SPAAC_48h.txt")
ieddac_MeTz <- read_lcms_data("data/reactivity_assay/4C-DMBO_iEDDAC_48h.txt")

# Function to create consistent LC-MS chromatogram plots
plot_chromatogram <- function(data, xlab = "", ylab = "Abs. 254 nm (a.u)", xlim = c(0.4, 3), show_xaxis = FALSE, col="#66c2a5
", lwd=2) {
  plot(data$Time.min, data$Value.mAU, type = "l", xlim = xlim, xlab = xlab, ylab = ylab, lwd=lwd, col=col,
       las = 1, axes = FALSE)
  if (show_xaxis) {
    axis(1, labels = c(1, 2, 3), at = c(1, 2, 3))
  } else {
    axis(1, labels = FALSE, at = c(1, 2, 3))
  }
  box()
}

### Figure 3, Panel B ###
quartz(width=260/30, height=150/30)
par(mfrow = c(3, 1), mar = c(4, 4, 1, 1))  # 3 stacked plots

plot_chromatogram(dmbo_photo, col="#66c2a5")
plot_chromatogram(dmbo_alkyne, col="#fc8d62")
plot_chromatogram(spaac_paz, xlab = "Time (min)", show_xaxis = TRUE, col="#8da0cb")

quartz.save(file="figures/pdf/SPAAC_MeOH.pdf", type='pdf')

### Figure 3, Panel C ###
quartz(width=260/30, height=150/30)
par(mfrow = c(3, 1), mar = c(4, 4, 1, 1))  # Reset for next panel

plot_chromatogram(dmbo_photo,  col="#66c2a5")
plot_chromatogram(dmbo_alkyne, col="#fc8d62")
plot_chromatogram(ieddac_MeTz, xlab = "Time (min)", show_xaxis = TRUE,  col="#8da0cb")

quartz.save(file="figures/pdf/iEDDA_MeOH.pdf", type='pdf')

