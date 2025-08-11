# Load and plot absorbance spectra for DMBO and DBCO derivatives

# --- Helper Function to Load Data ---
load_spectrum <- function(path) {
  read.table(path, skip = 12, sep = ',')
}

# --- Helper Function to Add Marker and Vertical Line at λ = 307 nm ---
mark_wavelength <- function(x, y, col) {
  # Draw a vertical dashed line from (x, 0) to (x, y)
  lines(c(x, x), c(0, y), lty = 2)

  # Add a colored point at (x, y)
  points(x, y, pch = 21, bg = col, cex = 2)
}

# --- Plot: 4C-DMBO Series ---
# Load data
dmbo_photo     <- load_spectrum('data/wavescan/4C-DMBO/4C_WS_A.csv')
dmbo_alkyne    <- load_spectrum('data/wavescan/4C-DMBO/4C_WS_B.csv')
dmbo_triazole  <- load_spectrum('data/wavescan/4C-DMBO/4C_WS_C.csv')

# Plot
quartz(width = 6.4, height = 5)
par(yaxs = 'i')
plot(dmbo_photo$V1, dmbo_photo$V2, type = 'l', lwd = 2,
     col = '#66c2a5', xlim = c(250, 400), ylim = c(0, 1),
     main = '4C-DMBO', xlab = 'Wavelength (nm)', ylab = 'Absorbance (AU)', las = 1)

# Add spectra lines
lines(dmbo_alkyne$V1, dmbo_alkyne$V2, col = '#fc8d62', lwd = 2)
lines(dmbo_triazole$V1, dmbo_triazole$V2, col = '#8da0cb', lwd = 2)

# Mark λ = 309 nm
λ <- 309
mark_wavelength(λ, dmbo_triazole$V2[dmbo_triazole$V1 == λ], '#8da0cb')
mark_wavelength(λ, dmbo_photo$V2[dmbo_photo$V1 == λ], '#66c2a5')
mark_wavelength(λ, dmbo_alkyne$V2[dmbo_alkyne$V1 == λ], '#fc8d62')

# Save figure
quartz.save(file = "figures/pdf/4C-DMBO-spectra.pdf", type = 'pdf')


# --- Plot: DBCO-Amine Series ---
# Load data
dbco_alkyne  <- load_spectrum('data/wavescan/DBCO-amine/DBCO_WS_D.csv')
dbco_triazole <- load_spectrum('data/wavescan/DBCO-amine/DBCO_WS_E.csv')

# Plot
quartz(width = 6.4, height = 5)
par(yaxs = 'i')
plot(dbco_alkyne$V1, dbco_alkyne$V2, type = 'l', lwd = 2,
     col = '#fc8d62', xlim = c(250, 400),
     main = 'ADIBO', xlab = 'Wavelength (nm)', ylab = 'Absorbance (AU)', las = 1)

# Add spectra line
lines(dbco_triazole$V1, dbco_triazole$V2, col = '#8da0cb', lwd = 2)

# Mark λ = 307 nm
mark_wavelength(λ, dbco_alkyne$V2[dbco_alkyne$V1 == λ], '#fc8d62')
mark_wavelength(λ, dbco_triazole$V2[dbco_triazole$V1 == λ], '#8da0cb')

# Save figure
quartz.save(file = "figures/pdf/DBCO-amine-spectra.pdf", type = 'pdf')
