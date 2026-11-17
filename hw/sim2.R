n <- 1000
phi <- 0.5

acvs <- phi**(0:(n-1)) * (1/(1 - phi**2))
  
Gamma <- toeplitz(acvs)

eigenvalues <- eigen(Gamma)$values
print(range(eigenvalues))
# [1] 0.4444454 3.9999215

spectral_density <- function(omega) {
  1 / (1 - 2 * phi * cos(omega) + phi**2)
}

curve(spectral_density, from = 0, to = 2 * pi)
print(range(curve(spectral_density, from = 0, to = 2 * pi)$y))
# [1] 0.4444444 4.0000000
