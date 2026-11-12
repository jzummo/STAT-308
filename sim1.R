folder_path <- "~/Downloads"

for (b in 1:10) {
  file_name <- file.path(folder_path, paste0("my_plot_", b, ".png"))
  png(filename = file_name, width = 6, height = 4, units = "in", res = 300)
  
  results = c()
  x <- c()
  for (n in 1:300) {
    x <- append(x, rnorm(1))
    
    hatGamma <- acf(x, lag.max = n - 1, type = "covariance", plot = F)$acf[, 1, 1]
    hatGamma <- toeplitz(hatGamma)
    
    Gamma = diag(n)
    
    results <- append(results, norm(hatGamma - Gamma, type = "2"))
  }
  plot(results)
  
  dev.off()
}

