Error_Landscape <- function(file_name, L = 10^seq(-3, 3, length.out = 100), isReg = TRUE, isSSRE = TRUE, isMatched = TRUE) {
  # Load data
  library(readxl)
  data <- read_excel(file_name, col_names = c("A", "B", "C", "D"))
  names(data) <- NULL
  data <- as.matrix(data)
  
  # 1. BOA_Condition
  source("BOA_Condition.R")
  BOA_Condition(file_name)
  
  
  # 2. CV_Inhibition
  source("CV_Inhibition.R")
  if (isReg) {
    lambda <- CV_Inhibition(file_name, L, isSSRE)
  } else {
    lambda <- 0
  }
  cat(sprintf('The regularization constant is %.3f.\n', lambda))
  
  # 3. Estimate Kic and Kiu
  Vmax <- data[1, 1]
  Km <- data[1, 2]
  IC50 <- data[1, 3]
  St_IC50 <- data[1, 4]
  St_setup <- data[2:nrow(data), 1]
  It_setup <- data[2:nrow(data), 2]
  V0 <- data[2:nrow(data), 3]
  
  X_setup <- cbind(St_setup, It_setup)
  C <- c(Vmax, Km)
  IC50s <- cbind(St_IC50, IC50)
  
  # Estimation
  K <- Fit_inhibition(X_setup, V0, C, IC50s, lambda, isSSRE)
  
  # Compute 95% confidence interval via bootstrapping
  bootSample <- 1000
  betaSample <- matrix(0, bootSample, 2)
  
  for (i in 1:bootSample) {
    # Resampling
    indices <- sample(1:length(V0), length(V0), replace = TRUE)
    X_boot <- X_setup[indices, ]
    V0_boot <- V0[indices]
    
    # Re-estimate parameter using bootstrap sample
    betaSample[i, ] <- Fit_inhibition(X_boot, V0_boot, C, IC50s, lambda)
  }
  
  # Compute confidence interval
  CI_lower <- apply(betaSample, 2, function(x) quantile(x, 0.025))
  CI_upper <- apply(betaSample, 2, function(x) quantile(x, 0.975))
  CI <- cbind(CI_lower, CI_upper)
  
  # Estimates of Kic and Kiu
  Estimate <- cbind(K, CI)
  cat(sprintf('Kic: %.4f, (%.4f, %.4f)\n', K[1], CI_lower[1], CI_upper[1]))
  cat(sprintf('Kiu: %.4f, (%.4f, %.4f)\n', K[2], CI_lower[2], CI_upper[2]))
  
  # Generate error landscape
  K_round <- round(log10(IC50))
  Kic_min <- K_round - 2
  Kic_max <- K_round + 2
  Kiu_min <- K_round - 2
  Kiu_max <- K_round + 2
  
  Kicrange <- 10^seq(Kic_min, Kic_max, length.out = 100)
  Kiurange <- 10^seq(Kiu_min, Kiu_max, length.out = 100)
  
  S1 <- numeric(length(Kicrange) * length(Kiurange))
  S2 <- numeric(length(Kicrange) * length(Kiurange))
  total_error <- numeric(length(Kicrange) * length(Kiurange))
  
  idx <- 1
  for (i in 1:length(Kicrange)) {
    for (j in 1:length(Kiurange)) {
      Kicr <- Kicrange[i]
      Kiur <- Kiurange[j]
      K <- c(Kicr, Kiur)
      total_error[idx] <- CV_loss(K, X_setup, V0, C, IC50s, lambda)
      S1[idx] <- Kicr
      S2[idx] <- Kiur
      idx <- idx + 1
    }
  }
  
  error_table <- data.frame(Kic = S1, Kiu = S2, Error = total_error)
  min_error <- min(total_error)
  
  # Generate heatmap
  library(ggplot2)
  library(scales)

  if (isMatched) {
    ggplot(error_table, aes(x = log10(Kic), y = log10(Kiu), fill = Error)) +
    geom_tile() +
    scale_fill_gradientn(colors = heat.colors(10), limits = c(min_error, 2 * min_error), oob = squish) +
    labs(x = expression(log[K[ic]]), y = expression(log[K[iu]]), fill = "Total error") +
    theme_minimal(base_size = 20)
  } else {
    ggplot(error_table, aes(x = log10(Kic), y = log10(Kiu), fill = Error)) +
    geom_tile() +
    scale_fill_gradientn(colors = heat.colors(10), limits = c(0, 0.05), oob = squish) +
    labs(x = expression(log[K[ic]]), y = expression(log[K[iu]]), fill = "Total error") +
    theme_minimal(base_size = 20)
  }
  
}
