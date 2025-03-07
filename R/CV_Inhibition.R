CV_Inhibition <- function(file_name, L, isSSRE = TRUE) {
  # Load data
  library(readxl)
  data <- read_excel(file_name, col_names = c("A", "B", "C", "D"))
  names(data) <- NULL
  data <- as.matrix(data)
  
  # Data check
  if (nrow(data) < 2 || ncol(data) < 4) {
    stop("Invalid Input")
  }
  
  # Load the formatted data
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
  
  # Cross-validation to select regularization constant
  cv_value <- numeric(length(L))
  
  for (i in seq_along(L)) {
    r <- L[i]
    cv_value[i] <- CV_error(X_setup, V0, C, IC50s, r, isSSRE)
  }
  
  best_r_idx <- which.min(cv_value)
  best_r <- L[best_r_idx]
  
  return(best_r)
}

# Inhibition model
Inhibition <- function(K, X, C) {
  v <- C[1] * X[, 1] / (C[2] * (1 + X[, 2] / K[1]) + X[, 1] * (1 + X[, 2] / K[2]))
  return(v)
}

# Cheng-Prusoff equation
Cheng_Prusoff <- function(K, X, C) {
  v <- (X + C) * K[1] * K[2] / (C * K[2] + X * K[1])
  return(v)
}

# Error Structure
Error_Structure <- function(X, Y, isSSRE) {
  if (isSSRE) {
    s <- Y
  } else {
    # User can assign his or her own model of standard deviation to "s" here.
    
    # Model 1
    # K1 <- 0.011
    # K2 <- 0.04
    # s <- K1 + K2 * Y
    
    # Model 2
    # K1 <- 0.036
    # K2 <- 0.008
    # s <- K1 + K2 * Y ^ 2
    
    # Model 3
    # K1 <- 0.05
    # K2 <- 0.889
    # s <- K1 * Y ^ K2
    
    # Model 4
    K1 <- 0.259
    K2 <- 0.529
    K3 <- 0.699
    K4 <- 0.561
    K5 <- 0.505
    s <- K1 * X[, 1] ^ K2 / (X[, 1] ^ K3 + K4 * X[, 2] ^ K5)
    
    # Model 5
    # K1 <- 0.308
    # K2 <- 0.063
    # K3 <- 0.273
    # K4 <- 0.095
    # s <- K1 * X[, 1] / (K2 + X[, 1] + K3 * X[, 1] ^ 2 + K4 * X[, 2])
    
    # Model 6
    # K1 <- 0.308
    # K2 <- 0.063
    # K3 <- 0.273
    # K4 <- 0.095
    # K5 <- 0.007
    # s <- (K1 * X[, 1] + K5 * X[, 2]) / (K2 + X[, 1] + K3 * X[, 1] ^ 2 + K4 * X[, 2])
  }
  
  return(s)
}

# Loss function with lambda
CV_loss <- function(K, X, Y, C, IC50s, lambda, isSSRE) {
  Y_predict <- Inhibition(K, X, C)
  loss <- mean(((Y - Y_predict) / Error_Structure(X, Y, isSSRE)) ^ 2) + 
    lambda * mean(((IC50s[, 2] - Cheng_Prusoff(K, IC50s[, 1], C[2])) / IC50s[, 2]) ^ 2)
  return(loss)
}

# Fitting
Fit_inhibition <- function(X, Y, C, IC50s, lambda, isSSRE) {
  K0 <- rep(max(IC50s[, 2]), 2)
  objFun <- function(K) CV_loss(K, X, Y, C, IC50s, lambda, isSSRE)
  
  params <- optim(K0, objFun, method = "Nelder-Mead", control = list(fnscale = 1))$par
  return(params)
}

# Test error
Test_error <- function(Xtrain, Ytrain, Xtest, Ytest, C, IC50s, lambda, isSSRE) {
  # Fit the model
  params <- Fit_inhibition(Xtrain, Ytrain, C, IC50s, lambda, isSSRE)
  
  # Predict
  Ypredict <- Inhibition(params, Xtest, C)
  
  # Calculate loss
  loss <- mean((Ytest - Ypredict)^2)
  
  return(loss)
}

# Cross-validation error
CV_error <- function(X, Y, C, IC50s, lambda, isSSRE) {
  #Leave-one-out cross-validation
  n <- nrow(X)
  loss_set <- numeric(n)
  
  for (i in 1:n) {
    # Split data
    trainIdx <- setdiff(1:n, i)
    testIdx <- i
    
    Xtrain <- X[trainIdx, , drop = FALSE]
    Ytrain <- Y[trainIdx]
    Xtest <- X[testIdx, , drop = FALSE]
    Ytest <- Y[testIdx]
    
    # Train & Test
    loss <- Test_error(Xtrain, Ytrain, Xtest, Ytest, C, IC50s, lambda, isSSRE)
    loss_set[i] <- loss
  }
  
  meanLoss <- mean(loss_set)
  return(meanLoss)
}