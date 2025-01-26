CV_Inhibition <- function(file_name, L) {
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
  Ka <- data[1, 2]
  Kb <- data[1, 3]
  a <- data[1, 4]
  IC50 <- data[1, 5]
  At_IC50 <- data[1, 6]
  Bt_IC50 <- data[1, 7]
  At_setup <- data[2:nrow(data), 1]
  Bt_setup <- data[2:nrow(data), 2]
  It_setup <- data[2:nrow(data), 3]
  V0 <- data[2:nrow(data), 4]
  
  X_setup <- cbind(St_setup, It_setup)
  C <- c(Vmax, Ka, Kb, a)
  IC50s <- cbind(At_IC50, Bt_IC50, IC50)
  
  # Cross-validation to select regularization constant
  cv_value <- numeric(length(L))
  
  for (i in seq_along(L)) {
    r <- L[i]
    cv_value[i] <- CV_error(X_setup, V0, C, IC50s, r)
  }
  
  best_r_idx <- which.min(cv_value)
  best_r <- L[best_r_idx]
  
  return(best_r)
}

# Inhibition model
Inhibition <- function(K, X, C) {
  v1 <- C[3] * (1 + X[, 3] / K[1]) + X[, 2] * (1 + X[, 3] / K[2])
  v2 <- X[, 2] + C[4] * C[3]
  v <- C[1] * X[, 1] / (C[4] * C[2] * v1 + X[, 1] * v2)
  return(v)
}

# Cheng-Prusoff equation
Cheng_Prusoff <- function(K, X, C) {
  r1 <- C[4] * C[2] * C[3] / (C[4] * C[2] * (C[3] + X[2]) + X[1] * (X[2] + C[4] + C[3]))
  r2 <- C[4] * C[2] * X[2] / (C[4] * C[2] * (C[3] + X[2]) + X[1] * (X[2] + C[4] + C[3]))
  v <- K[1]*K[2] / (r1 * K[2] + r2 * K[1])
  return(v)
}

# Loss function with lambda
CV_loss <- function(K, X, Y, C, IC50s, lambda) {
  Y_predict <- Inhibition(K, X, C)
  loss <- mean(((Y - Y_predict) / Y) ^ 2) + 
    lambda * mean(((IC50s[, 3] - Cheng_Prusoff(K, cbind(IC50s[, 1], IC50s[, 2]), C)) / IC50s[, 3]) ^ 2)
  return(loss)
}

# Fitting
Fit_inhibition <- function(X, Y, C, IC50s, lambda) {
  K0 <- rep(max(IC50s[, 2]), 2)
  objFun <- function(K) CV_loss(K, X, Y, C, IC50s, lambda)
  
  params <- optim(K0, objFun, method = "Nelder-Mead", control = list(fnscale = 1))$par
  return(params)
}

# Test error
Test_error <- function(Xtrain, Ytrain, Xtest, Ytest, C, IC50s, lambda) {
  # Fit the model
  params <- Fit_inhibition(Xtrain, Ytrain, C, IC50s, lambda)
  
  # Predict
  Ypredict <- Inhibition(params, Xtest, C)
  
  # Calculate loss
  loss <- mean((Ytest - Ypredict)^2)
  
  return(loss)
}

# Cross-validation error
CV_error <- function(X, Y, C, IC50s, lambda) {
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
    loss <- Test_error(Xtrain, Ytrain, Xtest, Ytest, C, IC50s, lambda)
    loss_set[i] <- loss
  }
  
  meanLoss <- mean(loss_set)
  return(meanLoss)
}