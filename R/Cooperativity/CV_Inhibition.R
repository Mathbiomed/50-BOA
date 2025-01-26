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
  Ks <- data[1, 2]
  IC50 <- data[1, 3]
  St_IC50 <- data[1, 4]
  St_setup <- data[2:nrow(data), 1]
  It_setup <- data[2:nrow(data), 2]
  V0 <- data[2:nrow(data), 3]
  a <- data[2, 4];
  
  X_setup <- cbind(St_setup, It_setup)
  C <- c(Vmax, Ks, a)
  IC50s <- cbind(St_IC50, IC50)
  
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
  v <- C[1] * (C[2] * X[, 1] * (1 + X[, 2] / K[2]) + 
                 (1 + X[, 2] / (K[2]^2 / K[1])) * X[, 1]^2 / C[3]) / 
    (C[2]^2 * (1 + X[, 2] / K[1]) + 
       2 * C[2] * X[, 1] * (1 + X[, 2] / K[2]) + 
       (1 + X[, 2] / (K[2]^2 / K[1])) * X[, 1]^2 / C[3])
  return(v)
}

# Cheng-Prusoff equation
Cheng_Prusoff <- function(K, X, C) {
  h1 <- function(S, a, x) {
    a * S * x + S^2
  }
  h2 <- function(S, a, x) {
    a * x^2 + 2 * a * S * x + S^2
  }
  numerator <- K[1] * (K[2] / K[1])^2 * (h1(X, C[2], C[1]) * h2(X, C[2], C[1]))
  denominator <- (h1(X, C[2], C[1]) * h2(X, C[2], C[1] * K[2] / K[1]) - 
                    2 * h1(X, C[2], C[1] * K[2] / K[1]) * h2(X, C[2], C[1]))
  v <- numerator / denominator 
  return(v)
}

# Loss function with lambda
CV_loss <- function(K, X, Y, C, IC50s, lambda) {
  Y_predict <- Inhibition(K, X, C)
  loss <- mean(((Y - Y_predict) / Y) ^ 2) + 
    lambda * mean(((IC50s[, 2] - Cheng_Prusoff(K, IC50s[, 1], cbind(C[2], C[3]))) / IC50s[, 2]) ^ 2)
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