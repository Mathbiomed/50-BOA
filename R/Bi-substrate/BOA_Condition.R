BOA_Condition <- function(file_name) {
  # Load the data
  library(readxl)
  data <- read_excel(file_name, col_names = c("A", "B", "C", "D"))
  names(data) <- NULL
  data <- as.matrix(data)
  
  # Data check
  if (nrow(data) < 2 || ncol(data) < 4) {
    stop("Invalid Input")
  }
  
  # Load the formatted data
  Ka <- data[1, 2]
  Kb <- data[1, 3]
  IC50 <- data[1, 4]
  At_setup <- data[2:nrow(data), 1]
  Bt_setup <- data[2:nrow(data), 2]
  It_setup <- data[2:nrow(data), 3]
  
  # 1. Check whether It >= IC50
  reject_1 <- ""
  check_1 <- FALSE
  check_It <- any(It_setup < IC50)
  
  if (check_It) {
    reject_1 <- "Inhibitor concentration < IC50"
    check_1 <- TRUE
  }
  
  # 2. Check whether At varies sufficiently
  reject_2 <- ""
  check_2 <- FALSE
  
  # Check whether St does not vary
  matrix_vary_A <- At_setup - At_setup[1]
  check_vary_A <- all(matrix_vary_A == 0)
  
  # Check whether At ranges sufficiently over 0.2Ka ~ 5Ka
  check_sufficient_A <- min(At_setup) > 0.2 * Ka || max(At_setup) < 5 * Ka
  
  if (check_vary_A || check_sufficient_A) {
    reject_2 <- "Substrate concentration should vary"
    check_2 <- TRUE
  }
  
  # 3. Check whether Bt varies sufficiently
  reject_3 <- ""
  check_3 <- FALSE
  
  # Check whether Bt does not vary
  matrix_vary_B <- Bt_setup - Bt_setup[1]
  check_vary_B <- all(matrix_vary == 0)
  
  # Check whether Bt ranges sufficiently over 0.2Kb ~ 5Kb
  check_sufficient_B <- min(Bt_setup) > 0.2 * Kb || max(Bt_setup) < 5 * Kb
  
  if (check_vary_B || check_sufficient_B) {
    reject_3 <- "Substrate concentration should vary"
    check_2 <- TRUE
  }
  
  # Output
  if (check_1 || check_2 || check_3) {
    cat("Estimation may be insufficient for precise results:", "\n")
    cat(reject_1, "\n")
    cat(reject_2, "\n")
    cat(reject_3, "\n")
  }
}