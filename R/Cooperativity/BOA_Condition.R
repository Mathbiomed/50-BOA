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
  Km <- data[1, 2]
  IC50 <- data[1, 3]
  St_setup <- data[2:nrow(data), 1]
  It_setup <- data[2:nrow(data), 2]
  
  # 1. Check whether It >= IC50
  reject_1 <- ""
  check_1 <- FALSE
  check_It <- any(It_setup < IC50)
  
  if (check_It) {
    reject_1 <- "Inhibitor concentration < IC50"
    check_1 <- TRUE
  }
  
  # 2. Check whether St varies sufficiently
  reject_2 <- ""
  check_2 <- FALSE
  
  # Check whether St does not vary
  matrix_vary <- St_setup - St_setup[1]
  check_vary <- all(matrix_vary == 0)
  
  # Check whether St ranges sufficiently over 0.2Km ~ 5Km
  check_sufficient <- min(St_setup) > 0.2 * Km || max(St_setup) < 5 * Km
  
  if (check_vary || check_sufficient) {
    reject_2 <- "Substrate concentration should vary"
    check_2 <- TRUE
  }
  
  # Output
  if (check_1 || check_2) {
    cat("Estimation may be insufficient for precise results:", "\n")
    cat(reject_1, "\n")
    cat(reject_2, "\n")
  }
}