#' Perform Sparse Subspace Non-negative Matrix Factorization using PyTorch
#'
#' @param X Input data matrix
#' @param k Number of components
#' @param Y Auxiliary matrix (optional)
#' @param A Initial A matrix (optional)
#' @param S Initial S matrix (optional)
#' @param N Number of iterations
#' @param lam Regularization parameter
#' @param mult Use multiplicative updates
#'
#' @return List containing the factorized matrices A and S
#' @export
#'
#' @examples
#' np <- import("numpy")
#' reticulate::use_condaenv("SSLAR-env")
#' # Example data
#' l <- as.integer(100)
#' h <- as.integer(50)
#' X <- np$random$rand(l, h)
#' k <- as.integer(10)
#' A <- np$random$rand(l, k)
#' S <- np$random$rand(k, h)
#'
#' # Perform SSNMF
#' nmf_result <- torch_ssNMF(X, k, Y = NULL, A = A, S = S, N = 100, lam = 100, mult = TRUE)
#'
#' # Print results
#' print("A matrix:")
#' print(nmf_result$A)
#' print("S matrix:")
#' print(nmf_result$S)
torch_ssNMF <- function(X, k, Y = NULL, A = NULL, S = NULL, N = 100, lam = 100, mult = TRUE) {
  Sys.setenv(PYTHONWARNINGS="ignore")
  Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
  Sys.setenv(CUDA_VISIBLE_DEVICES = "-1")

  ssnmf <- reticulate::import("ssnmf")
  torch <- reticulate::import("torch")

  device <- torch$device("cpu")
  # browser()
  # Convert input matrices to tensors
  Xt <- .getTensor(X, device)
  if (!is.null(Y)) {
    Yt <- .getTensor(Y, device)
  } else {
    Yt <- NULL
  }
  if (!is.null(A)) {
    At <- .getTensor(A, device)
  } else {
    At <- NULL
  }
  if (!is.null(S)) {
    St <- .getTensor(S, device)
  } else {
    St <- NULL
  }

  # Initialize SSNMF model
  nmf_mod <- ssnmf$SSNMF(Xt, as.integer(k), Y = Yt, lam = lam * torch$norm(Xt), A = At, S = St, modelNum = 6)

  if (mult) {
    # Perform multiplicative updates
    nmf_mod$mult(numiters = as.integer(N), saveerrs = TRUE)
  }

  # Convert results back to numpy arrays
  nmf_mod$A <- nmf_mod$A$numpy()
  nmf_mod$S <- nmf_mod$S$numpy()

  return(nmf_mod)
}

#' Convert a matrix to a PyTorch tensor
#'
#' @param m Matrix to convert
#' @param device Device to place the tensor on (e.g., "cpu" or "cuda")
#'
#' @return PyTorch tensor
#' @export
#'
#' @examples
#' Xt <- .getTensor(X, device)
.getTensor <- function(m, device) {
  torch <- reticulate::import("torch")
  copy <- reticulate::import("copy")
  np <- reticulate::import("numpy")
  mt <- suppressWarnings(torch$from_numpy(np$array(copy$deepcopy(m))))
  mt <- mt$type(torch$FloatTensor)
  mt <- mt$to(device)
  return(mt)
}

# 示例用法
# if (interactive()) {
#   np <- import("numpy")
#   reticulate::use_condaenv("SSLAR-env")
#   # 示例数据
#   l = as.integer(100)
#   h = as.integer(50)
#   X <- np$random$rand(l, h)
#   k <- as.integer(10)
#   A <- np$random$rand(l, k)
#   S <- np$random$rand(k, h)
#
#   # 执行 SSNMF
#   nmf_result <- torch_ssNMF(X, k, S, A, S, N = 100, lam = 100, mult = TRUE)
#
#   # 打印结果
#   # print("A matrix:")
#   # print(nmf_result$A)
#   # print("S matrix:")
#   # print(nmf_result$S)
# }
