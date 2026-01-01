# implementation of Ando, T. and J. Bai (2017). Clustering huge number of financial time series: A panel data approach with high-dimensional predictors and factor structures. Journal of the American Statistical Association 112 (519), 1182–1198. 3, 5, 6



update_common_components <- function(Y, Fs, Gammas, k_0){
  S <- length(Fs)
  Y_residual_study <- list()
  for(s in 1:S){
    Y_residual_study[[s]] <- Y[[s]] - tcrossprod(Fs[[s]], Gammas[[s]])
  }
  
  Y_residual_joint <- do.call(
    abind, c(Y_residual_study, along = 1)
  )
  
  n <- nrow(Y_residual_joint)
  
  s_Y_residual_joint <- safe_svd(Y_residual_joint)
  M_joint_new <- sqrt(n) * s_Y_residual_joint$u[,1:k_0]
  Lambda_c_new <- 1/sqrt(n) * s_Y_residual_joint$v[,1:k_0] %*% diag(s_Y_residual_joint$d[1:k_0])
  
  return(list(M_joint_new = M_joint_new, Lambda_c_new = Lambda_c_new))
}

update_specific_components <- function(Y_res_s, q_s){
  n_s <- nrow(Y_res_s)
  s_Y_res_s <- safe_svd(Y_res_s)
  F_new <- sqrt(n_s) * s_Y_res_s$u[,1:q_s]
  Gamma_new <- 1/sqrt(n_s) * s_Y_res_s$v[,1:q_s] %*% diag(s_Y_res_s$d[1:q_s])
  
  return(list(F_new = F_new, Gamma_new = Gamma_new))
  
}


compute_ando_bai_estimates <- function(Y, k_0, q_s, tol=0.01){
  
  ns <- sapply(Y, function(x) nrow(x))
  S <- length(Y)
  # inits
  Gammas <- list()
  Gammas_new <- list()
  Fs <- list()
  Fs_new <- list()
  Y_joint <-  do.call(abind, c(Y, along = 1))
  n <- nrow(Y_joint)
  s_Y <- safe_svd(Y_joint)
  M_joint_new <- sqrt(n) * s_Y$u[,1:k_0]
  Lambda_c <- 1/sqrt(n) *  s_Y$v[,1:k_0] %*% diag(s_Y$d[1:k_0])
  Ms_new <- compute_Ms(M_joint_new, ns, S)
  for(s in 1:S){
    Y_res_s <- Y[[s]] - Ms_new[[s]] %*% t(Lambda_c)
    Est_s <- update_specific_components(Y_res_s, q_s[s])
    Fs[[s]] <- Est_s$F_new; Gammas[[s]] <- Est_s$Gamma_new
  }

  # loop until convergence
  error <- 1
  #tol <- 0.001
  it <- 0
  while(error > tol){
    it <- it + 1
    # update common components
    Est_c <- update_common_components(Y, Fs, Gammas, k_0)
    M_joint_new <- Est_c$M_joint_new; Lambda_c_new <- Est_c$Lambda_c_new
    Ms_new <- compute_Ms(M_joint_new, ns, S)
    error_new <- sqrt(sum((tcrossprod(Lambda_c_new) - tcrossprod(Lambda_c))^2) / sum((tcrossprod(Lambda_c))^2))
    
    # update study-specific components
    for(s in 1:S){
      Y_res_s <- Y[[s]] - Ms_new[[s]] %*% t(Lambda_c_new)
      Est_s <- update_specific_components(Y_res_s, q_s[s])
      Fs_new[[s]] <- Est_s$F_new; Gammas_new[[s]] <- Est_s$Gamma_new
      error_new <- error_new + sqrt(sum((tcrossprod(Gammas_new[[s]]) - tcrossprod(Gammas[[s]]))^2) / sum((tcrossprod(Gammas[[s]]))^2))
    }
    error <- error_new
    Gammas <- Gammas_new
    Lambda_c <- Lambda_c_new
    Ms <- Ms_new 
    Fs <- Fs_new 
    
  }
  print(paste0('Optimizer has converged in ', it, ' iterations; relative error : ', error))
  
  return(list( 
    Gammas = Gammas_new,
    Lambda_c = Lambda_c_new,
    Ms = Ms_new,
    Fs = Fs_new,
    it = it)
    )
}


safe_svd <- function(x, nu = min(dim(x)), nv = min(dim(x)),
                     fallback = c("eigen",  "gesvd", "irlba"),
                     ...) {
  fallback <- match.arg(fallback)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  
  # 1) Try the default (fast) SVD first
  out <- tryCatch(
    svd(x, nu = nu, nv = nv, ...),
    error = function(e) e
  )
  if (!inherits(out, "error")) {
    out$method <- "svd_default"   # typically LAPACK dgesdd under the hood
    return(out)
  }
  
  # 2) Robust fallback
  if (fallback == "gesvd") {
    out2 <- tryCatch(
      La.svd(x, nu = nu, nv = nv, method = "gesvd"),
      error = function(e) e
    )
    if (!inherits(out2, "error")) {
      out2$method <- "gesvd"
      return(out2)
    }
    stop("safe_svd: svd() failed, and gesvd fallback failed too.\n",
         "svd() error:   ", conditionMessage(out), "\n",
         "gesvd error:   ", conditionMessage(out2))
  }
  
  if (fallback == "irlba") {
    if (!requireNamespace("irlba", quietly = TRUE)) {
      stop("safe_svd: fallback='irlba' but package 'irlba' is not installed.")
    }
    k <- min(nu, nv, min(dim(x)))
    out2 <- tryCatch(
      irlba::irlba(x, nu = min(nu, k), nv = min(nv, k)),
      error = function(e) e
    )
    if (!inherits(out2, "error")) {
      out2$method <- "irlba"
      return(out2)
    }
    stop("safe_svd: svd() failed, and irlba fallback failed too.\n",
         "svd() error:   ", conditionMessage(out), "\n",
         "irlba error:   ", conditionMessage(out2))
  }
  
  if (fallback == "eigen") {
    m <- nrow(x); n <- ncol(x)
    
    out2 <- tryCatch({
      if (m >= n) {
        ee <- eigen(crossprod(x), symmetric = TRUE)
        d  <- sqrt(pmax(ee$values, 0))
        ord <- order(d, decreasing = TRUE)
        d <- d[ord]; V <- ee$vectors[, ord, drop = FALSE]
        keep <- seq_len(min(nv, length(d)))
        d <- d[keep]; V <- V[, keep, drop = FALSE]
        U <- x %*% V
        U <- sweep(U, 2, d, "/")
        list(d = d,
             u = U[, seq_len(min(nu, ncol(U))), drop = FALSE],
             v = V)
      } else {
        ee <- eigen(tcrossprod(x), symmetric = TRUE)
        d  <- sqrt(pmax(ee$values, 0))
        ord <- order(d, decreasing = TRUE)
        d <- d[ord]; U <- ee$vectors[, ord, drop = FALSE]
        keep <- seq_len(min(nu, length(d)))
        d <- d[keep]; U <- U[, keep, drop = FALSE]
        V <- t(x) %*% U
        V <- sweep(V, 2, d, "/")
        list(d = d,
             u = U,
             v = V[, seq_len(min(nv, ncol(V))), drop = FALSE])
      }
    }, error = function(e) e)
    
    if (!inherits(out2, "error")) {
      out2$method <- "eigen"
      return(out2)
    }
    stop("safe_svd: svd() failed, and eigen fallback failed too.\n",
         "svd() error:   ", conditionMessage(out), "\n",
         "eigen error:   ", conditionMessage(out2))
  }
}

# Usage:
# res <- safe_svd(X, nu = 10, nv = 10, fallback = "gesvd")
# res <- safe_svd(X, nu = 50, nv = 50, fallback = "irlba")  # truncated, needs irlba
