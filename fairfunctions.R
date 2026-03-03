library(quantreg)
library(splines2)  
library(dfoptim)
library(Iso)   
library(isotone)
library(rqPen)
library(gmp)
library(CVXR)
library(mvtnorm)


# Mean squared error (MSE)
MSE<- function(YT, YF) {
  err <- sum((YT - YF)^2)/length(YT)
  return(err)
}

# Pinball loss for quantile regression
quantile_loss <- function(y, y_hat, tau) {
  loss <- ifelse(y < y_hat, (tau - 1) * (y - y_hat), tau * (y - y_hat))
  mean(loss)
}

# KS statistic
f_ks <- function(Y, S) {
  ts <- sort(unique(Y))
  
  # Get unique classes
  classes <- unique(S)
  num_classes <- length(classes)
  
  # Initialize the maximum KS statistic
  max_fai <- 0
  
  # Loop through all pairs of classes
  for (i in 1:(num_classes - 1)) {
    for (j in (i + 1):num_classes) {
      class_i <- classes[i]
      class_j <- classes[j]
      
      # Get the counts for each class
      nn_i <- sum(S == class_i)
      nn_j <- sum(S == class_j)
      
      # Compute KS statistic for this pair of classes
      fai <- 0
      for (t in ts) {
        tmp1 <- sum(Y[S == class_i] <= t) / nn_i
        tmp2 <- sum(Y[S == class_j] <= t) / nn_j
        fai <- max(fai, abs(tmp1 - tmp2))
      }
      
      # Update the maximum KS statistic
      max_fai <- max(max_fai, fai)
    }
  }
  
  return(max_fai)
}

# Poisson solver for PAVA method
poisson_solver = function(Y, w){ 
  
  return(sum(Y*w)/sum(w))  
  
}

# Poisson negative log-likelihood (up to an additive constant)
poissonloss = function(Y, expfhat){ 
  
  ans = expfhat - Y* log(expfhat)
  return( mean (ans) ) 
  
} 

# Huber loss 
huber_loss <- function(r, delta) {ifelse(abs(r) <= delta,
                                         0.5 * r^2,
                                         delta * abs(r) - 0.5 * delta^2)}

# Calculate the weighted Huber loss
L_obj <- function(m, y, w, M) sum(w * huber_loss(y - m, M))

# Huber solver for PAVA method
huber_solver <- function(y, w = NULL, M = delta) {
  if (is.null(w)) w <- rep(1, length(y))
  rng <- range(y) + c(-10*M, 10*M)       # wide bracket
  res <- optimize(L_obj, interval = rng, y = y, w = w, M = M)
  res$minimum
}

# Derivative of Huber loss
huber_psi <- function(r, delta) pmax(pmin(r, delta), -delta) 

# Huber regression via BFGS (unconstrained coefficients)
huber_fit_delta <- function(y, X, delta, intercept = TRUE) {
  X <- as.matrix(X)
  if (intercept && (ncol(X) == 0 || any(X[,1] != 1))) X <- cbind(1, X)
  p <- ncol(X)
  
  obj <- function(b) sum(huber_loss(y - drop(X %*% b), delta))
  grad <- function(b) {
    r <- y - drop(X %*% b)
    - crossprod(X, huber_psi(r, delta))
  }
  
  b0 <- rep(0, p)                  
  sol <- optim(b0, obj, grad, method = "BFGS", control = list(reltol = 1e-10))
  list(coef = sol$par, convergence = sol$convergence)
}

# Quantile regression with nonnegative spline coefficients via linear programming
fit_ispline_model_Quantile_LP <- function(rank_train, Y_train, tau = 0.75,
                                          knots = c(0.25, 0.5, 0.75), degree = 2) {
  rank_train <- as.numeric(rank_train)
  Y_train <- as.numeric(Y_train)
  if (length(rank_train) != length(Y_train)) stop("rank_train and Y_train have different length")
  if (anyNA(rank_train) || anyNA(Y_train)) stop("NA is found in the input")
  if (!(tau > 0 && tau < 1)) stop("tau must be in (0,1)")
  
  B <- as.matrix(iSpline(rank_train, knots = knots, Boundary.knots = c(0,1), degree = degree))
  X <- cbind(B, 1)  # intercept
  p_spline <- ncol(B)
  
  # R %*% beta >= 0
  R <- cbind(diag(p_spline), rep(0, p_spline))  
  r <- rep(0, p_spline)
  
  fit <- rq.fit.fnc(x = X, y = Y_train, tau = tau, R = R, r = r)
  list(coeff = fit$coefficients)
}

# Convex Poisson regression using CVXR (nonnegative spline coefficients)
fit_ispline_model_Poisson_Convex <- function(rank_train, Y_train,
                                             knots = c(0.25, 0.5, 0.75),
                                             degree = 2,
                                             solver = c("SCS", "ECOS")) {
  
  solver <- match.arg(solver)
  
  rank_train <- as.numeric(rank_train)
  Y_train <- as.numeric(Y_train)
  
  if (length(rank_train) != length(Y_train)) stop("rank_train and Y_train have different length")
  if (anyNA(rank_train) || anyNA(Y_train)) stop("there is NA in the input")
  if (any(Y_train < 0)) stop("Poisson requires nonnegative counts (Y_train >= 0)")
  
  B <- as.matrix(splines2::iSpline(rank_train, knots = knots,
                                   Boundary.knots = c(0, 1),
                                   degree = degree))
  X <- cbind(B, 1)
  
  n <- nrow(X)
  p_spline <- ncol(B)
  p <- ncol(X)
  
  beta <- CVXR::Variable(p)
  
  # linear predictor: eta = X %*% beta
  eta <- X %*% beta
  
  constraints <- list(
    # spline coefs >= 0; intercept (last entry) free
    beta[1:p_spline] >= 0   
  )
  
  # Poisson negative log-likelihood (up to constant):
  obj <- CVXR::Minimize(CVXR::sum_entries(exp(eta) - CVXR::multiply(Y_train, eta)))
  
  prob <- CVXR::Problem(obj, constraints)
  res <- CVXR::solve(prob, solver = solver)
  
  if (!(res$status %in% c("optimal", "optimal_inaccurate"))) {
    stop(paste("CVXR solve failed, status =", res$status))
  }
  
  coef_hat <- as.numeric(res$getValue(beta))
  list(coeff = coef_hat)
}

# Convex Huber regression using CVXR 
fit_ispline_model_Robust_Convex <- function(rank_train, Y_train,
                                            knots = c(0.25, 0.5, 0.75),
                                            degree = 2, delta = 1.0,
                                            solver = c("OSQP", "SCS")) {
  
  solver <- match.arg(solver)
  
  rank_train <- as.numeric(rank_train)
  Y_train <- as.numeric(Y_train)
  
  if (length(rank_train) != length(Y_train)) stop("rank_train and Y_train have different length")
  if (anyNA(rank_train) || anyNA(Y_train)) stop("there is NA in the input")
  if (delta <= 0) stop("delta must be greater than 0")
  
  
  B <- as.matrix(splines2::iSpline(rank_train, knots = knots, Boundary.knots = c(0,1), degree = degree))
  X <- cbind(B, 1)
  
  n <- nrow(X)
  p_spline <- ncol(B)
  p <- ncol(X)
  
  beta <- CVXR::Variable(p)
  r <- CVXR::Variable(n)
  
  constraints <- list(
    r == Y_train - X %*% beta,
    beta[1:p_spline] >= 0   # intercept (last entry) unconstrained
  )
  
  obj <- CVXR::Minimize(CVXR::sum_entries(CVXR::huber(r, delta)))
  prob <- CVXR::Problem(obj, constraints)
  
  res <- CVXR::solve(prob, solver = solver)
  
  if (!(res$status %in% c("optimal", "optimal_inaccurate"))) {
    stop(paste("CVXR solve failed, status =", res$status))
  }
  
  coef_hat <- as.numeric(res$getValue(beta))
  
  list(coeff = coef_hat)
}

# Predict using a fitted I-spline model at new rank values
ispline_predict <- function(rank_new, coeff, 
                            knots = c(0.25, 0.5, 0.75), degree = 2) {
  new_basis <- iSpline(rank_new, knots = knots,  Boundary.knots = c(0,1), degree = degree)
  as.vector(cbind(new_basis, 1) %*% coeff) # To do the prediction
}

# Fit a non-decreasing regression function using generalized PAVA under the pinball loss.
fit_pava_model_Quantile <- function(rank_train, Y_train, tau = NA, solverinput = weighted.fractile) {
  gpava_result <- gpava(rank_train, Y_train, 
                        solver = solverinput, ties = "secondary", p = tau)
  fitted_values <- gpava_result$x
  
  # Remove duplicates based on training ranks
  unique_idx <- !duplicated(rank_train)
  rank_unique <- rank_train[unique_idx]
  gpava_unique <- fitted_values[unique_idx]
  
  # Order the unique values
  ord <- order(rank_unique)
  rank_unique <- rank_unique[ord]
  gpava_unique <- gpava_unique[ord]
  
  
  # Determine knot positions: indices where the fitted values change
  knots_idx <- which(diff(gpava_unique) != 0) + 1
  y0 <- c(gpava_unique[1], gpava_unique[knots_idx])
  # Build a step function; here we use the index positions as knots.
  model_fun <- stepfun(knots_idx, y0)
  list(rank_unique = rank_unique, model = model_fun)
}

# Fit a non-decreasing regression function using generalized PAVA under the poisson loss.
fit_pava_model_Poisson <- function(rank_train, Y_train, solverinput = poisson_solver) {
  gpava_result <- gpava(rank_train, Y_train, 
                        solver = solverinput, ties = "secondary")
  fitted_values <- gpava_result$x
  
  # Remove duplicates based on training ranks
  unique_idx <- !duplicated(rank_train)
  rank_unique <- rank_train[unique_idx]
  gpava_unique <- fitted_values[unique_idx]
  
  # Order the unique values
  ord <- order(rank_unique)
  rank_unique <- rank_unique[ord]
  gpava_unique <- gpava_unique[ord]
  
  
  # Determine knot positions: indices where the fitted values change
  knots_idx <- which(diff(gpava_unique) != 0) + 1
  if(gpava_unique[1] == 0){
    gpava_unique[1] = unique(gpava_unique)[2]
    
  }
  y0 <- c(gpava_unique[1], gpava_unique[knots_idx])
  # Build a step function; here we use the index positions as knots.
  model_fun <- stepfun(knots_idx, y0)
  list(rank_unique = rank_unique, model = model_fun)
}

# Fit a non-decreasing regression function using generalized PAVA under the Huber loss.
fit_pava_model_Robust <- function(rank_train, Y_train, solverinput = huber_solver) {
  gpava_result <- gpava(rank_train, Y_train, 
                        solver = solverinput, ties = "secondary")
  fitted_values <- gpava_result$x
  
  # Remove duplicates based on training ranks
  unique_idx <- !duplicated(rank_train)
  rank_unique <- rank_train[unique_idx]
  gpava_unique <- fitted_values[unique_idx]
  
  # Order the unique values
  ord <- order(rank_unique)
  rank_unique <- rank_unique[ord]
  gpava_unique <- gpava_unique[ord]
  
  
  # Determine knot positions: indices where the fitted values change
  knots_idx <- which(diff(gpava_unique) != 0) + 1
  y0 <- c(gpava_unique[1], gpava_unique[knots_idx])
  # Build a step function; here we use the index positions as knots.
  model_fun <- stepfun(knots_idx, y0)
  list(rank_unique = rank_unique, model = model_fun)
}

# Predict from a fitted PAVA model (step function over the training ranks)
pava_predict <- function(pava_model, rank_new) {
  sapply(rank_new, function(x) {
    idx <- sum(x >= pava_model$rank_unique)
    pava_model$model(idx)
  })
}

# Helper: generate n_interior equally spaced interior knots in (0, 1)
make_knots <- function(n_interior){
  if(n_interior == 0) return(numeric(0))
  probs <- seq(0, 1, length.out = n_interior + 2)[-c(1, n_interior + 2)]
  return(probs)
}

# Read a CSV file created by appending multiple datasets with repeated headers
read_appended_csv <- function(path) {
  lines <- readLines(path)
  header <- lines[1]
  lines_clean <- c(header, lines[lines != header])  # drop repeated headers
  read.csv(text = paste(lines_clean, collapse = "\n"))
}
